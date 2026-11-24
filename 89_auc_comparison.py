# -*- coding: utf-8 -*-
# 89_auc_comparison.py - AUC vs LN_IC50 的度量对照（阻塞项2）
# 问题：审稿人会问"为什么不直接用 AUC/完整曲线规避删失"。
# 对 7 个正对照同时做 LN_IC50(capped, permutation MWW) 与 AUC(MWW) 检验，
# 对比方向、p 值与召回；并对全 panel 量化 AUC 在"WT 删失"单元上的压缩行为。
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, rankdata

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

print("加载 ...")
gdsc = pd.read_excel(GDSC2, usecols=["DRUG_NAME", "SANGER_MODEL_ID", "LN_IC50", "AUC", "MAX_CONC", "CANCER_TYPE"])
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False,
                  usecols=["Hugo_Symbol", "DepMap_ID", "Variant_Classification", "Protein_Change",
                           "isTCGAhotspot", "isCOSMIChotspot"])

s2d = {}
d2lineage = {}
for _, r in si.iterrows():
    sm_ = r.get('Sanger_Model_ID')
    if pd.notna(sm_):
        s2d[str(sm_).strip()] = r['DepMap_ID']
    d2lineage[r['DepMap_ID']] = r['lineage']

nonsilent = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
             'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)].copy()
def is_true(x):
    return str(x).strip().lower() in ('true', '1', 'yes', 't')
mut_ns['hot'] = mut_ns['isTCGAhotspot'].apply(is_true) | mut_ns['isCOSMIChotspot'].apply(is_true)
mut_genes_ns = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

def cells_with(gene, pattern):
    sub = mut_ns[(mut_ns['Hugo_Symbol'] == gene) &
                 (mut_ns['Protein_Change'].astype(str).str.contains(pattern, na=False, regex=True))]
    return set(sub['DepMap_ID'])

def cells_hotspot(gene):
    return set(mut_ns[(mut_ns['Hugo_Symbol'] == gene) & (mut_ns['hot'])]['DepMap_ID'])

gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage'] = gdsc['depmap'].map(d2lineage)

# (name, drug, gene, pattern_or_NONE_for_hotspot, lineage, direction)
controls = [
    ("BRAF V600E / Dabrafenib / skin", "Dabrafenib", "BRAF", "V600", "skin", "mut_sensitive"),
    ("BRAF V600E / Trametinib / skin", "Trametinib", "BRAF", "V600", "skin", "mut_sensitive"),
    ("TP53 / Nutlin-3a / pan", "Nutlin-3a", "TP53", None, None, "wt_sensitive"),
    ("EGFR L858R+ex19 / Afatinib / lung", "Afatinib", "EGFR", "L858|747_753|E746|ELREA", "lung", "mut_sensitive"),
    ("NRAS Q61 / Trametinib / pan", "Trametinib", "NRAS", "Q61", None, "mut_sensitive"),
    ("PIK3CA H1047R/E545 / Alpelisib / breast", "Alpelisib", "PIK3CA", "H1047|E545|E542", "breast", "mut_sensitive"),
    ("KRAS G12/13 / Trametinib / lung", "Trametinib", "KRAS", "G12|G13", "lung", "mut_sensitive"),
]

N_PERM = 20000
SEED = 20261122

def mww_perm(a, b, direction):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = len(a), len(b)
    ranks = rankdata(np.concatenate([a, b]))
    obs_S = ranks[:na].sum()
    alt = 'less' if direction == 'mut_sensitive' else 'greater'
    U, p_asym = mannwhitneyu(a, b, alternative=alt)
    delta = 2.0 * U / (na * nb) - 1.0
    rng = np.random.default_rng(SEED)
    cnt = 0
    for _ in range(N_PERM):
        perm = rng.permutation(na + nb)
        S = ranks[perm[:na]].sum()
        if direction == 'mut_sensitive':
            cnt += S <= obs_S
        else:
            cnt += S >= obs_S
    p_perm = (cnt + 1.0) / (N_PERM + 1.0)
    return delta, p_asym, p_perm

rows = []
print(f"\n{'对照':38s} {'nM':>4} {'nW':>4} | {'pi_WT':>6} | {'d_cap':>6} {'p_perm':>7} | {'d_AUC':>6} {'p_AUC':>7} | 一致")
for name, drug, gene, pattern, lineage, direction in controls:
    dm = [d for d in gdsc['DRUG_NAME'].unique() if drug.lower() in str(d).lower()]
    if not dm:
        print(f"{name}: drug 不在 GDSC2"); continue
    dsub = gdsc[gdsc['DRUG_NAME'] == dm[0]].copy()
    dsub = dsub[dsub['depmap'].notna() & dsub['LN_IC50'].notna() & dsub['AUC'].notna()]
    if lineage:
        dsub = dsub[dsub['lineage'] == lineage]
    if pattern is None:
        mcells = cells_hotspot(gene)
    else:
        mcells = cells_with(gene, pattern)
    wtcells = set(dsub['depmap']) - {c for c in dsub['depmap'] if gene in mut_genes_ns.get(c, set())}
    msub = dsub[dsub['depmap'].isin(set(dsub['depmap']) & mcells)]
    wsub = dsub[dsub['depmap'].isin(wtcells)]
    if len(msub) < 1 or len(wsub) < 1:
        print(f"{name}: 无格 n_mut={len(msub)} n_wt={len(wsub)}"); continue
    pi_wt = float(wsub['cens'].mean())
    m_cap = np.minimum(msub['LN_IC50'].values, np.log(msub['MAX_CONC'].values))
    w_cap = np.minimum(wsub['LN_IC50'].values, np.log(wsub['MAX_CONC'].values))
    dcap, _, pP = mww_perm(m_cap, w_cap, direction)
    d_auc, _, pA = mww_perm(msub['AUC'].values, wsub['AUC'].values, direction)
    ok_cap = ((dcap < 0) if direction == 'mut_sensitive' else (dcap > 0)) and pP < 0.05
    ok_auc = ((d_auc < 0) if direction == 'mut_sensitive' else (d_auc > 0)) and pA < 0.05
    agree = "一致" if ok_cap == ok_auc else "分歧"
    rows.append(dict(name=name, direction=direction, n_mut=len(msub), n_wt=len(wsub), pi_wt=pi_wt,
                     d_cap=dcap, p_perm=pP, d_auc=d_auc, p_auc=pA, agree=agree))
    print(f"{name:38s} {len(msub):4d} {len(wsub):4d} | {pi_wt:6.1%} | {dcap:+6.3f} {pP:7.4f} | {d_auc:+6.3f} {pA:7.4f} | {agree}")

print("\n=== 汇总 ===")
recall_cap = sum(1 for r in rows if ((r['d_cap'] < 0) if r['direction'] == 'mut_sensitive' else (r['d_cap'] > 0)) and r['p_perm'] < 0.05)
recall_auc = sum(1 for r in rows if ((r['d_auc'] < 0) if r['direction'] == 'mut_sensitive' else (r['d_auc'] > 0)) and r['p_auc'] < 0.05)
print(f"LN_IC50 capped 召回: {recall_cap}/{len(rows)}")
print(f"AUC 召回: {recall_auc}/{len(rows)}")
print(f"召回结论一致: {sum(1 for r in rows if r['agree']=='一致')}/{len(rows)}")

# AUC 压缩行为：按 (DRUG, CANCER_TYPE) 单元的删失率分层
print("\n=== AUC 在单元内的变异度（按删失率分层，全 panel）===")
gd2 = gdsc[gdsc['depmap'].notna() & gdsc['AUC'].notna() & gdsc['LN_IC50'].notna()].copy()
unit = gd2.groupby(["DRUG_NAME", "CANCER_TYPE"], as_index=False).apply(
    lambda g: pd.Series({"n": len(g), "pi": float(g['cens'].mean()),
                         "auc_std": float(g['AUC'].std()), "auc_mean": float(g['AUC'].mean())}),
    include_groups=False)
print(f"单元总数: {len(unit)}")
for lbl, lo, hi in [("pi=0 无删失", -0.001, 0.001),
                    ("0<pi<=0.5", 0.001, 0.5),
                    ("0.5<pi<1", 0.5, 0.999),
                    ("pi=1 全删失", 0.999, 1.001)]:
    u = unit[(unit['pi'] > lo) & (unit['pi'] <= hi)]
    if len(u):
        print(f"  {lbl:14s} 单元 {len(u):5d} | AUC均值 {u['auc_mean'].mean():7.3f} | AUC单元内SD均值 {u['auc_std'].mean():.3f} | n均值 {u['n'].mean():6.1f}")

out = r"D:\paper\论文5\89_auc_comparison_result.csv"
pd.DataFrame(rows).to_csv(out, index=False, encoding='utf-8-sig')
print(f"\n落盘: {out}")
