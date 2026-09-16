# -*- coding: utf-8 -*-
# 69_positive_controls.py - 全景推断：阳性对照召回率（3 估计量并列）
import pandas as pd, numpy as np, scipy.stats as st, statsmodels.api as sm

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

print("加载 ...")
gdsc = pd.read_excel(GDSC2, usecols=["DRUG_NAME","SANGER_MODEL_ID","LN_IC50","MAX_CONC"])
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False, usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification","Protein_Change","isTCGAhotspot","isCOSMIChotspot"])

s2d = {}
d2lineage = {}
for _, r in si.iterrows():
    sm_ = r.get('Sanger_Model_ID')
    if pd.notna(sm_):
        s2d[str(sm_).strip()] = r['DepMap_ID']
    d2lineage[r['DepMap_ID']] = r['lineage']

nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)].copy()
def is_true(x):
    return str(x).strip().lower() in ('true','1','yes','t')
mut_ns['hot'] = mut_ns['isTCGAhotspot'].apply(is_true) | mut_ns['isCOSMIChotspot'].apply(is_true)

mut_genes_ns = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

def cells_with(gene, pattern):
    sub = mut_ns[(mut_ns['Hugo_Symbol']==gene) & (mut_ns['Protein_Change'].astype(str).str.contains(pattern, na=False, regex=True))]
    return set(sub['DepMap_ID'])

def cells_hotspot(gene):
    sub = mut_ns[(mut_ns['Hugo_Symbol']==gene) & (mut_ns['hot'])]
    return set(sub['DepMap_ID'])

gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage'] = gdsc['depmap'].map(d2lineage)

def cliff_delta(a, b):
    # a 组相对 b 组：delta = 2*U/(na*nb) - 1，a 更小 -> delta < 0
    if len(a)==0 or len(b)==0:
        return np.nan
    U, _ = st.mannwhitneyu(a, b, alternative='two-sided')
    return 2.0*U/(len(a)*len(b)) - 1.0

# 阳性对照：name, drug, gene, pattern(位点) or hotspot, lineage, direction
# direction: 'mut_sensitive' (H1 mut 更敏感) / 'wt_sensitive' (H1 WT 更敏感，如 Nutlin/TP53)
controls = [
    ("BRAF V600E / Dabrafenib / skin", "Dabrafenib", "BRAF", "V600", "skin", "mut_sensitive"),
    ("BRAF V600E / Trametinib / skin", "Trametinib", "BRAF", "V600", "skin", "mut_sensitive"),
    ("TP53 / Nutlin-3a / pan", "Nutlin-3a", "TP53", "__hotspot__", None, "wt_sensitive"),
    ("EGFR L858R+ex19 / Afatinib / lung", "Afatinib", "EGFR", "L858|747_753|E746|ELREA", "lung", "mut_sensitive"),
    ("NRAS Q61 / Trametinib / skin", "Trametinib", "NRAS", "Q61", "skin", "mut_sensitive"),
    ("PIK3CA H1047R/E545 / Alpelisib / breast", "Alpelisib", "PIK3CA", "H1047|E545|E542", "breast", "mut_sensitive"),
    ("KRAS G12/13 / Trametinib / lung", "Trametinib", "KRAS", "G12|G13", "lung", "mut_sensitive"),
]

def one_test(mut_y, wt_y, direction):
    """对一组 mut/wt lnIC50，跑 3 个估计量，返回 (naive_result, clip_result, mww_result)
    每个 result = dict(beta, p, effect, recall)"""
    out = {}
    n_mut, n_wt = len(mut_y), len(wt_y)
    y = np.concatenate([mut_y, wt_y])
    g = np.array([1]*n_mut + [0]*n_wt)
    # direction: mut_sensitive -> mut 更小；wt_sensitive -> mut 更大
    sign = -1.0 if direction == 'mut_sensitive' else 1.0  # mut 系数期望的符号

    # naive OLS
    X = sm.add_constant(g)
    m = sm.OLS(y, X).fit()
    b = m.params[1]
    # 单侧 p：mut 系数 < 0（mut_sensitive）或 > 0（wt_sensitive）
    p_one = m.pvalues[1]/2 if sign*b < 0 else 1 - m.pvalues[1]/2
    recall_naive = (sign*b < 0) and (p_one < 0.05) and (abs(b) >= 0.5)
    out['naive'] = dict(beta=b, p=p_one, recall=recall_naive)

    # 封顶 OLS（lnIC50 截到 ln(MAX_CONC)）—— 这里用全局 c=ln(10) 简化，实际应每药 MAX_CONC
    # 用每药 MAX_CONC：由 gdsc 传入，此处简化用封顶到中位 MAX_CONC 的对数
    return out, y, g

# 完整实现：需要 MAX_CONC 封顶，重新组织
print("\n=== 阳性对照召回率（3 估计量）===\n")
print(f"{'对照':42s} {'n_mut':>5} {'n_WT':>5} | {'naive':>22s} {'封顶':>22s} {'MWW(秩)':>22s}")
rows = []
for name, drug, gene, pattern, lineage, direction in controls:
    dm = [d for d in gdsc['DRUG_NAME'].unique() if drug.lower() in str(d).lower()]
    if not dm:
        print(f"{name:42s}  drug {drug} 不在 GDSC2")
        continue
    dactual = dm[0]
    dsub = gdsc[gdsc['DRUG_NAME']==dactual].copy()
    dsub = dsub[dsub['depmap'].notna() & dsub['LN_IC50'].notna()]
    if lineage:
        dsub = dsub[dsub['lineage']==lineage]
    if len(dsub) == 0:
        print(f"{name:42s}  lineage={lineage} 无数据")
        continue
    # mut 格
    if pattern == '__hotspot__':
        mcells = cells_hotspot(gene)
    else:
        mcells = cells_with(gene, pattern)
    # WT 格：无该基因非沉默突变
    wtcells = set(dsub['depmap']) - set(c for c in dsub['depmap'] if gene in mut_genes_ns.get(c, set()))
    mut_cells_in = set(dsub['depmap']) & mcells
    wt_cells_in = set(dsub['depmap']) & wtcells

    mut_y = dsub[dsub['depmap'].isin(mut_cells_in)]['LN_IC50'].values
    wt_y = dsub[dsub['depmap'].isin(wt_cells_in)]['LN_IC50'].values
    n_mut, n_wt = len(mut_y), len(wt_y)
    if n_mut < 3 or n_wt < 3:
        print(f"{name:42s}  n_mut={n_mut} n_WT={n_wt} (不足)")
        continue

    sign = -1.0 if direction == 'mut_sensitive' else 1.0

    # naive OLS
    y = np.concatenate([mut_y, wt_y]); g = np.array([1]*n_mut + [0]*n_wt)
    m1 = sm.OLS(y, sm.add_constant(g)).fit()
    b1 = m1.params[1]
    dir_ok = (b1 < 0) if direction == 'mut_sensitive' else (b1 > 0)
    p1 = m1.pvalues[1]/2 if dir_ok else 1 - m1.pvalues[1]/2
    r1 = dir_ok and (p1 < 0.05) and (abs(b1) >= 0.5)

    # 封顶 OLS：截到该药 MAX_CONC
    c_ = np.log(dsub['MAX_CONC'].median())
    yc = np.clip(y, -np.inf, c_)
    m2 = sm.OLS(yc, sm.add_constant(g)).fit()
    b2 = m2.params[1]
    dir_ok2 = (b2 < 0) if direction == 'mut_sensitive' else (b2 > 0)
    p2 = m2.pvalues[1]/2 if dir_ok2 else 1 - m2.pvalues[1]/2
    r2 = dir_ok2 and (p2 < 0.05) and (abs(b2) >= 0.5)

    # MWW 秩：单侧（mut < wt 若 mut_sensitive；mut > wt 若 wt_sensitive）
    alt = 'less' if direction == 'mut_sensitive' else 'greater'
    U, p3 = st.mannwhitneyu(mut_y, wt_y, alternative=alt)
    d3 = 2.0*U/(n_mut*n_wt) - 1.0
    # 效应量门槛：|d3| >= 0.3
    r3 = (p3 < 0.05) and (abs(d3) >= 0.3)

    def fmt(b, p, r):
        return f"β={b:+.2f} p={p:.3f} {'✓' if r else '✗'}"
    print(f"{name:42s} {n_mut:5d} {n_wt:5d} | {fmt(b1,p1,r1):22s} {fmt(b2,p2,r2):22s} δ={d3:+.2f} p={p3:.3f} {'✓' if r3 else '✗'}")
    rows.append(dict(name=name, n_mut=n_mut, n_wt=n_wt,
                     naive=r1, clip=r2, mww=r3))

print("\n=== 召回率汇总 ===")
core = [r for r in rows if 'V600E' in r['name'] or 'Nutlin' in r['name'] or 'L858R' in r['name']]
for est in ['naive','clip','mww']:
    k_all = sum(1 for r in rows if r[est])
    k_core = sum(1 for r in core if r[est])
    print(f"  {est:6s}: 全部 {k_all}/{len(rows)}, 核心4 {k_core}/{len(core)}")
