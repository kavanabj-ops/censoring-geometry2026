# -*- coding: utf-8 -*-
# 70_negative_controls.py - 阴性对照：无关 drug x gene 配对的 FPR 校准（3 估计量）
import pandas as pd, numpy as np, scipy.stats as st, statsmodels.api as sm, random

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

print("加载 ...")
gdsc = pd.read_excel(GDSC2, usecols=["DRUG_NAME","PUTATIVE_TARGET","SANGER_MODEL_ID","LN_IC50","MAX_CONC"])
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False, usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification"])

s2d = {}
d2lineage = {}
for _, r in si.iterrows():
    sm_ = r.get('Sanger_Model_ID')
    if pd.notna(sm_):
        s2d[str(sm_).strip()] = r['DepMap_ID']
    d2lineage[r['DepMap_ID']] = r['lineage']

nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)]
mut_genes_ns = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage'] = gdsc['depmap'].map(d2lineage)

CYTO = ['crosslinker','alkylating','antimetabolite','microtubule','stabiliser','destabiliser','topoisomerase','nucleoside','nucleotide','dna synthesis','dna replication']
def classify(t):
    if pd.isna(t) or str(t).strip()=='':
        return 'unannotated'
    s = str(t).lower()
    if s.strip() in ('top1','top2') or 'top2a' in s or 'top2b' in s:
        return 'cytotoxic'
    if any(k in s for k in CYTO):
        return 'cytotoxic'
    return 'targeted'

drug_target = gdsc.groupby('DRUG_NAME')['PUTATIVE_TARGET'].first()
drug_class = {d: classify(t) for d, t in drug_target.items()}

import re
ALIAS = {'MEK1':'MAP2K1','MEK2':'MAP2K2','ERK1':'MAPK3','ERK2':'MAPK1','PI3K':'PIK3CA','PI3KALPHA':'PIK3CA','MTORC1':'MTOR','MTORC2':'MTOR','ABL':'ABL1','C-KIT':'KIT','KIT':'KIT','PDGFR':'PDGFRA','VEGFR':'KDR','PDK1':'PDPK1','CHK1':'CHEK1','CHK2':'CHEK2','BRAF':'BRAF','KRAS':'KRAS','NRAS':'NRAS','RAF1':'RAF1','CRAF':'RAF1','EGFR':'EGFR','ERBB2':'ERBB2','MET':'MET','ALK':'ALK','ROS1':'ROS1','RET':'RET','NTRK1':'NTRK1','NTRK2':'NTRK2','NTRK3':'NTRK3','FGFR1':'FGFR1','FGFR2':'FGFR2','FGFR3':'FGFR3','FGFR4':'FGFR4','AKT1':'AKT1','AKT2':'AKT2','AKT3':'AKT3','MDM2':'MDM2','MDM4':'MDM4','TP53':'TP53','JAK1':'JAK1','JAK2':'JAK2','JAK3':'JAK3','TYK2':'TYK2','IGF1R':'IGF1R','INSR':'INSR','FLT3':'FLT3','CSF1R':'CSF1R','IDH1':'IDH1','IDH2':'IDH2','EZH2':'EZH2','BRD4':'BRD4','DOT1L':'DOT1L','HDAC1':'HDAC1','HDAC6':'HDAC6','MCL1':'MCL1','BCL2':'BCL2','BCL6':'BCL6','XIAP':'XIAP','WEE1':'WEE1','ATM':'ATM','ATR':'ATR','PARP1':'PARP1','PARP2':'PARP2','TERT':'TERT','XPO1':'XPO1','SMO':'SMO','CDK2':'CDK2','CDK4':'CDK4','CDK6':'CDK6','CDK1':'CDK1','CDK7':'CDK7','CDK9':'CDK9','AURKA':'AURKA','AURKB':'AURKB','PLK1':'PLK1','PIM1':'PIM1','SIRT1':'SIRT1','BCL2L1':'BCL2L1','BCL2L2':'BCL2L2'}

def parse_genes(t):
    if pd.isna(t) or str(t).strip()=='':
        return []
    s = re.sub(r'\([^)]*\)', '', str(t))
    parts = re.split(r'[,;/\s]+', s)
    out = []
    for p in parts:
        p = p.strip().rstrip('.').strip('"').strip("'")
        if not p:
            continue
        m = re.fullmatch(r'([A-Z][A-Z0-9-]{0,11})', p.upper())
        if m:
            out.append(ALIAS.get(m.group(1), m.group(1)))
    seen = []
    for g in out:
        if g not in seen:
            seen.append(g)
    return seen

DRIVER_GENES = ['TP53','KRAS','PIK3CA','APC','BRAF','EGFR','PTEN','NRAS','CDKN2A','RB1',
'FBXW7','SMAD4','CTNNB1','ATM','ARID1A','KMT2D','KMT2C','ERBB2','MET','NF1',
'STK11','KEAP1','SMARCA4','VHL','IDH1','IDH2','FLT3','JAK2','NOTCH1','DNMT3A',
'TET2','GATA3','ESR1','KIT','PDGFRA','NFE2L2','RNF43','SOX9','CREBBP','EP300']

# 随机配对
random.seed(42)
targeted = [d for d, c in drug_class.items() if c == 'targeted']
pairs = []
for drug in targeted:
    tgt_genes = parse_genes(drug_target[drug])
    unrelated = [g for g in DRIVER_GENES if g not in tgt_genes]
    if not unrelated:
        continue
    dsub = gdsc[gdsc['DRUG_NAME']==drug]
    lineages = [l for l in dsub['lineage'].dropna().unique() if l is not None and not (isinstance(l,float) and pd.isna(l))]
    if not lineages:
        continue
    for gene in random.sample(unrelated, min(4, len(unrelated))):
        for lineage in random.sample(lineages, min(3, len(lineages))):
            pairs.append((drug, gene, lineage))

print(f"生成 {len(pairs)} 个无关 drug×gene×lineage 配对")

# 检验
def test_pair(drug, gene, lineage):
    dsub = gdsc[(gdsc['DRUG_NAME']==drug) & (gdsc['lineage']==lineage)].copy()
    dsub = dsub[dsub['depmap'].notna() & dsub['LN_IC50'].notna()]
    mcells = set(c for c, gs in mut_genes_ns.items() if gene in gs)
    cells = set(dsub['depmap'])
    cm = cells & mcells
    cw = cells - mcells
    mut_y = dsub[dsub['depmap'].isin(cm)]['LN_IC50'].values
    wt_y = dsub[dsub['depmap'].isin(cw)]['LN_IC50'].values
    n_mut, n_wt = len(mut_y), len(wt_y)
    if n_mut < 3 or n_wt < 5:
        return None
    y = np.concatenate([mut_y, wt_y]); g = np.array([1]*n_mut + [0]*n_wt)
    # naive OLS（双侧，看是否有假阳性方向）
    m1 = sm.OLS(y, sm.add_constant(g)).fit()
    b1 = m1.params[1]; p1 = m1.pvalues[1]
    # 封顶
    c_ = np.log(dsub['MAX_CONC'].median())
    yc = np.clip(y, -np.inf, c_)
    m2 = sm.OLS(yc, sm.add_constant(g)).fit()
    b2 = m2.params[1]; p2 = m2.pvalues[1]
    # MWW（双侧）
    _, p3 = st.mannwhitneyu(mut_y, wt_y, alternative='two-sided')
    # 判"假阳性"：双侧 p<0.05 且 |效应|≥门槛（不管方向，因为无关配对方向随机）
    f1 = (p1 < 0.05) and (abs(b1) >= 0.5)
    f2 = (p2 < 0.05) and (abs(b2) >= 0.5)
    # MWW 效应量
    U, _ = st.mannwhitneyu(mut_y, wt_y, alternative='two-sided')
    d3 = 2.0*U/(n_mut*n_wt) - 1.0
    f3 = (p3 < 0.05) and (abs(d3) >= 0.3)
    return (f1, f2, f3)

fpr = {'naive':0,'clip':0,'mww':0}
n_valid = 0
for drug, gene, lineage in pairs:
    r = test_pair(drug, gene, lineage)
    if r is None:
        continue
    n_valid += 1
    fpr['naive'] += r[0]
    fpr['clip'] += r[1]
    fpr['mww'] += r[2]

print(f"\n=== 阴性对照 FPR（无效配对，n_valid={n_valid}）===")
for k in ['naive','clip','mww']:
    print(f"  {k:6s}: FPR = {fpr[k]}/{n_valid} = {fpr[k]/max(n_valid,1):.1%}  (名义 5%)")
