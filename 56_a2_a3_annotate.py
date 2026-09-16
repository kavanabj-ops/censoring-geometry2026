# -*- coding: utf-8 -*-
"""
56_a2_a3_annotate.py — opus5 底线第2步：A2（8条耐药melanoma共变异）+ A3（Fisher功效+CRC背景）
A2: 8 条耐药 melanoma 系（BRAFi 下 IC50>max conc）的共变异注释，判断"双峰"能否被已知耐药机制解释。
A3: Fisher 精确功效曲线 + 6 条 CRC V600E 系的 MSI/CIMP/PIK3CA 背景。
"""
import re
import numpy as np
import pandas as pd
import scipy.stats as st

MUT = r'D:\paper\data\CCLE_mutations.csv'
SINFO = r'D:\paper\data\sample_info.csv'
GDSC2 = r'D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx'

def extract_codon(pc):
    m = re.search(r'p\.([A-Z])(\d+)([A-Z*])', str(pc))
    return m.group(1)+m.group(2)+m.group(3) if m else None

mut = pd.read_csv(MUT, low_memory=False)
nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins',
             'In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)].copy()
pc_col = 'Protein_Change' if 'Protein_Change' in mut.columns else 'HGVSp_Short'
mut_ns['codon'] = mut_ns[pc_col].apply(extract_codon)

si = pd.read_csv(SINFO)
s2d = {str(r['Sanger_Model_ID']).strip(): r['DepMap_ID'] for _, r in si.iterrows() if pd.notna(r.get('Sanger_Model_ID'))}
d2name = {}
for _, r in si.iterrows():
    dep = r['DepMap_ID']
    nm = r.get('model_name') or r.get('CellLineName') or r.get('cell_line_name')
    d2name[dep] = nm

gdsc = pd.read_excel(GDSC2)
sub = gdsc[gdsc['DRUG_NAME']=='Dabrafenib'].copy()
sub['depmap'] = sub['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)

# 找到 8 条耐药 melanoma 系（IC50 > ln(10)）
res_mel = sub[sub['LN_IC50'] > np.log(10)].copy()
# 需要 lineage=skin + BRAF V600E，用 CCLE 突变判定
braf_v600e = set(mut_ns[(mut_ns['Hugo_Symbol']=='BRAF') & (mut_ns['codon'].str.contains('V600', na=False))]['DepMap_ID'])

print("=== A2: BRAF V600E 耐药 melanoma 系（LN_IC50>ln10）共变异注释 ===")
res_mel_v600e = res_mel[res_mel['depmap'].isin(braf_v600e)]
print(f"耐药 melanoma V600E 系数量: {len(res_mel_v600e)}\n")

# 关注的耐药相关基因
genes = ['NRAS','KRAS','PTEN','CDKN2A','NF1','MAP2K1','MAP2K2','MITF','TP53','BRAF','EGFR','ERBB3','PIK3CA','AKT1']

for _, r in res_mel_v600e.iterrows():
    dep = r['depmap']
    nm = d2name.get(dep) or dep
    ln_ic50 = r['LN_IC50']
    b = mut_ns[(mut_ns['DepMap_ID']==dep)]
    # BRAF 突变详情
    braf_codons = [str(c) for c in b[b['Hugo_Symbol']=='BRAF']['codon'].tolist() if c]
    # 其他耐药基因共突变
    others = {}
    for g in ['NRAS','KRAS','PTEN','CDKN2A','NF1','MAP2K1','MAP2K2','TP53','PIK3CA','EGFR','ERBB3']:
        cc = [str(c) for c in b[b['Hugo_Symbol']==g]['codon'].tolist() if c]
        if cc:
            others[g] = cc
    print(f"  {nm:16s} LN_IC50={ln_ic50:.2f}")
    print(f"    BRAF={braf_codons}")
    if others:
        print(f"    共突变={others}")
    else:
        print(f"    共突变=无")
    print()

print("\n=== A2 对照：敏感 melanoma V600E 系（LN_IC50<-1，最敏感）的共突变 ===")
sens_mel = sub[(sub['LN_IC50'] < -1) & (sub['depmap'].isin(braf_v600e))]
for _, r in sens_mel.head(8).iterrows():
    dep = r['depmap']
    nm = d2name.get(dep) or dep
    b = mut_ns[mut_ns['DepMap_ID']==dep]
    others = {}
    for g in ['NRAS','PTEN','CDKN2A','NF1','MAP2K1','MAP2K2']:
        cc = [str(c) for c in b[b['Hugo_Symbol']==g]['codon'].tolist() if c]
        if cc:
            others[g] = cc
    print(f"  {nm:16s} LN_IC50={r['LN_IC50']:.2f}, 共突变={others if others else '无'}")

# ===== A3: Fisher 精确功效 + CRC 背景 =====
print("\n\n=== A3: Fisher 精确检验功效（CRC-mut n=6 vs skin-mut n=38）===")
print("模拟：给定真实耐药率差异，Fisher 单侧检验的检出率（功效）")
n1, n2 = 6, 38
for p1, p2 in [(0.5, 0.21), (0.7, 0.21), (0.9, 0.21), (0.5, 0.14), (0.7, 0.14)]:
    n_sig = 0
    n_iter = 5000
    for _ in range(n_iter):
        x1 = np.random.binomial(n1, p1)
        x2 = np.random.binomial(n2, p2)
        _, p = st.fisher_exact([[x1, n1-x1],[x2, n2-x2]], alternative='greater')
        if p < 0.05:
            n_sig += 1
    print(f"  CRC耐药率={p1:.0%} vs skin={p2:.0%}: 功效={n_sig/n_iter:.3f}")

print("\n=== A3: 6 条 CRC V600E 系背景（MSI/CIMP 相关基因 + PIK3CA）===")
# CRC V600E 系
crc_cells = ['HT-29','COLO-205','LS-411N','RKO','SW1417','SNU-C5']
# 用 CCLE 突变查 PIK3CA, RNF43, APC, TP53, MSI 相关
for nm in crc_cells:
    # 找 depmap
    dep = None
    for d, n in d2name.items():
        if n and nm.lower() in str(n).lower():
            dep = d
            break
    if dep is None:
        print(f"  {nm}: 未找到 DepMap ID")
        continue
    b = mut_ns[mut_ns['DepMap_ID']==dep]
    others = {}
    for g in ['PIK3CA','RNF43','APC','TP53','KRAS','BRAF','MLH1','MSH2','MSH6']:
        cc = [str(c) for c in b[b['Hugo_Symbol']==g]['codon'].tolist() if c]
        if cc:
            others[g] = cc
    print(f"  {nm:14s} ({dep}): {others if others else '无相关突变'}")
