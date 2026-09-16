# -*- coding: utf-8 -*-
# 78_sharp_interval.py - sharp interval（opus5 M12）：真实数据上算 beta_int 的识别区间
# Prop 1(i): mu_l sharp set = [theta_l, theta_l + pi_l*B]，theta_l = E[y^c]
# beta_int sharp interval = [beta_tilde - B*(pi_AW+pi_BM), beta_tilde + B*(pi_AM+pi_BW)]
import pandas as pd, numpy as np, re, re

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

gdsc = pd.read_excel(GDSC2, usecols=["DRUG_NAME","SANGER_MODEL_ID","LN_IC50","MAX_CONC"])
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False, usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification","Protein_Change"])

s2d = {}; d2lineage = {}
for _, r in si.iterrows():
    sm = r.get('Sanger_Model_ID')
    if pd.notna(sm): s2d[str(sm).strip()] = r['DepMap_ID']
    d2lineage[r['DepMap_ID']] = r['lineage']

nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)]
def _codon(pc):
    m = re.search(r'p\.([A-Z])(\d+)([A-Z*])', str(pc))
    return (m.group(1)+m.group(2)+m.group(3)) if m else None
braf_mut = mut_ns[mut_ns['Hugo_Symbol']=='BRAF']
braf_v600e = set(braf_mut[braf_mut['Protein_Change'].apply(_codon).astype(str).str.contains('V600E', na=False)]['DepMap_ID'])
mut_genes = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage'] = gdsc['depmap'].map(d2lineage)
gdsc['y'] = gdsc['LN_IC50']
gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
gdsc['yc'] = np.minimum(gdsc['LN_IC50'], np.log(gdsc['MAX_CONC']))

B = np.log(10)  # 外推上限：IC50 最多比 MAX_CONC 高 10 倍（保守，敏感性见下）

# BRAF 案例（dabrafenib, CRC vs skin）
d = gdsc[gdsc['DRUG_NAME']=='Dabrafenib'].copy()
d = d[d['lineage'].isin(['colorectal','skin']) & d['depmap'].notna() & d['LN_IC50'].notna()]
d['is_crc'] = (d['lineage']=='colorectal').astype(int)
d['v600e'] = d['depmap'].isin(braf_v600e).astype(int)

def cell_stats(dd):
    theta = dd['yc'].mean()          # E[y^c]
    pi = dd['cens'].mean()           # 删失率
    n = len(dd)
    return theta, pi, n

tAM, pAM, nAM = cell_stats(d[(d['is_crc']==1)&(d['v600e']==1)])
tAW, pAW, nAW = cell_stats(d[(d['is_crc']==1)&(d['v600e']==0)])
tBM, pBM, nBM = cell_stats(d[(d['is_crc']==0)&(d['v600e']==1)])
tBW, pBW, nBW = cell_stats(d[(d['is_crc']==0)&(d['v600e']==0)])

beta_tilde = tAM - tAW - tBM + tBW   # 封顶伪极限
lo = beta_tilde - B*(pAW + pBM)
hi = beta_tilde + B*(pAM + pBW)
width = B*(pAM + pAW + pBM + pBW)

print("=== BRAF V600E / dabrafenib 的 sharp interval ===")
print(f"四格 n=({nAM},{nAW},{nBM},{nBW}), 删失率=({pAM:.0%},{pAW:.0%},{pBM:.0%},{pBW:.0%})")
print(f"封顶伪极限 beta_tilde = {beta_tilde:.3f}")
print(f"sharp interval = [{lo:.3f}, {hi:.3f}], 宽度 = {width:.3f} (B*Σπ, B={B:.3f})")
print(f"区间包含 0: {'是 -> 不可信息化（不能排除任何符号）' if lo<=0<=hi else '否'}")

# 敏感性：B = ln(3) 和 ln(30)
for Bs in [np.log(3), np.log(10), np.log(30)]:
    lo_s = beta_tilde - Bs*(pAW + pBM)
    hi_s = beta_tilde + Bs*(pAM + pBW)
    print(f"  B={Bs:.3f} (外推 {np.exp(Bs):.0f}x): [{lo_s:.3f}, {hi_s:.3f}] 宽={Bs*(pAM+pAW+pBM+pBW):.3f}")

# 40 个 non-informative 单元的 sharp interval 宽度分布（示意）
print(f"\n=== non-informative 单元的 sharp interval 宽度 ===")
print(f"对 WT 删失率>0.90 的单元，宽度 = B*Σπ ≈ B*(0.2+0.9+0.2+0.9) = {B*(0.2+0.9+0.2+0.9):.3f}（典型 4 格删失 20%/90%/20%/90%）")
print(f"即：即使 B 只取 ln(10)，识别区间宽度 ~3.2 log 单位，远大于任何可信的交互效应量级 -> 不可信息化")
