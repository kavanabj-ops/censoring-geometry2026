# -*- coding: utf-8 -*-
# 64_a1b_survey.py - A1b 描述性普查：X% = 落入不可识别区(WT格删失率>c*)的 drug x gene x lineage 比例
import pandas as pd, numpy as np, re, json

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

C_STAR = 0.90
MIN_N = 5

# ---------- 1. 加载 ----------
print("加载 GDSC2 ...")
gdsc = pd.read_excel(GDSC2, usecols=["DRUG_NAME","DRUG_ID","PUTATIVE_TARGET","PATHWAY_NAME","SANGER_MODEL_ID","LN_IC50","MAX_CONC"])
si = pd.read_csv(SINFO)
print("加载 CCLE mutations ...")
mut = pd.read_csv(MUT, low_memory=False, usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification"])

# ---------- 2. 映射 ----------
s2d = {}
for _, r in si.iterrows():
    sm = r.get('Sanger_Model_ID')
    if pd.notna(sm):
        s2d[str(sm).strip()] = r['DepMap_ID']
d2lineage = {r['DepMap_ID']: r['lineage'] for _, r in si.iterrows()}

# ---------- 3. 非沉默突变 -> 每细胞系基因集合 ----------
nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)]
mut_genes = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()
print(f"非沉默突变: {len(mut_ns)} 行, {len(mut_genes)} 个细胞系")

# 癌症 driver 基因（细胞毒对照用，避免乘客基因 TTN/MUC16 无意义）
DRIVER_GENES = ['TP53','KRAS','PIK3CA','APC','BRAF','EGFR','PTEN','NRAS','CDKN2A','RB1',
'FBXW7','SMAD4','CTNNB1','ATM','ARID1A','KMT2D','KMT2C','ERBB2','MET','NF1',
'STK11','KEAP1','SMARCA4','VHL','IDH1','IDH2','FLT3','JAK2','NOTCH1','DNMT3A',
'TET2','GATA3','ESR1','KIT','PDGFRA','NFE2L2','RNF43','SOX9','CREBBP','EP300']
top_genes = DRIVER_GENES
print(f"driver 基因（细胞毒对照）: {len(top_genes)} 个")

# ---------- 4. 删失 ----------
gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage'] = gdsc['depmap'].map(d2lineage)

# ---------- 5. 靶向/细胞毒分类 ----------
CYTO = ['crosslinker','alkylating','antimetabolite','microtubule','stabiliser','destabiliser','topoisomerase','nucleoside','nucleotide','dna synthesis','dna replication']
def classify(t):
    if pd.isna(t) or str(t).strip()=='':
        return 'unannotated'
    s = str(t).lower()
    if 'top1' == s.strip() or 'top2' == s.strip() or 'top1, top2' in s or 'top2a' in s or 'top2b' in s:
        return 'cytotoxic'
    if any(k in s for k in CYTO):
        return 'cytotoxic'
    return 'targeted'

drug_target = gdsc.groupby('DRUG_NAME')['PUTATIVE_TARGET'].first()
drug_class = {d: classify(t) for d, t in drug_target.items()}
from collections import Counter
print("drug 分类:", Counter(drug_class.values()))

# ---------- 6. 基因解析 ----------
ALIAS = {
 'MEK1':'MAP2K1','MEK2':'MAP2K2','ERK1':'MAPK3','ERK2':'MAPK1',
 'PI3KALPHA':'PIK3CA','PI3KBETA':'PIK3CB','PI3K':'PIK3CA','P110ALPHA':'PIK3CA','PI3KALPHA_BETA':'PIK3CA',
 'MTORC1':'MTOR','MTORC2':'MTOR','FRAP1':'MTOR',
 'BCL-XL':'BCL2L1','BCL-W':'BCL2L2','BCL-B':'BCL2L10','BFL1':'BCL2A1','BCLXL':'BCL2L1',
 'PDGFR':'PDGFRA','VEGFR':'KDR','VEGFR2':'KDR',
 'IKK-1':'CHUK','IKK-2':'IKBKB','IKK1':'CHUK','IKK2':'IKBKB',
 'ABL':'ABL1','C-KIT':'KIT','KIT':'KIT','SRC':'SRC',
 'P38':'MAPK14','P38A':'MAPK14','JNK1':'MAPK8','JNK2':'MAPK9','JNK3':'MAPK10',
 'CHK1':'CHEK1','CHK2':'CHEK2','GSK3':'GSK3B','GSK3B':'GSK3B',
 'AURKB':'AURKB','AURKA':'AURKA','PLK1':'PLK1','PLK2':'PLK2','PLK3':'PLK3',
 'CDK4':'CDK4','CDK6':'CDK6','CDK2':'CDK2','CDK1':'CDK1','CDK9':'CDK9','CDK7':'CDK7',
 'EZH2':'EZH2','BRD4':'BRD4','BRD2':'BRD2','BRD3':'BRD3','DOT1L':'DOT1L',
 'HDAC1':'HDAC1','HDAC2':'HDAC2','HDAC3':'HDAC3','HDAC6':'HDAC6','SIRT1':'SIRT1',
 'MCL1':'MCL1','BCL2':'BCL2','BCL6':'BCL6','XIAP':'XIAP','CIAP':'BIRC2',
 'WEE1':'WEE1','WEE2':'WEE2','ATM':'ATM','ATR':'ATR','DNAPK':'PRKDC',
 'PARP1':'PARP1','PARP2':'PARP2','PARP':'PARP1','TERT':'TERT',
 'IDH1':'IDH1','IDH2':'IDH2','KDM1A':'KDM1A',
 'EGFR':'EGFR','ERBB2':'ERBB2','ERBB4':'ERBB4','MET':'MET','ALK':'ALK','ROS1':'ROS1',
 'NTRK1':'NTRK1','NTRK2':'NTRK2','NTRK3':'NTRK3','RET':'RET','FGFR1':'FGFR1','FGFR2':'FGFR2','FGFR3':'FGFR3','FGFR4':'FGFR4',
 'BRAF':'BRAF','KRAS':'KRAS','NRAS':'NRAS','HRAS':'HRAS','RAF1':'RAF1','CRAF':'RAF1',
 'AKT1':'AKT1','AKT2':'AKT2','AKT3':'AKT3','PDK1':'PDPK1',
 'MDM2':'MDM2','MDM4':'MDM4','TP53':'TP53','JAK1':'JAK1','JAK2':'JAK2','JAK3':'JAK3','TYK2':'TYK2',
 'IGF1R':'IGF1R','INSR':'INSR','FLT3':'FLT3','CSF1R':'CSF1R','EPHA2':'EPHA2','DDR1':'DDR1',
 'XPO1':'XPO1','SMO':'SMO','PTCH1':'PTCH1',
}

def parse_genes(t):
    if pd.isna(t) or str(t).strip()=='':
        return []
    s = str(t)
    s = re.sub(r'\([^)]*\)', '', s)
    parts = re.split(r'[,;/\s]+', s)
    out = []
    for p in parts:
        p = p.strip().rstrip('.').strip('"').strip("'")
        if not p:
            continue
        m = re.fullmatch(r'([A-Z][A-Z0-9-]{0,11})', p.upper())
        if m:
            g = m.group(1)
            out.append(ALIAS.get(g, g))
    seen = []
    for g in out:
        if g not in seen:
            seen.append(g)
    return seen

# ---------- 7. 主循环 ----------
def survey(drugs, genes_fn):
    recs = []
    for drug in drugs:
        genes = genes_fn(drug)
        if not genes:
            continue
        dsub = gdsc[gdsc['DRUG_NAME']==drug]
        for gene in genes:
            mcells = set(c for c, gs in mut_genes.items() if gene in gs)
            for lineage, lsub in dsub.groupby('lineage'):
                if lineage is None or (isinstance(lineage, float) and pd.isna(lineage)):
                    continue
                cells = set(lsub['depmap'].dropna())
                cm = cells & mcells
                cw = cells - mcells
                nm, nw = len(cm), len(cw)
                if nm < MIN_N or nw < MIN_N:
                    recs.append(dict(drug=drug, gene=gene, lineage=str(lineage), status='insufficient', n_mut=nm, n_wt=nw, r_mut=np.nan, r_wt=np.nan))
                    continue
                rm = lsub[lsub['depmap'].isin(cm)]['cens'].mean()
                rw = lsub[lsub['depmap'].isin(cw)]['cens'].mean()
                status = 'nonident' if rw > C_STAR else 'identifiable'
                recs.append(dict(drug=drug, gene=gene, lineage=str(lineage), status=status, n_mut=nm, n_wt=nw, r_mut=rm, r_wt=rw))
    return pd.DataFrame(recs, columns=['drug','gene','lineage','status','n_mut','n_wt','r_mut','r_wt'])

# 靶向药
targeted_drugs = [d for d, c in drug_class.items() if c == 'targeted']
print(f"\n靶向药 {len(targeted_drugs)} 个，普查中 ...")
df_t = survey(targeted_drugs, lambda d: parse_genes(drug_target[d]))

# 细胞毒药（用 top20 高频基因）
cyto_drugs = [d for d, c in drug_class.items() if c == 'cytotoxic']
print(f"细胞毒药 {len(cyto_drugs)} 个，普查中 ...")
df_c = survey(cyto_drugs, lambda d: top_genes)

# ---------- 8. 汇总 ----------
def summarize(df, label):
    print(f"\n=== {label} ===")
    print(f"总组合数: {len(df)}")
    suff = df[df['status']!='insufficient']
    ins = df[df['status']=='insufficient']
    print(f"sufficient: {len(suff)}, insufficient: {len(ins)} ({len(ins)/max(len(df),1):.1%})")
    for c in [0.80, C_STAR, 0.95]:
        nonid = suff[suff['r_wt'] > c]
        x = len(nonid)/max(len(suff),1)
        print(f"  c*={c:.2f}: 不可识别(WT删失>{c}) = {len(nonid)}/{len(suff)} = {x:.1%}")
    # 删失不平衡分布
    if len(suff):
        imb = suff['r_wt'] - suff['r_mut']
        print(f"  删失不平衡 r_wt-r_mut: 中位={imb.median():.2f}, 均值={imb.mean():.2f}, P75={imb.quantile(0.75):.2f}")
        print(f"  WT格删失率分布: 中位={suff['r_wt'].median():.2f}, P25={suff['r_wt'].quantile(0.25):.2f}, P75={suff['r_wt'].quantile(0.75):.2f}")
    return suff

suff_t = summarize(df_t, "靶向药 (drug x 靶点基因 x lineage)")
suff_c = summarize(df_c, "细胞毒药 (drug x 高频基因 x lineage)")

# 保存
out = {
    "C_STAR": C_STAR, "MIN_N": MIN_N,
    "targeted": {"n_total": len(df_t), "n_suff": int(len(df_t[df_t['status']!='insufficient'])),
                 "X_80": float((df_t[df_t['r_wt']>0.80]['r_wt']>0.80).sum()/max(len(df_t[df_t['status']!='insufficient']),1)),
                 "X_90": float((df_t[df_t['r_wt']>C_STAR]['r_wt']>C_STAR).sum()/max(len(df_t[df_t['status']!='insufficient']),1)),
                 "X_95": float((df_t[df_t['r_wt']>0.95]['r_wt']>0.95).sum()/max(len(df_t[df_t['status']!='insufficient']),1))},
    "cytotoxic": {"n_total": len(df_c), "n_suff": int(len(df_c[df_c['status']!='insufficient'])),
                  "X_80": float((df_c[df_c['r_wt']>0.80]['r_wt']>0.80).sum()/max(len(df_c[df_c['status']!='insufficient']),1)),
                  "X_90": float((df_c[df_c['r_wt']>C_STAR]['r_wt']>C_STAR).sum()/max(len(df_c[df_c['status']!='insufficient']),1)),
                  "X_95": float((df_c[df_c['r_wt']>0.95]['r_wt']>0.95).sum()/max(len(df_c[df_c['status']!='insufficient']),1))},
}
with open(r"D:\paper\02_analysis\results\a1b_survey.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n已保存 results/a1b_survey.json")
