# -*- coding: utf-8 -*-
# 68_external_gdsc1.py - 外部几何复现：GDSC2 vs GDSC1 的 X% 对比（删失几何非 artifact 验证）
import pandas as pd, numpy as np, re

GDSC2 = r"D:\paper\data\GDSC2_fitted_dose_response_27Oct23.xlsx"
GDSC1 = r"D:\paper\data\GDSC1_fitted_dose_response_27Oct23.xlsx"
SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"
C_STAR = 0.90
MIN_N = 5

print("加载共享数据 ...")
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False, usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification"])

s2d = {}
d2lineage = {}
for _, r in si.iterrows():
    sm = r.get('Sanger_Model_ID')
    if pd.notna(sm):
        s2d[str(sm).strip()] = r['DepMap_ID']
    d2lineage[r['DepMap_ID']] = r['lineage']

nonsilent = ['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns = mut[mut['Variant_Classification'].isin(nonsilent)]
mut_genes = mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

DRIVER_GENES = ['TP53','KRAS','PIK3CA','APC','BRAF','EGFR','PTEN','NRAS','CDKN2A','RB1',
'FBXW7','SMAD4','CTNNB1','ATM','ARID1A','KMT2D','KMT2C','ERBB2','MET','NF1',
'STK11','KEAP1','SMARCA4','VHL','IDH1','IDH2','FLT3','JAK2','NOTCH1','DNMT3A',
'TET2','GATA3','ESR1','KIT','PDGFRA','NFE2L2','RNF43','SOX9','CREBBP','EP300']

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

ALIAS = {
 'MEK1':'MAP2K1','MEK2':'MAP2K2','ERK1':'MAPK3','ERK2':'MAPK1',
 'PI3KALPHA':'PIK3CA','PI3KBETA':'PIK3CB','PI3K':'PIK3CA','P110ALPHA':'PIK3CA',
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

def survey(drugs, genes_fn, gdsc, min_n=MIN_N):
    insuff = 0; suff = 0; nonid = 0
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
                if nm < min_n or nw < min_n:
                    insuff += 1
                    continue
                suff += 1
                rw = lsub[lsub['depmap'].isin(cw)]['cens'].mean()
                if rw > C_STAR:
                    nonid += 1
    return insuff, suff, nonid

print("\n=== 外部几何复现：GDSC2 vs GDSC1 ===")
for label, path in [("GDSC2", GDSC2), ("GDSC1", GDSC1)]:
    print(f"\n--- {label} ---")
    gdsc = pd.read_excel(path, usecols=["DRUG_NAME","PUTATIVE_TARGET","SANGER_MODEL_ID","LN_IC50","MAX_CONC"])
    gdsc['cens'] = gdsc['LN_IC50'] > np.log(gdsc['MAX_CONC'])
    gdsc['depmap'] = gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
    gdsc['lineage'] = gdsc['depmap'].map(d2lineage)
    drug_target = gdsc.groupby('DRUG_NAME')['PUTATIVE_TARGET'].first()
    drug_class = {d: classify(t) for d, t in drug_target.items()}
    targeted = [d for d, c in drug_class.items() if c == 'targeted']
    cyto = [d for d, c in drug_class.items() if c == 'cytotoxic']
    print(f"  drug 数: 靶向 {len(targeted)}, 细胞毒 {len(cyto)}, 未注释 {sum(1 for c in drug_class.values() if c=='unannotated')}")
    it, st_, nid = survey(targeted, lambda d: parse_genes(drug_target[d]), gdsc)
    ic, sc, nidc = survey(cyto, lambda d: DRIVER_GENES, gdsc)
    print(f"  靶向: insuff={it/(it+st_):.1%}, suff={st_}, X%(WT删失>{C_STAR})={nid/max(st_,1):.1%}")
    print(f"  细胞毒: insuff={ic/(ic+sc):.1%}, suff={sc}, X%={nidc/max(sc,1):.1%}")
