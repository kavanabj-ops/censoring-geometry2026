# -*- coding: utf-8 -*-
# 91_published_associations_audit.py - 系统化"结论改变"审计（阻塞项3）
# 把文献/临床确立的 mutation-drug 敏感关联过一遍 gate（R0 样本量 + WT 删失），
# 并做 capped MWW 直接检验。报告：可检验比例、gate 失败比例（分层）、
# 以及"文献确立但数据上不可检验/方向不一致"的实例。
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, rankdata

SINFO = r"D:\paper\data\sample_info.csv"
MUT = r"D:\paper\data\CCLE_mutations.csv"

print("加载 ...")
gdsc = pd.read_pickle(r"D:\paper\data\gdsc2_cache.pkl")
si = pd.read_csv(SINFO)
mut = pd.read_csv(MUT, low_memory=False,
                  usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification","Protein_Change",
                           "isTCGAhotspot","isCOSMIChotspot"])

s2d={}; d2lineage={}
for _,r in si.iterrows():
    sm_=r.get('Sanger_Model_ID')
    if pd.notna(sm_): s2d[str(sm_).strip()]=r['DepMap_ID']
    d2lineage[r['DepMap_ID']]=r['lineage']

nonsilent=['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins',
           'In_Frame_Del','In_Frame_Ins','Splice_Site']
mut_ns=mut[mut['Variant_Classification'].isin(nonsilent)].copy()
def is_true(x): return str(x).strip().lower() in ('true','1','yes','t')
mut_ns['hot']=mut_ns['isTCGAhotspot'].apply(is_true)|mut_ns['isCOSMIChotspot'].apply(is_true)
mut_genes_ns=mut_ns.groupby('DepMap_ID')['Hugo_Symbol'].agg(set).to_dict()

def cells_with(gene,pattern):
    sub=mut_ns[(mut_ns['Hugo_Symbol']==gene)&(mut_ns['Protein_Change'].astype(str).str.contains(pattern,na=False,regex=True))]
    return set(sub['DepMap_ID'])
def cells_hotspot(gene):
    return set(mut_ns[(mut_ns['Hugo_Symbol']==gene)&(mut_ns['hot'])]['DepMap_ID'])

gdsc['cens']=gdsc['LN_IC50']>np.log(gdsc['MAX_CONC'])
gdsc['depmap']=gdsc['SANGER_MODEL_ID'].astype(str).str.strip().map(s2d)
gdsc['lineage']=gdsc['depmap'].map(d2lineage)

# 文献/临床确立的关联（方向 = mutant 敏感或 WT 敏感）
# (label, drug, gene, pattern/None=hotspot, lineage/None=pan, direction, 文献依据)
controls=[
 ("BRAF V600E / Dabrafenib / skin","Dabrafenib","BRAF","V600","skin","mut_sensitive","FDA-approved BRAF V600E melanoma"),
 ("BRAF V600E / PLX-4720 / skin","PLX-4720","BRAF","V600","skin","mut_sensitive","BRAF-selective tool compound"),
 ("BRAF V600E / Trametinib / skin","Trametinib","BRAF","V600","skin","mut_sensitive","MEKi active in BRAF-mutant melanoma"),
 ("BRAF V600E / Selumetinib / skin","Selumetinib","BRAF","V600","skin","mut_sensitive","MEKi active in BRAF-mutant melanoma"),
 ("EGFR L858R+ex19 / Afatinib / lung","Afatinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive","FDA-approved EGFR-mutant NSCLC"),
 ("EGFR L858R+ex19 / Erlotinib / lung","Erlotinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive","FDA-approved EGFR-mutant NSCLC"),
 ("EGFR L858R+ex19 / Gefitinib / lung","Gefitinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive","FDA-approved EGFR-mutant NSCLC"),
 ("EGFR L858R+ex19 / Osimertinib / lung","Osimertinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive","FDA-approved EGFR-mutant NSCLC"),
 ("NRAS Q61 / Trametinib / pan","Trametinib","NRAS","Q61",None,"mut_sensitive","NRAS-driven MEK dependence"),
 ("NRAS Q61 / Selumetinib / skin","Selumetinib","NRAS","Q61","skin","mut_sensitive","NRAS-mutant melanoma MEKi (Binimetinib)"),
 ("KRAS G12/13 / Trametinib / lung","Trametinib","KRAS","G12|G13","lung","mut_sensitive","KRAS-driven MEK dependence"),
 ("KRAS G12/13 / Selumetinib / lung","Selumetinib","KRAS","G12|G13","lung","mut_sensitive","KRAS-mutant NSCLC MEKi trials"),
 ("PIK3CA / Alpelisib / breast","Alpelisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive","FDA-approved PIK3CA-mutant breast"),
 ("PIK3CA / Taselisib / breast","Taselisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive","PI3K-selective inhibitor"),
 ("PIK3CA / Pictilisib / breast","Pictilisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive","pan-PI3K inhibitor"),
 ("TP53 / Nutlin-3a / pan","Nutlin-3a (-)","TP53",None,None,"wt_sensitive","MDM2 antagonist, TP53-wildtype active"),
]

N_PERM=20000
SEED=20261122
MIN_N=5
C_STAR=0.90

def mww_perm(a,b,direction):
    a=np.asarray(a,float); b=np.asarray(b,float)
    na,nb=len(a),len(b)
    ranks=rankdata(np.concatenate([a,b])); obs_S=ranks[:na].sum()
    alt='less' if direction=='mut_sensitive' else 'greater'
    U,p_asym=mannwhitneyu(a,b,alternative=alt)
    delta=2.0*U/(na*nb)-1.0
    rng=np.random.default_rng(SEED); cnt=0
    for _ in range(N_PERM):
        perm=rng.permutation(na+nb); S=ranks[perm[:na]].sum()
        if direction=='mut_sensitive': cnt+=S<=obs_S
        else: cnt+=S>=obs_S
    return delta,(cnt+1.0)/(N_PERM+1.0)

rows=[]
print(f"\n{'关联':42s} {'nM':>4} {'nW':>4} {'pi_WT':>6} | gate            | {'d_cap':>6} {'p_perm':>7} | 召回")
for label,drug,gene,pattern,lineage,direction,evidence in controls:
    dm=[d for d in gdsc['DRUG_NAME'].unique() if drug.lower() in str(d).lower()]
    if not dm:
        rows.append(dict(label=label,n_mut=0,n_wt=0,pi_wt=np.nan,gate="drug_missing",d_cap=np.nan,p_perm=np.nan,recall=None))
        print(f"{label:42s} [药不在 GDSC2]"); continue
    dsub=gdsc[gdsc['DRUG_NAME']==dm[0]].copy()
    dsub=dsub[dsub['depmap'].notna()&dsub['LN_IC50'].notna()&dsub['MAX_CONC'].notna()]
    if lineage: dsub=dsub[dsub['lineage']==lineage]
    if pattern is None: mcells=cells_hotspot(gene)
    else: mcells=cells_with(gene,pattern)
    wtcells=set(dsub['depmap'])-{c for c in dsub['depmap'] if gene in mut_genes_ns.get(c,set())}
    msub=dsub[dsub['depmap'].isin(set(dsub['depmap'])&mcells)]
    wsub=dsub[dsub['depmap'].isin(wtcells)]
    nm,nw=len(msub),len(wsub)
    if nm<1 or nw<1:
        rows.append(dict(label=label,n_mut=nm,n_wt=nw,pi_wt=np.nan,gate="no_data",d_cap=np.nan,p_perm=np.nan,recall=None))
        print(f"{label:42s} {nm:4d} {nw:4d} [无格]"); continue
    pi_wt=float(wsub['cens'].mean())
    flags=[]
    if nm<MIN_N or nw<MIN_N: flags.append("样本不足")
    if pi_wt>C_STAR: flags.append("WT高删失")
    gate="; ".join(flags) if flags else "PASS"
    m_cap=np.minimum(msub['LN_IC50'].values,np.log(msub['MAX_CONC'].values))
    w_cap=np.minimum(wsub['LN_IC50'].values,np.log(wsub['MAX_CONC'].values))
    dcap,pP=mww_perm(m_cap,w_cap,direction)
    ok=bool(((dcap<0) if direction=='mut_sensitive' else (dcap>0)) and pP<0.05)
    rows.append(dict(label=label,n_mut=nm,n_wt=nw,pi_wt=pi_wt,gate=gate,d_cap=dcap,p_perm=pP,recall=ok,evidence=evidence))
    print(f"{label:42s} {nm:4d} {nw:4d} {pi_wt:6.1%} | {gate:15s} | {dcap:+6.3f} {pP:7.4f} | {'✓' if ok else ('✗' if ok is not None else '-')}")

rdf=pd.DataFrame(rows)
print("\n=== 汇总 ===")
n_tot=len(rdf)
n_pass=sum(1 for g in rdf['gate'] if g=="PASS")
n_insuff=sum(1 for g in rdf['gate'] if isinstance(g,str) and "样本不足" in g)
n_cens=sum(1 for g in rdf['gate'] if isinstance(g,str) and "WT高删失" in g)
print(f"关联总数: {n_tot}")
print(f"gate PASS: {n_pass}/{n_tot}")
print(f"样本不足: {n_insuff} | WT高删失: {n_cens} | 数据缺失: {sum(1 for g in rdf['gate'] if g in ('no_data','drug_missing'))}")
rec=sum(1 for r in rows if r['recall'] is True)
print(f"直接检验召回: {rec}/{n_tot}")

out=r"D:\paper\论文5\91_published_associations_audit_result.csv"
rdf.to_csv(out,index=False,encoding='utf-8-sig')
print(f"\n落盘: {out}")
