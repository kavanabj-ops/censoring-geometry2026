# -*- coding: utf-8 -*-
# 92_prism_replication.py - PRISM (DepMap) 独立复现（阻塞项4）
# 在 PRISM Repurposing 数据上复现：16 个文献关联的 gate 审计 + 直接检验，
# 与 GDSC2 结果对比。PRISM 删失定义：ic50 > 10 uM（最高测试浓度）。
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, rankdata

PRISM = r"D:\paper\depmap portal\prism\secondary-screen-dose-response-curve-parameters.csv"
CELL  = r"D:\paper\depmap portal\prism\secondary-screen-cell-line-info.csv"
MUT   = r"D:\paper\data\CCLE_mutations.csv"
MAXDOSE = 10.0  # uM

print("加载 PRISM curve parameters（usecols）...")
df = pd.read_csv(PRISM, usecols=["name","depmap_id","ic50","auc"])
print("rows:", len(df), "| drugs:", df["name"].nunique(), "| cells:", df["depmap_id"].nunique())

print("加载 cell-line-info + mutations ...")
cli = pd.read_csv(CELL)
mut = pd.read_csv(MUT, low_memory=False,
                  usecols=["Hugo_Symbol","DepMap_ID","Variant_Classification","Protein_Change",
                           "isTCGAhotspot","isCOSMIChotspot"])
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

# lineage: PRISM primary_tissue 直接做 lineage
t2lineage={}
for _,r in cli.iterrows():
    t2lineage[r['depmap_id']]=str(r.get('primary_tissue','')).lower()

df=df[df['depmap_id'].notna()].copy()
df['ic50']=pd.to_numeric(df['ic50'],errors='coerce')
df['auc']=pd.to_numeric(df['auc'],errors='coerce')
df=df[(df['ic50']>0)&(df['auc'].notna())]
df['cens']=df['ic50']>MAXDOSE
df['ln_ic50_cap']=np.minimum(np.log10(df['ic50']), np.log10(MAXDOSE))
df['lineage']=df['depmap_id'].map(t2lineage)

# 16 关联（药名在 PRISM 里的大小写不敏感匹配）
controls=[
 ("BRAF V600E / Dabrafenib / skin","dabrafenib","BRAF","V600","skin","mut_sensitive"),
 ("BRAF V600E / PLX-4720 / skin","PLX-4720","BRAF","V600","skin","mut_sensitive"),
 ("BRAF V600E / Trametinib / skin","trametinib","BRAF","V600","skin","mut_sensitive"),
 ("BRAF V600E / Selumetinib / skin","selumetinib","BRAF","V600","skin","mut_sensitive"),
 ("EGFR L858R+ex19 / Afatinib / lung","afatinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive"),
 ("EGFR L858R+ex19 / Erlotinib / lung","erlotinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive"),
 ("EGFR L858R+ex19 / Gefitinib / lung","gefitinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive"),
 ("EGFR L858R+ex19 / Osimertinib / lung","osimertinib","EGFR","L858|747_753|E746|ELREA","lung","mut_sensitive"),
 ("NRAS Q61 / Trametinib / pan","trametinib","NRAS","Q61",None,"mut_sensitive"),
 ("NRAS Q61 / Selumetinib / skin","selumetinib","NRAS","Q61","skin","mut_sensitive"),
 ("KRAS G12/13 / Trametinib / lung","trametinib","KRAS","G12|G13","lung","mut_sensitive"),
 ("KRAS G12/13 / Selumetinib / lung","selumetinib","KRAS","G12|G13","lung","mut_sensitive"),
 ("PIK3CA / Alpelisib / breast","alpelisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive"),
 ("PIK3CA / Taselisib / breast","taselisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive"),
 ("PIK3CA / Pictilisib / breast","pictilisib","PIK3CA","H1047|E545|E542","breast","mut_sensitive"),
 ("TP53 / Nutlin-3a / pan","nutlin-3","TP53",None,None,"wt_sensitive"),
]

N_PERM=20000; SEED=20261122; MIN_N=5; C_STAR=0.90
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

print(f"\n{'关联':40s} {'nM':>4} {'nW':>4} {'pi_WT':>6} | gate        | {'d_cap':>6} {'p_perm':>7} | 召回 | PRISM药")
rows=[]
for label,drug, gene,pattern,lineage,direction in controls:
    dm=[d for d in df['name'].unique() if str(d).lower()==drug.lower()]
    if not dm:
        # 部分匹配
        dm=[d for d in df['name'].unique() if drug.lower() in str(d).lower()][:1]
    if not dm:
        rows.append(dict(label=label,drug=dm[0] if dm else None,n_mut=0,n_wt=0,pi_wt=np.nan,gate="drug_missing",d_cap=np.nan,p_perm=np.nan,recall=None))
        print(f"{label:40s} [药不在 PRISM]"); continue
    dsub=df[df['name']==dm[0]].copy()
    dsub=dsub[dsub['lineage'].notna()]
    if lineage: dsub=dsub[dsub['lineage'].astype(str).str.lower()==lineage]
    if pattern is None: mcells=cells_hotspot(gene)
    else: mcells=cells_with(gene,pattern)
    wtcells=set(dsub['depmap_id'])-{c for c in dsub['depmap_id'] if gene in mut_genes_ns.get(c,set())}
    msub=dsub[dsub['depmap_id'].isin(set(dsub['depmap_id'])&mcells)]
    wsub=dsub[dsub['depmap_id'].isin(wtcells)]
    nm,nw=len(msub),len(wsub)
    if nm<1 or nw<1:
        rows.append(dict(label=label,drug=dm[0],n_mut=nm,n_wt=nw,pi_wt=np.nan,gate="no_data",d_cap=np.nan,p_perm=np.nan,recall=None))
        print(f"{label:40s} {nm:4d} {nw:4d} [无格] drug={dm[0]}"); continue
    pi_wt=float(wsub['cens'].mean())
    flags=[]
    if nm<MIN_N or nw<MIN_N: flags.append("样本不足")
    if pi_wt>C_STAR: flags.append("WT高删失")
    gate="; ".join(flags) if flags else "PASS"
    dcap,pP=mww_perm(msub['ln_ic50_cap'].values,wsub['ln_ic50_cap'].values,direction)
    ok=bool(((dcap<0) if direction=='mut_sensitive' else (dcap>0)) and pP<0.05)
    rows.append(dict(label=label,drug=dm[0],n_mut=nm,n_wt=nw,pi_wt=pi_wt,gate=gate,d_cap=dcap,p_perm=pP,recall=ok))
    print(f"{label:40s} {nm:4d} {nw:4d} {pi_wt:6.1%} | {gate:11s} | {dcap:+6.3f} {pP:7.4f} | {'✓' if ok else '✗'} | {dm[0]}")

rdf=pd.DataFrame(rows)
print("\n=== PRISM 汇总 ===")
n_tot=len(rdf)
n_pass=sum(1 for g in rdf['gate'] if g=="PASS")
n_insuff=sum(1 for g in rdf['gate'] if isinstance(g,str) and "样本不足" in g)
n_cens=sum(1 for g in rdf['gate'] if isinstance(g,str) and "WT高删失" in g)
n_miss=sum(1 for g in rdf['gate'] if g in ('drug_missing','no_data'))
rec=sum(1 for r in rows if r['recall'] is True)
print(f"gate PASS: {n_pass}/{n_tot} | 样本不足: {n_insuff} | WT高删失: {n_cens} | 数据缺失: {n_miss}")
print(f"直接检验召回: {rec}/{n_tot}")
out=r"D:\paper\论文5\92_prism_replication_result.csv"
rdf.to_csv(out,index=False,encoding='utf-8-sig')
print(f"\n落盘: {out}")
