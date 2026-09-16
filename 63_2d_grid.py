# -*- coding: utf-8 -*-
"""
63_2d_grid.py — 二维网格 (删失率 p × n_min)：分离病因（opus5 round4 Q3 首要任务）

T1 有两个候选病因：删失几何 vs n_min=6 设计。
验收：balanced n=(40,40,40,40) 同删失率下失效仍在 ⇒ T1(删失几何)成立；
      若 n 增大后失效显著衰减 ⇒ 实为"小格×删失交互效应"，标题需改。

输出：对每个 (n设置 × WT删失率) 组合，报各估计量的 FPR / size-adj 功效 / σ̂比。
关键判读：n 从 6 增到 40 后，naive 的功效损失与封顶的 FPR 膨胀是否仍保持。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm
import time

rng = np.random.default_rng(777)

c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
rate_mut = 0.20
mu_mut = c + st.norm.ppf(rate_mut) * sigma_mut

def mu_for(rate, sigma):
    return c + st.norm.ppf(rate) * sigma

def gen_data(beta_int, ns, mu_wt):
    n1, n2, n3, n4 = ns
    y = np.concatenate([
        rng.normal(mu_mut + beta_int, sigma_mut, n1),
        rng.normal(mu_wt, sigma_wt, n2),
        rng.normal(mu_mut, sigma_mut, n3),
        rng.normal(mu_wt, sigma_wt, n4),
    ])
    mut = np.array([1]*n1 + [0]*n2 + [1]*n3 + [0]*n4)
    crc = np.array([1]*n1 + [1]*n2 + [0]*n3 + [0]*n4)
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o
def clip_obs(y, cens):
    o = y.copy(); o[cens] = c; return o
def aac_obs(y, cens=None):
    ic50 = np.exp(y); return 1.0/(1.0+(np.exp(c)/ic50)**0.5)

def fit_ols(y, mut, crc, cov='HC3'):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov)
    return m.params[3], m.pvalues[3], np.sqrt(m.mse_resid)

def fit_rank(y, mut, crc):
    a = y[(mut==1) & (crc==1)]; b = y[(mut==1) & (crc==0)]
    if len(a) == 0 or len(b) == 0:
        return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided')
    return p

N = 2000
MC_SE = np.sqrt(0.05*0.95/N)

n_settings = [
    ("n=(6,40,38,17)", (6, 40, 38, 17)),
    ("n=(20,20,20,20)", (20, 20, 20, 20)),
    ("n=(40,40,40,40)", (40, 40, 40, 40)),
]
rates_wt = [0.20, 0.50, 0.80, 0.90, 0.945]

print("=== 二维网格 (删失率 × n_min)：分离病因 ===")
print(f"mut 格删失率固定 {rate_mut:.0%}；扫 WT 格删失率 × 四格 n；N={N}，MC SE≈{MC_SE:.4f}\n")

results = {}
for nname, ns in n_settings:
    for rw in rates_wt:
        mu_wt = mu_for(rw, sigma_wt)
        row = {}
        for ename, obs in [("naive", naive_obs), ("封顶", clip_obs), ("AAC", aac_obs)]:
            ps0 = np.zeros(N); sr = []
            for i in range(N):
                y, mut, crc, cens = gen_data(0.0, ns, mu_wt)
                yo = obs(y, cens)
                b, p, sh = fit_ols(yo, mut, crc, 'HC3')
                ps0[i] = p
                if ename != "AAC":
                    sr.append(sh / sigma_mut)
            fpr = (ps0 < 0.05).mean()
            thresh = np.quantile(ps0, 0.05)
            ps1 = np.zeros(N)
            for i in range(N):
                y, mut, crc, cens = gen_data(1.8, ns, mu_wt)
                yo = obs(y, cens)
                b, p, _ = fit_ols(yo, mut, crc, 'HC3')
                ps1[i] = p
            pwr = (ps1 < thresh).mean()
            row[ename] = (fpr, pwr, np.mean(sr) if sr else np.nan)
        # mut-only 秩
        ps0 = np.array([fit_rank(*gen_data(0.0, ns, mu_wt)[:3]) for _ in range(N)])
        thresh = np.quantile(ps0, 0.05)
        ps1 = np.array([fit_rank(*gen_data(1.8, ns, mu_wt)[:3]) for _ in range(N)])
        row["mut秩"] = ((ps0 < 0.05).mean(), (ps1 < thresh).mean(), np.nan)
        results[(nname, rw)] = row
        print(f"{nname:18s} WT删失={rw:.3f}:")
        for ename in ["naive", "封顶", "AAC", "mut秩"]:
            fpr, pwr, sr = row[ename]
            srs = f"  σ̂比={sr:.2f}" if not np.isnan(sr) else ""
            print(f"    {ename:6s}: FPR={fpr:.3f}  功效={pwr:.3f}{srs}")

print("\n=== 病因判读 ===")
# 关键对照：固定 WT 删失 94.5%，对比 n 设置
rw_key = 0.945
print(f"WT 删失={rw_key} 下，n 设置对比：")
for ename in ["naive", "封顶", "mut秩"]:
    print(f"  {ename:6s}: ", end="")
    for nname, ns in n_settings:
        fpr, pwr, sr = results[(nname, rw_key)][ename]
        print(f"{nname.split('=')[1]:14s} FPR={fpr:.3f}/功效={pwr:.3f}  ", end="")
    print()
