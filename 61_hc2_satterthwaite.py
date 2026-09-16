# -*- coding: utf-8 -*-
"""
61_hc2_satterthwaite.py — HC2 + Satterthwaite 自由度（opus5 第四轮第1步）
测 AAC/naive/封顶 在 HC2 下的 FPR，以及 Welch-Satterthwaite 有效自由度。
目的：AAC+HC3 FPR=0.068 超标(+6.4 SE)，看 HC2/Satterthwaite 能否把 AAC 校准到 0.05。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm

rng = np.random.default_rng(2024)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
mu_mut = c + st.norm.ppf(0.20) * sigma_mut
mu_wt = c + st.norm.ppf(0.945) * sigma_wt

def gen_data(beta_int):
    y = np.concatenate([
        rng.normal(mu_mut + beta_int, sigma_mut, n_crc_mut),
        rng.normal(mu_wt, sigma_wt, n_crc_wt),
        rng.normal(mu_mut, sigma_mut, n_skin_mut),
        rng.normal(mu_wt, sigma_wt, n_skin_wt)])
    mut = np.array([1]*n_crc_mut + [0]*n_crc_wt + [1]*n_skin_mut + [0]*n_skin_wt)
    crc = np.array([1]*n_crc_mut + [1]*n_crc_wt + [0]*n_skin_mut + [0]*n_skin_wt)
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o
def clip_obs(y, cens):
    o = y.copy(); o[cens] = c; return o
def aac_obs(y, cens=None):
    ic50 = np.exp(y); return 1.0/(1.0+(np.exp(c)/ic50)**0.5)

def fit(y, mut, crc, cov_type):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov_type)
    return m.params[3], m.pvalues[3]

def fit_satterthwaite(y, mut, crc):
    """Welch-Satterthwaite 对交互对比：用四格组内 SD + 有效自由度"""
    groups = {}
    for (mi, ci) in [(1,1),(1,0),(0,1),(0,0)]:
        g = y[(mut==mi)&(crc==ci)]
        groups[(mi,ci)] = g
    # 交互对比：β_int = (μ_11 - μ_10) - (μ_01 - μ_00)
    s2 = {k: np.var(v, ddof=1) for k, v in groups.items()}
    n = {k: len(v) for k, v in groups.items()}
    b_int = (np.mean(groups[(1,1)]) - np.mean(groups[(1,0)])) - (np.mean(groups[(0,1)]) - np.mean(groups[(0,0)]))
    # SE（不等方差）
    se2 = sum(s2[k]/n[k] for k in groups)
    se = np.sqrt(se2)
    # Satterthwaite 自由度
    num = se2**2
    den = sum((s2[k]/n[k])**2 / (n[k]-1) for k in groups if n[k] > 1)
    df = num / den if den > 0 else 1.0
    t = b_int / se
    p = 2 * (1 - st.t.cdf(abs(t), df))
    return p

N = 8000
print("=== HC2 + Satterthwaite（决定 E1/E2 身份）===\n")

for name, obs in [("naive", naive_obs), ("封顶", clip_obs), ("AAC", aac_obs)]:
    for cov in ['HC0', 'HC2', 'HC3']:
        ps = np.zeros(N)
        for i in range(N):
            y, mut, crc, cens = gen_data(0.0)
            y_obs = obs(y, cens)
            _, p = fit(y_obs, mut, crc, cov)
            ps[i] = p
        fpr = (ps < 0.05).mean()
        se_mc = np.sqrt(0.05*0.95/N)
        nse = (fpr - 0.05) / se_mc
        print(f"  {name:6s} {cov:4s}: FPR={fpr:.3f} ({nse:+.1f} SE)")
    # Satterthwaite
    ps = np.zeros(N)
    for i in range(N):
        y, mut, crc, cens = gen_data(0.0)
        y_obs = obs(y, cens)
        ps[i] = fit_satterthwaite(y_obs, mut, crc)
    fpr = (ps < 0.05).mean()
    nse = (fpr - 0.05) / np.sqrt(0.05*0.95/N)
    print(f"  {name:6s} Satt: FPR={fpr:.3f} ({nse:+.1f} SE)")
    print()
