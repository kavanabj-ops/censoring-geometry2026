# -*- coding: utf-8 -*-
# 74_mc_se.py - 关键数字的 MC 不确定度（sigma_hat 比 + FPR 的 SE，回应 major 6）
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(2024)
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

def fit_ols(y, mut, crc, cov='HC3'):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov)
    return m.params[3], m.pvalues[3], np.sqrt(m.mse_resid)

N = 4000
n_settings = [
    ("n=(6,40,38,17)", (6, 40, 38, 17)),
    ("n=(20,20,20,20)", (20, 20, 20, 20)),
    ("n=(40,40,40,40)", (40, 40, 40, 40)),
]
mu_wt = mu_for(0.945, sigma_wt)  # WT 删失 94.5%

print("=== MC 不确定度（WT 删失 94.5%，N=4000）===")
for nname, ns in n_settings:
    for ename, obs in [("naive", naive_obs), ("封顶", clip_obs)]:
        sr_list = []
        fpr_list = []
        for i in range(N):
            y, mut, crc, cens = gen_data(0.0, ns, mu_wt)
            yo = obs(y, cens)
            b, p, sh = fit_ols(yo, mut, crc, 'HC3')
            sr_list.append(sh / sigma_mut)
            fpr_list.append(1.0 if p < 0.05 else 0.0)
        sr = np.array(sr_list)
        fpr = np.array(fpr_list)
        sr_mean, sr_se = sr.mean(), sr.std(ddof=1)/np.sqrt(N)
        fpr_mean, fpr_se = fpr.mean(), np.sqrt(fpr.mean()*(1-fpr.mean())/N)
        print(f"  {nname:18s} {ename:6s}: sigma_hat比={sr_mean:.3f} ± {sr_se:.3f}, FPR={fpr_mean:.3f} ± {fpr_se:.3f}")

print("\n判读：若三档 n 的 sigma_hat比 在各自 SE 内重叠，则 'constant across n' 成立")
