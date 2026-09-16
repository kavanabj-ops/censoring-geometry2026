# -*- coding: utf-8 -*-
"""
55_a1_censoring_sim.py — A1 删失几何模拟（校准版）
修正：删失率校准到实测（mut 14-21%, WT 94-95%）；AAC 用更平缓的 AUC 定义（避免耐药区饱和）。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm

rng = np.random.default_rng(42)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)  # 2.303

# 校准参数：使删失率匹配实测
mu_mut, sigma_mut = 0.93, 1.5   # mut 格删失率 ≈ 18%
mu_wt, sigma_wt = 3.58, 0.8     # WT 格删失率 ≈ 94.5%

def gen_data(beta_int):
    y = np.concatenate([
        rng.normal(mu_mut + beta_int, sigma_mut, n_crc_mut),
        rng.normal(mu_wt, sigma_wt, n_crc_wt),
        rng.normal(mu_mut, sigma_mut, n_skin_mut),
        rng.normal(mu_wt, sigma_wt, n_skin_wt),
    ])
    mut = np.array([1]*n_crc_mut + [0]*n_crc_wt + [1]*n_skin_mut + [0]*n_skin_wt)
    crc = np.array([1]*n_crc_mut + [1]*n_crc_wt + [0]*n_skin_mut + [0]*n_skin_wt)
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    y_naive = y.copy()
    y_naive[cens] = y[cens] + rng.normal(0, 2.0, cens.sum())
    return y_naive

def clip_obs(y, cens):
    y_clip = y.copy()
    y_clip[cens] = c
    return y_clip

def aac_obs(y, cens=None):
    """真实 AUC 形式：AUC = 1/(1+(c/IC50)^0.5)，耐药区仍有区分度（不饱和）"""
    ic50 = np.exp(y)
    return 1.0 / (1.0 + (np.exp(c)/ic50)**0.5)

def fit_ols(y, mut, crc):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit()
    return m.pvalues[3]

def run_sim(n_iter, beta_int, obs_fn):
    n_sig = 0
    for _ in range(n_iter):
        y, mut, crc, cens = gen_data(beta_int)
        y_obs = obs_fn(y, cens)
        if fit_ols(y_obs, mut, crc) < 0.05:
            n_sig += 1
    return n_sig / n_iter

N = 5000
print("=== A1 删失几何模拟（校准版）===\n")
y, mut, crc, cens = gen_data(0.0)
for name, cond in [("CRC-mut", (mut==1)&(crc==1)), ("CRC-WT", (mut==0)&(crc==1)),
                   ("skin-mut", (mut==1)&(crc==0)), ("skin-WT", (mut==0)&(crc==0))]:
    print(f"  {name}: n={cond.sum()}, 模拟删失率={cens[cond].mean():.1%} (实测 14.3/95.0/21.1/94.1)")

print("\n--- 零假设 β_int=0：经验 FPR（应≈0.05）---")
for name, obs_fn in [("naive OLS(外推)", naive_obs), ("封顶 OLS", clip_obs), ("AAC(有界)", aac_obs)]:
    print(f"  {name:16s}: FPR = {run_sim(N, 0.0, obs_fn):.3f}")

print("\n--- 备择 β_int=1.8：功效 ---")
for name, obs_fn in [("naive OLS(外推)", naive_obs), ("封顶 OLS", clip_obs), ("AAC(有界)", aac_obs)]:
    print(f"  {name:16s}: 功效 = {run_sim(N, 1.8, obs_fn):.3f}")
