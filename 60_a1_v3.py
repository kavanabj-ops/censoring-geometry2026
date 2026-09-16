# -*- coding: utf-8 -*-
"""
60_a1_v3.py — A1 修补（opus5 第三轮 must-fix a/c）
4×2 网格（估计量 naive/封顶/AAC/mut-only × SE 同方差/HC3）+ sanity check。
关键：验证 naive 在 HC3 下"保守"（σ̂膨胀→SE大）、HC0 下"反保守"（小样本高杠杆向下偏SE）。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm

rng = np.random.default_rng(2024)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
target_mut, target_wt = 0.20, 0.945
mu_mut = c + st.norm.ppf(target_mut) * sigma_mut
mu_wt = c + st.norm.ppf(target_wt) * sigma_wt

def gen_data(beta_int, mu_mut_=None, mu_wt_=None, sigma_mut_=None):
    mm = mu_mut if mu_mut_ is None else mu_mut_
    mw = mu_wt if mu_wt_ is None else mu_wt_
    smm = sigma_mut if sigma_mut_ is None else sigma_mut_
    y = np.concatenate([
        rng.normal(mm + beta_int, smm, n_crc_mut),
        rng.normal(mw, sigma_wt, n_crc_wt),
        rng.normal(mm, smm, n_skin_mut),
        rng.normal(mw, sigma_wt, n_skin_wt)])
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
    return m.params[3], m.pvalues[3], np.sqrt(m.mse_resid)

def fit_rank(y, mut, crc):
    a = y[(mut==1)&(crc==1)]; b = y[(mut==1)&(crc==0)]
    if len(a)==0 or len(b)==0: return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided'); return p

# 4×2 网格
grid = [
    ("naive", naive_obs, "ln"),
    ("封顶", clip_obs, "ln"),
    ("AAC", aac_obs, "aac"),
]

N = 6000
print("=== A1 v3：4×2 网格（估计量 × SE方法）===\n")

results = {}
for name, obs, scale in grid:
    for cov in ['HC0', 'HC3']:
        # 零假设
        ps0 = np.zeros(N); sr_list = []
        for i in range(N):
            y, mut, crc, cens = gen_data(0.0)
            y_obs = obs(y, cens)
            b, p, s_hat = fit(y_obs, mut, crc, cov)
            ps0[i] = p
            # σ_true：ln 尺度用 sigma_mut，aac 尺度用 aac 的真实 SD
            if scale == 'ln':
                sr_list.append(s_hat / sigma_mut)
            else:
                sr_list.append(s_hat)  # aac 尺度，后面单独处理
        fpr = (ps0 < 0.05).mean()
        thresh = np.quantile(ps0, 0.05)
        # 备择
        ps1 = np.zeros(N)
        for i in range(N):
            y, mut, crc, cens = gen_data(1.8)
            y_obs = obs(y, cens)
            b, p, _ = fit(y_obs, mut, crc, cov)
            ps1[i] = p
        pwr = (ps1 < thresh).mean()
        sr = np.mean(sr_list)
        results[f"{name}-{cov}"] = (fpr, pwr, sr)
        print(f"  {name:6s} {cov:4s}: FPR={fpr:.3f}, size-adj功效={pwr:.3f}, σ̂比={sr:.3f}" + (" (ln尺度)" if scale=='ln' else " (AAC尺度)"))

# mut-only 秩
ps0 = np.array([fit_rank(*gen_data(0.0)[:3]) for _ in range(N)])
ps1 = np.array([fit_rank(*gen_data(1.8)[:3]) for _ in range(N)])
fpr = (ps0 < 0.05).mean(); thresh = np.quantile(ps0, 0.05); pwr = (ps1 < thresh).mean()
print(f"  {'mut秩':6s} {'--':4s}: FPR={fpr:.3f}, size-adj功效={pwr:.3f}")
results['mut-rank'] = (fpr, pwr, np.nan)

print("\n=== Sanity check：完全无删失（c_mut=0, c_wt=0），所有臂应收敛 ===")
# 完全无删失：让 mut 和 WT 都远离删失点
mm0 = c - 5.0   # mut 格极敏感，删失率≈0
mw0 = c - 5.0   # WT 格也敏感（无删失）
for name, obs, scale in grid:
    # 零假设 FPR + σ̂ 比（应 FPR≈0.05，σ̂比≈1.0）
    ps = np.zeros(2000); sr_list = []
    for i in range(2000):
        y, mut, crc, cens = gen_data(0.0, mu_mut_=mm0, mu_wt_=mw0)
        y_obs = obs(y, cens)
        b, p, s_hat = fit(y_obs, mut, crc, 'HC0')
        ps[i] = p
        sr_list.append(s_hat / sigma_mut)
    print(f"  {name:6s} (无删失): FPR={ (ps<0.05).mean():.3f}, σ̂/σ_true={np.mean(sr_list):.3f} (应≈1.0)")
