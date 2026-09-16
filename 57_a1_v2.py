# -*- coding: utf-8 -*-
"""
57_a1_v2.py — A1 重跑（含 opus5 第二轮 must-fix）
修正：删失率精确校准（skin-mut 21.1%、CRC-WT 95%）；size-adjusted power；
β̂偏倚 + σ̂/σ_true 比值诊断；AAC 加 HC3 臂；mutant-only 秩检验。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm

rng = np.random.default_rng(42)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)  # 2.303

# 精确删失率校准：给定 σ，反推 μ 使删失率 = 目标
sigma_mut, sigma_wt = 1.5, 0.8
target_mut, target_wt = 0.20, 0.945  # mut 格取 20%（14.3 与 21.1 之间），WT 格 94.5%
mu_mut = c + st.norm.ppf(target_mut) * sigma_mut
mu_wt = c + st.norm.ppf(target_wt) * sigma_wt

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
    o = y.copy()
    o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum())
    return o

def clip_obs(y, cens):
    o = y.copy()
    o[cens] = c
    return o

def aac_obs(y, cens=None):
    ic50 = np.exp(y)
    return 1.0 / (1.0 + (np.exp(c)/ic50)**0.5)

def fit_ols(y, mut, crc, cov_type='HC0'):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov_type)
    return m.params[3], m.pvalues[3], np.sqrt(m.mse_resid)

def fit_mutonly_rank(y, mut, crc):
    """突变格内 Mann-Whitney：CRC-mut vs skin-mut"""
    a = y[(mut==1) & (crc==1)]
    b = y[(mut==1) & (crc==0)]
    if len(a) == 0 or len(b) == 0:
        return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided')
    return p

def run_batch(n_iter, beta_int, obs_fn, fit_fn, collect_diag=False):
    """返回 p 值数组；collect_diag=True 时额外返回 β̂偏倚和 σ̂比值"""
    ps = np.zeros(n_iter)
    biases, sratio = None, None
    if collect_diag:
        bias_list, sr_list = [], []
    for i in range(n_iter):
        y, mut, crc, cens = gen_data(beta_int)
        y_obs = obs_fn(y, cens)
        b, p, s_hat = fit_fn(y_obs, mut, crc)
        ps[i] = p
        if collect_diag:
            bias_list.append(b - beta_int)
            sr_list.append(s_hat / sigma_mut)
    if collect_diag:
        biases = np.mean(bias_list)
        sratio = np.mean(sr_list)
    return ps, biases, sratio

N = 8000
print(f"=== A1 v2（精确校准）===")
print(f"四格 n=(6,40,38,17)；μ_mut={mu_mut:.3f}(删失{target_mut:.0%}), μ_wt={mu_wt:.3f}(删失{target_wt:.0%})")
print(f"σ_mut={sigma_mut}, σ_wt={sigma_wt}\n")

estimators = [
    ("naive OLS", naive_obs, lambda y,m,cr: fit_ols(y,m,cr,'HC0')),
    ("封顶 OLS", clip_obs, lambda y,m,cr: fit_ols(y,m,cr,'HC0')),
    ("AAC(HC0)", aac_obs, lambda y,m,cr: fit_ols(y,m,cr,'HC0')),
    ("AAC(HC3)", aac_obs, lambda y,m,cr: fit_ols(y,m,cr,'HC3')),
]

# 零假设 p 值分布（用于 size-adjusted 阈值）+ 诊断
null_ps = {}
null_diag = {}
for name, obs, fit in estimators:
    ps, bias, sr = run_batch(N, 0.0, obs, fit, collect_diag=True)
    null_ps[name] = ps
    null_diag[name] = (bias, sr)

# mutant-only 秩检验零分布
null_rank = np.array([fit_mutonly_rank(*gen_data(0.0)[:3]) for _ in range(N)])

print("--- 零假设经验 FPR（名义 0.05，MC SE≈0.0024）---")
for name, obs, fit in estimators:
    fpr = (null_ps[name] < 0.05).mean()
    bias, sr = null_diag[name]
    print(f"  {name:10s}: FPR={fpr:.3f}, β̂偏倚={bias:+.3f}, σ̂/σ_true={sr:.3f}")
fpr_rank = (null_rank < 0.05).mean()
print(f"  {'mut-only秩':10s}: FPR={fpr_rank:.3f}")

# size-adjusted 阈值（各估计量零分布的经验 5% 分位）
print("\n--- size-adjusted 功效（β=1.8，用各自零分布 5% 分位作阈值）---")
for name, obs, fit in estimators:
    thresh = np.quantile(null_ps[name], 0.05)
    alt_ps, _, _ = run_batch(N, 1.8, obs, fit)
    power = (alt_ps < thresh).mean()
    raw_power = (alt_ps < 0.05).mean()
    print(f"  {name:10s}: size-adj功效={power:.3f}, 原始功效={raw_power:.3f}")
thresh_rank = np.quantile(null_rank, 0.05)
alt_rank = np.array([fit_mutonly_rank(*gen_data(1.8)[:3]) for _ in range(N)])
print(f"  {'mut-only秩':10s}: size-adj功效={(alt_rank < thresh_rank).mean():.3f}")

# 保存零分布供 calib 参考
np.save(r'D:\paper\02_analysis\results\a1_null_ps.npy', {k: v for k, v in null_ps.items()})
