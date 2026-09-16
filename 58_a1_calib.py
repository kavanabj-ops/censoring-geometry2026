# -*- coding: utf-8 -*-
"""
58_a1_calib.py — A1-calib 扫描（opus5 Q2 核心交付物）
扫描 WT 格删失率 c_wt 从 0 到 1.0，报各估计量的 FPR 与 size-adjusted power。
产出：FPR/power ~ c_wt 曲线 + 阈值 c*/c** + σ̂/σ_true 比值。
"""
import numpy as np
import scipy.stats as st
import statsmodels.api as sm

rng = np.random.default_rng(7)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
target_mut = 0.20
mu_mut = c + st.norm.ppf(target_mut) * sigma_mut

def mu_wt_for(cens_rate):
    return c + st.norm.ppf(cens_rate) * sigma_wt

def gen_data(beta_int, mu_wt):
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
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o
def clip_obs(y, cens):
    o = y.copy(); o[cens] = c; return o
def aac_obs(y, cens=None):
    ic50 = np.exp(y); return 1.0/(1.0+(np.exp(c)/ic50)**0.5)

def fit(y, mut, crc, cov_type='HC3'):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov_type)
    return m.params[3], m.pvalues[3], np.sqrt(m.mse_resid)

def fit_rank(y, mut, crc):
    a = y[(mut==1)&(crc==1)]; b = y[(mut==1)&(crc==0)]
    if len(a)==0 or len(b)==0: return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided'); return p

def run_null(n, beta_int, mu_wt, obs_fn, fit_fn, collect_diag=False):
    ps = np.zeros(n); sr = None; bias = None
    if collect_diag:
        sr_l, b_l = [], []
    for i in range(n):
        y, mut, crc, cens = gen_data(beta_int, mu_wt)
        y_obs = obs_fn(y, cens)
        b, p, s = fit_fn(y_obs, mut, crc)
        ps[i] = p
        if collect_diag:
            sr_l.append(s/sigma_mut); b_l.append(b - beta_int)
    if collect_diag:
        sr = np.mean(sr_l); bias = np.mean(b_l)
    return ps, sr, bias

estimators = [
    ("naive", naive_obs, lambda y,m,cr: fit(y,m,cr,'HC3')),
    ("封顶", clip_obs, lambda y,m,cr: fit(y,m,cr,'HC3')),
    ("AAC", aac_obs, lambda y,m,cr: fit(y,m,cr,'HC3')),
]

N = 2000
c_wt_range = np.arange(0.0, 1.01, 0.1)

print("c_wt | naive FPR/pwr | 封顶 FPR/pwr | AAC FPR/pwr | σ̂比值(封顶/naive)")
results = []
for c_wt in c_wt_range:
    mu_wt = mu_wt_for(c_wt)
    row = {}
    line = f"{c_wt:.1f} | "
    for name, obs, fit_fn in estimators:
        # 零假设 FPR
        ps0, sr0, _ = run_null(N, 0.0, mu_wt, obs, fit_fn)
        fpr = (ps0 < 0.05).mean()
        thresh = np.quantile(ps0, 0.05)
        # 备择 size-adjusted power
        ps1, _, _ = run_null(N, 1.8, mu_wt, obs, fit_fn)
        pwr = (ps1 < thresh).mean()
        row[name+'_fpr'] = fpr; row[name+'_pwr'] = pwr
        line += f"{name} {fpr:.2f}/{pwr:.2f} | "
    # σ̂ 比值诊断（封顶 vs naive 在零假设下）
    _, sr_clip, _ = run_null(N, 0.0, mu_wt, clip_obs, lambda y,m,cr: fit(y,m,cr,'HC3'), collect_diag=True)
    _, sr_naive, _ = run_null(N, 0.0, mu_wt, naive_obs, lambda y,m,cr: fit(y,m,cr,'HC3'), collect_diag=True)
    row['sr_clip'] = sr_clip; row['sr_naive'] = sr_naive
    line += f"σ̂:封顶{sr_clip:.2f}/naive{sr_naive:.2f}"
    print(line)
    results.append(row)

# 保存
import json
with open(r'D:\paper\02_analysis\results\a1_calib.json', 'w') as f:
    json.dump({'c_wt': c_wt_range.tolist(), 'results': results}, f)
print("\n已保存 a1_calib.json")
