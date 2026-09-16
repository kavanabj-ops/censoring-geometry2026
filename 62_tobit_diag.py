# -*- coding: utf-8 -*-
"""
62_tobit_diag.py — Tobit 三项诊断（opus5 Q2）
(a) 同方差 DGP + 同方差 Tobit：若 β̂偏倚消失 → "Tobit 失效由 σ 异质驱动"
对比之前的异方差 DGP 结果。
"""
import numpy as np
import scipy.stats as st
import scipy.optimize as opt

rng = np.random.default_rng(99)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)

def tobit_fit(y, mut, crc, c):
    X = np.column_stack([np.ones(len(y)), mut, crc, mut*crc])
    censored = (y >= c).astype(int)
    y_obs = np.minimum(y, c)
    def negloglik(params):
        beta = params[:-1]
        sigma = np.exp(params[-1])
        if sigma < 1e-4:
            return 1e10
        mu = X @ beta
        z = (y_obs - mu) / sigma
        zc = (c - mu) / sigma
        lik = np.where(censored == 1, st.norm.logsf(zc), st.norm.logpdf(z) - np.log(sigma))
        return -lik.sum()
    init = np.zeros(X.shape[1] + 1)
    init[-1] = np.log(max(np.std(y_obs), 0.5))
    try:
        res = opt.minimize(negloglik, init, method='BFGS', options={'maxiter': 1000})
    except Exception:
        return None, False, np.nan
    if not res.success:
        return None, False, np.nan
    beta = res.x[:-1]
    sigma = np.exp(res.x[-1])
    if sigma > 20 or sigma < 1e-3 or not np.all(np.isfinite(beta)):
        return None, False, sigma
    try:
        se = np.sqrt(np.diag(res.hess_inv))
    except Exception:
        return beta[3], False, sigma
    if not np.all(np.isfinite(se)) or se[3] <= 0:
        return beta[3], False, sigma
    z = beta[3] / se[3]
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return beta[3], True, p

def run(n_iter, sigma_mut, sigma_wt, beta_int=0.0):
    mu_mut = c + st.norm.ppf(0.20) * sigma_mut
    mu_wt = c + st.norm.ppf(0.945) * sigma_wt
    n_fail = 0
    bias_list = []
    fpr_count = 0
    n_valid = 0
    for i in range(n_iter):
        y = np.concatenate([
            rng.normal(mu_mut + beta_int, sigma_mut, n_crc_mut),
            rng.normal(mu_wt, sigma_wt, n_crc_wt),
            rng.normal(mu_mut, sigma_mut, n_skin_mut),
            rng.normal(mu_wt, sigma_wt, n_skin_wt)])
        mut = np.array([1]*n_crc_mut + [0]*n_crc_wt + [1]*n_skin_mut + [0]*n_skin_wt)
        crc = np.array([1]*n_crc_mut + [1]*n_crc_wt + [0]*n_skin_mut + [0]*n_skin_wt)
        b, conv, p = tobit_fit(y, mut, crc, c)
        if not conv:
            n_fail += 1
        else:
            bias_list.append(b - beta_int)
            n_valid += 1
            if p < 0.05:
                fpr_count += 1
    fpr = fpr_count / max(n_valid, 1)
    return n_fail/n_iter, np.mean(bias_list), fpr

N = 1000
print("=== Tobit 诊断 (a)：同方差 vs 异方差 DGP ===\n")

print("[同方差 DGP] σ_mut=σ_wt=1.0，删失率 20%/94.5%")
fail, bias, fpr = run(N, 1.0, 1.0)
print(f"  收敛失败率={fail:.1%}, β̂偏倚={bias:+.3f}, FPR={fpr:.3f}")

print("\n[异方差 DGP] σ_mut=1.5, σ_wt=0.8（之前的 59b）")
fail2, bias2, fpr2 = run(N, 1.5, 0.8)
print(f"  收敛失败率={fail2:.1%}, β̂偏倚={bias2:+.3f}, FPR={fpr2:.3f}")

print("\n[中异方差 DGP] σ_mut=1.3, σ_wt=1.0")
fail3, bias3, fpr3 = run(N, 1.3, 1.0)
print(f"  收敛失败率={fail3:.1%}, β̂偏倚={bias3:+.3f}, FPR={fpr3:.3f}")

print("\n结论：若同方差下 β̂偏倚≈0，则证实'偏倚由 σ 异质 + 边界非识别驱动，而非删失本身'")
