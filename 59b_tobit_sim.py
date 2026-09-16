# -*- coding: utf-8 -*-
"""
59b_tobit_sim.py — Tobit 操作特性（BFGS 快速版）
用 BFGS（自带逆 Hessian）+ N=1000，快。
"""
import numpy as np
import scipy.stats as st
import scipy.optimize as opt

rng = np.random.default_rng(99)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
mu_mut = c + st.norm.ppf(0.20) * sigma_mut
mu_wt = c + st.norm.ppf(0.945) * sigma_wt

def gen_data(beta_int, dist='normal'):
    if dist == 'normal':
        y = np.concatenate([
            rng.normal(mu_mut + beta_int, sigma_mut, n_crc_mut),
            rng.normal(mu_wt, sigma_wt, n_crc_wt),
            rng.normal(mu_mut, sigma_mut, n_skin_mut),
            rng.normal(mu_wt, sigma_wt, n_skin_wt)])
    else:
        y = np.concatenate([
            mu_mut + beta_int + rng.standard_t(3, n_crc_mut) * sigma_mut / 1.73,
            mu_wt + rng.standard_t(3, n_crc_wt) * sigma_wt / 1.73,
            mu_mut + rng.standard_t(3, n_skin_mut) * sigma_mut / 1.73,
            mu_wt + rng.standard_t(3, n_skin_wt) * sigma_wt / 1.73])
    mut = np.array([1]*n_crc_mut + [0]*n_crc_wt + [1]*n_skin_mut + [0]*n_skin_wt)
    crc = np.array([1]*n_crc_mut + [1]*n_crc_wt + [0]*n_skin_mut + [0]*n_skin_wt)
    return y, mut, crc

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
        return None, 1.0, False, np.nan
    if not res.success:
        return None, 1.0, False, np.nan
    beta = res.x[:-1]
    sigma = np.exp(res.x[-1])
    if sigma > 20 or sigma < 1e-3 or not np.all(np.isfinite(beta)):
        return None, 1.0, False, sigma
    # SE 从 BFGS 的逆 Hessian
    try:
        hess_inv = res.hess_inv
        se = np.sqrt(np.diag(hess_inv))
    except Exception:
        return beta[3], 1.0, False, sigma
    if not np.all(np.isfinite(se)) or se[3] <= 0:
        return beta[3], 1.0, False, sigma
    z = beta[3] / se[3]
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return beta[3], p, True, sigma

def run(n_iter, beta_int, dist='normal'):
    ps = np.zeros(n_iter)
    n_fail = 0
    bias_list = []
    for i in range(n_iter):
        y, mut, crc = gen_data(beta_int, dist)
        b, p, conv, sigma = tobit_fit(y, mut, crc, c)
        if not conv:
            n_fail += 1
            ps[i] = 1.0
            bias_list.append(np.nan)
        else:
            ps[i] = p
            bias_list.append(b - beta_int)
    return ps, n_fail, np.nanmean(bias_list)

N = 1000
print("=== Tobit 操作特性（BFGS，N=1000）===\n")
ps0, fail0, bias0 = run(N, 0.0, 'normal')
ps1, fail1, bias1 = run(N, 1.8, 'normal')
fpr = (ps0 < 0.05).mean()
thresh = np.quantile(ps0[ps0 < 1.0], 0.05)
pwr = (ps1 < thresh).mean()
print(f"[正确指定 高斯] FPR={fpr:.3f}, size-adj功效={pwr:.3f}")
print(f"  收敛失败率: 零假设{fail0/N:.1%}, 备择{fail1/N:.1%}, β̂偏倚={bias0:+.3f}")

ps0t, fail0t, bias0t = run(N, 0.0, 't')
fpr_t = (ps0t < 0.05).mean()
print(f"\n[误指定 t(df=3)] FPR={fpr_t:.3f}（正确指定 {fpr:.3f}）")
print(f"  收敛失败率: {fail0t/N:.1%}, β̂偏倚={bias0t:+.3f}")
