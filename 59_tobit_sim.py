# -*- coding: utf-8 -*-
"""
59_tobit_sim.py — Tobit 操作特性（opus5 第三轮 Q1：决定 E1-E5 层级）
报：正确指定下 FPR/size-adj功效、收敛失败率、误指定(t df=3)下 FPR、β̂偏倚。
"""
import numpy as np
import scipy.stats as st
import scipy.optimize as opt
import statsmodels.api as sm

rng = np.random.default_rng(99)

n_crc_mut, n_crc_wt, n_skin_mut, n_skin_wt = 6, 40, 38, 17
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
target_mut, target_wt = 0.20, 0.945
mu_mut = c + st.norm.ppf(target_mut) * sigma_mut
mu_wt = c + st.norm.ppf(target_wt) * sigma_wt

def gen_data(beta_int, dist='normal'):
    if dist == 'normal':
        y = np.concatenate([
            rng.normal(mu_mut + beta_int, sigma_mut, n_crc_mut),
            rng.normal(mu_wt, sigma_wt, n_crc_wt),
            rng.normal(mu_mut, sigma_mut, n_skin_mut),
            rng.normal(mu_wt, sigma_wt, n_skin_wt)])
    else:  # t(df=3)，重尾
        y = np.concatenate([
            mu_mut + beta_int + rng.standard_t(3, n_crc_mut) * sigma_mut / 1.73,
            mu_wt + rng.standard_t(3, n_crc_wt) * sigma_wt / 1.73,
            mu_mut + rng.standard_t(3, n_skin_mut) * sigma_mut / 1.73,
            mu_wt + rng.standard_t(3, n_skin_wt) * sigma_wt / 1.73])
    mut = np.array([1]*n_crc_mut + [0]*n_crc_wt + [1]*n_skin_mut + [0]*n_skin_wt)
    crc = np.array([1]*n_crc_mut + [1]*n_crc_wt + [0]*n_skin_mut + [0]*n_skin_wt)
    return y, mut, crc

def tobit_fit(y, mut, crc, c):
    """右删失 Tobit，返回 (β_int, p, 收敛, σ)"""
    X = np.column_stack([np.ones(len(y)), mut, crc, mut*crc])
    censored = (y >= c).astype(int)
    y_obs = np.minimum(y, c)

    def negloglik(params):
        beta = params[:-1]
        log_sigma = params[-1]
        sigma = np.exp(log_sigma)
        if sigma < 1e-4:
            return 1e10
        mu = X @ beta
        z = (y_obs - mu) / sigma
        zc = (c - mu) / sigma
        lik = np.where(censored == 1, st.norm.logsf(zc), st.norm.logpdf(z) - log_sigma)
        return -lik.sum()

    init = np.zeros(X.shape[1] + 1)
    init[-1] = np.log(max(np.std(y_obs), 0.5))
    try:
        res = opt.minimize(negloglik, init, method='Nelder-Mead',
                           options={'maxiter': 5000, 'xatol': 1e-5, 'fatol': 1e-5})
    except Exception:
        return None, 1.0, False, np.nan
    if not res.success:
        return None, 1.0, False, np.nan
    beta = res.x[:-1]
    sigma = np.exp(res.x[-1])
    # 收敛失败/σ爆掉
    if sigma > 20 or sigma < 1e-3 or not np.all(np.isfinite(beta)):
        return None, 1.0, False, sigma
    # Wald 检验 β_int（用数值 Hessian）
    eps = 1e-4
    H = np.zeros((len(res.x), len(res.x)))
    for i in range(len(res.x)):
        for j in range(len(res.x)):
            x1 = res.x.copy(); x1[i] += eps; x1[j] += eps
            x2 = res.x.copy(); x2[i] += eps; x2[j] -= eps
            x3 = res.x.copy(); x3[i] -= eps; x3[j] += eps
            x4 = res.x.copy(); x4[i] -= eps; x4[j] -= eps
            H[i,j] = (negloglik(x1)-negloglik(x2)-negloglik(x3)+negloglik(x4))/(4*eps**2)
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        return beta[3], 1.0, False, sigma
    if se[3] <= 0 or not np.isfinite(se[3]):
        return beta[3], 1.0, False, sigma
    z = beta[3] / se[3]
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return beta[3], p, True, sigma

def run(n_iter, beta_int, dist='normal'):
    ps = np.zeros(n_iter)
    n_conv_fail = 0
    bias_list = []
    for i in range(n_iter):
        y, mut, crc = gen_data(beta_int, dist)
        b, p, conv, sigma = tobit_fit(y, mut, crc, c)
        if not conv:
            n_conv_fail += 1
            ps[i] = 1.0  # 不收敛当作"无法拒绝"，保守
            bias_list.append(np.nan)
        else:
            ps[i] = p
            bias_list.append(b - beta_int)
    return ps, n_conv_fail, np.nanmean(bias_list)

N = 2000
print("=== Tobit 操作特性（N=2000）===\n")

# 正确指定：高斯
ps0, fail0, bias0 = run(N, 0.0, 'normal')
ps1, fail1, bias1 = run(N, 1.8, 'normal')
fpr = (ps0 < 0.05).mean()
thresh = np.quantile(ps0[ps0 < 1.0], 0.05)
pwr = (ps1 < thresh).mean()
print(f"[正确指定 高斯]")
print(f"  FPR={fpr:.3f}, size-adj功效={pwr:.3f}")
print(f"  收敛失败率: 零假设{fail0/N:.1%}, 备择{fail1/N:.1%}")
print(f"  β̂偏倚={bias0:+.3f}")

# 误指定：t(df=3)
ps0t, fail0t, bias0t = run(N, 0.0, 't')
fpr_t = (ps0t < 0.05).mean()
print(f"\n[误指定 t(df=3)，仍用高斯 Tobit]")
print(f"  FPR={fpr_t:.3f}（正确指定下是 {fpr:.3f}）")
print(f"  收敛失败率: {fail0t/N:.1%}")
print(f"  β̂偏倚={bias0t:+.3f}")
