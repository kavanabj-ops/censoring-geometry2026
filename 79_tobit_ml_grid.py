# -*- coding: utf-8 -*-
# 79_tobit_ml_grid.py - Tobit 偏倚随有效格数 m_l 变化（opus5 M11：bias -> 0 当 m_l 增大）
import numpy as np, scipy.stats as st, scipy.optimize as opt

rng = np.random.default_rng(42)
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
mu_mut = c + st.norm.ppf(0.20) * sigma_mut   # mut 格删失 20%
mu_wt = c + st.norm.ppf(0.945) * sigma_wt    # WT 格删失 94.5%
beta_int = 0.0  # null

def tobit_fit(y, mut, crc, c):
    X = np.column_stack([np.ones(len(y)), mut, crc, mut*crc])
    censored = (y >= c).astype(int)
    y_obs = np.minimum(y, c)
    def nll(params):
        beta = params[:-1]
        sigma = np.exp(params[-1])
        if sigma < 1e-4:
            return 1e10
        mu = X @ beta
        z = (y_obs - mu)/sigma
        zc = (c - mu)/sigma
        lik = np.where(censored==1, st.norm.logsf(zc), st.norm.logpdf(z) - np.log(sigma))
        return -lik.sum()
    init = np.zeros(X.shape[1]+1)
    init[-1] = np.log(max(np.std(y_obs), 0.5))
    try:
        res = opt.minimize(nll, init, method='BFGS', options={'maxiter':1000})
    except Exception:
        return None, False
    if not res.success:
        return None, False
    beta = res.x[:-1]
    sigma = np.exp(res.x[-1])
    if sigma > 20 or sigma < 1e-3 or not np.all(np.isfinite(beta)):
        return None, False
    return beta[3], True

# 扫 m_WT = n_WT*(1-pi_WT)，pi_WT=0.945 -> (1-pi)=0.055
pi_wt = 1 - st.norm.cdf((c - mu_wt)/sigma_wt)
print(f"WT 格删失率 = {pi_wt:.3f}, 未删失率 = {1-pi_wt:.3f}")

N = 500
print("\n=== Tobit 偏倚 vs 有效 WT 格数 m_WT ===")
print(f"{'m_WT':>6} {'n_WT':>6} | {'收敛率':>8} {'β̂偏倚':>10}")
for m_wt in [2, 3.5, 5, 10, 20]:
    n_wt = int(round(m_wt / (1-pi_wt)))
    n_mut = 40  # mut 格固定
    biases = []
    n_conv = 0
    for i in range(N):
        y = np.concatenate([
            rng.normal(mu_mut, sigma_mut, n_mut),   # CRC-mut
            rng.normal(mu_wt, sigma_wt, n_wt),       # CRC-WT
            rng.normal(mu_mut, sigma_mut, n_mut),   # skin-mut
            rng.normal(mu_wt, sigma_wt, n_wt),       # skin-WT
        ])
        mut = np.array([1]*n_mut + [0]*n_wt + [1]*n_mut + [0]*n_wt)
        crc = np.array([1]*n_mut + [1]*n_wt + [0]*n_mut + [0]*n_wt)
        b, conv = tobit_fit(y, mut, crc, c)
        if conv:
            n_conv += 1
            biases.append(b - beta_int)
    conv_rate = n_conv/N
    bias_mean = np.mean(biases) if biases else np.nan
    print(f"{m_wt:6.1f} {n_wt:6d} | {conv_rate:8.1%} {bias_mean:+10.3f}")

print("\n预期：m_WT 增大 -> 收敛率 -> 100%，β̂ 偏倚 -> 0（有限样本边界事件，非结构性失效）")
