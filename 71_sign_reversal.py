# -*- coding: utf-8 -*-
# 71_sign_reversal.py - 符号反转模拟：封顶 OLS 伪极限与真交互反号（opus5 major 7）
# 验证 Prop 1(iii) + Cor 1 符号翻转条件：delta_M delta_W > 0 且 Phi_W/Phi_M < |delta_M/delta_W| < 1
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(2024)

c = np.log(10)   # 2.303
sigma = 1.0

# 四格 mu（ln IC50 尺度）：delta_M = mu_AM - mu_BM = 1, delta_W = mu_AW - mu_BW = 3
mu_AM = c + sigma*st.norm.ppf(0.20)   # CRC-mut, 删失20%
mu_BM = mu_AM - 1.0                    # skin-mut  (delta_M = 1)
mu_AW = c + sigma*st.norm.ppf(0.90)   # CRC-WT, 删失90%
mu_BW = mu_AW - 3.0                    # skin-WT  (delta_W = 3)

beta_true = (mu_AM - mu_AW) - (mu_BM - mu_BW)   # = delta_M - delta_W = -2
print(f"真交互 beta_int = delta_M - delta_W = {mu_AM-mu_BM:.1f} - {mu_AW-mu_BW:.1f} = {beta_true:.2f}")

# 解析伪极限（封顶 E[y^c]）
def E_capped(mu):
    z = (c - mu)/sigma
    return c - sigma*(st.norm.pdf(z) + z*st.norm.cdf(z))
beta_tilde = (E_capped(mu_AM) - E_capped(mu_AW)) - (E_capped(mu_BM) - E_capped(mu_BW))
print(f"封顶 OLS 解析伪极限 beta_tilde = {beta_tilde:.3f}  (与真值反号: {np.sign(beta_tilde) != np.sign(beta_true)})")

n = (6, 40, 38, 17)
def gen():
    y = np.concatenate([
        rng.normal(mu_AM, sigma, n[0]),
        rng.normal(mu_AW, sigma, n[1]),
        rng.normal(mu_BM, sigma, n[2]),
        rng.normal(mu_BW, sigma, n[3]),
    ])
    mut = np.array([1]*n[0] + [0]*n[1] + [1]*n[2] + [0]*n[3])
    crc = np.array([1]*n[0] + [1]*n[1] + [0]*n[2] + [0]*n[3])
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o

N = 4000
X_cache = None
for name in ['true(无删失)', 'capped(封顶)', 'naive(外推)']:
    bs = np.zeros(N)
    for i in range(N):
        y, mut, crc, cens = gen()
        if name == 'true(无删失)':
            yo = y
        elif name == 'capped(封顶)':
            yo = np.minimum(y, c)
        else:
            yo = naive_obs(y, cens)
        X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
        m = sm.OLS(yo, X).fit()
        bs[i] = m.params[3]
    flip = (np.sign(bs) != np.sign(beta_true)).mean()
    print(f"{name:14s}: beta_hat 均值={np.mean(bs):+.3f}, 反号比例={flip:.1%}, 95%CI=[{np.quantile(bs,0.025):+.2f},{np.quantile(bs,0.975):+.2f}]")
