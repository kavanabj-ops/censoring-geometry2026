# -*- coding: utf-8 -*-
# 71c_sign_reversal_homoscedastic.py - 同方差符号反转（opus5 N1 选项 b）
# 同方差 sigma=1，WT 删失率不等（90%/78%），delta_M=0.2 < delta_W=0.5
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(2024)
c = np.log(10)
sigma = 1.0  # 同方差（与 Prop 1(iii)/2 一致）

mu_AM = c + st.norm.ppf(0.20)   # 1.461 CRC-mut 删失20%
mu_BM = mu_AM - 0.2             # 1.261 skin-mut (delta_M=0.2)
mu_AW = c + st.norm.ppf(0.90)   # 3.585 CRC-WT 删失90%
mu_BW = mu_AW - 0.5             # 3.085 skin-WT (delta_W=0.5, 删失78%)

beta_true = (mu_AM - mu_AW) - (mu_BM - mu_BW)
print(f"真交互 beta_int = delta_M - delta_W = 0.2 - 0.5 = {beta_true:.2f}")

# 检查删失率
for name, mu in [("AM",mu_AM),("BM",mu_BM),("AW",mu_AW),("BW",mu_BW)]:
    pi = 1 - st.norm.cdf(c - mu)
    print(f"  {name}: mu={mu:.3f} 删失率={pi:.1%}")

# 解析伪极限（封顶 E[y^c]）
def E_capped(mu):
    z = (c - mu)/sigma
    return c - sigma*(st.norm.pdf(z) + z*st.norm.cdf(z))
bt = (E_capped(mu_AM) - E_capped(mu_AW)) - (E_capped(mu_BM) - E_capped(mu_BW))
print(f"封顶 OLS 解析伪极限 beta_tilde = {bt:.3f} (反号: {np.sign(bt)!=np.sign(beta_true)})")

n = (40, 40, 40, 40)
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
for name in ['true(无删失)', 'capped(封顶)', 'naive(无偏外推)']:
    bs = np.zeros(N)
    for i in range(N):
        y, mut, crc, cens = gen()
        yo = y if name=='true(无删失)' else (np.minimum(y,c) if name=='capped(封顶)' else naive_obs(y,cens))
        X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
        m = sm.OLS(yo, X).fit()
        bs[i] = m.params[3]
    flip = (np.sign(bs) != np.sign(beta_true)).mean()
    print(f"{name:14s}: beta_hat 均值={np.mean(bs):+.3f}, 反号比例={flip:.1%}, 95%CI=[{np.quantile(bs,0.025):+.3f},{np.quantile(bs,0.975):+.3f}]")
