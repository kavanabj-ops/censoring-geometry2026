# -*- coding: utf-8 -*-
# 71b_sign_reversal.py - 符号反转（异方差 DGP）：让 WT 两格都高删失 且 delta_W > delta_M
# 删失率 pi_l = 1 - Phi((c-mu_l)/sigma_l)，用 sigma 异方差解耦 mu 差异与删失率
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(7)

c = np.log(10)  # 2.303

# 目标：mut 两格删失 20%（Phi_M=0.8），WT 两格删失 93%（Phi_W=0.07）
# mu = c + Phi^-1(pi)*sigma
def mu_for(pi, sigma):
    return c + st.norm.ppf(pi) * sigma

# delta_M = 0.5, delta_W = 3（都正，delta_W 更大 -> 真交互为负）
sig_AM = 1.5;  mu_AM = mu_for(0.20, sig_AM)             # 1.041
sig_BM = 2.09; mu_BM = mu_AM - 0.5                      # 0.541 (delta_M=0.5)
sig_AW = 2.5;  mu_AW = mu_for(0.93, sig_AW)             # 5.993
sig_BW = 0.47; mu_BW = mu_AW - 3.0                      # 2.993 (delta_W=3)

beta_true = (mu_AM - mu_AW) - (mu_BM - mu_BW)
print(f"真交互 beta_int = delta_M - delta_W = 0.5 - 3.0 = {beta_true:.2f}")

# 检查删失率
for name, mu, s in [("AM",mu_AM,sig_AM),("BM",mu_BM,sig_BM),("AW",mu_AW,sig_AW),("BW",mu_BW,sig_BW)]:
    pi = 1 - st.norm.cdf((c-mu)/s)
    print(f"  {name}: mu={mu:.3f} sigma={s:.2f} 删失率={pi:.1%}")

# 解析伪极限（封顶 E[y^c]，各格用各自 sigma）
def E_capped(mu, s):
    z = (c - mu)/s
    return c - s*(st.norm.pdf(z) + z*st.norm.cdf(z))
bt = (E_capped(mu_AM,sig_AM) - E_capped(mu_AW,sig_AW)) - (E_capped(mu_BM,sig_BM) - E_capped(mu_BW,sig_BW))
print(f"封顶 OLS 解析伪极限 beta_tilde = {bt:.3f}  (反号: {np.sign(bt)!=np.sign(beta_true)})")

# 模拟验证
n = (40, 40, 40, 40)
def gen():
    y = np.concatenate([
        rng.normal(mu_AM, sig_AM, n[0]),
        rng.normal(mu_AW, sig_AW, n[1]),
        rng.normal(mu_BM, sig_BM, n[2]),
        rng.normal(mu_BW, sig_BW, n[3]),
    ])
    mut = np.array([1]*n[0] + [0]*n[1] + [1]*n[2] + [0]*n[3])
    crc = np.array([1]*n[0] + [1]*n[1] + [0]*n[2] + [0]*n[3])
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o

N = 4000
for name in ['true(无删失)', 'capped(封顶)', 'naive(外推)']:
    bs = np.zeros(N)
    for i in range(N):
        y, mut, crc, cens = gen()
        yo = y if name=='true(无删失)' else (np.minimum(y,c) if name=='capped(封顶)' else naive_obs(y,cens))
        X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
        m = sm.OLS(yo, X).fit()
        bs[i] = m.params[3]
    flip = (np.sign(bs) != np.sign(beta_true)).mean()
    print(f"{name:14s}: beta_hat 均值={np.mean(bs):+.3f}, 反号比例={flip:.1%}, 95%CI=[{np.quantile(bs,0.025):+.2f},{np.quantile(bs,0.975):+.2f}]")
