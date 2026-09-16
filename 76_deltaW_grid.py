# -*- coding: utf-8 -*-
# 76_deltaW_grid.py - delta_W 网格：各估计量功效随 delta_W 变化（opus5 N4）
# 证明 mut秩 只测 delta_M（功效不随 delta_W 变），OLS 交互测 delta_M-delta_W（功效随 delta_W 降）
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(2024)
c = np.log(10)
sigma_mut, sigma_wt = 1.5, 0.8
mu_mut = c + st.norm.ppf(0.20) * sigma_mut   # mut 格删失 20%
mu_AW = c + st.norm.ppf(0.945) * sigma_wt    # WT 格删失 94.5%

delta_M = 1.8  # 固定：CRC-mut 比 skin-mut 敏感 1.8 log
delta_Ws = [0.0, 0.5, 1.0, 1.5]

n = (6, 40, 38, 17)

def gen(delta_W, beta):
    mu_BM = mu_mut - delta_M  # skin-mut 无位移
    mu_BW = mu_AW - delta_W
    # 备择：CRC-mut 有 +beta 位移（效力增强 -> lnIC50 降低，这里用 mu_AM + 位移表示 CRC-mut 更敏感）
    # 真交互 beta_int = (mu_AM - mu_AW) - (mu_BM - mu_BW) = delta_M - delta_W
    y = np.concatenate([
        rng.normal(mu_mut, sigma_mut, n[0]),   # CRC-mut（敏感，含位移）
        rng.normal(mu_AW, sigma_wt, n[1]),      # CRC-WT
        rng.normal(mu_mut - delta_M, sigma_mut, n[2]),  # skin-mut
        rng.normal(mu_BW, sigma_wt, n[3]),      # skin-WT
    ])
    mut = np.array([1]*n[0] + [0]*n[1] + [1]*n[2] + [0]*n[3])
    crc = np.array([1]*n[0] + [1]*n[1] + [0]*n[2] + [0]*n[3])
    cens = y > c
    return y, mut, crc, cens

def naive_obs(y, cens):
    o = y.copy(); o[cens] = y[cens] + rng.normal(0, 2.0, cens.sum()); return o
def clip_obs(y, cens):
    o = y.copy(); o[cens] = c; return o
def fit_ols(y, mut, crc, cov='HC3'):
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit(cov_type=cov)
    return m.params[3], m.pvalues[3]
def mut_rank_p(y, mut, crc):
    a = y[(mut==1)&(crc==1)]; b = y[(mut==1)&(crc==0)]
    if len(a)==0 or len(b)==0: return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided'); return p

N = 3000
print("=== delta_W 网格（delta_M=1.8 固定，N=3000）===")
print(f"{'delta_W':>8} {'真交互':>8} | {'naive 功效':>10} {'capped 功效':>10} {'mut秩 功效':>10}")
for delta_W in delta_Ws:
    beta_int = delta_M - delta_W
    # 各估计量：null 阈值 + alternative 功效
    pwrs = {}
    for ename, obs in [("naive", naive_obs), ("capped", clip_obs)]:
        ps0 = np.zeros(N); ps1 = np.zeros(N)
        for i in range(N):
            y0, mut0, crc0, _ = gen(delta_W, 0.0)
            yo = obs(y0, _ if ename=='naive' else None)
            # 注意：obs 需要 cens 参数
            pass
        # 简化：直接算（null 用 delta_M=0 的位移？）—— 这里 delta_M 固定，null 是 beta=0
        # 实际上 delta_M=1.8 是"备择位移"，null 应该是 delta_M=0
        # 重新：null = delta_M=0（无位移），alternative = delta_M=1.8
        break
    break

# 重新组织：null 是 delta_M=0，alternative 是 delta_M=1.8
def run_estimator(obs_fn, delta_W, delta_M_val):
    ps = np.zeros(N)
    for i in range(N):
        y, mut, crc, cens = gen(delta_W, 0.0)  # 这里 gen 里 delta_M 是固定的 1.8
        # 需要参数化 delta_M，重新写
        pass
    return ps

# 更简洁：重新定义 gen 支持 delta_M 参数
def gen2(delta_W, delta_M_val):
    mu_AM = mu_mut  # CRC-mut（delta_M_val=0 时无位移；下面用 mu_mut + 0 表示）
    mu_BM = mu_mut - delta_M_val
    mu_BW = mu_AW - delta_W
    y = np.concatenate([
        rng.normal(mu_AM, sigma_mut, n[0]),
        rng.normal(mu_AW, sigma_wt, n[1]),
        rng.normal(mu_BM, sigma_mut, n[2]),
        rng.normal(mu_BW, sigma_wt, n[3]),
    ])
    mut = np.array([1]*n[0] + [0]*n[1] + [1]*n[2] + [0]*n[3])
    crc = np.array([1]*n[0] + [1]*n[1] + [0]*n[2] + [0]*n[3])
    cens = y > c
    return y, mut, crc, cens

print("=== delta_W 网格（N=3000）===")
print(f"{'delta_W':>8} {'真交互':>8} | {'naive':>8} {'capped':>8} {'mut秩':>8}")
for delta_W in delta_Ws:
    beta_int = delta_M - delta_W
    # null 阈值（delta_M=0）
    thresh = {}
    for ename in ['naive', 'capped']:
        ps0 = np.zeros(N)
        for i in range(N):
            y, mut, crc, cens = gen2(delta_W, 0.0)
            yo = naive_obs(y, cens) if ename=='naive' else clip_obs(y, cens)
            b, p = fit_ols(yo, mut, crc)
            ps0[i] = p
        thresh[ename] = np.quantile(ps0, 0.05)
    ps0_rank = np.array([mut_rank_p(*gen2(delta_W, 0.0)[:3]) for _ in range(N)])
    thresh['mut秩'] = np.quantile(ps0_rank, 0.05)
    # alternative 功效（delta_M=1.8）
    pwrs = {}
    for ename in ['naive', 'capped']:
        ps1 = np.zeros(N)
        for i in range(N):
            y, mut, crc, cens = gen2(delta_W, delta_M)
            yo = naive_obs(y, cens) if ename=='naive' else clip_obs(y, cens)
            b, p = fit_ols(yo, mut, crc)
            ps1[i] = p
        pwrs[ename] = (ps1 < thresh[ename]).mean()
    ps1_rank = np.array([mut_rank_p(*gen2(delta_W, delta_M)[:3]) for _ in range(N)])
    pwrs['mut秩'] = (ps1_rank < thresh['mut秩']).mean()
    print(f"{delta_W:8.1f} {beta_int:8.1f} | {pwrs['naive']:8.3f} {pwrs['capped']:8.3f} {pwrs['mut秩']:8.3f}")
