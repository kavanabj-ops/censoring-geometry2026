# -*- coding: utf-8 -*-
# 72_rank_failure.py - 秩检验失效 DGP：交互由 WT 格驱动（delta_M=0, delta_W!=0）时 mut秩漏检
# 证明 R3（mut-only 秩）的适用边界：只测 delta_M，交互由 delta_W 驱动时失效
import numpy as np, scipy.stats as st, statsmodels.api as sm

rng = np.random.default_rng(11)

c = np.log(10)
sigma = 1.0

# DGP：delta_M = 0（mut 格无 lineage 差异），delta_W = 1（WT 格有差异）
# 真交互 beta_int = delta_M - delta_W = -1（由 WT 格驱动）
mu_AM = c + sigma*st.norm.ppf(0.20)   # 1.461 (CRC-mut)
mu_BM = mu_AM                          # delta_M = 0 (skin-mut = CRC-mut)
mu_AW = c + sigma*st.norm.ppf(0.90)   # 3.585 (CRC-WT, 删失90%)
mu_BW = mu_AW - 1.0                    # 2.585 (skin-WT, delta_W = 1, 删失39%)

beta_true = (mu_AM - mu_AW) - (mu_BM - mu_BW)
print(f"真交互 beta_int = delta_M - delta_W = 0 - 1 = {beta_true:.2f}")
print(f"mut秩检验的是 delta_M = {mu_AM-mu_BM:.2f}（=0，漏检交互）")

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

def mut_rank_p(y, mut, crc):
    a = y[(mut==1)&(crc==1)]; b = y[(mut==1)&(crc==0)]
    if len(a)==0 or len(b)==0:
        return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided')
    return p

N = 4000
# mut秩 功效
ps_rank = np.array([mut_rank_p(*gen()[:3]) for _ in range(N)])
pwr_rank = (ps_rank < 0.05).mean()

# 完整交互 OLS（封顶，因为 WT 格有删失）
ps_ols = np.zeros(N)
for i in range(N):
    y, mut, crc, cens = gen()
    yo = np.minimum(y, c)  # 封顶
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(yo, X).fit()
    ps_ols[i] = m.pvalues[3]
pwr_ols = (ps_ols < 0.05).mean()

# 完整交互 OLS（无删失，理论上限）
ps_true = np.zeros(N)
for i in range(N):
    y, mut, crc, cens = gen()
    X = sm.add_constant(np.column_stack([mut, crc, mut*crc]))
    m = sm.OLS(y, X).fit()
    ps_true[i] = m.pvalues[3]
pwr_true = (ps_true < 0.05).mean()

print(f"\n=== 秩检验失效 DGP（交互由 WT 格驱动）===")
print(f"mut-only 秩:   功效={pwr_rank:.3f}  (漏检，≈α=0.05)")
print(f"封顶 OLS 交互: 功效={pwr_ols:.3f}")
print(f"无删失 OLS 交互: 功效={pwr_true:.3f}  (理论上限)")
print(f"\n结论：mut秩只测 delta_M，交互由 delta_W 驱动时功效≈α，R3 适用边界=交互须由 mut 格驱动")
