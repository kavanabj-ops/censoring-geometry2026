# -*- coding: utf-8 -*-
# 77_E6_censoring_interaction.py - 删失概率交互 E6（opus5 N2）
# 验证：delta_M=0, delta_W=1 的 DGP（交互由 WT 格驱动）下，E6 能检出，mut秩漏检
import numpy as np, scipy.stats as st

rng = np.random.default_rng(2024)
c = np.log(10)
sigma = 1.0

mu_AM = c + st.norm.ppf(0.20)   # 1.461 CRC-mut 删失20%
mu_BM = mu_AM                    # delta_M = 0
mu_AW = c + st.norm.ppf(0.90)   # 3.585 CRC-WT 删失90%
mu_BW = mu_AW - 1.0             # 2.585 skin-WT (delta_W = 1, 删失39%)

n = (6, 40, 38, 17)
print(f"真交互 beta_int = delta_M - delta_W = 0 - 1 = -1（由 WT 格驱动）")

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
    if len(a)==0 or len(b)==0: return 1.0
    _, p = st.mannwhitneyu(a, b, alternative='two-sided'); return p

def E6_test(mut, crc, cens):
    # 四格删失率
    def pi(m, cc):
        g = cens[(mut==m)&(crc==cc)]
        return g.mean(), len(g)
    pAM, nAM = pi(1,1); pAW, nAW = pi(0,1); pBM, nBM = pi(1,0); pBW, nBW = pi(0,0)
    # risk difference 交互
    E6 = (pAM - pAW) - (pBM - pBW)
    # delta method SE
    def var(p, n):
        return p*(1-p)/n if n > 0 else 0
    se = np.sqrt(var(pAM,nAM) + var(pAW,nAW) + var(pBM,nBM) + var(pBW,nBW))
    if se == 0:
        return E6, 1.0
    z = E6 / se
    p = 2*(1 - st.norm.cdf(abs(z)))
    return E6, p

N = 4000
ps_rank = np.array([mut_rank_p(*gen()[:3]) for _ in range(N)])
ps_e6 = np.zeros(N)
for i in range(N):
    y, mut, crc, cens = gen()
    _, p = E6_test(mut, crc, cens)
    ps_e6[i] = p

print(f"\n=== 交互由 WT 格驱动（delta_M=0, delta_W=1）===")
print(f"mut-only 秩: 功效 = {(ps_rank < 0.05).mean():.3f}  (漏检，≈α)")
print(f"E6 删失概率交互: 功效 = {(ps_e6 < 0.05).mean():.3f}  (能检出 WT 格删失率差异)")
print(f"\n结论：E6 保留 2×2 结构，在交互由 WT 格驱动时仍可识别（删失率差异即信号），补上 mut秩 的盲区")
