# -*- coding: utf-8 -*-
# 83_figures_supplement.py - 补充关键计算图 Fig4/5/6
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG = r"D:\paper\02_analysis\figures"

# --- Figure 4: sigma-hat 双向偏倚 vs n（MC SE）---
fig, ax = plt.subplots(figsize=(7, 4.5))
n_labels = ["n=(6,40,38,17)", "n=(20×4)", "n=(40×4)"]
naive = [1.302, 1.283, 1.285]
naive_se = [0.002, 0.002, 0.001]
capped = [0.547, 0.587, 0.590]
capped_se = [0.001, 0.001, 0.001]
x = np.arange(3)
w = 0.35
ax.bar(x - w/2, naive, w, yerr=naive_se, capsize=4, color="#d62728", label="naive (inflation)")
ax.bar(x + w/2, capped, w, yerr=capped_se, capsize=4, color="#1f77b4", label="capped (deflation)")
ax.axhline(1.0, color="gray", ls="--", lw=1)
ax.set_xticks(x); ax.set_xticklabels(n_labels)
ax.set_ylabel(r"$\hat{\sigma}/\sigma_{\mathrm{true}}$")
ax.set_title("Bidirectional σ̂ bias is essentially n-invariant (WT censoring 94.5%)")
ax.legend()
ax.set_ylim(0.4, 1.5)
plt.tight_layout()
plt.savefig(FIG + r"\fig4_sigma_bias_bars.png", dpi=150)
plt.close()

# --- Figure 5: 功效 vs δ_W ---
fig, ax = plt.subplots(figsize=(7, 4.5))
dW = [0.0, 0.5, 1.0, 1.5]
naive_p = [0.431, 0.130, 0.010, 0.000]
capped_p = [0.662, 0.653, 0.592, 0.404]
rank_p = [0.736, 0.697, 0.744, 0.713]
ax.plot(dW, naive_p, "o-", color="#d62728", label="naive OLS")
ax.plot(dW, capped_p, "s-", color="#1f77b4", label="capped OLS")
ax.plot(dW, rank_p, "^-", color="#2ca02c", label="mutant-only rank")
ax.set_xlabel(r"$\delta_W$ (true WT-cell contrast, $\delta_M=1.8$ fixed)")
ax.set_ylabel("power")
ax.set_title("Estimators target different estimands across δ_W")
ax.legend()
ax.set_ylim(0, 0.9)
plt.tight_layout()
plt.savefig(FIG + r"\fig5_power_deltaW.png", dpi=150)
plt.close()

# --- Figure 6: sharp interval vs B ---
fig, ax = plt.subplots(figsize=(7, 3.8))
B_labels = ["B=ln3", "B=ln10", "B=ln30"]
intervals = [[0.911, 3.339], [-0.474, 4.615], [-1.737, 5.779]]
beta = 2.174
for i, (lo, hi) in enumerate(intervals):
    ax.plot([lo, hi], [i, i], lw=6, color="#9467bd", alpha=0.7)
    ax.plot(beta, i, "o", color="black", ms=6)
ax.axvline(0, color="red", ls="--", lw=1.2)
ax.set_yticks(range(3)); ax.set_yticklabels(B_labels)
ax.set_xlabel(r"$\beta_{\mathrm{int}}$ (interaction)")
ax.set_title("Sharp identified interval for β_int (dabrafenib, BRAF V600E) — contains 0 only for B ≥ ln6.6")
ax.text(beta, 2.5, r"$\tilde{\beta}=2.174$", ha="center")
plt.tight_layout()
plt.savefig(FIG + r"\fig6_sharp_interval.png", dpi=150)
plt.close()

print("Fig4/5/6 已生成")
