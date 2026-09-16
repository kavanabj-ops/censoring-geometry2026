# -*- coding: utf-8 -*-
# 75_figures.py - 3 张图：删失几何示意 + sigma_hat 比 vs n + CONSORT 瀑布图
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FIG = r"D:\paper\02_analysis\figures"

# ---------- 图 1: 删失几何示意（BRAF V600E 四格） ----------
fig, ax = plt.subplots(figsize=(6, 5))
cells = [
    ("CRC-mut", 6, 0.20, 0),
    ("CRC-WT", 40, 0.945, 1),
    ("skin-mut", 38, 0.20, 0),
    ("skin-WT", 17, 0.945, 1),
]
positions = [(0, 0), (1, 0), (0, 1), (1, 1)]  # (col, row)
for (label, n, pi, is_wt), (cx, cy) in zip(cells, positions):
    color = plt.cm.Reds(0.3 + 0.7*pi) if is_wt else plt.cm.Blues(0.3 + 0.7*pi)
    ax.add_patch(plt.Rectangle((cx, cy), 0.9, 0.9, facecolor=color, edgecolor='black'))
    txtc = 'white' if pi > 0.5 else 'black'
    ax.text(cx+0.45, cy+0.55, f"n={n}", ha='center', va='center', fontsize=11, fontweight='bold', color=txtc)
    ax.text(cx+0.45, cy+0.32, f"cens={pi:.0%}", ha='center', va='center', fontsize=10, color=txtc)
    ax.text(cx+0.45, cy+0.12, label, ha='center', va='center', fontsize=9, color='white')
ax.set_xlim(-0.2, 2.1); ax.set_ylim(-0.2, 2.1)
ax.set_aspect('equal'); ax.axis('off')
ax.set_title("Censoring geometry (dabrafenib, BRAF V600E)\nWT cells near-totally censored (94.5%), mutant cells informative (20%)")
plt.tight_layout(); plt.savefig(FIG + r"\fig1_censoring_geometry.png", dpi=150); plt.close()

# ---------- 图 2: sigma_hat 比 vs n（带 SE） ----------
fig, ax = plt.subplots(figsize=(6, 4))
labels = ["n=(6,40,38,17)", "n=(20x4)", "n=(40x4)"]
x = np.arange(3)
naive_mean = [1.302, 1.283, 1.285]; naive_se = [0.002, 0.002, 0.001]
capped_mean = [0.547, 0.587, 0.590]; capped_se = [0.001, 0.001, 0.001]
ax.errorbar(x, naive_mean, yerr=naive_se, fmt='o-', capsize=4, label='naive (inflation)', color='C0')
ax.errorbar(x, capped_mean, yerr=capped_se, fmt='s-', capsize=4, label='capped (deflation)', color='C1')
ax.axhline(1.0, ls='--', color='gray', lw=1, label='unbiased ($\\hat\\sigma/\\sigma=1$)')
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel('$\\hat\\sigma/\\sigma_{true}$')
ax.set_ylim(0.4, 1.45)
ax.legend(fontsize=9)
ax.set_title('$\\hat\\sigma$ bias is essentially invariant to n\n(naive 1.28-1.30; capped slightly more deflated at n=6)')
plt.tight_layout(); plt.savefig(FIG + r"\fig2_sigma_bias.png", dpi=150); plt.close()

# ---------- 图 3: CONSORT 瀑布图 ----------
fig, ax = plt.subplots(figsize=(6, 5))
steps = [
    ("Total drug x gene x lineage", 12737, 'C0'),
    ("No hotspot mutation", 12039, 'C3'),
    ("Candidates (n>=3)", 200, 'C0'),
    ("Insufficient (n<5)", 90, 'C3'),
    ("Sufficient (n>=5)", 110, 'C0'),
    ("Non-identifiable (WT cens>0.9)", 40, 'C1'),
]
y = np.arange(len(steps))[::-1]
vals = [s[1] for s in steps]
colors = [s[2] for s in steps]
ax.barh(y, vals, color=colors)
ax.set_yticks(y)
ax.set_yticklabels([s[0] for s in steps], fontsize=9)
ax.set_xscale('log')
ax.set_xlabel('count (log scale)')
ax.set_title('Three-level bottleneck (targeted agents)')
plt.tight_layout(); plt.savefig(FIG + r"\fig3_consort.png", dpi=150); plt.close()

print("3 张图已保存到", FIG)
import os
for f in sorted(os.listdir(FIG)):
    if f.startswith('fig'):
        print("  ", f, os.path.getsize(os.path.join(FIG, f)), "bytes")
