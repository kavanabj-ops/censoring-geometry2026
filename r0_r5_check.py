# -*- coding: utf-8 -*-
"""
r0_r5_check.py — 删失几何可识别性检查（R0–R5 守则的可执行版）

Companion script for:
  "When the wild type hits the ceiling: identifiability of mutation–lineage
   interactions in pharmacogenomic panels"

用法:
    from r0_r5_check import check_interaction
    report = check_interaction(n_mut_A, n_wt_A, n_mut_B, n_wt_B, wt_censoring_rate)

其中 A/B 是两条 lineage，mut/WT 是突变/野生型细胞组，n 是各组细胞系数，
wt_censoring_rate 是野生型组的删失率（IC50 超过最高测试浓度的比例）。

守则（R0–R5）:
  R0 两步门槛：先查逐格样本量，再查 WT 删失率（顺序不可反）
  R1 naive（外推）点估计不可靠，封顶值可用但须报删失率
  R2 封顶数据用秩/置换检验，模型标准误需 sandwich 校正
  R3 优先条件谱系秩对比（mut-only）或删失概率交互（E6）
  R4 强制报告删失矩阵（逐格 n + 删失率）
  R5 二值/尾部终点附精确功效
"""
import math


def check_interaction(n_mut_A, n_wt_A, n_mut_B, n_wt_B,
                      wt_censoring_rate, min_n=5, c_star=0.90,
                      delta_star=1.0, sigma=1.0, kappa=2.80):
    n = {"mut_A": n_mut_A, "wt_A": n_wt_A, "mut_B": n_mut_B, "wt_B": n_wt_B}

    report = {
        "input": n,
        "wt_censoring_rate": wt_censoring_rate,
        "rules": {},
        "recommendation": None,
    }

    # ---- R0 第一步：逐格样本量 ----
    insuff = [k for k, v in n.items() if v < min_n]
    if insuff:
        report["recommendation"] = (
            f"R0 (sample size) FAIL: cells {insuff} have n < {min_n}. "
            "The interaction is not estimable; do not run the 2x2 test."
        )
        report["rules"]["R0"] = "fail — insufficient sample size"
        return report

    # ---- R0 第二步：WT 删失率 ----
    if wt_censoring_rate > c_star:
        m_wt = min(n_wt_A, n_wt_B) * (1 - wt_censoring_rate)
        report["recommendation"] = (
            f"R0 (censoring) FAIL: WT censoring {wt_censoring_rate:.0%} > {c_star:.0%}. "
            f"The interaction is non-informative (effective WT cell size ~{m_wt:.1f}). "
            "Report the conditional-lineage (mutant-only) rank contrast instead."
        )
        report["rules"]["R0"] = "fail — WT censoring > threshold"
        report["rules"]["R3"] = "use mutant-only rank contrast or E6"
        return report

    # ---- 通过 R0：进入 R1–R5 建议 ----
    report["recommendation"] = (
        "R0 PASS: sample size and censoring are acceptable. "
        "Proceed with the interaction test, following R1–R5 below."
    )
    report["rules"]["R0"] = "pass"
    report["rules"]["R1"] = (
        "naive (extrapolated) point estimates are unreliable; "
        "use capped values and report the censoring fraction"
    )
    report["rules"]["R2"] = (
        "on capped data prefer rank/permutation tests; "
        "model-based SE requires a heteroscedasticity-robust (sandwich) correction"
    )
    report["rules"]["R3"] = (
        "prefer the conditional-lineage rank contrast (mutant-only) "
        "or the censoring-probability interaction (E6)"
    )
    report["rules"]["R4"] = "report the censoring matrix (per-cell n + censoring rate)"
    report["rules"]["R5"] = "attach exact power to binary/tail endpoints"

    # ---- 准入启发式（Corollary 1）：m_l = n_l(1 - pi_l) ----
    m = {k: v * (1 - wt_censoring_rate if "wt" in k else 0.80)
         for k, v in n.items()}  # mutant cells ~20% censoring
    sum_inv_m = sum(1.0 / max(v, 1e-6) for v in m.values())
    practical = sum_inv_m <= (delta_star / (kappa * sigma)) ** 2
    report["rules"]["Corollary1"] = (
        f"admission heuristic: sum(1/m_l) = {sum_inv_m:.4f} "
        f"{'<=' if practical else '>'} (delta*/kappa*sigma)^2 "
        f"= {(delta_star/(kappa*sigma))**2:.4f} -> "
        f"{'practical' if practical else 'not practical'} for interaction inference"
    )

    return report


if __name__ == "__main__":
    # 示例 1：BRAF V600E / dabrafenib（WT 删失 94.5%，应触发 R0-censoring）
    print("=== 示例 1: BRAF V600E / dabrafenib ===")
    r = check_interaction(6, 40, 38, 17, wt_censoring_rate=0.945)
    print(r["recommendation"])

    print("\n=== 示例 2: 良好几何（WT 删失 20%）===")
    r2 = check_interaction(40, 40, 40, 40, wt_censoring_rate=0.20)
    print(r2["recommendation"])
    for k, v in r2["rules"].items():
        print(f"  {k}: {v}")

    print("\n=== 示例 3: 样本量不足 ===")
    r3 = check_interaction(3, 40, 4, 30, wt_censoring_rate=0.20)
    print(r3["recommendation"])
