"""Core R0-R5 checklist logic."""
import math

__version__ = "1.0.0"


def check_interaction(n_mut_A, n_wt_A, n_mut_B, n_wt_B,
                      wt_censoring_rate, min_n=5, c_star=0.90,
                      delta_star=1.0, sigma=1.0, kappa=2.80):
    """Run the R0-R5 identifiability checklist for a 2x2 mutation x lineage contrast.

    Parameters
    ----------
    n_mut_A, n_wt_A, n_mut_B, n_wt_B : int
        Cell counts in the four strata (two lineages A/B x mutant/WT).
    wt_censoring_rate : float
        Wild-type censoring rate (fraction of WT IC50s at the assay ceiling).
    min_n : int
        Minimum per-cell count (R0, default 5).
    c_star : float
        WT-censoring threshold for non-informativeness (R0, default 0.90).
    delta_star : float
        Target resolution in log units (Corollary 1, default 1.0).
    sigma : float
        Residual SD (Corollary 1, default 1.0).
    kappa : float
        Normal-power constant (default 2.80).

    Returns
    -------
    dict with keys: input, wt_censoring_rate, rules, recommendation.
    """
    n = {"mut_A": n_mut_A, "wt_A": n_wt_A, "mut_B": n_mut_B, "wt_B": n_wt_B}

    report = {
        "input": n,
        "wt_censoring_rate": wt_censoring_rate,
        "rules": {},
        "recommendation": None,
    }

    # R0 step 1: per-cell sample size
    insuff = [k for k, v in n.items() if v < min_n]
    if insuff:
        report["recommendation"] = (
            f"R0 (sample size) FAIL: cells {insuff} have n < {min_n}. "
            "The interaction is not estimable; do not run the 2x2 test."
        )
        report["rules"]["R0"] = "fail - insufficient sample size"
        return report

    # R0 step 2: WT censoring rate
    if wt_censoring_rate > c_star:
        m_wt = min(n_wt_A, n_wt_B) * (1 - wt_censoring_rate)
        report["recommendation"] = (
            f"R0 (censoring) FAIL: WT censoring {wt_censoring_rate:.0%} > {c_star:.0%}. "
            f"The interaction is non-informative (effective WT cell size ~{m_wt:.1f}). "
            "Report the conditional-lineage (mutant-only) rank contrast instead."
        )
        report["rules"]["R0"] = "fail - WT censoring > threshold"
        report["rules"]["R3"] = "use mutant-only rank contrast or E6"
        return report

    report["recommendation"] = (
        "R0 PASS: sample size and censoring are acceptable. "
        "Proceed with the interaction test, following R1-R5 below."
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

    # Corollary 1 admission heuristic
    m = {k: v * (1 - wt_censoring_rate if "wt" in k else 0.80)
         for k, v in n.items()}
    sum_inv_m = sum(1.0 / max(v, 1e-6) for v in m.values())
    practical = sum_inv_m <= (delta_star / (kappa * sigma)) ** 2
    report["rules"]["Corollary1"] = (
        f"admission heuristic: sum(1/m_l) = {sum_inv_m:.4f} "
        f"{'<=' if practical else '>'} (delta*/kappa*sigma)^2 "
        f"= {(delta_star/(kappa*sigma))**2:.4f} -> "
        f"{'practical' if practical else 'not practical'} for interaction inference"
    )

    return report
