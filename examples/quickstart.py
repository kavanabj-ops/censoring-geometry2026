# -*- coding: utf-8 -*-
"""5-minute quickstart: run the R0-R5 gate on the paper's three headline examples.

Requires only the standard library (censoring_check.core uses math only).

Usage:
    python examples/quickstart.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from censoring_check import check_interaction


EXAMPLES = [
    # (label, n_mut_A, n_wt_A, n_mut_B, n_wt_B, wt_censoring_rate)
    ("BRAF V600E / dabrafenib (skin vs CRC) - WT 94.5% censored",
     6, 40, 38, 17, 0.945),
    ("BRAF V600E / trametinib (skin) - low censoring",
     40, 15, 30, 15, 0.13),
    ("Fictional well-powered panel - should pass",
     500, 500, 500, 500, 0.30),
]


def main():
    for label, n_mA, n_wA, n_mB, n_wB, pi_wt in EXAMPLES:
        print("=" * 70)
        print(label)
        r = check_interaction(n_mA, n_wA, n_mB, n_wB,
                              wt_censoring_rate=pi_wt)
        print(r["recommendation"])
        for k, v in r["rules"].items():
            print(f"  {k}: {v}")
    print("=" * 70)
    print("Done. See censoring_check/core.py for the full rule set (R0-R5).")


if __name__ == "__main__":
    main()
