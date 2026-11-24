# -*- coding: utf-8 -*-
"""Run the R0 gate on the bundled demo dataset (examples/demo_data_gdsc.csv).

Input format: CSV with columns cell_id, genotype, LN_IC50, MAX_CONC
  - LN_IC50: natural-log IC50
  - MAX_CONC: natural-log maximum tested concentration (right-censoring ceiling)

The script computes the wild-type censoring rate (fraction of WT cells with
LN_IC50 > MAX_CONC) and feeds the resulting 2x2 cell counts plus the censoring
rate into check_interaction.

Usage:
    python examples/run_gdsc_demo.py
"""
import csv, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from censoring_check import check_interaction

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_data_gdsc.csv")


def main():
    rows = list(csv.DictReader(open(DATA, encoding="utf-8")))
    mut = [r for r in rows if r["genotype"] == "mutant"]
    wt = [r for r in rows if r["genotype"] == "wildtype"]

    def censored(cells):
        return [r for r in cells if float(r["LN_IC50"]) > float(r["MAX_CONC"])]

    n_mut, n_wt = len(mut), len(wt)
    pi_wt = len(censored(wt)) / n_wt

    print(f"Demo dataset: {len(rows)} cells ({n_mut} mutant, {n_wt} wildtype)")
    print(f"Wild-type censoring rate: {len(censored(wt))}/{n_wt} = {pi_wt:.1%}")
    print()

    # Feed the counts into the R0-R5 gate (lineages A and B use the same demo cells)
    r = check_interaction(n_mut, n_wt, n_mut, n_wt, wt_censoring_rate=pi_wt)
    print(r["recommendation"])
    for k, v in r["rules"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
