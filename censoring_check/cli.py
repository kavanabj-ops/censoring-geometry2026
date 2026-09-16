"""Command-line interface for the R0-R5 checklist."""
import argparse
from .core import check_interaction, __version__


def main():
    p = argparse.ArgumentParser(
        prog="censoring-check",
        description="R0-R5 identifiability checklist for mutation x lineage interactions under censoring.",
    )
    p.add_argument("--mut-A", type=int, required=True, help="mutant cell count, lineage A")
    p.add_argument("--wt-A", type=int, required=True, help="WT cell count, lineage A")
    p.add_argument("--mut-B", type=int, required=True, help="mutant cell count, lineage B")
    p.add_argument("--wt-B", type=int, required=True, help="WT cell count, lineage B")
    p.add_argument("--wt-censoring", type=float, required=True, help="WT censoring rate (0-1)")
    p.add_argument("--min-n", type=int, default=5)
    p.add_argument("--c-star", type=float, default=0.90)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args()

    r = check_interaction(args.mut_A, args.wt_A, args.mut_B, args.wt_B,
                          args.wt_censoring, min_n=args.min_n, c_star=args.c_star)
    print(r["recommendation"])
    for k, v in r["rules"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
