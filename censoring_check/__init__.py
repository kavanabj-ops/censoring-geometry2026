"""R0-R5 identifiability checklist for mutation-lineage interaction tests under censoring.

Companion tool for "When the wild type hits the ceiling: identifiability of
mutation-lineage interactions in pharmacogenomic panels".

Usage:
    from censoring_check import check_interaction
    report = check_interaction(n_mut_A, n_wt_A, n_mut_B, n_wt_B, wt_censoring_rate)
"""
from .core import check_interaction, __version__

__all__ = ["check_interaction", "__version__"]
