"""Tests for the R0-R5 checklist (censoring_check.core)."""
import pytest

from censoring_check import check_interaction


def test_r0_sample_size_fail():
    # A1 example with one cell below min_n: fail before censoring check
    r = check_interaction(6, 40, 3, 17, wt_censoring_rate=0.945)
    assert r["rules"]["R0"].startswith("fail - insufficient sample size")
    assert "R0 (sample size) FAIL" in r["recommendation"]


def test_r0_censoring_fail():
    # BRAF V600E / dabrafenib geometry (WT censoring 94.5%) -> non-informative
    r = check_interaction(6, 40, 38, 17, wt_censoring_rate=0.945)
    assert r["rules"]["R0"].startswith("fail - WT censoring")
    assert "rank contrast" in r["recommendation"]
    assert "R3" in r["rules"]


def test_r0_pass_gives_rules_r1_r5():
    r = check_interaction(30, 40, 50, 60, wt_censoring_rate=0.30)
    assert r["rules"]["R0"] == "pass"
    for k in ["R1", "R2", "R3", "R4", "R5", "Corollary1"]:
        assert k in r["rules"]
    assert "R0 PASS" in r["recommendation"]


def test_corollary1_admission_flips_with_cell_size():
    # Large cells: practical; tiny cells: not practical
    big = check_interaction(500, 500, 500, 500, wt_censoring_rate=0.30)
    small = check_interaction(5, 5, 5, 5, wt_censoring_rate=0.30)
    assert "practical" in big["rules"]["Corollary1"]
    assert "not practical" in small["rules"]["Corollary1"]


def test_input_echoed():
    r = check_interaction(6, 40, 38, 17, wt_censoring_rate=0.945)
    assert r["input"] == {"mut_A": 6, "wt_A": 40, "mut_B": 38, "wt_B": 17}
    assert r["wt_censoring_rate"] == 0.945


def test_boundary_c_star_not_exceeded():
    # Exactly at threshold: censoring check uses strict >, so 0.90 passes
    r = check_interaction(20, 20, 20, 20, wt_censoring_rate=0.90)
    assert r["rules"]["R0"] == "pass"


def test_threshold_above_c_star_fails():
    r = check_interaction(20, 20, 20, 20, wt_censoring_rate=0.91)
    assert r["rules"]["R0"].startswith("fail - WT censoring")
