# censoring-check

Companion code for the manuscript:
*Censoring geometry limits mutation–lineage interaction testing in pharmacogenomic panels: a censoring-aware admission gate*

## The R0–R5 checklist tool (`censoring-check`)

An installable Python tool that encodes the paper's two-step gate (R0) and rules R1–R5 for deciding, before analysis, whether a mutation × lineage interaction is estimable under censoring.

### Install

```
pip install git+https://github.com/kavanabj-ops/censoring-geometry2026
```

Or use it without installing (Python 3.8+; the core check needs only the standard library):

```
from censoring_check import check_interaction
```

### Command line

```
censoring-check --mut-A 6 --wt-A 40 --mut-B 38 --wt-B 17 --wt-censoring 0.945
```

### 5-minute quickstart

```
python examples/quickstart.py
```

Runs the gate on the paper's three headline examples (BRAF V600E / dabrafenib; BRAF V600E / trametinib; a well-powered reference case).

## Input format

The gate takes five numbers for a 2×2 mutation × lineage contrast:

| Argument | Meaning |
|---|---|
| `--mut-A`, `--wt-A` | mutant and wild-type cell counts, lineage A |
| `--mut-B`, `--wt-B` | mutant and wild-type cell counts, lineage B |
| `--wt-censoring` | fraction of wild-type IC50 values recorded at the assay ceiling (0–1) |

Optional thresholds (defaults follow the paper):

| Argument | Default | Meaning |
|---|---|---|
| `--min-n` | 5 | minimum per-cell sample size (R0 step 1) |
| `--c-star` | 0.90 | WT-censoring rate above which the contrast is non-informative (R0 step 2) |

The Corollary 1 admission heuristic is reported on every pass: it combines per-cell effective sample size m_l = n_l (1 − censoring) into sum(1/m_l) and compares it against (Δ*/κσ)² with Δ* = 1 log unit, κ = 2.80, σ = 1.

## Rules (R0–R5)

- **R0** — two-step gate: per-cell sample size, then WT censoring. Fail either and do not run the 2×2 interaction test.
- **R1** — never trust naive (extrapolated) point estimates; use capped values and report the censoring fraction.
- **R2** — on capped data prefer rank/permutation tests; model-based SE requires a heteroscedasticity-robust correction.
- **R3** — prefer the conditional-lineage (mutant-only) rank contrast or the censoring-probability interaction (E6).
- **R4** — report the censoring matrix (per-cell n + censoring rate).
- **R5** — attach exact power to binary/tail endpoints.

## Example

```python
from censoring_check import check_interaction

# BRAF V600E / dabrafenib: WT censoring 94.5% -> non-informative
r = check_interaction(6, 40, 38, 17, wt_censoring_rate=0.945)
print(r["recommendation"])
# R0 (censoring) FAIL: WT censoring 94% > 90%. The interaction is non-informative
# (effective WT cell size ~0.9). Report the conditional-lineage (mutant-only)
# rank contrast instead.
```

## Tests

```
pip install pytest
pytest tests/ -v
```

CI runs the test suite and a CLI smoke test on Python 3.8/3.10/3.12 (`.github/workflows/tests.yml`).

## Analysis scripts (55–83)

| Script | Purpose |
|---|---|
| 63 | 2D grid (censoring rate × cell count), separates geometry from small-n effects |
| 64–65 | A1b census + n-threshold sensitivity |
| 67–68 | hotspot sub-analysis + GDSC1 external replication |
| 69–70 | positive controls + negative-control FPR calibration |
| 71c–72 | sign-reversal (homoscedastic) + rank-test failure DGP |
| 73 | CONSORT-style three-level bottleneck count |
| 74 | Monte-Carlo SE for σ̂ ratios |
| 75, 83 | figures |
| 76–79 | δ_W grid, E6, sharp interval, Tobit m_l grid |

## Data

- GDSC1/GDSC2 fitted dose-response: Genomics of Drug Sensitivity in Cancer (release 27 Oct 2023), https://www.cancerrxgene.org
- CCLE/DepMap: `CCLE_mutations.csv`, `sample_info.csv`, DepMap release 26Q1, https://depmap.org

## Environment

Python 3.8+; the core check needs only the standard library. Analysis scripts additionally use numpy, scipy, pandas, statsmodels.

## Citing

[Zenodo DOI placeholder — will be archived upon acceptance]
