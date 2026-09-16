# When the wild type hits the ceiling

Companion code for the manuscript:
*When the wild type hits the ceiling: identifiability of mutation–lineage interactions in pharmacogenomic panels*

## The R0–R5 checklist tool (`censoring-check`)

An installable Python tool that encodes the paper's two-step gate (R0) and rules R1–R5 for deciding, before analysis, whether a mutation × lineage interaction is estimable under censoring.

### Install

```
pip install git+https://github.com/kavanabj-ops/censoring-geometry
```

Or use it without installing (Python 3.8+, no dependencies beyond the standard library for the core check):

```
from censoring_check import check_interaction
```

### Command line

```
censoring-check --mut-A 6 --wt-A 40 --mut-B 38 --wt-B 17 --wt-censoring 0.945
```

### Example

```python
from censoring_check import check_interaction

# BRAF V600E / dabrafenib: WT censoring 94.5% -> non-informative
r = check_interaction(6, 40, 38, 17, wt_censoring_rate=0.945)
print(r["recommendation"])
# R0 (censoring) FAIL: ... report the conditional-lineage rank contrast instead.
```

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

Python 3.8+ with numpy, scipy, pandas, statsmodels.
