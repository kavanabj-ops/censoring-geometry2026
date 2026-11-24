# Independent installation & use test

This checklist is for a tester who is **not** the tool author. It takes ~5 minutes
and requires only Python 3.8+ (no other dependencies for the core check).

## 1. Install

```
pip install git+https://github.com/kavanabj-ops/censoring-geometry2026
```

Expected: pip completes without errors and reports `Successfully installed censoring-check-1.0.0`.

## 2. Run the CLI on the paper's headline example

```
censoring-check --mut-A 6 --wt-A 40 --mut-B 38 --wt-B 17 --wt-censoring 0.945
```

Expected output (verbatim):

```
R0 (censoring) FAIL: WT censoring 94% > 90%. The interaction is non-informative (effective WT cell size ~0.9). Report the conditional-lineage (mutant-only) rank contrast instead.
  R0: fail - WT censoring > threshold
  R3: use mutant-only rank contrast or E6
```

## 3. Run the bundled quickstart

Clone or download the repository, then:

```
python examples/quickstart.py
```

Expected: three blocks of output, starting with "BRAF V600E / dabrafenib (skin vs CRC) - WT 94.5% censored" and ending with "Done. See censoring_check/core.py for the full rule set (R0-R5)."

## 4. Run the bundled demo dataset

```
python examples/run_gdsc_demo.py
```

Expected: it reports "Demo dataset: 14 cells (6 mutant, 8 wildtype)", a wild-type
censoring rate of 100.0%, and the same R0 (censoring) FAIL message as step 2.

## 5. Run the test suite

```
pip install pytest
pytest tests/ -v
```

Expected: `7 passed`.

---

Record your result below and send it back:

- [ ] Steps 1-5 all completed successfully
- [ ] Any deviation (describe):

Tester name / date:
