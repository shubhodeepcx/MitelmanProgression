# Changelog

MitelmanProgression, Shubhodeep Chatterjee.

## 1.1.0 (2026-10-13)
- Model-free event timing from within-patient clone orders (Bradley-Terry, case bootstrap), compared
  with the MHN and CBN orders, plus a diagnostic for where the MHN orders disagree.
- 5-fold cross-validation for lambda and for every held-out comparison, with fold-level spread;
  100 MHN and 50 CBN bootstrap refits; 10 simulation-recovery data sets per model; 2,000 permutations.
- All thresholds, counts and seeds in `config/analysis.json`.
- Released case table, so steps 6-11 run from `results/tables/` without the raw download.
- `karyoprog` command-line tool, Makefile, Dockerfile, pinned versions, CI with lint and tests.
- Report: timing results, fold-level uncertainty, supplementary figures for every disease.

## 1.0.0 (2026-10-12)
- Co-occurrence baseline, CBN and MHN per disease with bootstrap, held-out likelihood, simulation
  recovery, prior-knowledge, clone-order and publication-bias checks; figures; IEEE report.
