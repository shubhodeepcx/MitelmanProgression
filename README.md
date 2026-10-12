# MitelmanProgression

**Author:** Shubhodeep Chatterjee

Per-disease models of the order in which chromosomal aberrations accumulate in
haematological malignancies, inferred from the karyotypes in the Mitelman Database
of Chromosome Aberrations and Gene Fusions in Cancer. The stemline/sideline structure
of each karyotype gives a within-patient partial order; cross-sectional event
patterns across thousands of cases give the rest.

Three models are fitted per disease and compared on the same data:

1. a co-occurrence / mutual-exclusivity baseline (Fisher's exact test and a
   permutation null that preserves each case's number of events, with
   Benjamini-Hochberg correction);
2. a conjunctive Bayesian network (CBN) with observation noise;
3. a Mutual Hazard Network (MHN), the headline model.

Every ordering claim is checked by bootstrap stability, held-out likelihood,
simulation recovery, agreement with established cytogenetic biology, agreement
with the within-patient clone orders, and refits on unselected series and on
recent publications (publication-bias sensitivity). The report is in `report/`.

## Run

    pip install -r requirements.txt     # python-docx only for the report
    bash run_all.sh

`run_all.sh` downloads the Mitelman bulk files if they are missing, rebuilds every
table, model, validation table and figure, and runs the tests. The event vocabulary
is frozen in `config/vocabulary_frozen.json`; set `REFREEZE=1` to rebuild it on
purpose. Each step can also be run on its own:

| Step | Script | Output |
|---|---|---|
| 1 | `scripts/01_build_clone_features.py` | per-clone features and case table (`data/interim/`) |
| 2 | `scripts/02_power_threshold.py` | minimum event count N per disease (`config/power_threshold.json`) |
| 3 | `scripts/03_freeze_vocabulary.py` | frozen vocabulary (only with `REFREEZE=1`) |
| 4 | `scripts/04_build_tables.py` | tidy event table, ordered pairs, case x event matrices |
| 5 | `scripts/05_descriptive.py` | selection-bias, era and data-quality tables |
| 6 | `scripts/06_cooccurrence.py` | co-occurrence baseline for every disease |
| 7 | `scripts/07_cbn.py` | CBN fits, bootstrap, held-out likelihood, simulation recovery |
| 8 | `scripts/08_mhn.py` | MHN fits, lambda by cross-validation, bootstrap, sensitivity, simulation recovery |
| 9 | `scripts/09_validation.py` | prior-knowledge, clone-order and stop-rule verdicts |
| 10 | `scripts/10_figures.py` | one figure per disease plus an overview |

The report is rebuilt from the results with `python report/build_report.py` (writes
`report/main.tex` and the DOCX); compile `main.tex` with `pdflatex` twice for the PDF.
`scripts/write_vocabulary_doc.py` regenerates `docs/event_vocabulary.md`;
`scripts/make_annotation_scaffold.py` drew the 50-karyotype annotation sample.

## Layout

- `karyoprog/io.py`: Mitelman loaders and disease definitions (morphology codes)
- `karyoprog/iscn.py`: clone tokenising, `idem` expansion, clone tree, ordered pairs
- `karyoprog/cn.py`: arm-level copy-number calls from CytoConverter segments
- `karyoprog/features.py`: atomic per-clone features (copy number, structural, ploidy)
- `karyoprog/dedup.py`: likely duplicate patients across papers
- `karyoprog/pipeline.py`: clone features, case table, tidy event table
- `karyoprog/cooccur.py`: co-occurrence baseline and count-preserving permutation null
- `karyoprog/cbn.py`: conjunctive Bayesian network with noise (EM, greedy BIC search)
- `karyoprog/mhn_ssr.py`, `karyoprog/mhnfit.py`: MHN likelihood with state-space restriction, fitting, sampling, order probabilities
- `karyoprog/common.py`: shared paths and loaders
- `config/`: candidate and frozen vocabularies, power threshold, prior-knowledge test set
- `results/tables/`, `results/models/`, `results/figures/`: outputs (see `docs/data_dictionary.md`)
- `tests/`: parser unit tests and the 50-karyotype hand-annotation harness
- `docs/`: frozen vocabulary and data dictionary
- `report/`: the IEEE-format report (LaTeX source, PDF, DOCX)
- `data/raw/`: Mitelman download (not tracked)

## Data

Mitelman F, Johansson B, Mertens F (Eds.). Mitelman Database of Chromosome Aberrations
and Gene Fusions in Cancer. https://mitelmandatabase.isb-cgc.org. The bulk files are
downloaded at run time and are not redistributed here; the derived tables in
`results/tables/` are released with this code.

## Scope

This is a research analysis of published karyotypes. It makes no clinical claims, and
the event orders it reports describe the published case mix, not individual patients.

## Citation and licence

See `CITATION.cff`. Code is released under the MIT licence (`LICENSE`),
copyright 2026 Shubhodeep Chatterjee.
