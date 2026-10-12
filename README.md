# karyotype-progression

Inferring the temporal order of chromosomal aberrations in haematological malignancies from
the Mitelman Database, using stemline/sideline structure as the ordering signal.

    pip install -r requirements.txt
    ./run_all.sh            # downloads Mitelman if needed, builds results/tables/, runs tests
    ./run_week1.sh          # the week-1 exploratory analyses (results/week1/)

Layout
- `karyoprog/io.py`: Mitelman loaders and disease definitions (MORPH codes)
- `karyoprog/iscn.py`: clone tokenising, `idem` expansion, clone tree, ordered pairs
- `karyoprog/cn.py`: arm-level copy-number calls from CytoConverter segments
- `karyoprog/features.py`: atomic per-clone features (copy number, structural, ploidy)
- `karyoprog/dedup.py`: likely duplicate patients across papers
- `karyoprog/pipeline.py`: clone features, case table, tidy event table
- `config/`: candidate and frozen event vocabularies, power-based threshold N
- `results/tables/`: released tables (see docs/data_dictionary.md)
- `karyoprog/events.py`, `karyoprog/build.py`: week-1 provisional versions, kept so results/week1 reproduces
- `scripts/`: week-1 analyses
- `tests/`: parser unit tests and the 50-karyotype hand-annotation harness
- `docs/`: week-1 and weeks 2-3 reports, frozen vocabulary, data dictionary
- `data/raw/`: Mitelman download (git-ignored)
