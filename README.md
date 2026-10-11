# karyotype-progression

Inferring the temporal order of chromosomal aberrations in haematological malignancies from
the Mitelman Database, using stemline/sideline structure as the ordering signal.

    pip install -r requirements.txt
    ./run_week1.sh          # downloads Mitelman if needed, writes results/week1/, runs tests

Layout
- `karyoprog/io.py`: Mitelman loaders and disease definitions (MORPH codes)
- `karyoprog/iscn.py`: clone tokenising, `idem` expansion, clone tree, ordered pairs
- `karyoprog/events.py`: provisional event vocabulary (not frozen)
- `karyoprog/build.py`: builds investigation, event and pair tables
- `scripts/`: week-1 analyses
- `tests/`: parser unit tests and the 50-karyotype hand-annotation harness
- `docs/week1_report.md`: week-1 findings and go/no-go
- `data/raw/`: Mitelman download (git-ignored)
