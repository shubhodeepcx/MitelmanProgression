#!/usr/bin/env bash
# MitelmanProgression, single entry point (Shubhodeep Chatterjee).
# Mitelman download -> clone features -> released tables -> baseline, CBN, MHN -> validation -> figures -> tests.
# The event vocabulary is frozen in config/vocabulary_frozen.json; set REFREEZE=1 to rebuild it deliberately.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
RAW="$HERE/data/raw"
if [ ! -f "$RAW/mitelman/CYTOGEN.TXT.DATA" ]; then
  mkdir -p "$RAW/mitelman"
  curl -sS -o "$RAW/mitelman_db.zip" https://storage.googleapis.com/mitelman-data-files/prod/mitelman_db.zip
  "$PY" -m zipfile -e "$RAW/mitelman_db.zip" "$RAW/mitelman"
fi
"$PY" "$HERE/scripts/01_build_clone_features.py"
"$PY" "$HERE/scripts/02_power_threshold.py"
if [ "${REFREEZE:-0}" = "1" ]; then
  "$PY" "$HERE/scripts/03_freeze_vocabulary.py"
fi
"$PY" "$HERE/scripts/04_build_tables.py"
"$PY" "$HERE/scripts/05_descriptive.py"
"$PY" "$HERE/scripts/06_cooccurrence.py"
"$PY" "$HERE/scripts/07_cbn.py"
"$PY" "$HERE/scripts/08_mhn.py"
"$PY" "$HERE/scripts/09_validation.py"
"$PY" "$HERE/scripts/10_figures.py"
"$PY" -m pytest -q "$HERE/tests"
