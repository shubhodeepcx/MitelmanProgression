#!/usr/bin/env bash
# MitelmanProgression, single entry point (Shubhodeep Chatterjee).
# Mitelman download -> clone features -> released tables -> baseline, CBN, MHN, clone timing
# -> validation -> figures -> tests. Settings: config/analysis.json.
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
cd "$HERE"
if [ "${REFREEZE:-0}" = "1" ]; then
  "$PY" -m karyoprog run --refreeze
else
  "$PY" -m karyoprog run
fi
"$PY" -m pytest -q "$HERE/tests"
