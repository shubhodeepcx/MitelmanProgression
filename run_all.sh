#!/usr/bin/env bash
# Single entry point: Mitelman download -> frozen vocabulary -> released tables -> tests.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RAW="$HERE/data/raw"
if [ ! -f "$RAW/mitelman/CYTOGEN.TXT.DATA" ]; then
  mkdir -p "$RAW/mitelman"
  curl -sS -o "$RAW/mitelman_db.zip" https://storage.googleapis.com/mitelman-data-files/prod/mitelman_db.zip
  python3 -m zipfile -e "$RAW/mitelman_db.zip" "$RAW/mitelman"
fi
python3 "$HERE/scripts/build_clone_features.py"
python3 "$HERE/scripts/power_threshold.py"
python3 "$HERE/scripts/freeze_vocabulary.py"   # re-freezing changes config/vocabulary_frozen.json; commit deliberately
python3 "$HERE/scripts/build_tables.py"
python3 -m pytest -q "$HERE/tests"
