#!/usr/bin/env bash
# Single entry point for the week-1 analyses. Re-downloads Mitelman if missing.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RAW="$HERE/data/raw"
if [ ! -f "$RAW/mitelman/CYTOGEN.TXT.DATA" ]; then
  mkdir -p "$RAW/mitelman"
  curl -sS -o "$RAW/mitelman_db.zip" https://storage.googleapis.com/mitelman-data-files/prod/mitelman_db.zip
  python3 -m zipfile -e "$RAW/mitelman_db.zip" "$RAW/mitelman"
fi
python3 "$HERE/scripts/week1_multiclone.py"
python3 "$HERE/scripts/week1_longitudinal.py"
python3 "$HERE/scripts/week1_cytoconverter_check.py"
python3 "$HERE/scripts/week1_aml_frequencies.py"
python3 "$HERE/scripts/make_annotation_scaffold.py"
python3 -m pytest -q "$HERE/tests"
