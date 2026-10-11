"""Assert the parser reproduces the hand annotation for every row marked status=done.

Format of expected_* columns (same as hand_curated_50_parser_output.tsv):
  events_by_clone:  "1:t(8;21); 2:-Y|t(8;21)"   (clone index, events sorted, '|'-joined)
  parent_by_clone:  "1:-; 2:1"                  (normal clones omitted)
  ordered_pairs:    "t(8;21)<-Y"                ('; '-joined, sorted)
"""
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from karyoprog.iscn import parse_investigation, relate, ordered_pairs
from karyoprog.events import events_of

ANN = Path(__file__).parent / "data" / "hand_curated_50.tsv"
rows = pd.read_csv(ANN, sep="\t", dtype=str, keep_default_na=False)
done = rows[rows.status == "done"]


def _set(s):
    return {x.strip() for x in s.split(";") if x.strip()}


@pytest.mark.skipif(done.empty, reason="no hand annotations completed yet")
@pytest.mark.parametrize("r", done.to_dict("records"), ids=lambda r: f"item{r['item']}")
def test_matches_hand_annotation(r):
    c = parse_investigation(r["karyotype"].split("/"))
    par, _ = relate(c)
    ev = {f"{x.index}:{'|'.join(sorted({e for a in x.abn for e in events_of(a)}))}" for x in c}
    assert ev == _set(r["expected_events_by_clone"])
    assert {f"{i}:{p if p is not None else '-'}" for i, p in par.items()} == _set(r["expected_parent_by_clone"])
    assert {f"{a}<{b}" for a, b in ordered_pairs(c, key=events_of)} == _set(r["expected_ordered_pairs"])


def test_annotation_set_has_50_rows():
    assert len(rows) == 50
