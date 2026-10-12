# Author: Shubhodeep Chatterjee
"""The parser must reproduce the hand annotation of the 50-karyotype set for every row with status=done.

Annotation format (tests/data/hand_curated_50.tsv):
  expected_features_by_clone  "1:t(8;21); 2:-Y|t(8;21)"  clone index : features from the ISCN text
                              ('|'-joined, any order; clones separated by '; '; normal clones omitted; "1:" = abnormal, no features)
  expected_parent_by_clone    "1:-; 2:1"                  inferred parent clone, '-' for founders
  expected_ordered_pairs      "t(8;21)<-Y; ..."           earlier<later pairs from parent->child edges
"""
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from karyoprog.iscn import parse_investigation, relate, ordered_pairs
from karyoprog.features import aberration_features, iscn_clone_features

ANN = Path(__file__).parent / "data" / "hand_curated_50.tsv"
rows = pd.read_csv(ANN, sep="\t", dtype=str, keep_default_na=False)
done = rows[rows.status == "done"]


def _bymap(s):
    out = {}
    for part in s.split("; "):
        part = part.strip()
        if part:
            k, _, v = part.partition(":")
            out[int(k)] = v
    return out


def _pairs(s):
    return {tuple(x.strip().split("<", 1)) for x in s.split("; ") if x.strip()}


@pytest.mark.parametrize("r", done.to_dict("records"), ids=lambda r: f"item{r['item']}")
def test_matches_hand_annotation(r):
    cl = parse_investigation(r["karyotype"].split("/"))
    got_f = {c.index: set(iscn_clone_features(c)) for c in cl if not c.is_normal}
    exp_f = {k: {x for x in v.split("|") if x} for k, v in _bymap(r["expected_features_by_clone"]).items()}
    assert got_f == exp_f
    par, _ = relate(cl)
    exp_p = {k: (None if v == "-" else int(v)) for k, v in _bymap(r["expected_parent_by_clone"]).items()}
    assert par == exp_p
    assert ordered_pairs(cl, key=aberration_features) == _pairs(r["expected_ordered_pairs"])


def test_annotation_set_complete():
    assert len(rows) == 50 and (rows.status == "done").all()
