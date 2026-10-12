# Author: Shubhodeep Chatterjee
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from karyoprog.features import clone_features, iscn_cn_features, structural_features, ploidy_features, STRUCTURAL
from karyoprog.pipeline import events_for


def test_structural():
    assert structural_features("der(11)t(9;11)(p22;q23)") == {"KMT2A"}
    assert structural_features("+der(22)t(9;22)(q34;q11)") == {"t(9;22)", "extraPh"}
    assert structural_features("t(3;21)(q26;q22)") == {"3q26"}
    assert structural_features("inv(3)(q21q26)") == {"inv(3)/t(3;3)"}
    assert structural_features("i(17)(q10)") == {"i(17q)"}
    assert structural_features("del(11)(q23)") == set()


def test_iscn_cn():
    assert iscn_cn_features("-7") == {"-7"}
    assert iscn_cn_features("del(5)(q13q33)") == {"5q-"}
    assert iscn_cn_features("i(17)(q10)") == {"17q+", "17p-"}
    assert iscn_cn_features("+i(17)(q10)") == {"17q+"}
    assert iscn_cn_features("t(8;21)(q22;q22)") == set()


def test_ploidy():
    assert ploidy_features("55", "57") == {"hyperdiploid51_67"}
    assert ploidy_features("92", "94") == {"near_tetraploid"}
    assert ploidy_features("46", "46") == set()


def test_complex_rule():
    assert "complex" in clone_features({"-5", "-7", "+8"})
    assert "complex" not in clone_features({"t(8;21)(q22;q22)", "-Y", "+8"})


def test_cytoconverter_calls_are_added():
    f = clone_features({"der(5)t(5;17)(q11;p11)"}, cn_calls=["5q-", "17p-"])
    assert {"5q-", "17p-"} <= f


def test_frozen_vocabulary_uses_known_features():
    vocab = json.loads((ROOT / "config/vocabulary_frozen.json").read_text())
    known = set(STRUCTURAL) | {"complex", "hyperdiploid51_67", "hypodiploid_lt44", "near_tetraploid",
                               "-X", "-Y", "+X", "+Y"}
    import re
    for d, spec in vocab["diseases"].items():
        assert 1 <= len(spec["events"]) <= vocab["max_events"]
        for e, anyof in spec["events"].items():
            for f in anyof:
                assert f in known or re.fullmatch(r"[+-]\d+|\d+[pq][+-]", f), (d, e, f)


def test_events_for():
    dv = {"5q loss": ["-5", "5q-"], "+8": ["+8"]}
    assert events_for("5q-|+8|complex", dv) == ["5q loss", "+8"]
    assert events_for("", dv) == []
