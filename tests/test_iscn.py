# Author: Shubhodeep Chatterjee
"""Parser unit cases. The 50-karyotype annotated set lives in test_hand_curated.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from karyoprog.iscn import split_top, normalise_token, parse_clone, parse_investigation, relate, ordered_pairs
from karyoprog.features import aberration_features


def test_split_top_respects_parentheses():
    assert split_top("46,XY,t(8;21)(q22;q22),der(1)t(1;3)(p36;q21)") == \
        ["46", "XY", "t(8;21)(q22;q22)", "der(1)t(1;3)(p36;q21)"]


def test_normalise_token():
    assert normalise_token("+8x2") == ("+8", 2, False, False)
    assert normalise_token("+der(20)del(20)(q11)x1-4")[0] == "+der(20)del(20)(q11)"
    assert normalise_token("+21c") == ("+21", 1, True, False)
    assert normalise_token("add(7)(q?22)")[0] == "add(7)(q22)"


def test_markers_and_split_sex_complement_ignored():
    assert parse_clone("48,XY,+8,+3-4mar,30-100dmin,inc").own == ["+8"]
    assert parse_clone("47,del(X)(p22p22),Y,+21c").own == ["del(X)(p22p22)"]


def test_idem_expansion_and_tree():
    cl = parse_investigation(["46,XY,t(8;21)(q22;q22)", "45,idem,-Y", "46,idem,-Y,+8"])
    assert cl[1].abn == {"t(8;21)(q22;q22)", "-Y"}
    parent, depth = relate(cl)
    assert parent == {1: None, 2: 1, 3: 2} and depth == {1: 0, 2: 1, 3: 2}
    assert ordered_pairs(cl, key=aberration_features) == {("t(8;21)", "-Y"), ("t(8;21)", "+8"), ("-Y", "+8")}


def test_idem_opposite_numerics():
    cl = parse_investigation(["45,XX,-21,del(5)(q13q33)", "46,idem,+21"])
    assert cl[1].abn == {"del(5)(q13q33)"}                       # plain idem: cancels
    cl = parse_investigation(["45,XX,-21,del(5)(q13q33)", "91,idemx2,+21"])
    assert cl[1].abn == {"-21", "del(5)(q13q33)"}                # idemx2: net loss remains


def test_no_cycle_when_sideline_loses_stemline_aberration():
    cl = parse_investigation(["47,XY,+8,del(9)(q13q22)", "46,idem,-8"])
    parent, _ = relate(cl)
    assert parent == {1: 2, 2: None}


def test_polyploid_idem_adds_no_events():
    cl = parse_investigation(["47,XY,+8", "94,idemx2"])
    assert cl[1].abn == {"+8"} and ordered_pairs(cl) == set()


def test_normal_clone_is_not_a_parent():
    cl = parse_investigation(["46,XX,inv(16)(p13q22)", "47,idem,+22", "46,XX"])
    assert relate(cl)[0] == {1: None, 2: 1}


def test_unrelated_clones_give_no_order():
    assert ordered_pairs(parse_investigation(["47,XY,+8", "45,X,-Y"])) == set()


def test_shared_late_event_not_used_by_default():
    cl = parse_investigation(["46,XY,del(5)(q13q33),-Y", "45,X,-Y,-7"])
    assert ordered_pairs(cl) == set()
    assert ("-Y", "-7") in ordered_pairs(cl, implied=True)
