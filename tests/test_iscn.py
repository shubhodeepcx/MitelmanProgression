"""Hand-written parser cases. The 50-karyotype curated set lives in test_hand_curated.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from karyoprog.iscn import split_top, normalise_token, parse_investigation, relate, ordered_pairs
from karyoprog.events import events_of, is_complex


def test_split_top_respects_parentheses():
    assert split_top("46,XY,t(8;21)(q22;q22),der(1)t(1;3)(p36;q21)") == \
        ["46", "XY", "t(8;21)(q22;q22)", "der(1)t(1;3)(p36;q21)"]


def test_normalise_token():
    assert normalise_token("+8x2") == ("+8", 2, False, False)
    assert normalise_token("+21c") == ("+21", 1, True, False)
    assert normalise_token("add(7)(q?22)")[0] == "add(7)(q22)"


def test_idem_expansion_and_tree():
    cl = parse_investigation(["46,XY,t(8;21)(q22;q22)", "45,idem,-Y", "46,idem,-Y,+8"])
    assert cl[1].abn == {"t(8;21)(q22;q22)", "-Y"}
    parent, depth = relate(cl)
    assert parent == {1: None, 2: 1, 3: 2} and depth == {1: 0, 2: 1, 3: 2}
    assert ordered_pairs(cl, key=events_of) == {("t(8;21)", "-Y"), ("t(8;21)", "+8"), ("-Y", "+8")}


def test_polyploid_idem_adds_no_events():
    cl = parse_investigation(["47,XY,+8", "94,idemx2"])
    assert cl[1].abn == {"+8"}
    assert ordered_pairs(cl) == set()


def test_constitutional_and_markers_ignored():
    cl = parse_investigation(["48,XY,+8,+21c,+mar"])
    assert cl[0].abn == {"+8"}


def test_normal_clone_is_not_a_parent():
    cl = parse_investigation(["46,XX,inv(16)(p13q22)", "47,idem,+22", "46,XX"])
    parent, _ = relate(cl)
    assert parent == {1: None, 2: 1}


def test_unrelated_clones_give_no_order():
    cl = parse_investigation(["47,XY,+8", "45,X,-Y"])
    assert ordered_pairs(cl) == set()


def test_shared_late_event_not_used_by_default():
    cl = parse_investigation(["46,XY,del(5)(q13q33),-Y", "45,X,-Y,-7"])
    assert ordered_pairs(cl) == set()
    assert ("-Y", "-7") in ordered_pairs(cl, implied=True)


def test_event_rules():
    assert events_of("der(11)t(9;11)(p22;q23)") == ["KMT2A(11q23)"]
    assert events_of("t(11;19)(q23;p13.1)") == ["KMT2A(11q23)"]
    assert events_of("t(8;21;12)(q22;q22;p13)") == ["t(8;21)"]
    assert events_of("t(16;16)(p13;q22)") == ["inv(16)"]
    assert events_of("del(5)(q13q33)") == ["-5/del(5q)"]
    assert events_of("del(11)(q23)") == []            # deletion is not a KMT2A rearrangement
    assert events_of("+8") == ["+8"] and events_of("+18") == []


def test_complex_excludes_recurrent_translocations():
    assert is_complex({"-5", "-7", "+8"})
    assert not is_complex({"t(15;17)(q24;q21)", "+8", "-7"})
