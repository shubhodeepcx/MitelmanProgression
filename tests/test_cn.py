# Author: Shubhodeep Chatterjee
import sys
from pathlib import Path
import pytest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from karyoprog.cn import call_clone, AUTOSOMES

# Toy arm map: every chromosome 0-100, centromere at 40 (acrocentric p arms ignored).
ARMS = {}
for c in AUTOSOMES:
    ARMS[(c, "p")] = (0, 40)
    ARMS[(c, "q")] = (40, 100)


def whole(c, v, n=1):
    return [(0, 100, v)] * n


def test_trisomy_and_monosomy():
    base, f = call_clone({"8": whole("8", 1), "7": whole("7", -1)}, ARMS)
    assert base == 2 and set(f) == {"+8", "-7"}


def test_arm_loss_threshold():
    _, f = call_clone({"5": [(60, 100, -1)]}, ARMS)     # 40/60 of 5q
    assert f == ["5q-"]
    _, f = call_clone({"5": [(80, 100, -1)]}, ARMS)     # 20/60 of 5q: below threshold
    assert f == []


def test_isochromosome_17q():
    _, f = call_clone({"17": [(0, 40, -1), (40, 100, 1)]}, ARMS)
    assert set(f) == {"17p-", "17q+"}


def test_tetraploid_baseline_removes_uniform_gain():
    segs = {c: whole(c, 1, 2) for c in AUTOSOMES}       # 92,XXYY: every autosome +2
    segs["8"] = whole("8", 1, 4)                         # plus two extra 8s
    base, f = call_clone(segs, ARMS)
    assert base == 4 and f == ["+8"]


def test_acrocentric_uses_q_arm_only():
    _, f = call_clone({"21": [(40, 100, 1)]}, ARMS)
    assert f == ["+21"]


DATA = ROOT / "data" / "raw" / "mitelman" / "CYTOCONVERTED.TXT.DATA"


@pytest.mark.skipif(not DATA.exists(), reason="Mitelman download not present")
def test_cytoconverter_clone_numbers_follow_kclone_order():
    """Case 10719:1 is 46,XY,t(8;21)/45,idem,-Y/46,idem,-Y,+8: clones 2 and 3 lose Y, only 3 gains 8."""
    from karyoprog import io
    cc = io.read_table("CYTOCONVERTED")
    g = cc[(cc.RefNo == "10719") & (cc.CaseNo == "1")]
    assert set(g[g.Chr == "chrY"].Clone) == {"2", "3"}
    assert set(g[g.Chr == "chr8"].Clone) == {"3"}
