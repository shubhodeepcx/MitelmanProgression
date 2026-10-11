"""Extra ordering signal from cases with several investigations (e.g. diagnosis and relapse).

Pools the abnormal clones of all investigations of a case (idem resolved within
each investigation) and counts cases that gain event-level ordering only this way.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog import io
from karyoprog.iscn import parse_investigation, ordered_pairs
from karyoprog.events import events_of

cl = io.load_clones(); cases = io.load_cases()
multi = cl.groupby(["RefNo", "CaseNo"]).InvNo.transform("nunique") > 1
rows = []
for (ref, case), g in cl[multi].groupby(["RefNo", "CaseNo"]):
    within, pooled = set(), []
    for inv, gi in g.groupby("InvNo"):
        c = parse_investigation(list(gi.clone))
        within |= ordered_pairs(c, key=events_of)
        pooled += c
    for i, c in enumerate(pooled):
        c.index = i + 1; c.uses_idem = False
    across = ordered_pairs(pooled, key=events_of)
    rows.append(dict(RefNo=ref, CaseNo=case, within=bool(within), pooled=bool(across)))
d = pd.DataFrame(rows).merge(cases[["RefNo", "CaseNo", "disease"]], on=["RefNo", "CaseNo"])
d["gain_only_longitudinal"] = d.pooled & ~d.within
t = d.groupby("disease")[["within", "pooled", "gain_only_longitudinal"]].sum()
t.insert(0, "multi_inv_cases", d.groupby("disease").size())
t.to_csv(ROOT / "results/week1/longitudinal_by_disease.csv")
print(t.sort_values("multi_inv_cases", ascending=False).to_string())
