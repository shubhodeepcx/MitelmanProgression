"""Week 1 go/no-go: how many cases per disease carry usable multi-clone order information?

Writes results/week1/{multiclone_by_disease.csv, event_table_provisional.tsv.gz,
ordered_pairs_provisional.tsv.gz, investigations.tsv.gz}.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog.build import build_tables

OUT = ROOT / "results" / "week1"
OUT.mkdir(parents=True, exist_ok=True)

cases, inv, ev, pairs = build_tables()
inv.to_csv(OUT / "investigations.tsv.gz", sep="\t", index=False)
ev.to_csv(OUT / "event_table_provisional.tsv.gz", sep="\t", index=False)
pairs.to_csv(OUT / "ordered_pairs_provisional.tsv.gz", sep="\t", index=False)

c = inv.groupby("case_id").agg(
    disease=("disease", "first"), Series=("Series", "first"),
    n_inv=("InvNo", "nunique"),
    multi_clone=("n_clones", lambda s: (s >= 2).any()),
    multi_abnormal=("n_abnormal", lambda s: (s >= 2).any()),
    related_edge=("n_edges", lambda s: (s > 0).any()),
    ab_order=("n_ab_pairs", lambda s: (s > 0).any()),
    ev_order=("n_ev_pairs", lambda s: (s > 0).any()),
    ev_order_implied=("n_ev_pairs_implied", lambda s: (s > 0).any()),
).reset_index()
c["longitudinal"] = c.n_inv >= 2

cols = ["multi_clone", "multi_abnormal", "related_edge", "ab_order", "ev_order",
        "ev_order_implied", "longitudinal"]
tab = c.groupby("disease")[cols].sum().astype(int)
tab.insert(0, "cases", c.groupby("disease").size())
tab["ev_order_unselected"] = c[c.Series == "U"].groupby("disease").ev_order.sum()
tab = tab.fillna(0).astype(int).sort_values("cases", ascending=False)
tab["go_ordering(ev_order>=200)"] = tab.ev_order >= 200
tab.to_csv(OUT / "multiclone_by_disease.csv")
print(tab.to_string())
print("\nall cases:", len(c), " with event-level order:", int(c.ev_order.sum()))
