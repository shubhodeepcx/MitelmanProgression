# Author: Shubhodeep Chatterjee
"""Step 4: the released tables, built with config/vocabulary_frozen.json.

results/tables/
  event_table.tsv.gz        one row per (case, investigation, abnormal clone, event): the tidy table
  ordered_pairs.tsv.gz      one row per (case, earlier event, later event) from parent->child clone edges
  matrix_<disease>.tsv.gz   case x event 0/1 matrix, first investigation, duplicates removed (model input)
  ordering_yield.csv        cases per disease with >=1 event-level ordered pair (go/no-go, frozen vocabulary)
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog.pipeline import load_vocabulary, tidy_event_table, INTERIM

OUT = ROOT / "results/tables"; OUT.mkdir(parents=True, exist_ok=True)
vocab = load_vocabulary()
cf = pd.read_csv(INTERIM / "clone_features.tsv.gz", sep="\t", keep_default_na=False,
                 dtype={"parent_clone": "Int64"})
cases = pd.read_csv(INTERIM / "cases.tsv.gz", sep="\t")
cases["is_duplicate"] = cases.is_duplicate.astype(bool)

ev = tidy_event_table(cf, cases, vocab)
ev.to_csv(OUT / "event_table.tsv.gz", sep="\t", index=False)

# ordered pairs from parent -> child edges, per investigation
evsets = ev.groupby(["case_id", "InvNo", "clone_id"]).event.apply(set)
pairs = []
for r in cf[cf.parent_clone.notna()].itertuples():
    child = evsets.get((r.case_id, r.InvNo, r.clone_id), set())
    par = evsets.get((r.case_id, r.InvNo, int(r.parent_clone)), set())
    for e in par:
        for l in child - par:
            pairs.append((r.case_id, r.InvNo, r.clone_id, e, l))
pairs = pd.DataFrame(pairs, columns=["case_id", "InvNo", "clone_id", "earlier", "later"]).drop_duplicates(
    ["case_id", "earlier", "later"])
meta = cases.set_index("case_id")[["disease", "Series", "is_duplicate"]]
pairs = pairs.join(meta, on="case_id")
pairs.to_csv(OUT / "ordered_pairs.tsv.gz", sep="\t", index=False)

# case x event matrices (first investigation, any clone, non-duplicates; all-zero rows kept)
# "first investigation" = the earliest investigation with an abnormal clone (as in 03), not the earliest
# one that happens to carry a vocabulary event
first_inv = cf.groupby("case_id").InvNo.min()
first = ev[ev.InvNo == ev.case_id.map(first_inv)]
for d, spec in vocab["diseases"].items():
    ids = cases[(cases.disease == d) & ~cases.is_duplicate].case_id
    m = (first[first.disease == d].assign(v=1).pivot_table(index="case_id", columns="event", values="v",
                                                          aggfunc="max", fill_value=0)
         .reindex(index=ids, columns=list(spec["events"]), fill_value=0).fillna(0).astype(int))
    m.to_csv(OUT / f"matrix_{d}.tsv.gz", sep="\t")

nd = pairs[~pairs.is_duplicate]
y = pd.DataFrame({
    "cases": cases[~cases.is_duplicate].groupby("disease").size(),
    "events_frozen": pd.Series({d: len(s["events"]) for d, s in vocab["diseases"].items()}),
    "cases_with_order": nd.groupby("disease").case_id.nunique(),
    "of_which_unselected": nd[nd.Series == "U"].groupby("disease").case_id.nunique(),
}).loc[list(vocab["diseases"])].fillna(0).astype(int)
y["meets_200_rule"] = y.cases_with_order >= 200
y.to_csv(OUT / "ordering_yield.csv")
print(y.to_string())
print("event rows:", len(ev), " ordered pairs:", len(pairs))
