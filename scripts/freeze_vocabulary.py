"""Freeze the per-disease event vocabulary.

For each disease: count non-duplicate cases whose first investigation has the event in
any abnormal clone; keep events with count >= N (config/power_threshold.json), at most
max_events by frequency; defining events are always kept. Writes
config/vocabulary_frozen.json and results/tables/vocabulary_selection.csv.
"""
import sys, json, hashlib
from pathlib import Path
from datetime import date
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog.pipeline import events_for, INTERIM

cand = json.loads((ROOT / "config/vocabulary_candidates.json").read_text())
power = json.loads((ROOT / "config/power_threshold.json").read_text())["diseases"]
cf = pd.read_csv(INTERIM / "clone_features.tsv.gz", sep="\t", keep_default_na=False)
cases = pd.read_csv(INTERIM / "cases.tsv.gz", sep="\t")
keep = cases[~cases.is_duplicate.astype(bool)].set_index("case_id").disease
cf = cf[cf.case_id.isin(keep.index)]
cf = cf[cf.InvNo == cf.groupby("case_id").InvNo.transform("min")]
cf["disease"] = cf.case_id.map(keep)

zip_md5 = hashlib.md5((ROOT / "data/raw/mitelman_db.zip").read_bytes()).hexdigest()
frozen = {"frozen_on": str(date.today()), "mitelman_zip_md5": zip_md5,
          "rule": "count >= N (power_threshold.json), top max_events by count, defining always kept",
          "max_events": cand["max_events"], "diseases": {}}
rows = []
for d, spec in cand["diseases"].items():
    sub = cf[cf.disease == d]
    n = keep.eq(d).sum()
    N = power[d]["N"]
    counts = {e: sub[sub.features.map(lambda f: e in events_for(f, {e: anyof}))].case_id.nunique()
              for e, anyof in spec["events"].items()}
    ranked = sorted(counts, key=lambda e: -counts[e])
    chosen = [e for e in ranked if e in spec["defining"]]
    for e in ranked:
        if len(chosen) >= cand["max_events"]:
            break
        if e not in chosen and counts[e] >= N:
            chosen.append(e)
    for e in ranked:
        rows.append(dict(disease=d, event=e, cases=counts[e], pct=round(100 * counts[e] / n, 2),
                         N=N, selected=e in chosen, defining=e in spec["defining"]))
    frozen["diseases"][d] = {"n_cases": int(n), "N": N,
                             "events": {e: spec["events"][e] for e in ranked if e in chosen}}
    print(f"{d}: n={n} N={N} -> {len(chosen)} events: {', '.join(e for e in ranked if e in chosen)}")
(ROOT / "config/vocabulary_frozen.json").write_text(json.dumps(frozen, indent=2))
out = ROOT / "results/tables"; out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out / "vocabulary_selection.csv", index=False)
