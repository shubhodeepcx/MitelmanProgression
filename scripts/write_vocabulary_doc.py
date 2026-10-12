# Author: Shubhodeep Chatterjee
"""Regenerate docs/event_vocabulary.md from config/vocabulary_frozen.json and
results/tables/vocabulary_selection.csv (run after 03; it does not re-freeze anything)."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd

frozen = json.loads((ROOT / "config/vocabulary_frozen.json").read_text())
sel = pd.read_csv(ROOT / "results/tables/vocabulary_selection.csv")
L = ["# Frozen event vocabulary (v1)", "",
     "Shubhodeep Chatterjee, MitelmanProgression.", "",
     f"Frozen on {frozen['frozen_on']} against the Mitelman download with md5 `{frozen['mitelman_zip_md5']}`.",
     "Rule: an event enters a disease's vocabulary if it is present (any abnormal clone, first investigation, "
     "duplicates removed) in at least N cases, where N comes from the power calculation in "
     "`config/power_threshold.json`; at most 20 events, ranked by count; defining events are always kept. "
     "Candidates were listed in `config/vocabulary_candidates.json` before any ordering result was looked at.", "",
     "An event is present in a clone when any of its atomic features is (see `karyoprog/features.py`). "
     "Copy-number features are the union of what the ISCN text states and CytoConverter arm-level calls "
     "(at least 50% of an arm, at least 90% of a chromosome for whole gains and losses, relative to the clone's "
     "own baseline ploidy).", "",
     "Count-derived events (`complex`, `high hyperdiploid`, `hypodiploid`) are kept in the co-occurrence "
     "baseline but excluded from CBN and MHN, because they are functions of the other events.", ""]
for d, spec in frozen["diseases"].items():
    s = sel[sel.disease == d]
    L += [f"## {d} (n = {spec['n_cases']:,}, N = {spec['N']})", "",
          "| Event | Atomic features | Cases | % |", "|---|---|---|---|"]
    for e, feats in spec["events"].items():
        r = s[s.event == e].iloc[0]
        L.append(f"| {e} | {', '.join(feats)} | {int(r.cases):,} | {r.pct} |")
    rest = s[~s.selected.astype(bool)]
    if len(rest):
        L += ["", "Below N or beyond the cap: " + ", ".join(f"{r.event} ({int(r.cases)})" for r in rest.itertuples())]
    L.append("")
L += ["## Known limitations", "",
      "- `complex` follows ELN 2022 (at least 3 aberrations, no recurrent balanced rearrangement). It is partly a "
      "function of other events (for example 5q loss, -7, 17p loss), so it is used only in the co-occurrence baseline. "
      "It is not a candidate in ALL or MM, where hyperdiploidy is the relevant count-derived class.",
      "- In ALL, `high hyperdiploid` and the individual trisomies describe one process. Ordering between the "
      "trisomies mostly reflects how the hyperdiploid clone accumulated gains, not separate steps.",
      "- In CML, `17p loss` includes i(17q), so `17p loss` and `i(17q)` are nested by definition; an edge between "
      "them is not biological evidence.",
      "- t(12;21) (ETV6-RUNX1) is cryptic by karyotyping, so Mitelman understates it.",
      "- APL ends up with 2 events and MM has fewer than 200 cases with a clone order, so both are "
      "co-occurrence-only under the stop rule.", ""]
(ROOT / "docs/event_vocabulary.md").write_text("\n".join(L))
print("wrote docs/event_vocabulary.md")
