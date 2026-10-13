# Author: Shubhodeep Chatterjee
"""Model-free event timing per ordering disease from within-patient clone orders (Bradley-Terry,
karyoprog/timing.py), with case-bootstrap intervals. Count-derived events are left out.

results/tables/clone_timing_<disease>.csv   earliness score, 90% interval, rank and rank interval
results/models/clone_timing_<disease>.npz   scores and bootstrap scores (for validation and figures)
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from karyoprog.common import TABLES, MODELS, CONFIG, DERIVED_EVENTS, ordering_diseases
from karyoprog.pipeline import load_vocabulary
from karyoprog.timing import clone_timing

cfg = CONFIG["clone_timing"]
vocab = load_vocabulary()
pairs = pd.read_csv(TABLES / "ordered_pairs.tsv.gz", sep="\t")
pairs = pairs[~pairs.is_duplicate.astype(bool)]
for d in ordering_diseases():
    events = [e for e in vocab["diseases"][d]["events"] if e not in DERIVED_EVENTS]
    tab, keep, s, S = clone_timing(pairs[pairs.disease == d], events, cfg["n_bootstrap"],
                                   cfg["min_pairs_per_event"], cfg["ridge"], cfg["seed"])
    tab.to_csv(TABLES / f"clone_timing_{d}.csv", index=False)
    np.savez_compressed(MODELS / f"clone_timing_{d}.npz", events=np.array(keep), s=s, boot=S)
    print(d, " > ".join(f"{r.event} [{r.rank_ci05}-{r.rank_ci95}]" for r in tab.itertuples()), flush=True)
