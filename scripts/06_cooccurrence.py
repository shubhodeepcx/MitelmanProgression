# Author: Shubhodeep Chatterjee
"""Baseline: pairwise co-occurrence / mutual exclusivity for every disease (incl. APL, MM).

Writes results/tables/cooccurrence_<disease>.csv.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
from karyoprog.common import TABLES, CONFIG, ALL_DISEASES, load_matrix
from karyoprog.cooccur import pairwise

for d in ALL_DISEASES:
    m, _ = load_matrix(d)
    cfg = CONFIG["cooccurrence"]
    df = pairwise(m.values, list(m.columns), n_perm=cfg["n_permutations"], rng=np.random.default_rng(cfg["seed"]), q=CONFIG["thresholds"]["fdr_q"])
    df.to_csv(TABLES / f"cooccurrence_{d}.csv", index=False)
    sig = df[df.significant]
    print(f"{d}: {len(df)} pairs, {len(sig)} significant "
          f"({(sig.direction == 'co-occurring').sum()} co-occurring, {(sig.direction == 'exclusive').sum()} exclusive)")
