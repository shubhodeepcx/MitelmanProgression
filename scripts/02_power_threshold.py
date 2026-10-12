# Author: Shubhodeep Chatterjee
"""Choose the per-disease minimum event count N before any modelling.

Criterion: with at most 20 events (190 pairs, Bonferroni alpha = 0.05/190, two-sided),
an event seen in N cases must give >= 80% power, against a partner event of 10%
frequency, to detect BOTH
  - co-occurrence with odds ratio 3, and
  - mutual exclusivity with odds ratio 0.2,
using the exact conditional (Fisher) test. Cohort size n = non-duplicate cases.
N is the larger of the two minima. Written to config/power_threshold.json.
"""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, nchypergeom_fisher

ALPHA = 0.05 / 190
POWER = 0.8
PARTNER_F = 0.10
OR_CO, OR_EX = 3.0, 0.2


def power(n, k, m, odds):
    """Power of the two-sided exact test at ALPHA for an event in k of n cases, partner in m."""
    x = np.arange(0, min(k, m) + 1)
    p_hi = hypergeom.sf(x - 1, n, m, k)        # P(X >= x) under the null
    p_lo = hypergeom.cdf(x, n, m, k)           # P(X <= x)
    reject = (p_hi <= ALPHA / 2) | (p_lo <= ALPHA / 2)
    pmf = nchypergeom_fisher.pmf(x, n, m, k, odds)
    return pmf[reject].sum()


def min_count(n, odds):
    m = int(round(PARTNER_F * n))
    for k in range(5, n):
        if power(n, k, m, odds) >= POWER:
            return k
    return None


cases = pd.read_csv(ROOT / "data/interim/cases.tsv.gz", sep="\t")
cases = cases[~cases.is_duplicate.astype(bool)]
out = {"criterion": __doc__.strip(), "diseases": {}}
for d in ["AML", "APL", "MDS", "CML", "ALL", "CLL", "MM"]:
    n = int((cases.disease == d).sum())
    k_co, k_ex = min_count(n, OR_CO), min_count(n, OR_EX)
    N = max(k_co, k_ex)
    out["diseases"][d] = {"n_cases": n, "N_cooccurrence_OR3": k_co, "N_exclusivity_OR0.2": k_ex,
                          "N": N, "N_as_fraction": round(N / n, 4)}
    print(d, out["diseases"][d])
(ROOT / "config/power_threshold.json").write_text(json.dumps(out, indent=2))
