# Author: Shubhodeep Chatterjee
"""Model-free event timing from within-patient clone orders (Bradley-Terry).

Each parent->child clone edge says that every event in the parent came before every event the
child adds. Pooling these pairs over the cases of one disease, the Bradley-Terry model
    P(a before b) = 1 / (1 + exp(-(s_a - s_b)))
gives every event an earliness score s (higher = earlier), identified by centring and a small
ridge penalty. Uncertainty comes from resampling cases, so a case that contributes many pairs
is resampled as a unit. Bradley and Terry (1952); this is the standard paired-comparison model.
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize


def pair_counts(pairs, events):
    """W[i, j] = number of cases with events[i] before events[j]."""
    idx = {e: i for i, e in enumerate(events)}
    W = np.zeros((len(events), len(events)))
    for (a, b), n in pairs.groupby(["earlier", "later"]).size().items():
        if a in idx and b in idx:
            W[idx[a], idx[b]] += n
    return W


def fit_bt(W, ridge=0.01):
    """Earliness scores s (centred) maximising the penalised Bradley-Terry likelihood."""
    n = len(W)
    N = W + W.T

    def f(s):
        d = s[:, None] - s[None, :]
        ll = (W * -np.logaddexp(0, -d)).sum()
        p = 1 / (1 + np.exp(-d))                       # p_ij = P(i before j)
        g = (W - N * p).sum(1)                         # d ll / d s_i
        return -(ll - ridge * (s ** 2).sum()), -(g - 2 * ridge * s)
    r = minimize(f, np.zeros(n), jac=True, method="L-BFGS-B")
    s = r.x - r.x.mean()
    return s


def clone_timing(pairs, events, n_boot=500, min_pairs=10, ridge=0.01, seed=0):
    """Scores with case-bootstrap intervals for events seen in >= min_pairs ordered pairs.

    pairs: rows of (case_id, earlier, later) for one disease (duplicates already removed).
    Returns (table, events_used, s, S_boot)."""
    W = pair_counts(pairs, events)
    keep = [e for i, e in enumerate(events) if W[i].sum() + W[:, i].sum() >= min_pairs]
    W = pair_counts(pairs, keep)
    s = fit_bt(W, ridge)
    rng = np.random.default_rng(seed)
    by_case = {c: g for c, g in pairs[pairs.earlier.isin(keep) & pairs.later.isin(keep)].groupby("case_id")}
    case_ids = np.array(list(by_case))
    Wc = {c: pair_counts(g, keep) for c, g in by_case.items()}
    S = np.empty((n_boot, len(keep)))
    for b in range(n_boot):
        draw = rng.choice(case_ids, len(case_ids), replace=True)
        Wb = np.sum([Wc[c] for c in draw], axis=0)
        S[b] = fit_bt(Wb, ridge)
    lo, hi = np.percentile(S, [5, 95], axis=0)
    rank_boot = (-S).argsort(1).argsort(1) + 1
    tab = pd.DataFrame(dict(event=keep, earliness=s, ci05=lo, ci95=hi,
                            rank=(-s).argsort().argsort() + 1,
                            rank_ci05=np.percentile(rank_boot, 5, axis=0).astype(int),
                            rank_ci95=np.percentile(rank_boot, 95, axis=0).astype(int),
                            pairs_involving=(W.sum(1) + W.sum(0)).astype(int)))
    return tab.sort_values("earliness", ascending=False).reset_index(drop=True), keep, s, S


def p_before(s):
    """Matrix of P(i before j) implied by Bradley-Terry scores."""
    return 1 / (1 + np.exp(-(s[:, None] - s[None, :])))
