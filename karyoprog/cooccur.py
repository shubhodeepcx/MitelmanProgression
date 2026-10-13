# Author: Shubhodeep Chatterjee
"""Pairwise co-occurrence and mutual exclusivity: the baseline floor.

For each event pair: Fisher's exact test with Benjamini-Hochberg correction, and
a permutation null that keeps every case's number of events fixed and the event
frequencies matched approximately in expectation. Null cases draw their k events without
replacement with weights calibrated by iterative scaling (Gumbel top-k), which
removes the "busy karyotype" effect: complex cases carry many events and would
otherwise make everything look co-occurring.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.empty_like(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


def _null_weights(X, iters=50):
    k = X.sum(1)
    target = X.sum(0).astype(float)
    w = np.maximum(target, 0.5) / target.sum()
    for _ in range(iters):
        # expected inclusion approximated by k * w / sum(w) per row, capped at 1
        incl = np.minimum(1.0, k[:, None] * w[None, :] / w.sum())
        exp = incl.sum(0)
        w = w * np.where(exp > 0, target / exp, 1.0)
    return w


def permutation_null(X, n_perm=1000, rng=None, chunk=50):
    """Yields n_perm null co-occurrence count matrices (n x n)."""
    rng = rng or np.random.default_rng(0)
    X = np.asarray(X, dtype=np.int8)
    k = X.sum(1)
    logw = np.log(_null_weights(X))
    n = X.shape[1]
    done = 0
    while done < n_perm:
        b = min(chunk, n_perm - done)
        g = rng.gumbel(size=(b, X.shape[0], n)) + logw
        rank = np.argsort(-g, axis=2).argsort(axis=2)            # rank of each event in each row
        Y = (rank < k[None, :, None]).astype(np.float32)
        for y in Y:
            yield y.T @ y
        done += b


def pairwise(X, names, n_perm=1000, rng=None, q=0.05):
    X = np.asarray(X, dtype=np.int64)
    N, n = X.shape
    C = X.T @ X
    nulls = np.array(list(permutation_null(X, n_perm, rng)))
    mean = nulls.mean(0)
    rows = []
    for a in range(n):
        for b in range(a + 1, n):
            n11 = C[a, b]; na = C[a, a]; nb = C[b, b]
            n10, n01 = na - n11, nb - n11
            n00 = N - n11 - n10 - n01
            OR = ((n11 + .5) * (n00 + .5)) / ((n10 + .5) * (n01 + .5))
            p = fisher_exact([[n11, n10], [n01, n00]])[1]
            dev = abs(n11 - mean[a, b])
            p_perm = (1 + (np.abs(nulls[:, a, b] - mean[a, b]) >= dev - 1e-9).sum()) / (n_perm + 1)
            rows.append(dict(event_a=names[a], event_b=names[b], n_a=int(na), n_b=int(nb), n_both=int(n11),
                             expected_null=round(float(mean[a, b]), 2), odds_ratio=OR, p_fisher=p, p_perm=p_perm))
    df = pd.DataFrame(rows)
    df["q_fisher"] = bh(df.p_fisher)
    df["q_perm"] = bh(df.p_perm)
    df["direction"] = np.where(df.n_both > df.expected_null, "co-occurring", "exclusive")
    # significant only if both tests reject and agree on the direction (Fisher's is relative to
    # independence, the permutation test's to the count-preserving null)
    same_dir = (df.odds_ratio > 1) == (df.direction == "co-occurring")
    df["significant"] = (df.q_fisher < q) & (df.q_perm < q) & same_dir
    return df
