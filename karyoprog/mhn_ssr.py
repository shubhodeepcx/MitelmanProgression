# Author: Shubhodeep Chatterjee
"""Mutual Hazard Network likelihood and gradient with state-space restriction.

Model (Schill, Solbrig and Spang, Bioinformatics 2020): events accumulate as a
continuous-time Markov chain; the rate of event i in state s is
    lambda_i(s) = Theta_ii * prod_{j in s} Theta_ij,
and a tumour is observed at an Exp(1) time, so p = (I - Q)^{-1} e_empty.
For an observed state x only states that are subsets of x carry probability
towards x, so the marginal likelihood is a triangular solve on the 2^|x|
subsets of x, with the diagonal of Q still summing over all n events.

This re-implements the computation of the reference `mhn` package
(mhn.ssr.state_space_restriction) in vectorised numpy, batching samples with
the same number of events. It agrees with the package to ~1e-10 wherever the
package returns finite values (tests/test_mhn.py); in our build the package
returned NaN for every state with 5 or 6 events, which is why it is not
called directly. Parameters are log(Theta).
"""
import numpy as np


class StateData:
    """Unique observed states grouped by event count, with multiplicities."""

    def __init__(self, X):
        X = np.asarray(X, dtype=np.int8)
        self.n = X.shape[1]
        self.N = X.shape[0]
        uniq, counts = np.unique(X, axis=0, return_counts=True)
        self.groups = []
        k_of = uniq.sum(1)
        for k in np.unique(k_of):
            sel = k_of == k
            U = uniq[sel]
            active = np.array([np.flatnonzero(u) for u in U], dtype=np.int64).reshape(len(U), k)
            nb = 1 << k
            bits = ((np.arange(nb)[:, None] >> np.arange(k)[None, :]) & 1).astype(np.float64) if k else np.zeros((1, 0))
            self.groups.append(dict(k=int(k), active=active, counts=counts[sel].astype(float), bits=bits))


def _group_terms(log_theta, g, want_grad):
    k, A, bits = g["k"], g["active"], g["bits"]
    n = log_theta.shape[0]
    m = A.shape[0]
    nb = 1 << k
    diag = np.diag(log_theta)
    # L[s, b, i] = log rate of event i in subset b of sample s
    if k:
        thA = log_theta[:, A].transpose(1, 2, 0)            # (m, k, n): theta[i, a_l]
        L = diag[None, None, :] + np.einsum("bl,sli->sbi", bits, thA)
    else:
        L = np.broadcast_to(diag, (m, 1, n)).copy()
    R = np.exp(L)
    # zero rates of events already present in b
    present = np.zeros((m, nb, n), dtype=bool)
    for l in range(k):
        present[np.arange(m)[:, None], np.flatnonzero(bits[:, l])[None, :], A[:, l][:, None]] = True
    R[present] = 0.0
    d = 1.0 + R.sum(2)                                       # (m, nb) diagonal of I - Q
    # forward solve p = (I-Q)^{-1} e_0, layer by popcount
    pop = bits.sum(1).astype(int) if k else np.zeros(1, int)
    p = np.zeros((m, nb))
    p[:, 0] = 1.0 / d[:, 0]
    rows = np.arange(m)
    for layer in range(1, k + 1):
        idx = np.flatnonzero(pop == layer)
        acc = np.zeros((m, len(idx)))
        for l in range(k):
            has = bits[idx, l] == 1
            if not has.any():
                continue
            b = idx[has]
            prev = b ^ (1 << l)
            acc[:, has] += R[rows[:, None], prev[None, :], A[:, l][:, None]] * p[:, prev]
        p[:, idx] = acc / d[:, idx]
    px = p[:, nb - 1]
    if not want_grad:
        return px, None
    # backward solve q = (I-Q)^{-T} e_full
    q = np.zeros((m, nb))
    q[:, nb - 1] = 1.0 / d[:, nb - 1]
    for layer in range(k - 1, -1, -1):
        idx = np.flatnonzero(pop == layer)
        acc = np.zeros((m, len(idx)))
        for l in range(k):
            lacks = bits[idx, l] == 0
            if not lacks.any():
                continue
            b = idx[lacks]
            nxt = b | (1 << l)
            acc[:, lacks] += R[rows[:, None], b[None, :], A[:, l][:, None]] * q[:, nxt]
        q[:, idx] = acc / d[:, idx]
    # W[s, b, i] = rate_i(b) p(b) (q(b+i)[i in A] - q(b))
    qnext = np.zeros((m, nb, n))
    for l in range(k):
        lacks = np.flatnonzero(bits[:, l] == 0)
        qnext[rows[:, None], lacks[None, :], A[:, l][:, None]] = q[:, lacks | (1 << l)]
    W = R * p[:, :, None] * (qnext - q[:, :, None])
    W /= px[:, None, None]
    G = np.zeros((m, n, n))
    Wsum = W.sum(1)                                          # (m, n)
    G[:, np.arange(n), np.arange(n)] = Wsum
    if k:
        off = np.einsum("sbi,bl->sil", W, bits)              # (m, n, k): sum_b W[b,i] [a_l in b]
        for l in range(k):
            G[rows, :, A[:, l]] += off[:, :, l]
    return px, G


def score_and_gradient(log_theta, data):
    """Mean log-likelihood and its gradient (n x n) over all samples in `data`."""
    n = data.n
    total, grad = 0.0, np.zeros((n, n))
    for g in data.groups:
        px, G = _group_terms(log_theta, g, True)
        total += (g["counts"] * np.log(px)).sum()
        grad += np.einsum("s,sij->ij", g["counts"], G)
    return total / data.N, grad / data.N


def score(log_theta, data):
    total = 0.0
    for g in data.groups:
        px, _ = _group_terms(log_theta, g, False)
        total += (g["counts"] * np.log(px)).sum()
    return total / data.N
