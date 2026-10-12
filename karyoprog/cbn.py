# Author: Shubhodeep Chatterjee
"""Conjunctive Bayesian network with observation noise: the discrete-time CBN of
Beerenwinkel, Eriksson and Sturmfels (Bernoulli, 2007) with the uniform per-bit
error model used by H-CBN (Gerstung et al., Bioinformatics, 2009).

A poset P on the events says which events must precede which: an event can only
be present if all its parents are. Hidden genotypes compatible with P follow
    P(g) = prod over events i whose parents are all in g of theta_i^g_i (1-theta_i)^(1-g_i),
and each observed bit is flipped independently with probability epsilon.
theta and epsilon are fitted by EM; the poset is chosen by greedy hill-climbing
on BIC with add / remove / reverse moves, starting from the empty poset.
The reference implementations (ct-cbn, h-cbn) are C programs that could not be
obtained in the analysis environment; this follows the published model, with a
discrete (Bernoulli) rather than continuous-time parameterisation.
Hidden states are enumerated, so use at most ~10-12 events.
"""
import itertools
import numpy as np


def _all_genotypes(n):
    return ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(np.int8)


def is_acyclic(P):
    n = len(P)
    indeg = P.sum(0).astype(int).copy()
    stack = [i for i in range(n) if indeg[i] == 0]
    seen = 0
    while stack:
        v = stack.pop(); seen += 1
        for w in np.nonzero(P[v])[0]:
            indeg[w] -= 1
            if indeg[w] == 0:
                stack.append(w)
    return seen == n


def transitive_reduction(P):
    n = len(P)
    C = P.astype(bool).copy()
    for k in range(n):
        C |= C[:, [k]] & C[[k], :]
    R = P.astype(bool).copy()
    for a, b in zip(*np.nonzero(R)):
        if any(C[a, c] and C[c, b] for c in range(n) if c not in (a, b)):
            R[a, b] = False
    return R.astype(np.int8)


class Data:
    """Unique observations with counts, plus Hamming distances to every genotype."""

    def __init__(self, X):
        X = np.asarray(X, dtype=np.int8)
        self.uniq, self.counts = np.unique(X, axis=0, return_counts=True)
        self.N, self.n = X.shape
        self.allg = _all_genotypes(self.n)
        self.D_all = (self.uniq[:, None, :] != self.allg[None, :, :]).sum(2).astype(np.int16)


def _compatible_mask(P, allg):
    ok = np.ones(len(allg), bool)
    for a, b in zip(*np.nonzero(P)):
        ok &= ~((allg[:, b] == 1) & (allg[:, a] == 0))
    return ok


def em(P, data, theta=None, eps=0.05, iters=200, tol=1e-7):
    """Returns (loglik, theta, eps)."""
    mask = _compatible_mask(P, data.allg)
    G = data.allg[mask]
    D = data.D_all[:, mask]
    n = data.n
    pa_ok = np.ones((len(G), n), bool)                     # parents of i all present in g
    for a, b in zip(*np.nonzero(P)):
        pa_ok[:, b] &= G[:, a] == 1
    G1 = G == 1
    theta = np.full(n, 0.2) if theta is None else theta.copy()
    ll_old = -np.inf
    for _ in range(iters):
        th = np.clip(theta, 1e-6, 1 - 1e-6)
        logpg = (pa_ok * np.where(G1, np.log(th), np.log(1 - th))).sum(1)
        e = min(max(eps, 1e-6), 0.499)
        L = D * np.log(e / (1 - e)) + n * np.log(1 - e) + logpg[None, :]
        mx = L.max(1, keepdims=True)
        W = np.exp(L - mx)
        tot = W.sum(1, keepdims=True)
        ll = float((data.counts * (np.log(tot[:, 0]) + mx[:, 0])).sum())
        W *= data.counts[:, None] / tot
        wg = W.sum(0)                                        # expected count of each hidden genotype
        num = (wg[:, None] * (G1 & pa_ok)).sum(0)
        den = (wg[:, None] * pa_ok).sum(0)
        theta = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        eps = float((W * D).sum() / (data.N * n))
        if ll - ll_old < tol * abs(ll):
            break
        ll_old = ll
    return ll, theta, eps


def _bic(ll, P, data):
    return -2 * ll + (P.sum() + data.n + 1) * np.log(data.N)


def fit(X, max_steps=60, verbose=False, candidate_iters=15):
    """Greedy BIC search. Candidates are scored with a short warm-started EM; the
    best few are refitted to convergence. Returns dict(P, P_reduced, theta, eps, loglik, bic)."""
    data = Data(X)
    n = data.n
    P = np.zeros((n, n), np.int8)
    ll, th, ep = em(P, data)
    bic = _bic(ll, P, data)
    for step in range(max_steps):
        cands = []
        for a, b in itertools.permutations(range(n), 2):
            for move in ("add", "remove", "reverse"):
                Q = P.copy()
                if move == "add" and not P[a, b] and not P[b, a]:
                    Q[a, b] = 1
                elif move == "remove" and P[a, b]:
                    Q[a, b] = 0
                elif move == "reverse" and P[a, b]:
                    Q[a, b] = 0; Q[b, a] = 1
                else:
                    continue
                if is_acyclic(Q):
                    l2, _, _ = em(Q, data, th, ep, iters=candidate_iters)
                    cands.append((_bic(l2, Q, data), Q))
        improved = False
        for _, Q in sorted(cands, key=lambda c: c[0])[:3]:
            l3, t3, e3 = em(Q, data, th, ep)
            b3 = _bic(l3, Q, data)
            if b3 < bic - 1e-6:
                bic, P, ll, th, ep = b3, Q, l3, t3, e3
                improved = True
                break
        if not improved:
            break
        if verbose:
            print(step, round(bic, 1), int(P.sum()))
    return dict(P=P, P_reduced=transitive_reduction(P), theta=th, eps=ep, loglik=ll, bic=bic)


def loglik(model, X):
    """Mean log-likelihood per case of X under a fitted CBN (parameters held fixed)."""
    data = Data(X)
    return em(model["P"], data, model["theta"], model["eps"], iters=1)[0] / data.N


def sample(model, n_samples, rng):
    P, th, eps = model["P"], model["theta"], model["eps"]
    n = len(th)
    order = []
    indeg = P.sum(0).astype(int).copy(); stack = [i for i in range(n) if indeg[i] == 0]
    while stack:
        v = stack.pop(); order.append(v)
        for w in np.nonzero(P[v])[0]:
            indeg[w] -= 1
            if indeg[w] == 0:
                stack.append(w)
    G = np.zeros((n_samples, n), np.int8)
    for i in order:
        parents = np.nonzero(P[:, i])[0]
        allowed = G[:, parents].all(1) if len(parents) else np.ones(n_samples, bool)
        G[:, i] = allowed & (rng.random(n_samples) < th[i])
    flip = rng.random(G.shape) < eps
    return np.where(flip, 1 - G, G).astype(np.int8)
