# Author: Shubhodeep Chatterjee
"""Checks on the model code: MHN likelihood and gradient, CBN recovery, BH, permutation null."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from karyoprog import mhn_ssr, mhnfit, cbn, cooccur


def test_mhn_diagonal_theta_matches_closed_form():
    # with no interactions, P(state x) = int_0^inf e^-t prod_{i in x}(1 - e^{-r_i t}) prod_{i not in x} e^{-r_i t} dt
    from scipy.integrate import quad
    r = np.array([0.6, 0.25, 0.1])
    T = np.diag(np.log(r))
    states = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 1], [1, 1, 1]])
    for x in states:
        f = lambda t: np.exp(-t) * np.prod(np.where(x == 1, 1 - np.exp(-r * t), np.exp(-r * t)))
        exact = np.log(quad(f, 0, np.inf)[0])
        assert np.isclose(mhnfit.mean_loglik(T, x[None, :]), exact, atol=1e-7)


def test_mhn_gradient_matches_finite_differences():
    rng = np.random.default_rng(1)
    X = (rng.random((300, 4)) < 0.3).astype(int)
    data = mhn_ssr.StateData(X)
    T = rng.normal(0, 0.5, (4, 4))
    _, g = mhn_ssr.score_and_gradient(T, data)
    h = 1e-6
    for i, j in [(0, 0), (1, 2), (3, 0), (2, 3)]:
        Tp, Tm = T.copy(), T.copy()
        Tp[i, j] += h; Tm[i, j] -= h
        fd = (mhn_ssr.score(Tp, data) - mhn_ssr.score(Tm, data)) / (2 * h)
        assert np.isclose(g[i, j], fd, rtol=1e-4, atol=1e-7)


def test_mhn_recovers_a_strong_promotion():
    rng = np.random.default_rng(2)
    T = np.diag(np.log([0.5, 0.05, 0.3]))
    T[1, 0] = 2.5                       # event 0 strongly promotes event 1
    X = mhnfit.sample(T, 3000, rng)
    That = mhnfit.fit(X, lam=1 / len(X))
    assert That[1, 0] > 1.0
    P, _ = mhnfit.order_probabilities(That, 3000, np.random.default_rng(3))
    assert P[0, 1] > 0.7


def test_cbn_recovers_a_chain():
    rng = np.random.default_rng(4)
    model = dict(P=np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], np.int8),
                 theta=np.array([0.6, 0.5, 0.5]), eps=0.02)
    X = cbn.sample(model, 3000, rng)
    fit = cbn.fit(X)
    assert fit["P_reduced"][0, 1] == 1 and fit["P_reduced"][1, 2] == 1
    assert fit["P_reduced"].sum() == 2
    assert abs(fit["eps"] - 0.02) < 0.02


def test_bh_matches_definition():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    q = cooccur.bh(p)
    assert np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.2])


def test_permutation_null_keeps_event_total():
    rng = np.random.default_rng(5)
    X = (rng.random((200, 6)) < 0.25).astype(int)
    for M in cooccur.permutation_null(X, n_perm=5, rng=rng):
        assert np.isclose(np.trace(M), X.sum())        # total events preserved
