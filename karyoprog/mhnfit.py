# Author: Shubhodeep Chatterjee
"""Fitting and using Mutual Hazard Networks (Schill et al., Bioinformatics 2020).

Likelihood and gradient: karyoprog.mhn_ssr (state-space restriction, checked
against the reference `mhn` package). Optimiser: L-BFGS-B on the mean
log-likelihood minus lambda times the smooth L1 penalty on off-diagonal
entries used by the reference package. Default lambda = 1/N as in the package;
scripts/08_mhn.py chooses lambda by held-out likelihood.

Theta convention: log_theta[i, i] is the log base rate of event i;
log_theta[i, j] (i != j) is the log multiplicative effect of event j, once
present, on the rate of event i. Positive = j promotes i, negative = j inhibits i.
"""
import numpy as np
from scipy.optimize import minimize
from . import mhn_ssr

EPS = 1e-5


def l1(th):
    t = th.copy(); np.fill_diagonal(t, 0)
    return np.sqrt(t ** 2 + EPS).sum()


def l1_grad(th):
    t = th.copy(); np.fill_diagonal(t, 0)
    return t / np.sqrt(t ** 2 + EPS)


def independence_theta(X):
    """Independence MHN start: Theta_ii = empirical odds, off-diagonals 0 (log scale)."""
    f = np.clip(np.asarray(X).mean(0), 1e-6, 1 - 1e-6)
    return np.diag(np.log(f / (1 - f)))


def fit(X, lam=None, maxit=3000, init=None):
    """Fit an MHN to a 0/1 matrix X (cases x events). Returns log_theta (n x n)."""
    X = np.asarray(X)
    data = mhn_ssr.StateData(X)
    n = X.shape[1]
    lam = 1.0 / X.shape[0] if lam is None else lam

    def f(x):
        th = x.reshape(n, n)
        s, g = mhn_ssr.score_and_gradient(th, data)
        return -(s - lam * l1(th)), -(g - lam * l1_grad(th)).ravel()

    x0 = (independence_theta(X) if init is None else init).ravel()
    r = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": maxit, "gtol": 1e-6})
    if not np.isfinite(r.fun):
        raise RuntimeError("MHN optimisation produced a non-finite score")
    return r.x.reshape(n, n)


def mean_loglik(log_theta, X):
    """Mean log-likelihood per case of X under the MHN."""
    return mhn_ssr.score(log_theta, mhn_ssr.StateData(X))


def independence_loglik(Xtrain, Xtest):
    """Mean log-likelihood of Xtest under independent Bernoulli events fitted on Xtrain."""
    p = np.clip(np.asarray(Xtrain).mean(0), 1e-6, 1 - 1e-6)
    Xtest = np.asarray(Xtest)
    return float((Xtest * np.log(p) + (1 - Xtest) * np.log(1 - p)).sum(1).mean())


def _trajectory(th, rng):
    """One accumulation trajectory up to an Exp(1) observation time; returns event order."""
    n = th.shape[0]
    base = np.diag(th)
    state = np.zeros(n, dtype=bool)
    order = []
    t_obs, t = rng.exponential(1.0), 0.0
    while True:
        off = np.where(state[None, :], th, 1.0)
        np.fill_diagonal(off, 1.0)
        rates = np.where(state, 0.0, base * off.prod(axis=1))   # rate_i = theta_ii * prod_{j present} theta_ij
        total = rates.sum()
        if total <= 0:
            break
        t += rng.exponential(1.0 / total)
        if t > t_obs:
            break
        i = rng.choice(n, p=rates / total)
        state[i] = True
        order.append(i)
    return order


def sample(log_theta, n_samples, rng):
    """Draw observed 0/1 states from the MHN."""
    th = np.exp(log_theta)
    out = np.zeros((n_samples, th.shape[0]), dtype=np.int32)
    for k in range(n_samples):
        out[k, _trajectory(th, rng)] = 1
    return out


def order_probabilities(log_theta, n_traj=20000, rng=None):
    """P(event i occurs before event j | both occur before observation), by simulation.

    Returns (P, both) where P[i, j] = P(i before j | both observed) and both[i, j] is the
    fraction of trajectories in which both were observed.
    """
    rng = rng or np.random.default_rng(0)
    th = np.exp(log_theta)
    n = th.shape[0]
    before = np.zeros((n, n))
    both = np.zeros((n, n))
    for _ in range(n_traj):
        order = _trajectory(th, rng)
        for a_pos, a in enumerate(order):
            for b in order[a_pos + 1:]:
                before[a, b] += 1
                both[a, b] += 1
                both[b, a] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(both > 0, before / both, np.nan)
    return P, both / n_traj


def edge_recovery(T_true, T_hat, strong=0.5):
    """Agreement between two MHNs on their edges (|log theta| >= strong, same sign)."""
    off = ~np.eye(len(T_true), dtype=bool)
    true = off & (np.abs(T_true) >= strong)
    hat = off & (np.abs(T_hat) >= strong)
    tp = (true & hat & (np.sign(T_true) == np.sign(T_hat))).sum()
    return dict(n_true=int(true.sum()), n_called=int(hat.sum()),
                recall=float(tp / max(true.sum(), 1)), precision=float(tp / max(hat.sum(), 1)),
                corr_offdiag=float(np.corrcoef(T_true[off], T_hat[off])[0, 1]))
