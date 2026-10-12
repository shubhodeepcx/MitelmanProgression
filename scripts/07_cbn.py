# Author: Shubhodeep Chatterjee
"""Conjunctive Bayesian networks per ordering disease, on its 10 most frequent events
(count-derived events such as 'complex' excluded, see karyoprog.common).

For each disease: fit (greedy BIC), 2-fold held-out log-likelihood against an
independence model and an MHN on the same events, 20 bootstrap refits (edge
stability), and a simulation-recovery check (3 datasets drawn from the fitted
CBN at the real sample size, refitted). Writes results/models/cbn_<disease>.json.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from multiprocessing import Pool
import numpy as np
from karyoprog import cbn, mhnfit
from karyoprog.common import MODELS, ordering_diseases, load_matrix, save_json

N_EVENTS, N_BOOT, N_SIM = 10, 20, 3


def closure(P):
    C = np.asarray(P, bool).copy()
    for k in range(len(C)):
        C |= C[:, [k]] & C[[k], :]
    return C


def _boot(args):
    X, seed = args
    rng = np.random.default_rng(seed)
    return cbn.fit(X[rng.integers(0, len(X), len(X))])["P"]


def _sim(args):
    model, N, seed = args
    rng = np.random.default_rng(seed)
    S = cbn.sample(model, N, rng)
    return cbn.fit(S)["P"]


def _cv(args):
    Xa, Xb = args
    M = cbn.fit(Xa)
    T = mhnfit.fit(Xa)
    return cbn.loglik(M, Xb), mhnfit.mean_loglik(T, Xb), mhnfit.independence_loglik(Xa, Xb)


if __name__ == "__main__":
    MODELS.mkdir(parents=True, exist_ok=True)
    for d in ordering_diseases():
        m, _ = load_matrix(d, for_model=True)
        events = list(m.sum().sort_values(ascending=False).index[:N_EVENTS])
        X = m[events].values.astype(np.int8)
        model = cbn.fit(X)
        rng = np.random.default_rng(7)
        perm = rng.permutation(len(X)); half = len(X) // 2
        A, B = X[perm[:half]], X[perm[half:]]
        with Pool(4) as pool:
            boots = pool.map(_boot, [(X, 100 + b) for b in range(N_BOOT)])
            sims = pool.map(_sim, [(model, len(X), 200 + s) for s in range(N_SIM)])
            cv = pool.map(_cv, [(A, B), (B, A)])
        Ctrue = closure(model["P"])
        stab = np.mean([closure(P) for P in boots], axis=0)
        rec = []
        for P in sims:
            Cs = closure(P)
            tp = (Cs & Ctrue).sum()
            rec.append(dict(precision=tp / max(Cs.sum(), 1), recall=tp / max(Ctrue.sum(), 1)))
        cvm = np.mean(cv, axis=0)
        out = dict(disease=d, events=events, n_cases=len(X), P=model["P"], P_reduced=model["P_reduced"],
                   theta=model["theta"], eps=model["eps"], loglik=model["loglik"], bic=model["bic"],
                   closure_bootstrap_frequency=stab, simulation_recovery=rec,
                   heldout_loglik=dict(cbn=cvm[0], mhn_same_events=cvm[1], independence=cvm[2]))
        save_json(out, MODELS / f"cbn_{d}.json")
        edges = [f"{events[a]}->{events[b]} ({stab[a, b]:.2f})" for a, b in zip(*np.nonzero(model["P_reduced"]))]
        print(f"{d}: eps={model['eps']:.3f} edges={len(edges)} heldout cbn={cvm[0]:.3f} mhn={cvm[1]:.3f} indep={cvm[2]:.3f}")
        print("   ", "; ".join(edges))
        print("    sim recovery", [{k: round(v, 2) for k, v in r.items()} for r in rec])
