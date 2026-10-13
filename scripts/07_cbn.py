# Author: Shubhodeep Chatterjee
"""Conjunctive Bayesian networks per ordering disease, on its most frequent model events
(count-derived events such as 'complex' excluded, see karyoprog.common).

For each disease: fit (greedy BIC); k-fold held-out log-likelihood against an independence
model and an MHN on the same events; bootstrap refits (stability of each implied order);
simulation recovery (data sets drawn from the fitted CBN at the real sample size, refitted).
Counts, folds and seeds come from config/analysis.json. Writes results/models/cbn_<disease>.json.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from multiprocessing import Pool
import numpy as np
from karyoprog import cbn, mhnfit
from karyoprog.common import MODELS, CONFIG, ordering_diseases, load_matrix, save_json, kfold_indices

CFG = CONFIG["cbn"]


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
    S = cbn.sample(model, N, np.random.default_rng(seed))
    return cbn.fit(S)["P"]


def _cv(args):
    Xa, Xb = args
    M = cbn.fit(Xa)
    T = mhnfit.fit(Xa)
    return cbn.loglik(M, Xb), mhnfit.mean_loglik(T, Xb), mhnfit.independence_loglik(Xa, Xb)


if __name__ == "__main__":
    MODELS.mkdir(parents=True, exist_ok=True)
    diseases = sys.argv[1:] or ordering_diseases()
    for d in diseases:
        m, _ = load_matrix(d, for_model=True)
        events = list(m.sum().sort_values(ascending=False).index[:CFG["n_events"]])
        X = m[events].values.astype(np.int8)
        model = cbn.fit(X)
        folds = kfold_indices(len(X), CFG["cv_folds"], CFG["seed"])
        s0 = CFG["seed"] * 1000
        with Pool(CONFIG["processes"]) as pool:
            boots = pool.map(_boot, [(X, s0 + b) for b in range(CFG["n_bootstrap"])])
            sims = pool.map(_sim, [(model, len(X), s0 + 500 + s) for s in range(CFG["n_simulations"])])
            cv = np.array(pool.map(_cv, [(X[tr], X[te]) for tr, te in folds]))
        Ctrue = closure(model["P"])
        stab = np.mean([closure(P) for P in boots], axis=0)
        rec = []
        for P in sims:
            Cs = closure(P)
            tp = (Cs & Ctrue).sum()
            rec.append(dict(precision=tp / max(Cs.sum(), 1), recall=tp / max(Ctrue.sum(), 1)))
        cvm, cvs = cv.mean(0), cv.std(0, ddof=1)
        out = dict(disease=d, events=events, n_cases=len(X), P=model["P"], P_reduced=model["P_reduced"],
                   theta=model["theta"], eps=model["eps"], loglik=model["loglik"], bic=model["bic"],
                   closure_bootstrap_frequency=stab, simulation_recovery=rec,
                   heldout_loglik=dict(cbn=cvm[0], mhn_same_events=cvm[1], independence=cvm[2]),
                   heldout_loglik_sd=dict(cbn=cvs[0], mhn_same_events=cvs[1], independence=cvs[2]),
                   heldout_folds=cv, settings=CFG)
        save_json(out, MODELS / f"cbn_{d}.json")
        edges = [f"{events[a]}->{events[b]} ({stab[a, b]:.2f})" for a, b in zip(*np.nonzero(model["P_reduced"]))]
        print(f"{d}: eps={model['eps']:.3f} edges={len(edges)} heldout cbn={cvm[0]:.3f} mhn={cvm[1]:.3f} indep={cvm[2]:.3f}")
        print("   ", "; ".join(edges))
        print("    sim recovery", np.mean([r["precision"] for r in rec]).round(2), np.mean([r["recall"] for r in rec]).round(2), flush=True)
