# Author: Shubhodeep Chatterjee
"""Mutual Hazard Networks per ordering disease (the headline model), with validation fits.

For each disease (count-derived events excluded):
  1. lambda = c/N chosen by k-fold held-out log-likelihood over the c grid; the independence
     model is scored on the same folds.
  2. final fit on all cases; event-order probabilities by simulation.
  3. bootstrap refits (edge and order stability), warm-started from the final fit.
  4. sensitivity refits: Unselected series only (publication bias), cases published from 2000 on (era).
  5. simulation recovery: data sets drawn from the fitted MHN at the real sample size, refitted
     from the independence start.
Grid, counts, folds and seeds come from config/analysis.json.
Writes results/models/mhn_<disease>.npz and mhn_<disease>_summary.json.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from multiprocessing import Pool
import numpy as np
from karyoprog import mhnfit
from karyoprog.common import MODELS, CONFIG, THRESH, ordering_diseases, load_matrix, save_json, kfold_indices

CFG = CONFIG["mhn"]


def _fit(job):
    X, lam, init = job
    return mhnfit.fit(X, lam=lam, init=init)


def _heldout(job):
    Xa, Xb, c = job
    return mhnfit.mean_loglik(mhnfit.fit(Xa, lam=c / len(Xa)), Xb)


def _orders(job):
    T, n, seed = job
    return mhnfit.order_probabilities(T, n, np.random.default_rng(seed))[0]


if __name__ == "__main__":
    MODELS.mkdir(parents=True, exist_ok=True)
    diseases = sys.argv[1:] or ordering_diseases()
    grid = CFG["lambda_c_grid"]
    for d in diseases:
        m, cov = load_matrix(d, for_model=True)
        X = m.values.astype(np.int8)
        N = len(X)
        folds = kfold_indices(N, CFG["cv_folds"], CFG["seed"])
        s0 = CFG["seed"] * 1000
        with Pool(CONFIG["processes"]) as pool:
            jobs = [(X[tr], X[te], c) for c in grid for tr, te in folds]
            hl = np.array(pool.map(_heldout, jobs)).reshape(len(grid), len(folds))
            cv = {c: float(hl[i].mean()) for i, c in enumerate(grid)}
            cv_sd = {c: float(hl[i].std(ddof=1)) for i, c in enumerate(grid)}
            ind = np.array([mhnfit.independence_loglik(X[tr], X[te]) for tr, te in folds])
            c_best = max(cv, key=cv.get)
            lam = c_best / N
            T = mhnfit.fit(X, lam=lam)
            boot_idx = [np.random.default_rng(s0 + b).integers(0, N, N) for b in range(CFG["n_bootstrap"])]
            U = X[(cov.Series == "U").values]
            E = X[(cov.Year >= 2000).values]
            sims = [mhnfit.sample(T, N, np.random.default_rng(s0 + 5000 + s)) for s in range(CFG["n_simulations"])]
            jobs = ([(X[i], lam, T) for i in boot_idx] + [(U, c_best / len(U), T), (E, c_best / len(E), T)]
                    + [(S, lam, None) for S in sims])
            fits = pool.map(_fit, jobs)
            nb = CFG["n_bootstrap"]
            boots = np.array(fits[:nb]); T_U, T_E = fits[nb], fits[nb + 1]; T_sims = fits[nb + 2:]
            P_boot = np.array(pool.map(_orders, [(t, CFG["bootstrap_order_trajectories"], s0 + 9000 + b)
                                                 for b, t in enumerate(boots)]))
        P, both = mhnfit.order_probabilities(T, CFG["order_trajectories"], np.random.default_rng(s0 + 1))
        np.savez_compressed(MODELS / f"mhn_{d}.npz", events=np.array(m.columns), theta=T, boot=boots, theta_unselected=T_U,
                            theta_post2000=T_E, order_P=P, order_both=both, order_P_boot=P_boot)
        off = ~np.eye(len(T), dtype=bool)
        strong = THRESH["edge_abs_log_theta"]
        summary = dict(disease=d, n_cases=N, n_events=len(m.columns), events=list(m.columns), lambda_c=c_best, lambda_=lam,
                       heldout_loglik=dict(mhn_by_c=cv, mhn_by_c_sd=cv_sd, independence=float(ind.mean()),
                                           independence_sd=float(ind.std(ddof=1)), mhn_folds=hl[grid.index(c_best)],
                                           independence_folds=ind),
                       n_unselected=int(len(U)), n_post2000=int(len(E)),
                       simulation_recovery=[mhnfit.edge_recovery(T, t, strong) for t in T_sims],
                       sensitivity=dict(unselected=mhnfit.edge_recovery(T, T_U, strong),
                                        post2000=mhnfit.edge_recovery(T, T_E, strong)),
                       n_strong_edges=int((off & (np.abs(T) >= strong)).sum()), settings=CFG)
        save_json(summary, MODELS / f"mhn_{d}_summary.json")
        sim = {k: round(float(np.mean([r[k] for r in summary["simulation_recovery"]])), 2) for k in ("precision", "recall", "corr_offdiag")}
        print(d, "c*", c_best, "heldout", {k: round(v, 4) for k, v in cv.items()}, "indep", round(float(ind.mean()), 4))
        print("   sim", sim, " sens U r", round(summary["sensitivity"]["unselected"]["corr_offdiag"], 2),
              "post2000 r", round(summary["sensitivity"]["post2000"]["corr_offdiag"], 2), flush=True)
