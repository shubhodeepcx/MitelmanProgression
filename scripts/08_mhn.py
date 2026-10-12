# Author: Shubhodeep Chatterjee
"""Mutual Hazard Networks per ordering disease (the headline model), with validation fits.

For each disease (count-derived events excluded):
  1. lambda chosen by 2-fold held-out log-likelihood over c/N, c in {0.1, 0.3, 1, 3};
     the independence model is scored on the same folds.
  2. final fit on all cases; event-order probabilities by simulation.
  3. 25 bootstrap refits (edge and order stability).
  4. sensitivity refits: Unselected series only (publication bias), cases published
     from 2000 on (era).
  5. simulation recovery: 3 datasets drawn from the fitted MHN at the real sample
     size, refitted from the independence start.
Writes results/models/mhn_<disease>.npz and mhn_<disease>_summary.json.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from multiprocessing import Pool
import numpy as np
from karyoprog import mhnfit
from karyoprog.common import MODELS, ordering_diseases, load_matrix, save_json

C_GRID = [0.1, 0.3, 1.0, 3.0]
N_BOOT, N_SIM = 25, 3
STRONG = 0.5            # |log theta| threshold for an interaction to count as an edge


def _fit(job):
    X, lam, init = job
    return mhnfit.fit(X, lam=lam, init=init)


def _heldout(job):
    Xa, Xb, lam = job
    T = mhnfit.fit(Xa, lam=lam)
    return mhnfit.mean_loglik(T, Xb)


def edge_recovery(T_true, T_hat):
    off = ~np.eye(len(T_true), dtype=bool)
    true = off & (np.abs(T_true) >= STRONG)
    hat = off & (np.abs(T_hat) >= STRONG)
    same = np.sign(T_true) == np.sign(T_hat)
    tp = (true & hat & same).sum()
    return dict(n_true=int(true.sum()), n_called=int(hat.sum()),
                recall=float(tp / max(true.sum(), 1)), precision=float(tp / max(hat.sum(), 1)),
                corr_offdiag=float(np.corrcoef(T_true[off], T_hat[off])[0, 1]))


if __name__ == "__main__":
    MODELS.mkdir(parents=True, exist_ok=True)
    diseases = sys.argv[1:] or ordering_diseases()
    for d in diseases:
        m, cov = load_matrix(d, for_model=True)
        X = m.values.astype(np.int8)
        N = len(X)
        rng = np.random.default_rng(11)
        perm = rng.permutation(N); h = N // 2
        A, B = X[perm[:h]], X[perm[h:]]
        with Pool(4) as pool:
            jobs = [(A, B, c / h) for c in C_GRID] + [(B, A, c / (N - h)) for c in C_GRID]
            hl = pool.map(_heldout, jobs)
            cv = {c: float(np.mean([hl[i], hl[i + len(C_GRID)]])) for i, c in enumerate(C_GRID)}
            indep = float(np.mean([mhnfit.independence_loglik(A, B), mhnfit.independence_loglik(B, A)]))
            c_best = max(cv, key=cv.get)
            lam = c_best / N
            T = mhnfit.fit(X, lam=lam)
            boot_idx = [np.random.default_rng(1000 + b).integers(0, N, N) for b in range(N_BOOT)]
            U = X[(cov.Series == "U").values]
            E = X[(cov.Year >= 2000).values]
            sims = [mhnfit.sample(T, N, np.random.default_rng(2000 + s)) for s in range(N_SIM)]
            jobs = ([(X[i], lam, T) for i in boot_idx] + [(U, c_best / len(U), T), (E, c_best / len(E), T)]
                    + [(S, lam, None) for S in sims])
            fits = pool.map(_fit, jobs)
        boots = np.array(fits[:N_BOOT]); T_U, T_E = fits[N_BOOT], fits[N_BOOT + 1]; T_sims = fits[N_BOOT + 2:]
        P, both = mhnfit.order_probabilities(T, 20000, np.random.default_rng(5))
        P_boot = np.array([mhnfit.order_probabilities(t, 4000, np.random.default_rng(6 + b))[0] for b, t in enumerate(boots)])
        np.savez_compressed(MODELS / f"mhn_{d}.npz", events=np.array(m.columns), theta=T, boot=boots, theta_unselected=T_U,
                            theta_post2000=T_E, order_P=P, order_both=both, order_P_boot=P_boot)
        off = ~np.eye(len(T), dtype=bool)
        summary = dict(disease=d, n_cases=N, n_events=len(m.columns), events=list(m.columns), lambda_c=c_best, lambda_=lam,
                       heldout_loglik=dict(mhn_by_c=cv, independence=indep),
                       n_unselected=int(len(U)), n_post2000=int(len(E)),
                       simulation_recovery=[edge_recovery(T, t) for t in T_sims],
                       sensitivity=dict(unselected=edge_recovery(T, T_U), post2000=edge_recovery(T, T_E)),
                       n_strong_edges=int((off & (np.abs(T) >= STRONG)).sum()))
        save_json(summary, MODELS / f"mhn_{d}_summary.json")
        print(d, "c*", c_best, "heldout", {k: round(v, 4) for k, v in cv.items()}, "indep", round(indep, 4))
        print("   sim", [{k: round(v, 2) for k, v in r.items()} for r in summary["simulation_recovery"]])
        print("   sens U", {k: round(v, 2) for k, v in summary["sensitivity"]["unselected"].items()},
              "post2000", {k: round(v, 2) for k, v in summary["sensitivity"]["post2000"].items()}, flush=True)
