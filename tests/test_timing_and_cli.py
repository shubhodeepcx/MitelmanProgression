# Author: Shubhodeep Chatterjee
"""Bradley-Terry clone timing, k-fold splits, the analysis config and the command-line entry point."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from karyoprog.timing import fit_bt, clone_timing, p_before, pair_counts
from karyoprog.common import kfold_indices, CONFIG, THRESH


def test_bt_gradient_matches_finite_differences():
    from scipy.optimize import approx_fprime
    rng = np.random.default_rng(0)
    W = rng.integers(0, 20, (4, 4)).astype(float); np.fill_diagonal(W, 0)
    s = fit_bt(W, ridge=0.01)
    # at the optimum the penalised gradient is ~0
    ll = lambda x: (W * -np.logaddexp(0, -(x[:, None] - x[None, :]))).sum() - 0.01 * (x ** 2).sum()
    g = approx_fprime(s, ll, 1e-6)
    assert np.allclose(g, 0, atol=1e-3)


def test_bt_recovers_a_known_order():
    rng = np.random.default_rng(1)
    true = np.array([2.0, 1.0, 0.0, -1.5])
    P = p_before(true)
    rows = []
    for case in range(400):
        i, j = rng.choice(4, 2, replace=False)
        a, b = (i, j) if rng.random() < P[i, j] else (j, i)
        rows.append((f"c{case}", f"e{a}", f"e{b}"))
    pairs = pd.DataFrame(rows, columns=["case_id", "earlier", "later"])
    tab, keep, s, S = clone_timing(pairs, [f"e{k}" for k in range(4)], n_boot=50, min_pairs=1, seed=2)
    assert list(tab.event) == ["e0", "e1", "e2", "e3"]
    assert np.corrcoef(s, true - true.mean())[0, 1] > 0.95


def test_pair_counts():
    pairs = pd.DataFrame({"case_id": ["a", "b", "b"], "earlier": ["x", "x", "y"], "later": ["y", "y", "z"]})
    W = pair_counts(pairs, ["x", "y", "z"])
    assert W[0, 1] == 2 and W[1, 2] == 1 and W.sum() == 3


def test_kfold_partitions_every_case_once():
    folds = kfold_indices(103, 5, seed=0)
    test = np.concatenate([te for _, te in folds])
    assert sorted(test) == list(range(103))
    for tr, te in folds:
        assert not set(tr) & set(te) and len(tr) + len(te) == 103


def test_analysis_config_is_complete():
    cfg = json.loads((ROOT / "config" / "analysis.json").read_text())
    for k in ("edge_abs_log_theta", "stable_sign_agreement", "sim_recovery_precision", "sim_recovery_recall",
              "supported_share_of_stable_edges", "fdr_q", "min_cases_with_clone_order"):
        assert k in cfg["thresholds"]
    assert THRESH == cfg["thresholds"] and CONFIG["mhn"]["cv_folds"] >= 2


def test_cli_lists_steps():
    out = subprocess.run([sys.executable, "-m", "karyoprog", "steps"], cwd=ROOT, capture_output=True, text=True, check=True)
    assert "08_mhn.py" in out.stdout and "09_clone_timing.py" in out.stdout and len(out.stdout.splitlines()) == 11
