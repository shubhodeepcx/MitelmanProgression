# Author: Shubhodeep Chatterjee
"""Validation summary across models (run after 06, 07 and 08).

results/tables/
  validation_prior_knowledge.csv  each relationship in config/prior_knowledge.json checked against the
                                  co-occurrence baseline, the MHN (with bootstrap support), the CBN and the
                                  within-patient clone orders
  validation_clone_order.csv      per disease: do the MHN's simulated orders, and the CBN's implied orders,
                                  agree with the directions seen in parent->child clone pairs (pairs observed
                                  in >= 10 cases)?
  mhn_edges_<disease>.csv         every MHN interaction with its bootstrap sign stability
  validation_summary.csv          one row per disease: held-out likelihoods, simulation recovery,
                                  stability, sensitivity, and the stop-rule verdicts

Pre-specified thresholds (set before the fits): an MHN interaction is an edge if |log theta| >= 0.5;
an edge or an order is stable if it keeps its sign in >= 80% of bootstrap refits; simulation recovery is
adequate if mean precision >= 0.8 and mean recall >= 0.5; a disease's ordering claims are "supported"
when >= 50% of its edges are stable, otherwise "underdetermined beyond the stable subset".
"""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from karyoprog.common import ROOT as R, TABLES, MODELS, ALL_DISEASES, ordering_diseases

STRONG, STABLE = 0.5, 0.8
prior = json.loads((R / "config/prior_knowledge.json").read_text())
pairs = pd.read_csv(TABLES / "ordered_pairs.tsv.gz", sep="\t")
pairs = pairs[~pairs.is_duplicate.astype(bool)]
yield_ = pd.read_csv(TABLES / "ordering_yield.csv", index_col=0)
order_d = ordering_diseases()


def clone_counts(d, a, b):
    p = pairs[pairs.disease == d]
    return int(((p.earlier == a) & (p.later == b)).sum()), int(((p.earlier == b) & (p.later == a)).sum())


def baseline(d, a, b):
    c = pd.read_csv(TABLES / f"cooccurrence_{d}.csv")
    r = c[((c.event_a == a) & (c.event_b == b)) | ((c.event_a == b) & (c.event_b == a))]
    if r.empty:
        return None
    r = r.iloc[0]
    return r.direction if r.significant else "not significant"


pk_rows, co_rows, sum_rows = [], [], []
for d in ALL_DISEASES:
    mhn = np.load(MODELS / f"mhn_{d}.npz", allow_pickle=True) if d in order_d else None
    cbnm = json.loads((MODELS / f"cbn_{d}.json").read_text()) if d in order_d else None
    if mhn is not None:
        ev = list(mhn["events"]); T = mhn["theta"]; B = mhn["boot"]; P = mhn["order_P"]; PB = mhn["order_P_boot"]
        idx = {e: i for i, e in enumerate(ev)}
        # edge table
        er = []
        for i in range(len(ev)):
            for j in range(len(ev)):
                if i != j:
                    er.append(dict(source=ev[j], target=ev[i], log_theta=T[i, j],
                                   sign_stability=float(np.mean(np.sign(B[:, i, j]) == np.sign(T[i, j]))) if T[i, j] != 0 else np.nan,
                                   strong=abs(T[i, j]) >= STRONG))
        er = pd.DataFrame(er)
        er.to_csv(TABLES / f"mhn_edges_{d}.csv", index=False)
    for c in prior.get(d, []):
        a, b, kind = c["a"], c["b"], c["kind"]
        row = dict(disease=d, kind=kind, a=a, b=b, why=c["why"], baseline=baseline(d, a, b))
        fwd, rev = clone_counts(d, a, b)
        row.update(clone_a_before_b=fwd, clone_b_before_a=rev)
        if kind in ("exclusive", "cooccur"):
            row["baseline_pass"] = row["baseline"] == ("exclusive" if kind == "exclusive" else "co-occurring")
        if kind == "before":
            row["clone_pass"] = (fwd > rev) if fwd + rev >= 5 else None
        if mhn is not None and a in idx and b in idx:
            i, j = idx[a], idx[b]
            if kind == "before":
                row["mhn_P_a_before_b"] = round(float(P[i, j]), 3)
                row["mhn_order_support"] = round(float(np.mean(PB[:, i, j] > 0.5)), 2)
                row["mhn_pass"] = bool(P[i, j] > 0.5 and np.mean(PB[:, i, j] > 0.5) >= STABLE)
                C = np.asarray(cbnm["closure_bootstrap_frequency"]) if cbnm and a in cbnm["events"] and b in cbnm["events"] else None
                if C is not None:
                    ca, cb = cbnm["events"].index(a), cbnm["events"].index(b)
                    row["cbn_edge_support"] = round(float(C[ca, cb]), 2)
                    row["cbn_pass"] = bool(C[ca, cb] >= 0.5)
            else:
                s = T[i, j] + T[j, i]
                sb = B[:, i, j] + B[:, j, i]
                want = -1 if kind == "exclusive" else 1
                row["mhn_theta_sum"] = round(float(s), 2)
                row["mhn_sign_support"] = round(float(np.mean(np.sign(sb) == want)), 2)
                row["mhn_pass"] = bool(np.sign(s) == want and np.mean(np.sign(sb) == want) >= STABLE)
        pk_rows.append(row)
    if mhn is not None:
        p = pairs[pairs.disease == d]
        cnt = p.groupby(["earlier", "later"]).size()
        obs, pred, wts = [], [], []
        for (a, b), n_ab in cnt.items():
            if a >= b or a not in idx or b not in idx:
                continue
            n_ba = int(cnt.get((b, a), 0))
            tot = n_ab + n_ba
            if tot < 10 or np.isnan(P[idx[a], idx[b]]):
                continue
            obs.append(n_ab / tot); pred.append(P[idx[a], idx[b]]); wts.append(tot)
        # also count pairs seen only in the b<a direction
        for (b, a), n_ba in cnt.items():
            if a >= b or a not in idx or b not in idx or (a, b) in cnt.index:
                continue
            if n_ba >= 10 and not np.isnan(P[idx[a], idx[b]]):
                obs.append(0.0); pred.append(P[idx[a], idx[b]]); wts.append(n_ba)
        obs, pred, wts = map(np.array, (obs, pred, wts))
        agree = (np.sign(obs - 0.5) == np.sign(pred - 0.5)) & (obs != 0.5)
        # the CBN's implied orders (transitive closure, bootstrap frequency >= 0.5) against the same clone pairs
        cev = cbnm["events"]; Cb = np.asarray(cbnm["closure_bootstrap_frequency"])
        cb_n = cb_ok = 0
        for a_i, b_i in zip(*np.nonzero(Cb >= 0.5)):
            fwd, rev = int(cnt.get((cev[a_i], cev[b_i]), 0)), int(cnt.get((cev[b_i], cev[a_i]), 0))
            if fwd + rev >= 10:
                cb_n += 1; cb_ok += fwd > rev
        co_rows.append(dict(disease=d, pairs_compared=len(obs), cases_in_pairs=int(wts.sum()),
                            direction_agreement=float(agree.mean()) if len(obs) else np.nan,
                            weighted_agreement=float((agree * wts).sum() / wts.sum()) if len(obs) else np.nan,
                            spearman=float(spearmanr(obs, pred)[0]) if len(obs) > 2 else np.nan,
                            cbn_orders_compared=cb_n, cbn_agreement=cb_ok / cb_n if cb_n else np.nan))
        s = json.loads((MODELS / f"mhn_{d}_summary.json").read_text())
        strong = er[er.strong]
        stable = strong[strong.sign_stability >= STABLE]
        sim = pd.DataFrame(s["simulation_recovery"]).mean()
        frac_stable = len(stable) / max(len(strong), 1)
        sum_rows.append(dict(
            disease=d, cases=s["n_cases"], model_events=s["n_events"], cases_with_clone_order=int(yield_.loc[d, "cases_with_order"]),
            lambda_c=s["lambda_c"], heldout_mhn=max(s["heldout_loglik"]["mhn_by_c"].values()),
            heldout_independence=s["heldout_loglik"]["independence"],
            cbn_events=len(cbnm["events"]), heldout_cbn=cbnm["heldout_loglik"]["cbn"],
            heldout_mhn_cbn_events=cbnm["heldout_loglik"]["mhn_same_events"], heldout_indep_cbn_events=cbnm["heldout_loglik"]["independence"],
            cbn_eps=cbnm["eps"], cbn_edges=int(np.asarray(cbnm["P_reduced"]).sum()),
            cbn_sim_precision=float(np.mean([r["precision"] for r in cbnm["simulation_recovery"]])),
            cbn_sim_recall=float(np.mean([r["recall"] for r in cbnm["simulation_recovery"]])),
            mhn_strong_edges=len(strong), mhn_stable_edges=len(stable), mhn_frac_stable=frac_stable,
            mhn_sim_precision=sim["precision"], mhn_sim_recall=sim["recall"], mhn_sim_corr=sim["corr_offdiag"],
            unselected_agreement_recall=s["sensitivity"]["unselected"]["recall"], unselected_corr=s["sensitivity"]["unselected"]["corr_offdiag"],
            post2000_agreement_recall=s["sensitivity"]["post2000"]["recall"], post2000_corr=s["sensitivity"]["post2000"]["corr_offdiag"],
            verdict_sim="adequate" if (sim["precision"] >= 0.8 and sim["recall"] >= 0.5) else "inadequate",
            verdict_stability="supported" if frac_stable >= 0.5 else "underdetermined beyond the stable subset"))
    else:
        sum_rows.append(dict(disease=d, cases=int(yield_.loc[d, "cases"]), cases_with_clone_order=int(yield_.loc[d, "cases_with_order"]),
                             verdict_sim="not fitted (co-occurrence only: fewer than 200 cases with clone order)"))

pk = pd.DataFrame(pk_rows); pk.to_csv(TABLES / "validation_prior_knowledge.csv", index=False)
co = pd.DataFrame(co_rows); co.to_csv(TABLES / "validation_clone_order.csv", index=False)
sm = pd.DataFrame(sum_rows); sm.to_csv(TABLES / "validation_summary.csv", index=False)
pd.set_option("display.width", 250)
print(pk[["disease", "kind", "a", "b", "baseline", "clone_a_before_b", "clone_b_before_a"] + [c for c in pk.columns if c.endswith("pass")]].to_string())
print(co.round(3).to_string())
print(sm.round(3).T.to_string())
