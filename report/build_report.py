# Author: Shubhodeep Chatterjee
"""Build the IEEE-format report (report/main.tex, report/MitelmanProgression_report.docx)
from the results tables. Compile main.tex with pdflatex (twice) for the PDF.
"""
import sys, json, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
from render import to_latex, to_docx
from karyoprog.common import TABLES, MODELS, FIGURES, ALL_DISEASES, DERIVED_EVENTS, ordering_diseases

order_d = ordering_diseases()
yld = pd.read_csv(TABLES / "ordering_yield.csv", index_col=0)
dq = pd.read_csv(TABLES / "data_quality.csv", index_col=0)
power = json.loads((ROOT / "config/power_threshold.json").read_text())["diseases"]
frozen = json.loads((ROOT / "config/vocabulary_frozen.json").read_text())
vs = pd.read_csv(TABLES / "validation_summary.csv").set_index("disease")
pk = pd.read_csv(TABLES / "validation_prior_knowledge.csv")
co = pd.read_csv(TABLES / "validation_clone_order.csv").set_index("disease")
fs = pd.read_csv(TABLES / "frequencies_by_series.csv")
cbn = {d: json.loads((MODELS / f"cbn_{d}.json").read_text()) for d in order_d}
mhn = {d: json.loads((MODELS / f"mhn_{d}_summary.json").read_text()) for d in order_d}
cooc = {d: pd.read_csv(TABLES / f"cooccurrence_{d}.csv") for d in ALL_DISEASES}
edges = {d: pd.read_csv(TABLES / f"mhn_edges_{d}.csv") for d in order_d}
f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:.3f}"
n_total = int(yld.cases.sum())
n_dup = int(dq.duplicates_flagged.sum())
n_order = int(yld.cases_with_order.sum())

# ------------------------------------------------------------------ tables
t_data = ("table", "tab:data", "Cohorts, vocabulary and ordering yield (non-duplicate cases). N: minimum event count from the power calculation. Clone order: cases with at least one parent-to-child event pair, of which Unsel. are from unselected series.",
          ["Disease", "Cases", "N", "Events", "Clone order", "Unsel.", "Models"],
          [[d, f"{int(yld.loc[d, 'cases']):,}", power[d]["N"], int(yld.loc[d, "events_frozen"]),
            f"{int(yld.loc[d, 'cases_with_order']):,}", f"{int(yld.loc[d, 'of_which_unselected']):,}",
            "all three" if d in order_d else "baseline"] for d in ALL_DISEASES],
          "lrrrrrl", False)

rows = []
for d in ALL_DISEASES:
    c = cooc[d]; s = c[c.significant]
    rows.append([d, len(c), len(s), int((s.direction == "co-occurring").sum()), int((s.direction == "exclusive").sum())])
t_cooc = ("table", "tab:cooc", "Pairwise baseline: pairs significant under both Fisher's test and the count-preserving null (BH q < 0.05, same direction).",
          ["Disease", "Pairs", "Significant", "Co-occurring", "Exclusive"], rows, "lrrrr", False)

rows = []
for d in order_d:
    v = vs.loc[d]
    rows.append([d, int(v.model_events), f2(v.lambda_c), f3(v.heldout_mhn), f3(v.heldout_independence),
                 int(v.mhn_strong_edges), int(v.mhn_stable_edges), f2(v.mhn_frac_stable),
                 f"{v.mhn_sim_precision:.2f}/{v.mhn_sim_recall:.2f}", f2(v.unselected_corr), f2(v.post2000_corr)])
t_mhn = ("table", "tab:mhn", "MHN fits and validation. Held-out values are mean log-likelihood per case (2-fold). Stable: bootstrap sign agreement ≥ 0.8. Sim P/R: simulation-recovery precision and recall. U and 2000+: correlation of off-diagonal log θ between the full fit and the refits on unselected series and on cases published from 2000.",
         ["Disease", "Events", "c", "MHN", "Indep.", "Edges", "Stable", "Share", "Sim P/R", "U r", "2000+ r"],
         rows, "lrrrrrrrrrr", True)

rows = []
for d in order_d:
    v = vs.loc[d]; m = cbn[d]
    ev = m["events"]; R = np.asarray(m["P_reduced"]); C = np.asarray(m["closure_bootstrap_frequency"])
    es = "; ".join(f"{ev[a]} → {ev[b]} ({C[a, b]:.2f})" for a, b in zip(*np.nonzero(R))) or "none"
    rows.append([d, f3(v.cbn_eps), f3(v.heldout_cbn), f3(v.heldout_mhn_cbn_events), f3(v.heldout_indep_cbn_events),
                 f"{v.cbn_sim_precision:.2f}/{v.cbn_sim_recall:.2f}", es])
t_cbn = ("table", "tab:cbn", "CBN on the 10 most frequent model events per disease: noise ε, held-out log-likelihood per case against an MHN and an independence model on the same events, simulation recovery, and the cover relations of the fitted poset with their bootstrap frequency.",
         ["Disease", "ε", "CBN", "MHN", "Indep.", "Sim P/R", "Poset cover relations (bootstrap)"],
         rows, "lrrrrrp{8.2cm}", True)


def yesno(x):
    if pd.isna(x):
        return "n/a"
    return "yes" if (x is True or str(x) == "True") else "no"


rows = []
for r in pk.itertuples():
    rel = {"before": f"{r.a} before {r.b}", "exclusive": f"{r.a} / {r.b} exclusive", "cooccur": f"{r.a} + {r.b} co-occur"}[r.kind]
    clone = f"{int(r.clone_a_before_b)}:{int(r.clone_b_before_a)}" if r.kind == "before" else "-"
    rows.append([r.disease, rel, yesno(getattr(r, "baseline_pass", np.nan)) if r.kind != "before" else "-",
                 yesno(getattr(r, "mhn_pass", np.nan)), yesno(getattr(r, "cbn_pass", np.nan)) if r.kind == "before" else "-",
                 clone, yesno(getattr(r, "clone_pass", np.nan)) if r.kind == "before" else "-"])
t_pk = ("table", "tab:pk", "Prior-knowledge test set, written before the fits. Clone a:b counts parent-to-child clone pairs ordered a→b and b→a (non-duplicate cases); clone check needs ≥ 5 pairs. n/a: an event outside that model's vocabulary.",
        ["Disease", "Relationship", "Baseline", "MHN", "CBN", "Clone a:b", "Clone"], rows, "llccccc", True)

rows = [[d, int(co.loc[d, "pairs_compared"]), f"{int(co.loc[d, 'cases_in_pairs']):,}", f2(co.loc[d, "direction_agreement"]),
         f2(co.loc[d, "weighted_agreement"]), f2(co.loc[d, "spearman"]) if pd.notna(co.loc[d, "spearman"]) else "n/a"]
        for d in order_d]
t_clone = ("table", "tab:clone", "Out-of-sample check: MHN order probabilities (fitted on cross-sectional presence/absence) against the direction of event pairs in parent-to-child clone edges (pairs seen in ≥ 10 cases).",
           ["Disease", "Pairs", "Clone pairs", "Agree", "Weighted", "Spearman"], rows, "lrrrrr", False)

from scipy.stats import binomtest
pairs = pd.read_csv(TABLES / "ordered_pairs.tsv.gz", sep="\t")
pairs = pairs[~pairs.is_duplicate.astype(bool)]
rows = []
for d in order_d:
    cnt = pairs[pairs.disease == d].groupby(["earlier", "later"]).size()
    seen, k = set(), 0
    for (a, b), n_ab in cnt.sort_values(ascending=False).items():
        if frozenset((a, b)) in seen or {a, b} & DERIVED_EVENTS:
            continue
        seen.add(frozenset((a, b)))
        n_ba = int(cnt.get((b, a), 0))
        if n_ab < n_ba:
            continue
        pv = binomtest(int(n_ab), int(n_ab + n_ba), 0.5).pvalue
        rows.append([d if k == 0 else "", f"{a} → {b}", f"{int(n_ab)}:{n_ba}", f"{pv:.0e}" if pv < 0.001 else f"{pv:.3f}"])
        k += 1
        if k == 4:
            break
t_cloneorders = ("table", "tab:cloneorders", "Model-free orders: the event pairs most often split between a parent clone and its child (non-duplicate cases), with a two-sided sign test.",
                 ["Disease", "Order", "a→b : b→a", "p"], rows, "lllr", False)

# ------------------------------------------------------------------ text
exec((HERE / "report_text.py").read_text())     # defines META, BLOCKS, REFS from the variables above

out_tex = HERE / "main.tex"
(HERE / "figures").mkdir(exist_ok=True)
for b in BLOCKS:
    if b[0] == "figure":
        for ext in (".pdf", ".png"):
            shutil.copy(Path(b[2]).with_suffix(ext), HERE / "figures" / (Path(b[2]).stem + ext))
to_latex(META, BLOCKS, REFS, out_tex)
to_docx(META, BLOCKS, REFS, HERE / "MitelmanProgression_report.docx")
print("wrote", out_tex, "and the DOCX")
