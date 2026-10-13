# Author: Shubhodeep Chatterjee
"""One figure per disease, with uncertainty (run after 06-09).

Ordering diseases (>= 200 cases with clone order): four panels
  a  MHN log-theta matrix (row = affected event, column = acting event); a dot marks
     interactions that keep their sign in >= 80% of bootstrap refits
     (MHN log theta is relative to the association that a shared observation time alone
     induces, so a mildly positive raw odds ratio can carry a negative log theta)
  b  progression map: the 25 strongest stable promoting MHN edges (log theta >= 0.5) width = log theta, label = bootstrap sign stability
  c  co-occurrence baseline: log odds ratio, * = significant after BH and permutation null
  d  for the event pairs with the most within-patient clone evidence: the MHN's P(a before b |
     both), with bootstrap 10-90% interval, next to the share of clone pairs ordered a -> b
Co-occurrence-only diseases (APL, MM): panel c alone.
Writes results/figures/<disease>.png and .pdf, and results/figures/overview.png.
"""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import networkx as nx
from karyoprog.common import TABLES, MODELS, FIGURES, ALL_DISEASES, ordering_diseases

STRONG, STABLE = 0.5, 0.8
MAX_EDGES = 25          # progression map shows the strongest stable edges only
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "font.family": "DejaVu Sans", "savefig.dpi": 200})
pairs = pd.read_csv(TABLES / "ordered_pairs.tsv.gz", sep="\t")
pairs = pairs[~pairs.is_duplicate.astype(bool)]


def theta_panel(ax, ev, T, B):
    n = len(ev)
    M = T.copy(); np.fill_diagonal(M, np.nan)
    lim = max(np.nanmax(np.abs(M)), 0.5)
    im = ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
    stab = np.mean(np.sign(B) == np.sign(T)[None], axis=0)
    for i in range(n):
        for j in range(n):
            if i != j and abs(T[i, j]) >= STRONG and stab[i, j] >= STABLE:
                ax.plot(j, i, "o", ms=2.5, color="black")
    ax.set_xticks(range(n), ev, rotation=90); ax.set_yticks(range(n), ev)
    ax.set_xlabel("acting event"); ax.set_ylabel("affected event")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log theta")
    return stab


def graph_panel(ax, ev, T, stab, base_rate, max_edges=MAX_EDGES):
    G = nx.DiGraph()
    G.add_nodes_from(ev)
    cand = [(abs(T[i, j]), i, j) for i in range(len(ev)) for j in range(len(ev))
            if i != j and T[i, j] >= STRONG and stab[i, j] >= STABLE]
    for _, i, j in sorted(cand, reverse=True)[:max_edges]:
        G.add_edge(ev[j], ev[i], w=T[i, j], s=stab[i, j])
    iso = [v for v in G if G.degree(v) == 0]
    G.remove_nodes_from(iso)
    if G.number_of_edges() == 0:
        ax.text(0.5, 0.5, "no stable promoting edges", ha="center", va="center"); ax.axis("off"); return
    pos = nx.circular_layout(G)                      # circular: no overlapping nodes
    size = [250 + 1000 * base_rate[v] for v in G]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=size, node_color="#e8e8e8", edgecolors="#555")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
    for u, v, d in G.edges(data=True):
        nx.draw_networkx_edges(G, pos, [(u, v)], ax=ax, width=0.6 + 0.8 * min(abs(d["w"]), 4),
                               edge_color="#1b7837" if d["w"] > 0 else "#b2182b", arrowsize=10,
                               connectionstyle="arc3,rad=0.12", node_size=size)
    nx.draw_networkx_edge_labels(G, pos, {(u, v): f"{d['s']:.2f}" for u, v, d in G.edges(data=True)},
                                 ax=ax, font_size=6, label_pos=0.4, connectionstyle="arc3,rad=0.12")
    ax.axis("off"); ax.margins(0.15)
    ax.text(1.0, -0.02, "promotions only (log theta >= 0.5); inhibitions are in panel a",
            transform=ax.transAxes, fontsize=6, ha="right", va="top")
    if iso:
        ax.text(0.0, -0.06, "not shown: " + ", ".join(iso), transform=ax.transAxes, fontsize=6, va="top", wrap=True)


def cooc_panel(ax, d):
    c = pd.read_csv(TABLES / f"cooccurrence_{d}.csv")
    ev = list(dict.fromkeys(list(c.event_a) + list(c.event_b)))
    idx = {e: i for i, e in enumerate(ev)}
    M = np.full((len(ev), len(ev)), np.nan)
    for r in c.itertuples():
        orr = np.clip(r.odds_ratio, 1e-2, 1e2) if np.isfinite(r.odds_ratio) and r.odds_ratio > 0 else (1e-2 if r.odds_ratio == 0 else 1e2)
        M[idx[r.event_a], idx[r.event_b]] = M[idx[r.event_b], idx[r.event_a]] = np.log(orr)
        if r.significant:
            ax.text(idx[r.event_b], idx[r.event_a], "*", ha="center", va="center", fontsize=6)
            ax.text(idx[r.event_a], idx[r.event_b], "*", ha="center", va="center", fontsize=6)
    lim = max(np.nanmax(np.abs(M)) if np.isfinite(M).any() else 1, 0.5)
    im = ax.imshow(M, cmap="PuOr_r", norm=TwoSlopeNorm(0, -lim, lim))
    ax.set_xticks(range(len(ev)), ev, rotation=90); ax.set_yticks(range(len(ev)), ev)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log odds ratio")


def order_panel(ax, d, ev, P, PB, n_show=12):
    """Event pairs with the most within-patient clone evidence: MHN order probability (bootstrap
    10-90%) against the share of parent->child clone pairs in that direction."""
    p = pairs[pairs.disease == d]
    cnt = p.groupby(["earlier", "later"]).size()
    idx = {e: i for i, e in enumerate(ev)}
    rows, seen = [], set()
    for (a, b), n_ab in cnt.sort_values(ascending=False).items():
        if a not in idx or b not in idx or frozenset((a, b)) in seen:
            continue
        seen.add(frozenset((a, b)))
        n_ba = int(cnt.get((b, a), 0))
        if n_ab < n_ba or n_ab + n_ba < 10 or not np.isfinite(P[idx[a], idx[b]]):
            continue
        lo, hi = np.nanpercentile(PB[:, idx[a], idx[b]], [10, 90])
        rows.append((f"{a} -> {b}", P[idx[a], idx[b]], lo, hi, int(n_ab), n_ba))
        if len(rows) == n_show:
            break
    if not rows:
        ax.text(0.5, 0.5, "no pair with >= 10 clone pairs", ha="center", va="center"); ax.axis("off"); return
    rows = rows[::-1]
    y = np.arange(len(rows))
    ax.errorbar([r[1] for r in rows], y, xerr=[[max(r[1] - r[2], 0) for r in rows], [max(r[3] - r[1], 0) for r in rows]],
                fmt="o", color="#2166ac", ms=4, capsize=2, label="MHN P(a before b) (bootstrap 10-90%)")
    ax.plot([r[4] / (r[4] + r[5]) for r in rows], y, "D", color="#f4a582", ms=4, mec="#b35806",
            label="share of clone pairs a -> b")
    for k, r in enumerate(rows):
        ax.text(1.02, k, f"{r[4]}:{r[5]}", va="center", fontsize=6, transform=ax.get_yaxis_transform())
    ax.axvline(0.5, color="grey", lw=0.6, ls="--")
    ax.set_yticks(y, [r[0] for r in rows]); ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("probability that a precedes b")
    ax.legend(loc="lower left", fontsize=6, frameon=True, framealpha=0.9)


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    order_d = ordering_diseases()
    for d in ALL_DISEASES:
        if d in order_d:
            z = np.load(MODELS / f"mhn_{d}.npz", allow_pickle=True)
            ev, T, B, P, PB = list(z["events"]), z["theta"], z["boot"], z["order_P"], z["order_P_boot"]
            m = pd.read_csv(TABLES / f"matrix_{d}.tsv.gz", sep="\t", index_col=0)
            base = {e: float(m[e].mean()) for e in ev}
            fig, axs = plt.subplots(2, 2, figsize=(11, 10))
            stab = theta_panel(axs[0, 0], ev, T, B); axs[0, 0].set_title("a  MHN interactions (dot: stable edge)")
            graph_panel(axs[0, 1], ev, T, stab, base); axs[0, 1].set_title(f"b  Strongest stable promotions (up to {MAX_EDGES}; label: bootstrap stability)")
            cooc_panel(axs[1, 0], d); axs[1, 0].set_title("c  Co-occurrence baseline (*: significant)")
            order_panel(axs[1, 1], d, ev, P, PB); axs[1, 1].set_title("d  MHN order vs within-patient clone order (right: a->b : b->a)")
            s = json.loads((MODELS / f"mhn_{d}_summary.json").read_text())
            fig.suptitle(f"{d}: {s['n_cases']:,} cases, {len(ev)} events in the MHN", fontsize=11, fontweight="bold")
        else:
            fig, ax = plt.subplots(figsize=(6, 5.5))
            cooc_panel(ax, d)
            ax.set_title(f"{d}: co-occurrence only (fewer than 200 cases with clone order)\n* significant after BH and permutation null")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(FIGURES / f"{d}.{ext}", bbox_inches="tight")
        plt.close(fig)
        print("wrote", d, flush=True)
    # overview: held-out log-likelihood gain over independence, and stability
    sm = pd.read_csv(TABLES / "validation_summary.csv").dropna(subset=["heldout_mhn"])
    fig, axs = plt.subplots(1, 2, figsize=(8, 3))
    x = np.arange(len(sm))
    axs[0].bar(x - 0.2, sm.heldout_mhn - sm.heldout_independence, 0.4, label="MHN", color="#2166ac")
    axs[0].bar(x + 0.2, sm.heldout_cbn - sm.heldout_indep_cbn_events, 0.4, label="CBN (10 events)", color="#f4a582")
    axs[0].set_xticks(x, sm.disease); axs[0].set_ylabel("held-out log-lik gain\nover independence (per case)")
    axs[0].legend(frameon=False, fontsize=7); axs[0].axhline(0, color="black", lw=0.5)
    axs[1].bar(x, sm.mhn_frac_stable, color="#4d4d4d"); axs[1].axhline(0.5, color="red", lw=0.8, ls="--")
    axs[1].set_xticks(x, sm.disease); axs[1].set_ylim(0, 1); axs[1].set_ylabel("share of MHN edges stable")
    fig.tight_layout(); fig.savefig(FIGURES / "overview.png", bbox_inches="tight"); fig.savefig(FIGURES / "overview.pdf", bbox_inches="tight")
    print("wrote overview")
