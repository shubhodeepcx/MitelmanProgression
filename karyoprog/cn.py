"""Arm-level copy-number calls per clone from CytoConverter segments.

CytoConverter (shipped in the Mitelman download as CYTOCONVERTED) reports each
clone's Gain/Loss segments on hg38 relative to a diploid genome, one row per
extra or missing copy. A tetraploid clone therefore shows a gain of every
chromosome. We rebuild an integer copy-number profile per clone, take the
clone's own baseline ploidy as the most common whole-autosome copy number, and
call arm and whole-chromosome changes relative to that baseline.

Clone numbers in CYTOCONVERTED are positions in CYTOVAL.CytoConvInput, which
lists clones in the same order as KCLONE (checked in tests/test_cn.py).
"""
from collections import Counter
import numpy as np
import pandas as pd
from . import io

ACROCENTRIC = {"13", "14", "15", "21", "22"}
AUTOSOMES = [str(i) for i in range(1, 23)]
ARM_MIN = 0.5      # fraction of an arm that must be changed to call it
WHOLE_MIN = 0.9    # fraction of a chromosome (q arm only for acrocentrics)


def arm_table(data_dir=io.DEFAULT_DIR):
    cb = io.read_table("CYTOBANDS_HG38", data_dir)
    cb["chrom"] = cb.chromosome.str.replace("chr", "", regex=False)
    cb["arm"] = np.where(cb.cytoband_name.str.contains("p"), "p", "q")
    cb[["hg38_start", "hg38_stop"]] = cb[["hg38_start", "hg38_stop"]].astype(int)
    arms = cb.groupby(["chrom", "arm"]).agg(start=("hg38_start", "min"), end=("hg38_stop", "max")).reset_index()
    return {(r.chrom, r.arm): (r.start, r.end) for r in arms.itertuples()}


def _profile(segs, lo, hi):
    """Piecewise-constant copy-number delta on [lo, hi) from (start, end, +1/-1) segments.
    Returns list of (length, delta)."""
    pts = {lo, hi}
    for s, e, _ in segs:
        if e > lo and s < hi:
            pts.add(max(s, lo)); pts.add(min(e, hi))
    pts = sorted(pts)
    out = []
    for a, b in zip(pts, pts[1:]):
        d = sum(v for s, e, v in segs if s <= a and e >= b)
        out.append((b - a, d))
    return out


def _frac(profile, pred):
    tot = sum(l for l, _ in profile)
    return sum(l for l, d in profile if pred(d)) / tot if tot else 0.0


def call_clone(segs_by_chrom, arms):
    """segs_by_chrom: {chrom: [(start, end, +1|-1), ...]} for one clone.
    Returns (baseline_ploidy, sorted list of atomic CN features)."""
    prof, whole_cn = {}, []
    for c in AUTOSOMES:
        segs = segs_by_chrom.get(c, [])
        for arm in ("p", "q"):
            if c in ACROCENTRIC and arm == "p":
                continue
            lo, hi = arms[(c, arm)]
            prof[(c, arm)] = _profile(segs, lo, hi)
        q = prof[(c, "q")]
        p = prof.get((c, "p"), [])
        both = p + q
        # modal integer copy number of the chromosome (length-weighted)
        w = Counter()
        for l, d in both:
            w[2 + d] += l
        whole_cn.append(w.most_common(1)[0][0])
    base = Counter(whole_cn).most_common(1)[0][0]
    rel = base - 2   # delta that corresponds to "no change" for this clone
    feats = []
    for c in AUTOSOMES:
        armsc = ["q"] if c in ACROCENTRIC else ["p", "q"]
        allp = sum((prof[(c, a)] for a in armsc), [])
        if _frac(allp, lambda d: d < rel) >= WHOLE_MIN:
            feats.append(f"-{c}")
            continue
        if _frac(allp, lambda d: d > rel) >= WHOLE_MIN:
            feats.append(f"+{c}")
            continue
        for a in armsc:
            if _frac(prof[(c, a)], lambda d: d < rel) >= ARM_MIN:
                feats.append(f"{c}{a}-")
            if _frac(prof[(c, a)], lambda d: d > rel) >= ARM_MIN:
                feats.append(f"{c}{a}+")
    return base, feats


def call_all(data_dir=io.DEFAULT_DIR):
    """DataFrame: RefNo, CaseNo, InvNo, CloneNo, cc_ploidy, cn_features ('|'-joined)."""
    arms = arm_table(data_dir)
    cc = io.read_table("CYTOCONVERTED", data_dir)
    cc = cc[cc.Chr.isin([f"chr{c}" for c in AUTOSOMES])]
    cc["chrom"] = cc.Chr.str[3:]
    cc["v"] = np.where(cc.Type == "Gain", 1, -1)
    cc[["Start", "End"]] = cc[["Start", "End"]].astype(int)
    rows = []
    for (ref, case, inv, clone), g in cc.groupby(["RefNo", "CaseNo", "InvNo", "Clone"], sort=False):
        sb = {}
        for ch, s, e, v in zip(g.chrom, g.Start, g.End, g.v):
            sb.setdefault(ch, []).append((s, e, v))
        base, feats = call_clone(sb, arms)
        rows.append((ref, case, int(inv), int(clone), base, "|".join(feats)))
    return pd.DataFrame(rows, columns=["RefNo", "CaseNo", "InvNo", "CloneNo", "cc_ploidy", "cn_features"])
