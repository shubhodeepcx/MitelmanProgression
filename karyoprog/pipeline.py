# Author: Shubhodeep Chatterjee
"""Mitelman download -> clone feature table -> frozen vocabulary -> tidy event table.

Entry points used by scripts/01_build_clone_features.py and scripts/04_build_tables.py:
  clone_feature_table()  one row per abnormal clone with atomic features, tree position
  case_table()           one row per case with final disease group and duplicate flag
  tidy_event_table()     the released table: case_id, disease, event, clone_id, inferred_depth, ...
"""
from pathlib import Path
import json
import pandas as pd
from . import io, cn
from .iscn import parse_investigation, relate
from .features import clone_features
from .dedup import flag_duplicates

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"


def cn_calls(data_dir=io.DEFAULT_DIR, cache=True):
    path = INTERIM / "cn_calls.tsv.gz"
    if cache and path.exists():
        return pd.read_csv(path, sep="\t", dtype={"RefNo": str, "CaseNo": str}, keep_default_na=False)
    d = cn.call_all(data_dir)
    INTERIM.mkdir(parents=True, exist_ok=True)
    d.to_csv(path, sep="\t", index=False)
    return d


def clone_feature_table(data_dir=io.DEFAULT_DIR):
    clones = io.load_clones(data_dir)
    cv = io.read_table("CYTOVAL", data_dir)
    cv["InvNo"] = cv.InvNo.astype(int)
    # CytoConverter clone numbers are positions in CytoConvInput. Use its calls only when the
    # syntax checker accepted the investigation AND kept the same number of clones as KCLONE;
    # otherwise positions may be shifted and we fall back to the ISCN text alone.
    n_kclone = clones.groupby(["RefNo", "CaseNo", "InvNo"]).size().to_dict()
    n_input = {(r, c, i): len([x for x in s.split("/") if x.strip()])
               for r, c, i, s in zip(cv.RefNo, cv.CaseNo, cv.InvNo, cv.CytoConvInput)}
    valid = {(r, c, i) for r, c, i, m in zip(cv.RefNo, cv.CaseNo, cv.InvNo, cv.Message)
             if m != "Invalid" and n_input.get((r, c, i)) == n_kclone.get((r, c, i))}
    calls = cn_calls(data_dir)
    callmap = {(r, c, int(i), int(k)): (f.split("|") if f else [])
               for r, c, i, k, f in zip(calls.RefNo, calls.CaseNo, calls.InvNo, calls.CloneNo, calls.cn_features)}
    rows = []
    for (ref, case, inv), g in clones.groupby(["RefNo", "CaseNo", "InvNo"], sort=False):
        cl = parse_investigation(list(g.clone))
        parent, depth = relate(cl)
        is_valid = (ref, case, inv) in valid
        for c, mn, mx in zip(cl, g.ChromoMin, g.ChromoMax):
            if c.is_normal:
                continue
            calls_c = callmap.get((ref, case, inv, c.index), []) if is_valid else None
            f = clone_features(c.abn, calls_c, (mn, mx))
            rows.append(dict(RefNo=ref, CaseNo=case, InvNo=inv, clone_id=c.index,
                             parent_clone=parent[c.index], inferred_depth=depth[c.index],
                             n_aberrations=len(c.abn), cn_source="iscn+cytoconverter" if is_valid else "iscn",
                             features="|".join(sorted(f)), clone=c.text))
    out = pd.DataFrame(rows)
    out.insert(0, "case_id", out.RefNo + ":" + out.CaseNo)
    return out


def case_table(cf, data_dir=io.DEFAULT_DIR):
    """Final disease group per case (APL = t(15;17) in a myeloid case, or APL morphology) and duplicate flags."""
    cases = io.load_cases(data_dir)
    cases["case_id"] = cases.RefNo + ":" + cases.CaseNo
    first = cf[cf.InvNo == cf.groupby("case_id").InvNo.transform("min")]
    t1517 = set(first[first.features.str.contains(r"(?:^|\|)t\(15;17\)(?:\||$)", regex=True)].case_id)
    myeloid_aml = cases.Morph.str.startswith("11")
    cases.loc[myeloid_aml & cases.case_id.isin(t1517), "disease"] = "APL"
    clones = io.load_clones(data_dir)
    clones["case_id"] = clones.RefNo + ":" + clones.CaseNo
    kfull = clones.groupby("case_id").clone.apply(lambda s: "/".join(s)).rename("karyotype_full")
    cases = cases.merge(kfull, on="case_id", how="left")
    cases["karyotype_full"] = cases.karyotype_full.fillna("")
    cases = cases.merge(flag_duplicates(cases), on="case_id", how="left")
    return cases


def load_vocabulary(path=ROOT / "config" / "vocabulary_frozen.json"):
    return json.loads(Path(path).read_text())


def events_for(features, disease_vocab):
    fs = set(features.split("|")) if features else set()
    return [e for e, anyof in disease_vocab.items() if fs & set(anyof)]


def tidy_event_table(cf, cases, vocab):
    """case_id, disease, InvNo, clone_id, parent_clone, inferred_depth, event, + case covariates."""
    meta = cases.set_index("case_id")[["disease", "Series", "Year", "is_duplicate"]]
    cf = cf.join(meta, on="case_id")
    rows = []
    for r in cf[cf.disease.isin(vocab["diseases"].keys())].itertuples():
        dv = vocab["diseases"][r.disease]["events"]
        for e in events_for(r.features, dv):
            rows.append((r.case_id, r.disease, r.InvNo, r.clone_id, r.parent_clone, r.inferred_depth,
                         e, r.Series, r.Year, r.is_duplicate))
    return pd.DataFrame(rows, columns=["case_id", "disease", "InvNo", "clone_id", "parent_clone",
                                       "inferred_depth", "event", "Series", "Year", "is_duplicate"])
