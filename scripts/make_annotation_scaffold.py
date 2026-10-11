"""Draw the 50-karyotype hand-annotation set (stratified, seed fixed) and write
tests/data/hand_curated_50.tsv with empty `expected_*` columns for a human to fill,
plus tests/data/hand_curated_50_parser_output.tsv (what the parser says today).

Fill the expected columns *without* looking at the parser output file.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog import io
from karyoprog.iscn import parse_investigation, relate, ordered_pairs
from karyoprog.events import events_of

SEED = 20261011
inv = pd.read_csv(ROOT / "results/week1/investigations.tsv.gz", sep="\t")
cl = io.load_clones(); cl["case_id"] = cl.RefNo + ":" + cl.CaseNo
texts = cl.groupby(["case_id", "InvNo"]).clone.apply(list)
inv["karyotype"] = [ "/".join(texts[(c, i)]) for c, i in zip(inv.case_id, inv.InvNo)]
k = inv.karyotype
strata = {  # name: (mask, n)
    "AML single clone": ((inv.disease == "AML") & (inv.n_clones == 1) & ~inv.complex, 6),
    "AML complex single clone": ((inv.disease == "AML") & (inv.n_clones == 1) & inv.complex, 4),
    "AML idem sideline": ((inv.disease == "AML") & inv.uses_idem & (inv.n_clones == 2), 8),
    "AML idem, >=3 clones": ((inv.disease == "AML") & inv.uses_idem & (inv.n_clones >= 3), 6),
    "AML related, no idem": ((inv.disease == "AML") & ~inv.uses_idem & (inv.n_edges > 0), 5),
    "AML unrelated clones": ((inv.disease == "AML") & (inv.n_abnormal >= 2) & (inv.n_edges == 0), 4),
    "AML polyploid idemx2": ((inv.disease == "AML") & k.str.contains("idemx"), 3),
    "constitutional c": (inv.disease.isin(["AML", "MDS", "ALL"]) & k.str.contains(r"\dc[,/x]|\dc$", regex=True), 3),
    "uncertain ?": ((inv.disease == "AML") & k.str.contains(r"\?") & (inv.n_clones >= 2), 3),
    "KMT2A variant": ((inv.disease == "AML") & k.str.contains(r";11\)\([^)]*;q23|t\(11;[^)]*\)\(q23", regex=True), 3),
    "MDS multi-clone": ((inv.disease == "MDS") & (inv.n_edges > 0), 3),
    "CML multi-clone": ((inv.disease == "CML") & (inv.n_edges > 0), 2),
}
picked, rows = set(), []
for name, (mask, n) in strata.items():
    pool = inv[mask & ~inv.case_id.isin(picked)]
    s = pool.sample(n=min(n, len(pool)), random_state=SEED)
    picked |= set(s.case_id)
    for _, r in s.iterrows():
        rows.append(dict(item=len(rows) + 1, stratum=name, case_id=r.case_id, InvNo=r.InvNo,
                         disease=r.disease, karyotype=r.karyotype))
blank = pd.DataFrame(rows)
for col in ["expected_events_by_clone", "expected_parent_by_clone", "expected_ordered_pairs",
            "annotator", "status", "notes"]:
    blank[col] = ""
blank["status"] = "todo"
out = ROOT / "tests/data"; out.mkdir(parents=True, exist_ok=True)
blank.to_csv(out / "hand_curated_50.tsv", sep="\t", index=False)

def fmt_parser(kary):
    c = parse_investigation(kary.split("/"))
    par, _ = relate(c)
    ev = "; ".join(f"{x.index}:{'|'.join(sorted({e for a in x.abn for e in events_of(a)}))}" for x in c)
    pa = "; ".join(f"{i}:{p if p is not None else '-'}" for i, p in par.items())
    pr = "; ".join(sorted(f"{a}<{b}" for a, b in ordered_pairs(c, key=events_of)))
    return ev, pa, pr
po = blank[["item", "case_id", "InvNo", "karyotype"]].copy()
po[["parser_events_by_clone", "parser_parent_by_clone", "parser_ordered_pairs"]] = [fmt_parser(x) for x in po.karyotype]
po.to_csv(out / "hand_curated_50_parser_output.tsv", sep="\t", index=False)
print(blank.stratum.value_counts().to_string(), "\ntotal", len(blank))
