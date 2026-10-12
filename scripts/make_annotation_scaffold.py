# Author: Shubhodeep Chatterjee
"""Draw the 50-karyotype hand-annotation set (stratified, fixed seed).

Writes tests/data/hand_curated_50.tsv with empty expected_* columns, and refuses to
overwrite it once it exists so annotations cannot be lost. Run after
scripts/01_build_clone_features.py.
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd

OUT = ROOT / "tests/data/hand_curated_50.tsv"
if OUT.exists():
    sys.exit(f"{OUT} exists; delete it deliberately to redraw the sample.")
SEED = 20261011
cf = pd.read_csv(ROOT / "data/interim/clone_features.tsv.gz", sep="\t", keep_default_na=False)
cases = pd.read_csv(ROOT / "data/interim/cases.tsv.gz", sep="\t").set_index("case_id")
inv = cf.groupby(["case_id", "InvNo"]).agg(n_abn=("clone_id", "size"), edges=("parent_clone", lambda s: (s != "").sum()),
                                           complex=("features", lambda s: s.str.contains("complex").any()),
                                           karyotype=("clone", "/".join)).reset_index()
inv["disease"] = inv.case_id.map(cases.disease)
k = inv.karyotype
strata = {
    "AML single clone": ((inv.disease == "AML") & (inv.n_abn == 1) & ~inv.complex, 6),
    "AML complex single clone": ((inv.disease == "AML") & (inv.n_abn == 1) & inv.complex, 4),
    "AML idem sideline": ((inv.disease == "AML") & k.str.contains("idem") & (inv.n_abn == 2), 8),
    "AML idem, >=3 clones": ((inv.disease == "AML") & k.str.contains("idem") & (inv.n_abn >= 3), 6),
    "AML related, no idem": ((inv.disease == "AML") & ~k.str.contains("idem") & (inv.edges > 0), 5),
    "AML unrelated clones": ((inv.disease == "AML") & (inv.n_abn >= 2) & (inv.edges == 0), 4),
    "AML polyploid idemx2": ((inv.disease == "AML") & k.str.contains("idemx"), 3),
    "constitutional c": (inv.disease.isin(["AML", "MDS", "ALL"]) & k.str.contains(r"\dc[,/x]|\dc$", regex=True), 3),
    "uncertain ?": ((inv.disease == "AML") & k.str.contains(r"\?") & (inv.n_abn >= 2), 3),
    "KMT2A variant": ((inv.disease == "AML") & k.str.contains(r";11\)\([^)]*;q23|t\(11;[^)]*\)\(q23", regex=True), 3),
    "MDS multi-clone": ((inv.disease == "MDS") & (inv.edges > 0), 3),
    "CML multi-clone": ((inv.disease == "CML") & (inv.edges > 0), 2),
}
picked, rows = set(), []
for name, (mask, n) in strata.items():
    s = inv[mask & ~inv.case_id.isin(picked)].sample(n=n, random_state=SEED)
    picked |= set(s.case_id)
    rows += [dict(item=len(rows) + 1, stratum=name, case_id=r.case_id, InvNo=r.InvNo, disease=r.disease,
                  karyotype=r.karyotype) for r in s.itertuples()]
df = pd.DataFrame(rows)
for col in ["expected_features_by_clone", "expected_parent_by_clone", "expected_ordered_pairs", "annotator", "notes"]:
    df[col] = ""
df["status"] = "todo"
df.to_csv(OUT, sep="\t", index=False)
print(df.stratum.value_counts().to_string())
