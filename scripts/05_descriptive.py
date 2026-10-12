# Author: Shubhodeep Chatterjee
"""Descriptive tables: selection bias (Selected vs Unselected series), era, data quality.

results/tables/
  frequencies_by_series.csv   per disease and frozen event: % in all, Unselected, Selected, S/U ratio
  frequencies_by_era.csv      per disease and event: % in Unselected cases by publication era
  data_quality.csv            per disease: cases, duplicates, CytoConverter-usable investigations,
                              cases with several investigations
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog.common import TABLES, INTERIM, ALL_DISEASES, load_matrix

rows, era_rows = [], []
for d in ALL_DISEASES:
    m, cov = load_matrix(d)
    era = pd.cut(cov.Year, [0, 1989, 1999, 2009, 3000], labels=["<1990", "1990-1999", "2000-2009", ">=2010"])
    for e in m.columns:
        x = m[e]
        u, s = x[cov.Series == "U"], x[cov.Series == "S"]
        rows.append(dict(disease=d, event=e, n_all=len(x), pct_all=100 * x.mean(), n_unselected=len(u),
                         pct_unselected=100 * u.mean(), n_selected=len(s), pct_selected=100 * s.mean(),
                         selected_to_unselected=(s.mean() / u.mean()) if u.mean() > 0 else float("nan")))
        for lab in era.cat.categories:
            sel = (era == lab) & (cov.Series == "U")
            era_rows.append(dict(disease=d, event=e, era=lab, n=int(sel.sum()), pct=100 * x[sel].mean() if sel.any() else float("nan")))
pd.DataFrame(rows).round(3).to_csv(TABLES / "frequencies_by_series.csv", index=False)
pd.DataFrame(era_rows).round(3).to_csv(TABLES / "frequencies_by_era.csv", index=False)

cases = pd.read_csv(INTERIM / "cases.tsv.gz", sep="\t")
cf = pd.read_csv(INTERIM / "clone_features.tsv.gz", sep="\t", keep_default_na=False, usecols=["case_id", "InvNo", "cn_source"])
inv = cf.drop_duplicates(["case_id", "InvNo"]).merge(cases[["case_id", "disease"]], on="case_id")
q = pd.DataFrame({
    "cases": cases.groupby("disease").size(),
    "duplicates_flagged": cases.groupby("disease").is_duplicate.sum(),
    "investigations_with_abnormal_clone": inv.groupby("disease").size(),
    "pct_with_cytoconverter": 100 * inv.assign(c=inv.cn_source.eq("iscn+cytoconverter")).groupby("disease").c.mean(),
    "cases_with_several_investigations": inv.groupby(["disease", "case_id"]).size().gt(1).groupby("disease").sum(),
}).loc[ALL_DISEASES]
q.round(1).to_csv(TABLES / "data_quality.csv")
print(q.round(1).to_string())
s = pd.DataFrame(rows)
print(s[s.disease == "AML"].sort_values("pct_all", ascending=False)[["event", "pct_unselected", "pct_selected", "selected_to_unselected"]].round(2).head(12).to_string(index=False))
