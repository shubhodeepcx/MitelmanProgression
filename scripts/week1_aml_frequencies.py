"""AML (non-APL) and APL+AML aberration frequencies at diagnosis (first investigation),
split by Mitelman's Selected/Unselected series flag and publication era.

Mitelman holds abnormal karyotypes only, so frequencies are 'among cytogenetically
abnormal cases'. To compare with a consecutive series, divide that series' rate by
its abnormal fraction (1 - normal-karyotype fraction).
"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog.events import RULES

OUT = ROOT / "results" / "week1"
ev = pd.read_csv(OUT / "event_table_provisional.tsv.gz", sep="\t")
inv = pd.read_csv(OUT / "investigations.tsv.gz", sep="\t")

def freq_table(diseases, label):
    base = inv[(inv.InvNo == 1) & inv.disease.isin(diseases)].copy()
    base["era"] = pd.cut(base.Year, [0, 1989, 2004, 3000], labels=["<1990", "1990-2004", ">=2005"])
    e = ev[(ev.InvNo == 1) & ev.disease.isin(diseases)]
    present = e.groupby(["case_id", "event"]).size().unstack(fill_value=0).gt(0)
    present = present.reindex(base.case_id).fillna(False).astype(bool)
    present["complex(>=3, ELN)"] = base.set_index("case_id").complex.astype(bool)
    present["any abnormal clone"] = base.set_index("case_id").n_abnormal.gt(0)
    g = {"all": base.case_id, "Unselected": base[base.Series == "U"].case_id,
         "Selected": base[base.Series == "S"].case_id}
    for era in ["<1990", "1990-2004", ">=2005"]:
        g[f"U {era}"] = base[(base.Series == "U") & (base.era == era)].case_id
    out = pd.DataFrame({k: present.loc[v].mean() * 100 for k, v in g.items()}).round(1)
    out.loc["n cases"] = [len(v) for v in g.values()]
    out["S/U ratio"] = (out.Selected / out.Unselected).round(2)
    out.to_csv(OUT / f"frequencies_{label}.csv")
    print(f"== {label}\n", out.sort_values("all", ascending=False).to_string())

freq_table(["AML"], "AML_nonAPL")
freq_table(["AML", "APL"], "AML_incl_APL")
