"""How much of Mitelman did CytoConverter (shipped in the download) convert, and do its clones line up with KCLONE?"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import pandas as pd
from karyoprog import io

cc = io.read_table("CYTOCONVERTED"); log = io.read_table("CYTOCONVERTED_LOG"); cv = io.read_table("CYTOVAL")
kc = io.load_clones(); cases = io.load_cases()
k = ["RefNo", "CaseNo", "InvNo"]
for d in (cc, log, cv): d["InvNo"] = d.InvNo.astype(int)
n_inv = kc.groupby(k).CloneNo.nunique().rename("kclone_n")
cc_n = cc.groupby(k).Clone.nunique().rename("cc_n")
err = log[log.Type == "Error"].groupby(k).size().rename("errors")
inval = cv.set_index(k).Message.eq("Invalid").rename("invalid")
df = pd.concat([n_inv, cc_n, err, inval], axis=1).reset_index().merge(cases[["RefNo", "CaseNo", "disease"]], on=["RefNo", "CaseNo"], how="left")
df["converted"] = df.cc_n.notna()  # has >=1 Gain/Loss segment; balanced-only karyotypes have none
df["invalid"] = df.invalid.fillna(False).astype(bool)
df["clone_count_match"] = df.cc_n == df.kclone_n
print("investigations:", len(df))
print("message types:", log.Type.value_counts().to_dict())
print("CYTOVAL Message:", cv.Message.value_counts().head(5).to_dict())
summ = df.groupby("disease").agg(inv=("converted", "size"), iscn_invalid=("invalid", "mean"), has_cn_segments=("converted", "mean"),
                                 clone_count_match_given_converted=("clone_count_match", lambda s: s[df.loc[s.index, "converted"]].mean()))
summ.loc["ALL_DISEASES"] = [len(df), df.invalid.mean(), df.converted.mean(), df[df.converted].clone_count_match.mean()]
print(summ.round(3).to_string())
summ.round(4).to_csv(ROOT / "results/week1/cytoconverter_coverage.csv")
print("segment types:", cc.Type.value_counts().to_dict())
