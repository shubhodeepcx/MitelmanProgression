# Author: Shubhodeep Chatterjee
"""Loaders for the Mitelman bulk download (mitelman_db.zip, tab-separated *.TXT.DATA)."""
from pathlib import Path
import pandas as pd

DEFAULT_DIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "mitelman"
URL = "https://storage.googleapis.com/mitelman-data-files/prod/mitelman_db.zip"


def read_table(name, data_dir=DEFAULT_DIR):
    return pd.read_csv(Path(data_dir) / f"{name}.TXT.DATA", sep="\t", dtype=str,
                       keep_default_na=False, quoting=3, encoding="latin-1")


# Disease entities, defined on Mitelman MORPH codes (see KODER table).
DISEASES = {
    "AML": lambda m: m.startswith("11") and m != "1106",   # non-APL AML, all FAB types
    "APL": lambda m: m == "1106",
    "MDS": lambda m: m.startswith("15"),
    "MDS/MPN": lambda m: m.startswith("13"),
    "MPN": lambda m: m.startswith("14"),
    "CML": lambda m: m.startswith("12"),
    "ALL": lambda m: m.startswith("16"),
    "MM": lambda m: m == "1702",
    "CLL": lambda m: m == "1802",
    "DLBCL": lambda m: m == "1820",
    "FL": lambda m: m == "1810",
    "MCL": lambda m: m == "1808",
    "BL": lambda m: m == "1822",
}


def disease_of(morph):
    morph = (morph or "").strip()
    for name, pred in DISEASES.items():
        if pred(morph):
            return name
    return None


def load_cases(data_dir=DEFAULT_DIR):
    """One row per case: RefNo, CaseNo, disease, Series (S/U), Year, Morph, Topo."""
    cg = read_table("CYTOGEN", data_dir)
    ref = read_table("REF", data_dir)[["RefNo", "Year", "Pubmed"]]
    cg = cg.merge(ref, on="RefNo", how="left")
    cg["Morph"] = cg.Morph.str.strip()
    cg["disease"] = cg.Morph.map(disease_of)
    cg["Year"] = pd.to_numeric(cg.Year, errors="coerce")
    return cg[["RefNo", "CaseNo", "disease", "Morph", "Topo", "Series", "Year",
               "Pubmed", "Sex", "Age", "PrevTreat"]]


def load_clones(data_dir=DEFAULT_DIR):
    """One row per clone: RefNo, CaseNo, InvNo, CloneNo, clone (full ISCN text)."""
    kc = read_table("KCLONE", data_dir)
    kc["clone"] = kc.CloneLong.where(kc.CloneLong.str.strip() != "", kc.CloneShort).str.strip()
    kc["InvNo"] = kc.InvNo.astype(int)
    kc["CloneNo"] = kc.CloneNo.astype(int)
    return kc[["RefNo", "CaseNo", "InvNo", "CloneNo", "clone", "ChromoMin", "ChromoMax"]].sort_values(
        ["RefNo", "CaseNo", "InvNo", "CloneNo"]).reset_index(drop=True)
