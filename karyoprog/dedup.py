"""Flag cases that are probably the same patient reported in more than one paper.

Mitelman has no usable patient identifier (PatientID is deprecated), so a
collision is called when two cases from *different* references share disease
group, sex, age (both present) and the exact full karyotype, and that
karyotype is rare enough to be informative (seen in at most `max_string_freq`
cases database-wide). Common simple karyotypes such as 46,XX,t(15;17)(q24;q21)
are excluded because identical strings there are expected by chance.
The earliest reference in each group is kept; the rest get is_duplicate=True.
"""
import pandas as pd


def flag_duplicates(case_tbl, max_string_freq=3):
    """case_tbl: one row per case with case_id, RefNo, disease, Sex, Age, Year, karyotype_full."""
    d = case_tbl.copy()
    d["Age"] = pd.to_numeric(d.Age, errors="coerce")
    d["Sex"] = d.Sex.str.strip()
    freq = d.karyotype_full.map(d.karyotype_full.value_counts())
    ok = d.Age.notna() & d.Sex.isin(["M", "F"]) & (d.karyotype_full.str.len() > 0) & (freq <= max_string_freq)
    d["dup_group"] = pd.NA
    d["is_duplicate"] = False
    key = ["disease", "Sex", "Age", "karyotype_full"]
    cand = d[ok]
    grp = cand.groupby(key)
    for gi, (_, g) in enumerate(grp):
        if g.RefNo.nunique() < 2:
            continue
        g = g.sort_values(["Year", "RefNo"])
        d.loc[g.index, "dup_group"] = gi
        d.loc[g.index[1:], "is_duplicate"] = True
    return d[["case_id", "dup_group", "is_duplicate"]]
