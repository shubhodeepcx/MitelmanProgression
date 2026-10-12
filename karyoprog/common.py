# Author: Shubhodeep Chatterjee
"""Shared paths, disease lists and loaders for the analysis scripts."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
MODELS = ROOT / "results" / "models"
FIGURES = ROOT / "results" / "figures"
INTERIM = ROOT / "data" / "interim"

ALL_DISEASES = ["AML", "APL", "MDS", "CML", "ALL", "CLL", "MM"]


def ordering_diseases():
    """Diseases that pass the 200-case stop rule (results/tables/ordering_yield.csv)."""
    y = pd.read_csv(TABLES / "ordering_yield.csv", index_col=0)
    return [d for d in y.index if bool(y.loc[d, "meets_200_rule"])]


# Events defined by counting other aberrations. They stay in the vocabulary and the
# co-occurrence baseline but are left out of progression models, where they would
# duplicate the information in their components and invite spurious precedence
# (e.g. "complex before 5q loss").
DERIVED_EVENTS = {"complex", "high hyperdiploid", "hypodiploid"}


def load_matrix(disease, for_model=False):
    """Case x event 0/1 matrix plus aligned case covariates (Series, Year).
    for_model=True drops the count-derived events in DERIVED_EVENTS."""
    m = pd.read_csv(TABLES / f"matrix_{disease}.tsv.gz", sep="\t", index_col=0)
    if for_model:
        m = m[[c for c in m.columns if c not in DERIVED_EVENTS]]
    cases = pd.read_csv(INTERIM / "cases.tsv.gz", sep="\t").set_index("case_id")
    cov = cases.loc[m.index, ["Series", "Year"]]
    return m, cov


def save_json(obj, path):
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(obj, indent=2, default=conv))
