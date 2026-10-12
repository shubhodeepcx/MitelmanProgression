"""Step 1: clone feature table and case table (writes data/interim/, results/tables/)."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from karyoprog.pipeline import clone_feature_table, case_table, INTERIM

cf = clone_feature_table()
cases = case_table(cf)
INTERIM.mkdir(parents=True, exist_ok=True)
cf.to_csv(INTERIM / "clone_features.tsv.gz", sep="\t", index=False)
cases.drop(columns=["karyotype_full"]).to_csv(INTERIM / "cases.tsv.gz", sep="\t", index=False)
print("abnormal clones:", len(cf), " cases:", len(cases))
print("CN source:", cf.cn_source.value_counts().to_dict())
d = cases[cases.dup_group.notna()]
print("duplicate groups:", d.dup_group.nunique(), " cases flagged as duplicates:", int(cases.is_duplicate.sum()))
print(cases.groupby("disease").is_duplicate.sum().astype(int).to_dict())
