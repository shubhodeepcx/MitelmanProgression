# Released tables (results/tables/)

Built by `./run_all.sh` from the Mitelman download using the frozen vocabulary in
`config/vocabulary_frozen.json`. `case_id` is `RefNo:CaseNo` from Mitelman.

## event_table.tsv.gz (the tidy event table)
One row per case, investigation, abnormal clone and event.

| Column | Meaning |
|---|---|
| case_id | Mitelman reference and case number |
| disease | AML (non-APL), APL (t(15;17) in a myeloid case or APL morphology), MDS, CML, ALL, CLL, MM |
| InvNo | investigation number (separate samples of the same patient, e.g. diagnosis and relapse) |
| clone_id | clone position in the karyotype string (1 = stemline) |
| parent_clone | inferred parent clone (largest abnormal clone whose aberrations are a proper subset); empty for founders |
| inferred_depth | number of ancestors in the inferred clone tree (0 = founder) |
| event | event name from the frozen vocabulary for that disease |
| Series | U = unselected consecutive series, S = selected for an unusual feature |
| Year | publication year |
| is_duplicate | likely the same patient as a case in an earlier paper (see karyoprog/dedup.py); exclude for modelling |

## ordered_pairs.tsv.gz
One row per case and (earlier, later) event pair, from parent-to-child clone edges: every
event in the parent precedes every event the child adds. Columns: case_id, InvNo,
clone_id (the child), earlier, later, disease, Series, is_duplicate.

## matrix_<disease>.tsv.gz
Model input for co-occurrence, CBN and MHN: one row per non-duplicate case, one 0/1
column per frozen event, from the first investigation, an event counting if present in
any abnormal clone. Cases with none of the events are kept as all-zero rows.

## ordering_yield.csv
Per disease: non-duplicate cases, frozen event count, cases with at least one ordered pair,
how many of those are in unselected series, and whether the 200-case stop rule is met.

## vocabulary_selection.csv
Every candidate event with its count, percentage, N and whether it was selected.
