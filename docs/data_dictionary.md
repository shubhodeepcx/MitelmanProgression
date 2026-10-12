# Released tables and model outputs

Shubhodeep Chatterjee, MitelmanProgression.

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
| parent_clone | inferred parent clone (the largest abnormal clone whose aberrations are a proper subset; ties go to the earlier-listed clone); empty for founders |
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
column per frozen event, from the case's first investigation with an abnormal clone, an
event counting if present in any abnormal clone. Cases with none of the events are kept
as all-zero rows. CBN and MHN drop the count-derived columns (complex, high hyperdiploid,
hypodiploid).

## ordering_yield.csv
Per disease: non-duplicate cases, frozen event count, cases with at least one ordered pair,
how many of those are in unselected series, and whether the 200-case stop rule is met.

## vocabulary_selection.csv
Every candidate event with its count, percentage, N and whether it was selected.

## frequencies_by_series.csv, frequencies_by_era.csv, data_quality.csv
Selection bias and data quality (scripts/05_descriptive.py): each frozen event's percentage
in all, Unselected (U) and Selected (S) series with the S/U ratio; percentages by
publication era; and per disease the cases, flagged duplicates, investigations with an
abnormal clone, the share CytoConverter could convert, and cases with several investigations.

## cooccurrence_<disease>.csv
Pairwise baseline (scripts/06_cooccurrence.py). Columns: event_a, event_b, n_a, n_b, n_both,
expected_null (mean co-occurrence under the count-preserving permutation null), odds_ratio
(with a 0.5 continuity correction), p_fisher, p_perm (1,000 permutations), q_fisher and
q_perm (Benjamini-Hochberg), direction (relative to the permutation null), significant
(both q < 0.05 and both tests agree on the direction).

## mhn_edges_<disease>.csv
Every MHN interaction (scripts/09_validation.py): source (acting event), target (affected
event), log_theta, sign_stability (share of 25 bootstrap refits with the same sign), strong
(|log theta| >= 0.5).

## validation_prior_knowledge.csv, validation_clone_order.csv, validation_summary.csv
Validation (scripts/09_validation.py): each relationship in config/prior_knowledge.json
checked against the baseline, the MHN, the CBN and the within-patient clone orders; how
often the MHN's simulated order agrees with the direction seen in parent-to-child clone
pairs; and one row per disease with held-out log-likelihoods, simulation recovery,
bootstrap stability, sensitivity refits and the stop-rule verdicts.

# Model outputs (results/models/)

- `cbn_<disease>.json`: events, poset P and its transitive reduction, theta, epsilon,
  log-likelihood, BIC, bootstrap frequency of each implied order (closure), simulation
  recovery and held-out log-likelihoods (CBN, MHN on the same events, independence).
- `mhn_<disease>.npz`: events, log theta, 25 bootstrap log thetas, the Unselected-only and
  post-2000 refits, order probabilities P(i before j | both) and their bootstrap copies.
- `mhn_<disease>_summary.json`: lambda choice, held-out log-likelihoods, simulation recovery
  and sensitivity agreement.

# Figures (results/figures/)

One figure per disease (`<disease>.png` and `.pdf`) and `overview.png`; see
scripts/10_figures.py for the panels.

