# Week 1 report: data access, multi-clone yield, go/no-go

Data: Mitelman Database bulk download `mitelman_db.zip` (release files dated 2026-07-16),
fetched from `storage.googleapis.com/mitelman-data-files/prod/`. 81,148 cases,
84,519 investigations, 105,867 clones. Everything below is reproducible with `./run_week1.sh`.

## 1. Decision: ordering is a go for AML; marginal for MDS, ALL and CML; co-occurrence only for the rest

The stop rule is "fewer than ~200 cases per disease with usable multi-clone information".
"Usable" here means a case has at least one stemline-to-sideline (or parent-to-subclone)
relationship in which the parent and the subclone both carry an event from the provisional
vocabulary and the subclone adds one, so at least one (earlier, later) event pair is implied.

| Disease | Cases | ≥2 clones | Related clones (any aberration) | **Usable event-level order** | of which Unselected series | Extra from multi-investigation cases |
|---|---|---|---|---|---|---|
| AML (non-APL) | 21,061 | 4,110 | 3,249 | **1,219** | 466 | +123 |
| MDS | 5,741 | 1,171 | 841 | **276** | 83 | +32 |
| ALL | 12,935 | 2,104 | 1,726 | **270** | 65 | +42 |
| CML | 4,854 | 753 | 525 | **217** | 127 | +60 |
| APL | 2,134 | 267 | 215 | 98 | 41 | +7 |
| MM | 2,100 | 326 | 206 | 36 | 23 | +3 |
| CLL | 3,474 | 808 | 529 | 16 | 8 | +3 |
| DLBCL / FL / MCL / BL / MPN / MDS-MPN | | | | 3 to 22 each | | |

Source: `results/week1/multiclone_by_disease.csv`, `results/week1/longitudinal_by_disease.csv`.

Caveats on these numbers:
- The provisional vocabulary is myeloid-centred. ALL, CLL and MM counts are **underestimates**
  until each gets its own events (ALL: hyperdiploidy, t(12;21), KMT2A, iAMP21, del(9p);
  CLL: del(13q), +12, del(11q), del(17p); MM: t(11;14), t(4;14), 1q gain, del(13), hyperdiploidy).
  The "related clones (any aberration)" column is the ceiling.
- Recommendation: run AML as the primary ordering analysis. Treat MDS and CML as secondary
  ordering analyses, re-checking the count after the vocabulary freeze. APL, MM, CLL and the
  lymphomas go co-occurrence only unless the disease-specific vocabularies lift them past 200.
- APL defined by morphology code misses APL coded as "AML NOS": 1.1% of non-APL AML cases carry
  t(15;17). Weeks 2-3 should define APL as t(15;17)-positive regardless of code.

## 2. First-pass ordering signal agrees with known biology

Cases supporting each direction (parent clone carries the first event, subclone adds the second):

| Disease | Pair | Forward | Reverse |
|---|---|---|---|
| AML | t(8;21) before -Y | 76 | 2 |
| AML | t(8;21) before del(9q) | 49 | 1 |
| AML | inv(16) before +22 | 69 | 2 |
| AML | inv(16) before +8 | 44 | 1 |
| AML | KMT2A(11q23) before +8 | 61 | 4 |
| AML | -5/del(5q) before +8 | 72 | 8 |
| AML | -5/del(5q) before -17/17p- | 72 | 4 |
| AML | -5/del(5q) before -7 | 54 | 15 |
| AML | +8 before +21 | 36 | 19 |
| APL | t(15;17) before +8 | 50 | 3 |
| CML | t(9;22) before +8 | 123 | 10 |
| MDS | -5/del(5q) before +8 | 38 | 7 |

Founding translocations come before secondary gains and losses, and loss of 5q precedes +8 and
17p loss, in line with established biology. This is a raw count, not a model fit; it is a
sanity floor, not a finding.

## 3. CytoGPS

cytogps.org, CRAN and Bioconductor are blocked by this environment's network policy, and R is not
installed, so RCytoGPS could not be run here. The Mitelman download already ships
**CytoConverter** output (same group as CytoGPS; per-clone hg38 Gain/Loss segments, table
`CYTOCONVERTED`). Coverage: 7.1% of investigations are flagged invalid by Mitelman's syntax
checker; 75.6% have at least one copy-number segment (balanced-only karyotypes have none, which
is why APL is 36%). CytoConverter has no fusion/"F" dimension, so balanced rearrangements are
taken from the ISCN text by `karyoprog/events.py`. Source: `results/week1/cytoconverter_coverage.csv`.

To get true CytoGPS LGF vectors as a cross-check, RCytoGPS has to run on a machine with R and
internet access (e.g. your laptop).

## 4. AML aberration frequencies and the publication-bias handle

Mitelman stores abnormal karyotypes only, so frequencies are "among abnormal cases". Mitelman's
own `Series` flag splits cases into Selected (reported for an unusual feature) and Unselected
(consecutive series). The Unselected subset is a built-in unbiased comparator.

Non-APL AML, first investigation (`results/week1/frequencies_AML_nonAPL.csv`):

| Event | Unselected % (n=8,700) | Selected % (n=12,353) | S/U |
|---|---|---|---|
| complex (≥3, ELN) | 26.8 | 26.5 | 0.99 |
| +8 | 18.0 | 14.7 | 0.82 |
| -5/del(5q) | 14.8 | 14.0 | 0.95 |
| -7 | 12.8 | 13.5 | 1.05 |
| t(8;21) | 8.5 | 11.1 | 1.31 |
| KMT2A (11q23) | 6.8 | 9.6 | 1.41 |
| inv(16)/t(16;16) | 5.2 | 6.8 | 1.31 |
| inv(3)/t(3;3) | 1.7 | 4.3 | 2.53 |
| t(9;22) | 1.2 | 2.1 | 1.75 |

- Publication bias is measurable and has a direction: rare balanced rearrangements are
  over-represented 1.3 to 2.5 fold in Selected reports, while +8 is under-represented.
- Era effect: in Unselected cases, inv(16) rises from 3.2% (before 1990) to 5.7% (1990-2004) and
  KMT2A from 4.0% to 7.8%, consistent with improved detection of subtle rearrangements. Year
  stratification is needed, as the brief anticipated.
- Literature check (**not yet done against the source tables**): consecutive adult AML series
  are commonly quoted at roughly 7% t(8;21), 5% inv(16), 10-12% complex, with ~40-45% normal
  karyotypes; that's from memory, not checked. Rescaled to abnormal cases that is about 12%,
  9% and 18-20%. Mitelman Unselected is lower for t(8;21) and inv(16) and higher for complex
  karyotype, which fits an older and more secondary-AML-heavy population. The publisher and
  PubMed Central sites are blocked from here, so the exact Grimwade 2010 (MRC, n=5,876) or
  Byrd 2002 (CALGB 8461) tables still need to be pulled in.

## 5. Parser and hand-annotation set

- `karyoprog/iscn.py` splits each clone, expands `idem` against the stemline (Mitelman uses
  `idem` only; no `sl`/`sdl`), drops constitutional (`c`) changes and markers, builds a clone
  tree by set containment, and emits (earlier, later) pairs from parent-to-child edges only.
  Clones that overlap without containment (e.g. sharing only -Y) are not used by default.
- Tidy table: `results/week1/event_table_provisional.tsv.gz` (case_id, InvNo, clone_id, event,
  inferred_depth, parent_clone, disease, Series, Year), plus `ordered_pairs_provisional.tsv.gz`.
- 50-karyotype set drawn over 12 strata (fixed seed): `tests/data/hand_curated_50.tsv`. The
  `expected_*` columns are **empty and need a human annotator**; fill them without opening
  `hand_curated_50_parser_output.tsv`, set `status=done`, and `tests/test_hand_curated.py`
  asserts the parser matches. 11 hand-written unit tests pass today.

## 6. Open items for weeks 2-3

1. Hand-annotate the 50 (human task).
2. Freeze the vocabulary per disease; add CytoConverter-derived arm-level CN events so losses
   through unbalanced translocations (e.g. der(5)t(5;17)) are counted.
3. De-duplication: the PatientID field is deprecated, so collisions must be matched on
   karyotype, age, sex and reference overlap.
4. Pull the published consecutive-series frequency tables for the external comparison.
5. Power calculation for the per-disease event-frequency threshold N, before modelling.
