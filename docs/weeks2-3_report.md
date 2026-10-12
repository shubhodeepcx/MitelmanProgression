# Weeks 2–3: frozen vocabulary, parser, tidy table

Everything here reproduces with `./run_all.sh`.

## What changed since week 1
- **Copy number now uses CytoConverter.** `karyoprog/cn.py` rebuilds each clone's copy-number profile from the CytoConverter segments that ship with Mitelman, takes the clone's own baseline ploidy (so a tetraploid clone isn't read as "gain of everything"), and calls arm and whole-chromosome changes. Those calls are unioned with what the ISCN text states. Against the ISCN text, CytoConverter reproduces +8 in 95% of clones, +21 in 95% and -7 in 87%. Most disagreements are real: for example, +8 written alongside add(8)(p11)x2 is net 8p loss and 8q gain. Interstitial deletions such as del(7)(q22q36) often cover less than half the arm, which is why the ISCN text is kept in the union.
- **The power calculation fixed N before selection.** An event must appear in enough cases to give 80% power (exact Fisher, Bonferroni over 190 pairs) to detect both co-occurrence at OR 3 and exclusivity at OR 0.2 against a 10% partner event. That gives N ≈ 200 cases in every disease (`config/power_threshold.json`).
- **The vocabulary is frozen per disease.** Candidates were written down before any ordering was looked at (`config/vocabulary_candidates.json`), then filtered by N and capped at 20. Result: AML 20 events, ALL 20, MM 16, MDS 13, CML 9, CLL 6, APL 2. Full list with counts: [event_vocabulary.md](event_vocabulary.md).
- **APL is redefined** as t(15;17) in any myeloid case plus APL morphology, because 1.1% of "non-APL" AML carried t(15;17).
- **De-duplication:** 237 groups, 250 cases flagged as likely repeat reports of the same patient. A collision means the same rare karyotype (≤3 occurrences database-wide), sex, age and disease in different papers. Most are in AML (100), ALL (93) and CML (35). Flagged cases are excluded from the matrices and counts.
- **Tests:** 24 unit tests pass (parser, copy-number caller, feature rules, frozen-vocabulary consistency, CytoConverter clone alignment on a real case). The 50-karyotype hand-annotation test is still skipped until the expected columns are filled in.

## Go/no-go with the frozen vocabulary
| Disease | Cases | Events | Cases with ordering | of which unselected | ≥200 |
|---|---|---|---|---|---|
| AML | 20,723 | 20 | 1,471 | 569 | yes |
| APL | 2,371 | 2 | 61 | 29 | no |
| MDS | 5,734 | 13 | 408 | 114 | yes |
| CML | 4,819 | 9 | 239 | 145 | yes |
| ALL | 12,842 | 20 | 610 | 133 | yes |
| CLL | 3,469 | 6 | 227 | 113 | yes |
| MM | 2,096 | 16 | 58 | 37 | no |

The disease-specific vocabularies changed the picture from week 1: ALL rises from 270 to 610 and CLL from 16 to 227, so both now pass. APL (61) and MM (58) stay co-occurrence-only. In ALL, much of the ordering comes from the hyperdiploid clone gaining trisomies, which week 4 will need to treat as one process (see the limitations in event_vocabulary.md).

## Still open
1. Hand-annotate the 50 karyotypes in `tests/data/hand_curated_50.tsv` (human task).
2. Week 4: simulation recovery at these sample sizes and sparsities, using `matrix_<disease>.tsv.gz` and `ordered_pairs.tsv.gz` as the template.
