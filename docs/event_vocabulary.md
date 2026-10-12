# Frozen event vocabulary (v1)

Shubhodeep Chatterjee, MitelmanProgression.

Frozen on 2026-10-12 against the Mitelman download with md5 `434cf3759ba633e4d9d99e2bb8208283`.
Rule: an event enters a disease's vocabulary if it is present (any abnormal clone, first investigation, duplicates removed) in at least N cases, where N comes from the power calculation in `config/power_threshold.json`; at most 20 events, ranked by count; defining events are always kept. Candidates were listed in `config/vocabulary_candidates.json` before any ordering result was looked at.

An event is present in a clone when any of its atomic features is (see `karyoprog/features.py`). Copy-number features are the union of what the ISCN text states and CytoConverter arm-level calls (at least 50% of an arm, at least 90% of a chromosome for whole gains and losses, relative to the clone's own baseline ploidy).

Count-derived events (`complex`, `high hyperdiploid`, `hypodiploid`) are kept in the co-occurrence baseline but excluded from CBN and MHN, because they are functions of the other events.

## AML (n = 20,753, N = 201)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| complex | complex | 5,605 | 27.01 |
| 5q loss | -5, 5q- | 3,508 | 16.9 |
| +8 | +8 | 3,377 | 16.27 |
| -7 | -7 | 2,783 | 13.41 |
| 17p loss | -17, 17p-, i(17q) | 2,081 | 10.03 |
| t(8;21) | t(8;21) | 2,070 | 9.97 |
| KMT2A | KMT2A | 1,762 | 8.49 |
| 7q- | 7q- | 1,602 | 7.72 |
| +21 | +21 | 1,443 | 6.95 |
| -Y | -Y | 1,322 | 6.37 |
| inv(16) | inv(16) | 1,262 | 6.08 |
| 18 loss | -18, 18q- | 1,111 | 5.35 |
| +22 | +22 | 873 | 4.21 |
| +13 | +13 | 843 | 4.06 |
| +19 | +19 | 770 | 3.71 |
| 9q- | 9q- | 745 | 3.59 |
| 12p loss | 12p- | 742 | 3.58 |
| +11 | +11 | 738 | 3.56 |
| inv(3)/t(3;3) | inv(3)/t(3;3) | 677 | 3.26 |
| -X | -X | 650 | 3.13 |

Below N or beyond the cap: 20q- (538), +4 (488), 3q26 other (402), 3p- (377), t(9;22) (345), t(6;9) (185)

## APL (n = 2,341, N = 196)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| t(15;17) | t(15;17) | 1,920 | 82.02 |
| +8 | +8 | 279 | 11.92 |

Below N or beyond the cap: 17p loss (54), 7q- (49), 9q- (48), +21 (41), -Y (39), 5q loss (34), -7 (22), 6q- (16), 12p loss (13), +22 (10), -X (7), 1q gain (5)

## MDS (n = 5,734, N = 200)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| 5q loss | -5, 5q- | 2,104 | 36.69 |
| complex | complex | 1,823 | 31.79 |
| -7 | -7 | 1,104 | 19.25 |
| +8 | +8 | 951 | 16.59 |
| 7q- | 7q- | 757 | 13.2 |
| 17p loss | -17, 17p-, i(17q) | 612 | 10.67 |
| 20q- | 20q- | 550 | 9.59 |
| 1q gain | 1q+ | 371 | 6.47 |
| 18 loss | -18, 18q- | 371 | 6.47 |
| 13q loss | -13, 13q- | 316 | 5.51 |
| -Y | -Y | 289 | 5.04 |
| 12p loss | 12p- | 264 | 4.6 |
| +21 | +21 | 254 | 4.43 |

Below N or beyond the cap: 11q- (191), +19 (189), 3p- (140), +22 (130), 9q- (110), 3q26 other (100), inv(3)/t(3;3) (88), -X (82), idic(X)(q13) (23)

## CML (n = 4,819, N = 199)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| t(9;22) | t(9;22) | 3,517 | 72.98 |
| +8 | +8 | 1,108 | 22.99 |
| extra Ph | extraPh | 734 | 15.23 |
| 17p loss | -17, 17p- | 593 | 12.31 |
| i(17q) | i(17q) | 416 | 8.63 |
| +19 | +19 | 306 | 6.35 |
| -Y | -Y | 294 | 6.1 |
| +21 | +21 | 203 | 4.21 |
| 3q26 | 3q26, inv(3)/t(3;3) | 203 | 4.21 |

Below N or beyond the cap: -7 (169), complex (130), +22 (125), +10 (94), +6 (94), +17 (89), 9p loss (81), near-tetraploid (26)

## ALL (n = 12,842, N = 201)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| +21 | +21 | 2,960 | 23.05 |
| 9p loss | -9, 9p- | 2,106 | 16.4 |
| high hyperdiploid | hyperdiploid51_67 | 2,050 | 15.96 |
| +X | +X | 1,887 | 14.69 |
| +14 | +14 | 1,585 | 12.34 |
| +6 | +6 | 1,562 | 12.16 |
| +10 | +10 | 1,496 | 11.65 |
| +18 | +18 | 1,452 | 11.31 |
| +4 | +4 | 1,362 | 10.61 |
| t(9;22) | t(9;22) | 1,347 | 10.49 |
| +17 | +17 | 1,112 | 8.66 |
| KMT2A | KMT2A | 1,040 | 8.1 |
| 13q loss | -13, 13q- | 943 | 7.34 |
| 6q- | 6q- | 931 | 7.25 |
| -7 | -7 | 900 | 7.01 |
| 17p loss | -17, 17p-, i(17q) | 817 | 6.36 |
| 1q gain | 1q+ | 765 | 5.96 |
| 12p loss | 12p- | 734 | 5.72 |
| hypodiploid | hypodiploid_lt44 | 551 | 4.29 |
| t(1;19) | t(1;19) | 545 | 4.24 |

Below N or beyond the cap: t(12;21) (531), extra Ph (202), t(8;14) (91)

## CLL (n = 3,469, N = 198)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| complex | complex | 1,068 | 30.79 |
| +12 | +12 | 978 | 28.19 |
| 13q loss | -13, 13q- | 633 | 18.25 |
| 17p loss | -17, 17p-, i(17q) | 404 | 11.65 |
| 11q- | 11q- | 392 | 11.3 |
| 6q- | 6q- | 331 | 9.54 |

Below N or beyond the cap: t(14;18) (147), 14q- (124), t(11;14) (101), +18 (101), 8p- (84), +19 (83), +3 (56), 2p gain (39)

## MM (n = 2,096, N = 195)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| 13 loss | -13, 13q- | 849 | 40.51 |
| +9 | +9 | 693 | 33.06 |
| +15 | +15 | 675 | 32.2 |
| +19 | +19 | 595 | 28.39 |
| +3 | +3 | 533 | 25.43 |
| 1q gain | 1q+ | 521 | 24.86 |
| +11 | +11 | 520 | 24.81 |
| +7 | +7 | 519 | 24.76 |
| +5 | +5 | 511 | 24.38 |
| 1p- | 1p- | 475 | 22.66 |
| +21 | +21 | 439 | 20.94 |
| -X | -X | 381 | 18.18 |
| 17p loss | -17, 17p- | 315 | 15.03 |
| -Y | -Y | 303 | 14.46 |
| 6q- | 6q- | 301 | 14.36 |
| t(11;14) | t(11;14) | 208 | 9.92 |

Below N or beyond the cap: 8p- (172), 16q- (158), 14q- (77), t(4;14) (31), t(14;16) (29)

## Known limitations

- `complex` follows ELN 2022 (at least 3 aberrations, no recurrent balanced rearrangement). It is partly a function of other events (for example 5q loss, -7, 17p loss), so it is used only in the co-occurrence baseline. It is not a candidate in ALL or MM, where hyperdiploidy is the relevant count-derived class.
- In ALL, `high hyperdiploid` and the individual trisomies describe one process. Ordering between the trisomies mostly reflects how the hyperdiploid clone accumulated gains, not separate steps.
- In CML, `17p loss` includes i(17q), so `17p loss` and `i(17q)` are nested by definition; an edge between them is not biological evidence.
- t(12;21) (ETV6-RUNX1) is cryptic by karyotyping, so Mitelman understates it.
- APL ends up with 2 events and MM has fewer than 200 cases with a clone order, so both are co-occurrence-only under the stop rule.
