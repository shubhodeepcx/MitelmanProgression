# Frozen event vocabulary (v1)

Frozen on 2026-10-12 against the Mitelman download with md5 `434cf3759ba633e4d9d99e2bb8208283`.
Rule: an event enters a disease's vocabulary if it is present (any abnormal clone, first investigation, duplicates removed) in at least N cases, where N comes from the power calculation in `config/power_threshold.json`; at most 20 events, ranked by count; defining events are always kept. Candidates were listed in `config/vocabulary_candidates.json` before any ordering result was looked at.

An event is present in a clone when any of its atomic features is (see `karyoprog/features.py`). Copy-number features are the union of what the ISCN text states and CytoConverter arm-level calls (≥50% of an arm, ≥90% of a chromosome for whole gains/losses, relative to the clone's own baseline ploidy).

## AML (n = 20,723, N = 201)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| complex | complex | 5,570 | 26.88 |
| 5q loss | -5, 5q- | 3,495 | 16.87 |
| +8 | +8 | 3,369 | 16.26 |
| -7 | -7 | 2,775 | 13.39 |
| t(8;21) | t(8;21) | 2,113 | 10.2 |
| 17p loss | -17, 17p-, i(17q) | 2,071 | 9.99 |
| KMT2A | KMT2A | 1,760 | 8.49 |
| 7q- | 7q- | 1,602 | 7.73 |
| +21 | +21 | 1,444 | 6.97 |
| -Y | -Y | 1,321 | 6.37 |
| inv(16) | inv(16) | 1,297 | 6.26 |
| 18 loss | -18, 18q- | 1,106 | 5.34 |
| +22 | +22 | 874 | 4.22 |
| +13 | +13 | 845 | 4.08 |
| +19 | +19 | 771 | 3.72 |
| 9q- | 9q- | 744 | 3.59 |
| 12p loss | 12p- | 744 | 3.59 |
| +11 | +11 | 739 | 3.57 |
| inv(3)/t(3;3) | inv(3)/t(3;3) | 677 | 3.27 |
| -X | -X | 650 | 3.14 |

Below N or beyond the cap: 20q- (540), +4 (489), 3q26 other (402), 3p- (379), t(9;22) (366), t(6;9) (211)

## APL (n = 2,371, N = 196)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| t(15;17) | t(15;17) | 1,970 | 83.09 |
| +8 | +8 | 287 | 12.1 |

Below N or beyond the cap: 17p loss (70), 7q- (56), 5q loss (51), 9q- (49), +21 (45), -Y (40), -7 (30), 12p loss (18), 6q- (16), +22 (12), -X (7), 1q gain (5)

## MDS (n = 5,734, N = 200)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| 5q loss | -5, 5q- | 2,107 | 36.75 |
| complex | complex | 1,822 | 31.78 |
| -7 | -7 | 1,104 | 19.25 |
| +8 | +8 | 951 | 16.59 |
| 7q- | 7q- | 761 | 13.27 |
| 17p loss | -17, 17p-, i(17q) | 615 | 10.73 |
| 20q- | 20q- | 553 | 9.64 |
| 1q gain | 1q+ | 377 | 6.57 |
| 18 loss | -18, 18q- | 372 | 6.49 |
| 13q loss | -13, 13q- | 317 | 5.53 |
| -Y | -Y | 289 | 5.04 |
| 12p loss | 12p- | 264 | 4.6 |
| +21 | +21 | 257 | 4.48 |

Below N or beyond the cap: 11q- (191), +19 (189), 3p- (142), +22 (130), 9q- (110), 3q26 other (100), inv(3)/t(3;3) (88), -X (82), idic(X)(q13) (23)

## CML (n = 4,819, N = 199)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| t(9;22) | t(9;22) | 3,552 | 73.71 |
| +8 | +8 | 1,108 | 22.99 |
| extra Ph | extraPh | 734 | 15.23 |
| 17p loss | -17, 17p- | 593 | 12.31 |
| i(17q) | i(17q) | 416 | 8.63 |
| +19 | +19 | 306 | 6.35 |
| -Y | -Y | 294 | 6.1 |
| +21 | +21 | 203 | 4.21 |
| 3q26 | 3q26, inv(3)/t(3;3) | 203 | 4.21 |

Below N or beyond the cap: -7 (169), +22 (125), complex (116), +10 (94), +6 (94), +17 (89), 9p loss (81), near-tetraploid (26)

## ALL (n = 12,842, N = 201)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| +21 | +21 | 2,961 | 23.06 |
| 9p loss | -9, 9p- | 2,120 | 16.51 |
| high hyperdiploid | hyperdiploid51_67 | 2,050 | 15.96 |
| +X | +X | 1,888 | 14.7 |
| +14 | +14 | 1,590 | 12.38 |
| +6 | +6 | 1,563 | 12.17 |
| +10 | +10 | 1,496 | 11.65 |
| +18 | +18 | 1,455 | 11.33 |
| t(9;22) | t(9;22) | 1,366 | 10.64 |
| +4 | +4 | 1,364 | 10.62 |
| +17 | +17 | 1,113 | 8.67 |
| KMT2A | KMT2A | 1,040 | 8.1 |
| 13q loss | -13, 13q- | 945 | 7.36 |
| 6q- | 6q- | 931 | 7.25 |
| -7 | -7 | 902 | 7.02 |
| 17p loss | -17, 17p-, i(17q) | 822 | 6.4 |
| 1q gain | 1q+ | 766 | 5.96 |
| 12p loss | 12p- | 737 | 5.74 |
| t(1;19) | t(1;19) | 618 | 4.81 |
| hypodiploid | hypodiploid_lt44 | 552 | 4.3 |

Below N or beyond the cap: t(12;21) (540), t(8;14) (243), extra Ph (201)

## CLL (n = 3,469, N = 198)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| complex | complex | 1,061 | 30.59 |
| +12 | +12 | 978 | 28.19 |
| 13q loss | -13, 13q- | 633 | 18.25 |
| 17p loss | -17, 17p-, i(17q) | 406 | 11.7 |
| 11q- | 11q- | 393 | 11.33 |
| 6q- | 6q- | 331 | 9.54 |

Below N or beyond the cap: t(14;18) (158), 14q- (124), t(11;14) (120), +18 (101), 8p- (84), +19 (83), +3 (56), 2p gain (39)

## MM (n = 2,096, N = 195)

| Event | Atomic features | Cases | % |
|---|---|---|---|
| 13 loss | -13, 13q- | 849 | 40.51 |
| +9 | +9 | 695 | 33.16 |
| +15 | +15 | 675 | 32.2 |
| +19 | +19 | 595 | 28.39 |
| +3 | +3 | 534 | 25.48 |
| 1q gain | 1q+ | 528 | 25.19 |
| +7 | +7 | 520 | 24.81 |
| +11 | +11 | 520 | 24.81 |
| +5 | +5 | 512 | 24.43 |
| 1p- | 1p- | 483 | 23.04 |
| +21 | +21 | 441 | 21.04 |
| -X | -X | 381 | 18.18 |
| 17p loss | -17, 17p- | 316 | 15.08 |
| -Y | -Y | 303 | 14.46 |
| 6q- | 6q- | 302 | 14.41 |
| t(11;14) | t(11;14) | 231 | 11.02 |

Below N or beyond the cap: 8p- (175), 16q- (162), 14q- (79), t(4;14) (39), t(14;16) (33)

## Known limitations

- `complex` follows ELN 2022 (≥3 aberrations, no recurrent balanced rearrangement). It is partly a function of other events (e.g. 5q loss, -7, 17p loss), so models must treat it as such. It is not a candidate in ALL or MM, where hyperdiploidy alone meets the count.
- In ALL, `high hyperdiploid` and the individual trisomies (+4, +6, +10, +14, +17, +18, +21, +X) describe one process. Ordering between them mostly reflects how the hyperdiploid clone accumulated gains, not separate steps.
- t(12;21) (ETV6-RUNX1) is missing from ALL because it is cryptic by karyotyping; Mitelman frequencies understate it.
- APL ends up with 2 events and MM's ordering yield is under 200 cases, so both are co-occurrence-only under the stop rule.
