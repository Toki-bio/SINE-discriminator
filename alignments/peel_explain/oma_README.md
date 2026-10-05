# Over-split groups: oma

Peel settings: Min size 2, at least 2 diagnostic columns, refine 1 (the defaults in ViewAlign v237).
Source alignment: `CURATE__oma__subfam600_23seeds.aln.fa`. A group of his is listed when its chunks end up in more than one peel group.

How to read each file: rows are `<piece>|<chunk name>`, ordered piece by piece. Above each piece, `DIAG_<piece>` shows the characters
that made the loop peel that piece: at those columns the piece was (at least 90%) fixed for that character, and at most 2% of the
sequences still in the pool at that moment had it. Everything else in the DIAG row is a gap.

## Summary

| his group | chunks | pieces | mean distance inside his group | how it was split |
|---|---|---|---|---|
| SINE27 | 37 | 2 | 0.070 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| grp128 | 8 | 3 | 0.361 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE22 | 67 | 2 + 1 unassigned | 0.209 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| long12 | 12 | 5 | 0.083 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); a piece is defined mostly by GAP columns (shared truncation / indel) |
| SINE31 | 15 | 2 | 0.278 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); a piece is defined mostly by GAP columns (shared truncation / indel) |
| ancient78 | 78 | 4 + 70 unassigned | 0.572 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE25 | 31 | 3 | 0.066 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE19 | 10 | 5 | 0.444 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE26 | 9 | 2 | 0.032 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE24 | 24 | 2 | 0.054 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE21 | 12 | 6 | 0.119 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE18 | 8 | 3 + 1 unassigned | 0.129 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE32 | 8 | 2 | 0.026 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE2 | 4 | 2 | 0.069 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE5 | 17 | 3 | 0.055 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |

## SINE27 (37 chunks) -> 2 pieces

File: `peel_explain/oma__SINE27__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 36 | 0 | step 44, level 1 | 137 | 22 | 0.061 | 0.233 |
| B | 1 | 44 | step 45, level 1 | 101 | 31 | NaN | 0.233 |

Pattern sharing: B carries 41% of A's 22 diagnostic characters; A carries 33% of B's 31 diagnostic characters.

## grp128 (8 chunks) -> 3 pieces

File: `peel_explain/oma__grp128__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 3 | 0 | leftover of a refined group | - | 0 | 0.462 | 0.385 |
| B | 3 | 0 | step 47, level 1 | 8 | 7 | 0.175 | 0.359 |
| C | 2 | 0 | step 48, level 1 | 5 | 7 | 0.319 | 0.380 |

Pattern sharing: A carries 0% of B's 7 diagnostic characters; C carries 0% of B's 7 diagnostic characters; A carries 0% of C's 7 diagnostic characters; B carries 24% of C's 7 diagnostic characters.

## SINE22 (67 chunks) -> 2 pieces + 1 unassigned

File: `peel_explain/oma__SINE22__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 64 | 0 | step 57, level 1 | 66 | 12 | 0.189 | 0.360 |
| B | 2 | 0 | step 58, level 1 | 2 | 0 | 0.348 | 0.360 |

Pattern sharing: B carries 0% of A's 12 diagnostic characters.

## long12 (12 chunks) -> 5 pieces

File: `peel_explain/oma__long12__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 5 | 0 | step 36, level 1 | 7 | 9 | 0.009 | 0.081 |
| B | 2 | 0 | step 3, level 0 | 589 | 101 | 0.157 | 0.109 |
| C | 2 | 0 | step 2, level 0 | 591 | 113 | 0.080 | 0.098 |
| D | 2 | 0 | step 37, level 1 | 2 | 0 | 0.025 | 0.064 |
| E | 1 | 1 | leftover of a refined group | - | 0 | NaN | 0.178 |

Pattern sharing: B carries 50% of A's 9 diagnostic characters; C carries 17% of A's 9 diagnostic characters; D carries 0% of A's 9 diagnostic characters; E carries 0% of A's 9 diagnostic characters; A carries 99% of B's 101 diagnostic characters; C carries 95% of B's 101 diagnostic characters; D carries 100% of B's 101 diagnostic characters; E carries 75% of B's 101 diagnostic characters; A carries 91% of C's 113 diagnostic characters; B carries 88% of C's 113 diagnostic characters; D carries 97% of C's 113 diagnostic characters; E carries 82% of C's 113 diagnostic characters.

## SINE31 (15 chunks) -> 2 pieces

File: `peel_explain/oma__SINE31__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 14 | 0 | step 40, level 1 | 16 | 10 | 0.271 | 0.323 |
| B | 1 | 1 | leftover of a refined group | - | 0 | NaN | 0.323 |

Pattern sharing: B carries 0% of A's 10 diagnostic characters.

## ancient78 (78 chunks) -> 4 pieces + 70 unassigned

File: `peel_explain/oma__ancient78__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 2 | 0 | step 33, level 0 | 88 | 3 | 0.310 | 0.566 |
| B | 2 | 0 | step 32, level 0 | 90 | 2 | 0.333 | 0.604 |
| C | 2 | 0 | step 34, level 0 | 86 | 2 | 0.413 | 0.628 |
| D | 2 | 0 | step 31, level 0 | 92 | 3 | 0.444 | 0.652 |

Pattern sharing: B carries 33% of A's 3 diagnostic characters; C carries 0% of A's 3 diagnostic characters; D carries 0% of A's 3 diagnostic characters; A carries 25% of B's 2 diagnostic characters; C carries 0% of B's 2 diagnostic characters; D carries 0% of B's 2 diagnostic characters; A carries 0% of C's 2 diagnostic characters; B carries 0% of C's 2 diagnostic characters; D carries 25% of C's 2 diagnostic characters; A carries 0% of D's 3 diagnostic characters; B carries 0% of D's 3 diagnostic characters; C carries 0% of D's 3 diagnostic characters.

## SINE25 (31 chunks) -> 3 pieces

File: `peel_explain/oma__SINE25__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 28 | 0 | step 60, level 1 | 30 | 16 | 0.044 | 0.162 |
| B | 2 | 0 | step 59, level 1 | 32 | 9 | 0.000 | 0.180 |
| C | 1 | 1 | leftover of a refined group | - | 0 | NaN | 0.138 |

Pattern sharing: B carries 69% of A's 16 diagnostic characters; C carries 0% of A's 16 diagnostic characters; A carries 0% of B's 9 diagnostic characters; C carries 0% of B's 9 diagnostic characters.

## SINE19 (10 chunks) -> 5 pieces

File: `peel_explain/oma__SINE19__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 2 | 0 | step 16, level 0 | 542 | 34 | 0.413 | 0.480 |
| B | 2 | 0 | step 12, level 0 | 550 | 41 | 0.196 | 0.471 |
| C | 2 | 0 | step 15, level 0 | 544 | 34 | 0.368 | 0.452 |
| D | 2 | 0 | step 14, level 0 | 546 | 36 | 0.268 | 0.438 |
| E | 2 | 0 | step 29, level 0 | 126 | 24 | 0.287 | 0.464 |

Pattern sharing: B carries 93% of A's 34 diagnostic characters; C carries 60% of A's 34 diagnostic characters; D carries 65% of A's 34 diagnostic characters; E carries 22% of A's 34 diagnostic characters; A carries 77% of B's 41 diagnostic characters; C carries 52% of B's 41 diagnostic characters; D carries 55% of B's 41 diagnostic characters; E carries 12% of B's 41 diagnostic characters; A carries 51% of C's 34 diagnostic characters; B carries 53% of C's 34 diagnostic characters; D carries 51% of C's 34 diagnostic characters; E carries 32% of C's 34 diagnostic characters; A carries 51% of D's 36 diagnostic characters; B carries 54% of D's 36 diagnostic characters; C carries 56% of D's 36 diagnostic characters; E carries 33% of D's 36 diagnostic characters; A carries 48% of E's 24 diagnostic characters; B carries 46% of E's 24 diagnostic characters; C carries 75% of E's 24 diagnostic characters; D carries 83% of E's 24 diagnostic characters.

## SINE26 (9 chunks) -> 2 pieces

File: `peel_explain/oma__SINE26__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 7 | 0 | step 39, level 1 | 7 | 0 | 0.024 | 0.045 |
| B | 2 | 0 | step 38, level 1 | 9 | 4 | 0.026 | 0.045 |

Pattern sharing: A carries 0% of B's 4 diagnostic characters.

## SINE24 (24 chunks) -> 2 pieces

File: `peel_explain/oma__SINE24__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 18 | 0 | leftover of a refined group | - | 0 | 0.052 | 0.062 |
| B | 6 | 0 | step 53, level 1 | 24 | 4 | 0.014 | 0.062 |

Pattern sharing: A carries 0% of B's 4 diagnostic characters.

## SINE21 (12 chunks) -> 6 pieces

File: `peel_explain/oma__SINE21__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 2 | 0 | step 13, level 0 | 548 | 37 | 0.171 | 0.167 |
| B | 2 | 0 | step 11, level 0 | 552 | 48 | 0.085 | 0.110 |
| C | 2 | 0 | step 10, level 0 | 554 | 41 | 0.014 | 0.126 |
| D | 2 | 0 | step 7, level 0 | 560 | 44 | 0.055 | 0.108 |
| E | 2 | 0 | step 4, level 0 | 587 | 41 | 0.000 | 0.107 |
| F | 2 | 0 | step 8, level 0 | 558 | 51 | 0.101 | 0.125 |

Pattern sharing: B carries 97% of A's 37 diagnostic characters; C carries 95% of A's 37 diagnostic characters; D carries 99% of A's 37 diagnostic characters; E carries 100% of A's 37 diagnostic characters; F carries 99% of A's 37 diagnostic characters; A carries 80% of B's 48 diagnostic characters; C carries 95% of B's 48 diagnostic characters; D carries 96% of B's 48 diagnostic characters; E carries 100% of B's 48 diagnostic characters; F carries 95% of B's 48 diagnostic characters; A carries 83% of C's 41 diagnostic characters; B carries 96% of C's 41 diagnostic characters; D carries 93% of C's 41 diagnostic characters; E carries 95% of C's 41 diagnostic characters; F carries 91% of C's 41 diagnostic characters; A carries 68% of D's 44 diagnostic characters; B carries 84% of D's 44 diagnostic characters; C carries 77% of D's 44 diagnostic characters; E carries 100% of D's 44 diagnostic characters; F carries 92% of D's 44 diagnostic characters; A carries 49% of E's 41 diagnostic characters; B carries 72% of E's 41 diagnostic characters; C carries 57% of E's 41 diagnostic characters; D carries 88% of E's 41 diagnostic characters; F carries 95% of E's 41 diagnostic characters; A carries 61% of F's 51 diagnostic characters; B carries 80% of F's 51 diagnostic characters; C carries 69% of F's 51 diagnostic characters; D carries 88% of F's 51 diagnostic characters; E carries 100% of F's 51 diagnostic characters.

## SINE18 (8 chunks) -> 3 pieces + 1 unassigned

File: `peel_explain/oma__SINE18__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 3 | 0 | step 54, level 1 | 7 | 10 | 0.039 | 0.112 |
| B | 2 | 0 | step 55, level 1 | 4 | 2 | 0.020 | 0.085 |
| C | 2 | 0 | step 56, level 1 | 2 | 0 | 0.072 | 0.097 |

Pattern sharing: B carries 0% of A's 10 diagnostic characters; C carries 0% of A's 10 diagnostic characters; A carries 67% of B's 2 diagnostic characters; C carries 0% of B's 2 diagnostic characters.

## SINE32 (8 chunks) -> 2 pieces

File: `peel_explain/oma__SINE32__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 6 | 0 | leftover of a refined group | - | 0 | 0.022 | 0.033 |
| B | 2 | 0 | step 52, level 1 | 8 | 3 | 0.000 | 0.033 |

Pattern sharing: A carries 0% of B's 3 diagnostic characters.

## SINE2 (4 chunks) -> 2 pieces

File: `peel_explain/oma__SINE2__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 2 | 0 | step 17, level 0 | 540 | 33 | 0.105 | 0.074 |
| B | 2 | 0 | step 9, level 0 | 556 | 41 | 0.012 | 0.074 |

Pattern sharing: B carries 100% of A's 33 diagnostic characters; A carries 87% of B's 41 diagnostic characters.

## SINE5 (17 chunks) -> 3 pieces

File: `peel_explain/oma__SINE5__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 12 | 0 | step 41, level 1 | 17 | 9 | 0.032 | 0.083 |
| B | 3 | 0 | step 43, level 1 | 3 | 0 | 0.046 | 0.075 |
| C | 2 | 0 | step 42, level 1 | 5 | 3 | 0.016 | 0.079 |

Pattern sharing: B carries 0% of A's 9 diagnostic characters; C carries 0% of A's 9 diagnostic characters; A carries 39% of C's 3 diagnostic characters; B carries 0% of C's 3 diagnostic characters.
