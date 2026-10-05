# Over-split groups: oma

Peel settings: Min size 3, at least 2 diagnostic columns, refine 1, second pass 4 columns (indel x2) (the defaults in ViewAlign v238).
Source alignment: `CURATE__oma__subfam600_23seeds.aln.fa`. A group of his is listed when its chunks end up in more than one peel group.

How to read each file: rows are `<piece>|<chunk name>`, ordered piece by piece. Above each piece, `DIAG_<piece>` shows the characters
that made the loop peel that piece: at those columns the piece was (at least 90%) fixed for that character, and at most 2% of the
sequences still in the pool at that moment had it. Everything else in the DIAG row is a gap.

## Summary

| his group | chunks | pieces | mean distance inside his group | how it was split |
|---|---|---|---|---|
| SINE27 | 37 | 2 | 0.070 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| grp128 | 8 | 3 | 0.361 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE22 | 67 | 2 | 0.209 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| long12 | 12 | 2 | 0.083 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE25 | 31 | 2 | 0.066 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE19 | 10 | 2 | 0.444 | peeled separately from the start: each piece had its own pattern against the whole pool |
| SINE24 | 24 | 2 + 1 unassigned | 0.054 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE21 | 12 | 3 | 0.119 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); a piece is defined mostly by GAP columns (shared truncation / indel) |
| SINE18 | 8 | 2 + 1 unassigned | 0.129 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| SINE5 | 17 | 2 | 0.055 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns) |
| g1 | 4 | 2 | 0.198 | peeled separately from the start: each piece had its own pattern against the whole pool; part of it was absorbed into a group made mostly of OTHER chunks |

## SINE27 (37 chunks) -> 2 pieces

File: `peel_explain/oma__SINE27__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 36 | 0 | step 26, level 1 | 87 | 25 | 0.061 | 0.233 |
| B | 1 | 4 | leftover of a refined group | - | 0 | NaN | 0.233 |

Pattern sharing: B carries 48% of A's 27 diagnostic characters.

## grp128 (8 chunks) -> 3 pieces

File: `peel_explain/oma__grp128__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 3 | 2 | leftover of a refined group | - | 0 | 0.444 | 0.383 |
| B | 3 | 0 | step 28, level 1 | 8 | 4 | 0.175 | 0.359 |
| C | 2 | 1 | step 21, level 0 | 100 | 7 | 0.404 | 0.378 |

Pattern sharing: A carries 0% of B's 4 diagnostic characters; C carries 25% of B's 4 diagnostic characters; A carries 76% of C's 7 diagnostic characters; B carries 100% of C's 7 diagnostic characters.

## SINE22 (67 chunks) -> 2 pieces

File: `peel_explain/oma__SINE22__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 64 | 0 | step 37, level 1 | 67 | 6 | 0.189 | 0.410 |
| B | 3 | 0 | leftover of a refined group | - | 0 | 0.500 | 0.410 |

Pattern sharing: B carries 0% of A's 6 diagnostic characters.

## long12 (12 chunks) -> 2 pieces

File: `peel_explain/oma__long12__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 9 | 0 | step 1, level 0 | 598 | 121 | 0.036 | 0.136 |
| B | 3 | 0 | step 2, level 0 | 589 | 81 | 0.161 | 0.136 |

Pattern sharing: B carries 74% of A's 121 diagnostic characters; A carries 97% of B's 81 diagnostic characters.

## SINE25 (31 chunks) -> 2 pieces

File: `peel_explain/oma__SINE25__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 27 | 0 | step 36, level 1 | 31 | 4 | 0.041 | 0.143 |
| B | 4 | 0 | leftover of a refined group | - | 0 | 0.156 | 0.143 |

Pattern sharing: B carries 8% of A's 3 diagnostic characters.

## SINE19 (10 chunks) -> 2 pieces

File: `peel_explain/oma__SINE19__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 6 | 0 | step 20, level 0 | 106 | 11 | 0.379 | 0.506 |
| B | 4 | 0 | step 8, level 0 | 546 | 28 | 0.358 | 0.506 |

Pattern sharing: B carries 80% of A's 11 diagnostic characters; A carries 48% of B's 28 diagnostic characters.

## SINE24 (24 chunks) -> 2 pieces + 1 unassigned

File: `peel_explain/oma__SINE24__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 16 | 0 | step 34, level 1 | 23 | 4 | 0.025 | 0.052 |
| B | 7 | 0 | step 35, level 1 | 7 | 0 | 0.018 | 0.052 |

Pattern sharing: B carries 0% of A's 3 diagnostic characters.

## SINE21 (12 chunks) -> 3 pieces

File: `peel_explain/oma__SINE21__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 6 | 0 | step 3, level 0 | 586 | 41 | 0.067 | 0.141 |
| B | 3 | 0 | leftover of a refined group | - | 0 | 0.181 | 0.141 |
| C | 3 | 0 | step 23, level 1 | 6 | 4 | 0.036 | 0.130 |

Pattern sharing: B carries 58% of A's 41 diagnostic characters; C carries 59% of A's 41 diagnostic characters; A carries 25% of C's 4 diagnostic characters; B carries 0% of C's 4 diagnostic characters.

## SINE18 (8 chunks) -> 2 pieces + 1 unassigned

File: `peel_explain/oma__SINE18__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 4 | 0 | step 33, level 1 | 4 | 0 | 0.055 | 0.112 |
| B | 3 | 0 | step 32, level 1 | 7 | 11 | 0.039 | 0.112 |

Pattern sharing: A carries 0% of B's 10 diagnostic characters.

## SINE5 (17 chunks) -> 2 pieces

File: `peel_explain/oma__SINE5__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 12 | 0 | step 24, level 1 | 18 | 6 | 0.032 | 0.083 |
| B | 5 | 1 | leftover of a refined group | - | 0 | 0.042 | 0.083 |

Pattern sharing: B carries 0% of A's 4 diagnostic characters.

## g1 (4 chunks) -> 2 pieces

File: `peel_explain/oma__g1__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 3 | 0 | step 13, level 0 | 250 | 24 | 0.055 | 0.342 |
| B | 1 | 2 | step 21, level 0 | 100 | 7 | NaN | 0.342 |

Pattern sharing: B carries 75% of A's 24 diagnostic characters; A carries 86% of B's 7 diagnostic characters.
