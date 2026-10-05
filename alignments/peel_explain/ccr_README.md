# Over-split groups: ccr

Peel settings: Min size 2, at least 2 diagnostic columns, refine 1 (the defaults in ViewAlign v237).
Source alignment: `CURATE__ccr__subfam608_g1-g7_accr.aln.fa`. A group of his is listed when its chunks end up in more than one peel group.

How to read each file: rows are `<piece>|<chunk name>`, ordered piece by piece. Above each piece, `DIAG_<piece>` shows the characters
that made the loop peel that piece: at those columns the piece was (at least 90%) fixed for that character, and at most 2% of the
sequences still in the pool at that moment had it. Everything else in the DIAG row is a gap.

## Summary

| his group | chunks | pieces | mean distance inside his group | how it was split |
|---|---|---|---|---|
| g3 | 137 | 3 | 0.041 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g4 | 40 | 2 | 0.014 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g6 | 61 | 2 | 0.028 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g5 | 5 | 3 | 0.011 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g2 | 69 | 2 | 0.056 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| new13 | 13 | 6 | 0.036 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |

## g3 (137 chunks) -> 3 pieces

File: `peel_explain/ccr__g3__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 75 | 2 | step 18, level 1 | 563 | 4 | 0.020 | 0.058 |
| B | 32 | 39 | step 19, level 1 | 486 | 4 | 0.016 | 0.056 |
| C | 30 | 381 | step 21, level 1 | 413 | 8 | 0.023 | 0.053 |

Pattern sharing: B carries 0% of A's 4 diagnostic characters; C carries 6% of A's 4 diagnostic characters; A carries 26% of B's 4 diagnostic characters; C carries 6% of B's 4 diagnostic characters; A carries 88% of C's 8 diagnostic characters; B carries 100% of C's 8 diagnostic characters.

## g4 (40 chunks) -> 2 pieces

File: `peel_explain/ccr__g4__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 39 | 32 | step 19, level 1 | 486 | 4 | 0.013 | 0.030 |
| B | 1 | 410 | step 21, level 1 | 413 | 8 | NaN | 0.030 |

Pattern sharing: B carries 100% of A's 4 diagnostic characters; A carries 100% of B's 8 diagnostic characters.

## g6 (61 chunks) -> 2 pieces

File: `peel_explain/ccr__g6__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 60 | 351 | step 21, level 1 | 413 | 8 | 0.027 | 0.048 |
| B | 1 | 76 | step 18, level 1 | 563 | 4 | NaN | 0.048 |

Pattern sharing: B carries 88% of A's 8 diagnostic characters; A carries 1% of B's 4 diagnostic characters.

## g5 (5 chunks) -> 3 pieces

File: `peel_explain/ccr__g5__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 2 | 0 | step 3, level 0 | 596 | 9 | 0.000 | 0.012 |
| B | 2 | 0 | step 2, level 0 | 598 | 10 | 0.000 | 0.010 |
| C | 1 | 76 | step 18, level 1 | 563 | 4 | NaN | 0.023 |

Pattern sharing: B carries 100% of A's 9 diagnostic characters; C carries 100% of A's 9 diagnostic characters; A carries 90% of B's 10 diagnostic characters; C carries 100% of B's 10 diagnostic characters; A carries 75% of C's 4 diagnostic characters; B carries 75% of C's 4 diagnostic characters.

## g2 (69 chunks) -> 2 pieces

File: `peel_explain/ccr__g2__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 67 | 344 | step 21, level 1 | 413 | 8 | 0.054 | 0.091 |
| B | 2 | 0 | step 1, level 0 | 600 | 11 | 0.035 | 0.091 |

Pattern sharing: B carries 94% of A's 8 diagnostic characters; A carries 0% of B's 11 diagnostic characters.

## new13 (13 chunks) -> 6 pieces

File: `peel_explain/ccr__new13__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 4 | 0 | step 10, level 1 | 581 | 6 | 0.006 | 0.039 |
| B | 2 | 0 | step 12, level 1 | 575 | 7 | 0.036 | 0.039 |
| C | 2 | 0 | step 11, level 1 | 577 | 6 | 0.010 | 0.035 |
| D | 2 | 0 | step 13, level 1 | 573 | 6 | 0.015 | 0.039 |
| E | 2 | 0 | step 14, level 1 | 571 | 6 | 0.021 | 0.038 |
| F | 1 | 410 | step 21, level 1 | 413 | 8 | NaN | 0.051 |

Pattern sharing: B carries 100% of A's 6 diagnostic characters; C carries 83% of A's 6 diagnostic characters; D carries 83% of A's 6 diagnostic characters; E carries 83% of A's 6 diagnostic characters; F carries 83% of A's 6 diagnostic characters; A carries 100% of B's 7 diagnostic characters; C carries 86% of B's 7 diagnostic characters; D carries 86% of B's 7 diagnostic characters; E carries 86% of B's 7 diagnostic characters; F carries 86% of B's 7 diagnostic characters; A carries 100% of C's 6 diagnostic characters; B carries 100% of C's 6 diagnostic characters; D carries 100% of C's 6 diagnostic characters; E carries 100% of C's 6 diagnostic characters; F carries 100% of C's 6 diagnostic characters; A carries 100% of D's 6 diagnostic characters; B carries 100% of D's 6 diagnostic characters; C carries 100% of D's 6 diagnostic characters; E carries 100% of D's 6 diagnostic characters; F carries 100% of D's 6 diagnostic characters; A carries 100% of E's 6 diagnostic characters; B carries 100% of E's 6 diagnostic characters; C carries 100% of E's 6 diagnostic characters; D carries 100% of E's 6 diagnostic characters; F carries 100% of E's 6 diagnostic characters; A carries 100% of F's 8 diagnostic characters; B carries 100% of F's 8 diagnostic characters; C carries 100% of F's 8 diagnostic characters; D carries 100% of F's 8 diagnostic characters; E carries 100% of F's 8 diagnostic characters.
