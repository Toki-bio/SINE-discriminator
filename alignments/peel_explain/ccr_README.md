# Over-split groups: ccr

Peel settings: Min size 3, at least 2 diagnostic columns, refine 1, second pass 4 columns (indel x2) (the defaults in ViewAlign v238).
Source alignment: `CURATE__ccr__subfam608_g1-g7_accr.aln.fa`. A group of his is listed when its chunks end up in more than one peel group.

How to read each file: rows are `<piece>|<chunk name>`, ordered piece by piece. Above each piece, `DIAG_<piece>` shows the characters
that made the loop peel that piece: at those columns the piece was (at least 90%) fixed for that character, and at most 2% of the
sequences still in the pool at that moment had it. Everything else in the DIAG row is a gap.

## Summary

| his group | chunks | pieces | mean distance inside his group | how it was split |
|---|---|---|---|---|
| g3 | 137 | 3 | 0.041 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g6 | 61 | 2 | 0.028 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |
| g2 | 69 | 3 | 0.056 | split inside a group the loop had first peeled as one (the refine step found a sub-group with its own columns); part of it was absorbed into a group made mostly of OTHER chunks |

## g3 (137 chunks) -> 3 pieces

File: `peel_explain/ccr__g3__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 75 | 1 | step 6, level 0 | 571 | 4 | 0.020 | 0.058 |
| B | 33 | 40 | step 7, level 0 | 495 | 4 | 0.017 | 0.056 |
| C | 29 | 71 | step 13, level 1 | 323 | 4 | 0.020 | 0.053 |

Pattern sharing: B carries 1% of A's 4 diagnostic characters; C carries 5% of A's 4 diagnostic characters; A carries 26% of B's 4 diagnostic characters; C carries 3% of B's 4 diagnostic characters; A carries 100% of C's 4 diagnostic characters; B carries 99% of C's 4 diagnostic characters.

## g6 (61 chunks) -> 2 pieces

File: `peel_explain/ccr__g6__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 60 | 40 | step 13, level 1 | 323 | 4 | 0.027 | 0.048 |
| B | 1 | 75 | step 6, level 0 | 571 | 4 | NaN | 0.048 |

Pattern sharing: B carries 100% of A's 4 diagnostic characters; A carries 1% of B's 4 diagnostic characters.

## g2 (69 chunks) -> 3 pieces

File: `peel_explain/ccr__g2__oversplit.aln.fa`

| piece | his chunks | other chunks in it | peeled at | pool then | diag columns | within | to other pieces |
|---|---|---|---|---|---|---|---|
| A | 64 | 22 | step 11, level 0 | 412 | 2 | 0.055 | 0.063 |
| B | 3 | 0 | step 8, level 0 | 422 | 3 | 0.014 | 0.058 |
| C | 2 | 98 | step 13, level 1 | 323 | 4 | 0.031 | 0.073 |

Pattern sharing: B carries 100% of A's 2 diagnostic characters; C carries 25% of A's 2 diagnostic characters; A carries 1% of B's 3 diagnostic characters; C carries 0% of B's 3 diagnostic characters; A carries 49% of C's 4 diagnostic characters; B carries 50% of C's 4 diagnostic characters.
