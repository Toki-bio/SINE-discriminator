# Findings — what his judgements have established

Separate from `HANDOFF.md` (working history) and `PLAN.md` (what to do next).
This file holds only **results that change the tool**, newest first. Updated
after every exchange.

*(Build-1 findings from 2026-08-31 moved to `FINDINGS_old_2026-08.md`.)*

---

## 2026-09-01 — Review batch 2 opened; two failures he caught

### R7. `CONSENSUS_OVEREXTENDED` has a false negative

He judged `NEGCHIM__ccr__g1_180seqs`: *"overextended consensus and no degapped
right flank"*. The tool flags `FLANKS_UNMEASURED, RECOVERABLE_CORE` — **the
overextension flag does not fire**, although it fires correctly on `g3` and `g5`
which he confirmed.

So the flag is right when it fires (4 of 4 in batch 1) but misses cases. Its
condition — a core window beating the whole span by ≥0.15 with the remainder at
background — is too strict. Needs re-fitting against `g1` as a positive.

### R8. I asked about an alignment whose answer he had already given

> "it feels like you havent tried to learn anything from my previous answers
> because it perfectly matches what i have already answered"

He had judged **eight** NEGCHIM sets, five of them *"consensus too long but good
sine"*. I put a ninth in the next batch. He answered it identically.

Fixed by `already_answered.py`: before asking, check whether he has answered
enough alignments of the same construction, consistently, and if so predict
rather than ask.

Grouping by exact label was too literal — his nine NEGCHIM answers use six
different labels and the top one is only 44%. Grouped by **substance** (is it a
SINE at all, ignoring which extra property is named) it is **89% "a SINE"**.
That is the level at which his answers agree.

Constructions now predictable without asking:

| construction | judged | agreement | answer |
|---|---|---|---|
| NEGSPLICE | 20 | 100% | a SINE |
| NEGSAT | 3 | 100% | not a SINE |
| NEGCHIM | 9 | 89% | a SINE |
| ERI | 4 | 75% | a SINE |
| HUM | 26 | 62% | a SINE |

### R9. The flank trim needed two more fixes

His *"no degapped right flank"* on `g1_180seqs` was still true after the first
fix. Two causes:

1. A 75th-percentile cap is not enough when flank lengths are skewed — right
   flank median 20 bases but p75 58, leaving the panel 53% gaps. Now the cap is
   chosen as the widest width whose panel is at most 25% gaps.
2. A constant floor of 25 columns cannot help when copies carry ~10 bases. The
   floor now follows the copies' own median.

| alignment | right-flank gaps: raw → first fix → now |
|---|---|
| `NEGCHIM__ccr__g1_180seqs` | 0.78 → 0.53 → **0.36** |
| `NEGCHIM__ccr__g3_71seqs` | 0.85 → 0.55 → **0.42** |
| `MIXED10__ccr__g4_77seqs` | 0.29 → 0.25 → **0.25** (left 0.80 → **0.23**) |

The residue is intrinsic: those copies genuinely carry only 10–20 bases of right
flank, so any panel showing them has gaps for the shorter ones.

---

## 2026-09-01 — Review batch 1: 16 judgements, six results

Recorded verbatim in `calls.tsv` (now 87 judgements).

### R1. The mixture test is wrong 3 times out of 4 — must be re-fitted

| alignment | what he said |
|---|---|
| `NEGCHIM__saq__s5_5seqs` | "not mixture" |
| `NEGCHIM__saq__s8_225seqs` | "not mixture" |
| `NEGCHIM__ccr__g5_7seqs` | "i dont see mixture" |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | "mixture of 2 subfamilies" — the only true one |

`HETEROGENEOUS_SELECTION` fires on all four. **75 % false positive rate.**

Its threshold (relative group length difference ≥ 0.12) was fitted to three
examples and was already flagged as overfit in `HANDOFF.md` §53a. This confirms
it. And the one true positive is a *subfamily* mixture — which by his own rule
does not affect the sine/not-sine verdict anyway.

**Action:** re-fit against these four; reconsider whether a subfamily mixture
should trigger it at all.

### R2. The over-long-consensus test is right 4 of 4 — but it is a NOTE, not a fault

| alignment | what he said |
|---|---|
| `NEGCHIM__ccr__g3_71seqs` | "yes consensus is too long, but its a good sine" |
| `NEGCHIM__ccr__g6_58seqs` | same |
| `NEGCHIM__dmo__d1_16seqs` | same |
| `NEGCHIM__ccr__g5_7seqs` | "consensus too long in 3' end" |

`CONSENSUS_OVEREXTENDED` is correct every time, **but he calls all of them good
SINEs**. It must never lower the verdict — it is an instruction to trim, not
evidence against the family. He also localises it further than the tool does
("3' end"): the flag should say which end.

### R3. The copy-count threshold is bounded, and 30 is too high

| alignment | n | what he said |
|---|---|---|
| `HUM__hum__AluYb9` | 21 | "no problems good sine" |
| `TIMB__timb__SINE_46` | 11 | "not enough sequences but left flank is probably bad or incomplete" |

`INSUFFICIENT_COPIES` fires below 30, so it wrongly fires on `AluYb9`.
**The line sits between 11 and 21.** Needs a judgement near n≈15 to place it.

### R4. Display bug — flank width set by one outlier copy *(FIXED)*

> "flanks are badly degapped and it hinders my estimates"

The flanks **were** justified — every copy sits hard against the element,
gap-to-element 0 for all. The fault was the column **width**, set by the single
longest copy while the median copy has far fewer bases.

`NEGCHIM__ccr__g3_71seqs`: right flank 125 columns, median copy has 10 bases in
it, panel 85 % gaps. One outlier padding seventy copies with whitespace.

Fixed by `trim_flanks.py` — cap each flank at the 75th percentile of the copies'
own flank lengths:

| alignment | flank gaps before | after |
|---|---|---|
| `NEGCHIM__ccr__g3_71seqs` | 0.56 | **0.16** |
| `NEGSAT__eri__r02` | 0.68 | **0.26** |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | 0.67 | **0.07** |
| `NEGCHIM__saq__s8_225seqs` | 0.45 | **0.23** |

All review alignments republished as `.trim.aln.fa`.

### R5. Some published alignments are too poor to judge — confirmed

> "you provided it unaligned; too long and no flanks"

Technically aligned (all rows equal length) but **59–70 % gaps, smeared to 4–6×
the element length**, against 33 % and 2.3× for a good one. Any alignment past
roughly 3× its consensus length should not be shown for judgement without
saying so.

### R6. Batch size

> "actually i asked by 6 not to overload me"

**Six alignments per review batch.** Batch 1 had 28.

---

## 2026-09-01 — Tests 1–3

### T1. Three existing measurements carry what the score discards

- **`edge_drop`** (how far copy-to-copy similarity falls at the element edge):
  clean 0.52, grey 0.15, badly-presented 0.18, wrong-consensus 0.13 — all of
  which score 89–100 and are identical by score.
- **`pair_bg`** explains `ERI e1-4`: flanks 0.57 similar against a normal 0.25.
  The copies share flanking sequence — a different statement from "no element".
- **`cov_min` = 0.000** marks every badly-presented alignment.

### T2. The composition-aware 3′ edge rule WORKS — the plan's go/no-go, passed

Identity alone fails: it collapses in the poly-A tail because A-run length
varies between copies (AluY 0.93 in the body → 0.27 in the last 90 bp, while
A+T rises 0.26 → 0.63). An identity-only rule amputates the tail from every SINE.

Watching composition as well as identity, the tail correctly kept:

| what he said | tail rescued |
|---|---|
| SINE | **53.9 bp** |
| SINE_MINOR_NOTE | 42.5 bp |
| MIXTURE_SPLIT_FIRST | 56.0 bp |
| NOT_SINE | **10.6 bp** |
| ANCIENT_BADLY_PRESENTED | 9.3 bp |
| SINE_CONSENSUS_WRONG | **0.0 bp** |

Real SINEs carry ~50 bp of A/T-rich tail; things he called not-SINEs carry
almost none. So it is both a boundary rule and a discriminating measurement.

### T3. `rowsplit_dident` earns its place; local instability does not

Of 10 comparable pairs, all separate. **`rowsplit_dident`** — how much two
groups of *copies* differ in identity — is the best separator for both mixture
comparisons, a property nothing measured before (everything else is per-column;
he was describing groups of rows).

**`dip_depth`** (local instability) separates nothing: 0.55 for clean SINEs vs
0.61 for the "island of instability" case. Dropped.

### T3-bug. A zero-variance artifact, caught before reporting it as success

The first Test 3 run reported separations of ~1e11 and "0 of 120 pairs not told
apart". Nine of his sixteen properties have a single example, whose variance is
zero, so the pooled spread collapses to the epsilon. Now requires 3 examples per
side — leaving **only 5 of 16 properties comparable at all**. That is why Test 4
exists.

---

## Standing conclusions

- **He names a property and attaches an action; the tool emits a number.** Four
  of his properties all score exactly 100.0. No threshold tuning fixes this.
- **Human Dfam is not ground truth.** His judgements and the curated `tim/` set
  are.
- **Two thresholds are fitted to almost nothing** and must be re-fitted: the
  mixture test (now measured at 75 % false positive) and the over-long-consensus
  tail bound.
