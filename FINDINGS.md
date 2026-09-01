# Findings — what his judgements have established

Separate from `HANDOFF.md` (working history) and `PLAN.md` (what to do next).
This file holds only **results that change the tool**, newest first. Updated
after every exchange.

*(Build-1 findings from 2026-08-31 moved to `FINDINGS_old_2026-08.md`.)*

---

## 2026-09-02 — R12. The mosaic he means, caught in the wild — and a measurement for it

He said months of examples earlier: *"What i mean with mosaic looks truly
different and maybe we can get back to it when we'll catch it in the wild."*
None of the six synthetic MOSAIC* constructions was it — he called all six SINEs.

**`aca_SINE_1` (starfish) is it.** His words: *"sine_1 is bad and it includes
mosaicism - inconsistent mixture of conserved and discordant columns/spots"*.
The tool had already rejected it at 26.3, but for the wrong reason (weak overall
identity), not for being a mosaic.

### What does NOT measure it

- **Column jumpiness** — `aca_SINE_1` 0.261 vs a legitimate `pom_SINE_0` at
  0.239. No separation.
- **Block structure** (runs test against chance): `aca_SINE_1` blockiness 1.17,
  while the legitimate `aca_SINE_0` is the second-blockiest at 1.30. Backwards.
- **Mean identity** — just says "less conserved", not "mosaic".

### What does: per-copy variation across regions (`patch2d`)

For each copy, its concordance with the consensus in each 20 bp window; then how
much that varies **within** a copy. A clean family: every copy is uniformly
concordant. A mosaic: a copy is concordant in some regions and discordant in
others, and which regions differ between copies.

| candidate | his judgement | patch2d |
|---|---|---|
| `pom_SINE_0` | strongest — 86 % TSD, perfect B-box | **0.061** |
| `aca_SINE_0` | "sine_0 looks legit" | **0.097** |
| `aca_SINE_1` | **"bad ... includes mosaicism"** | **0.143** |
| `pom_SINE_1` | "possible sine too" — only 35 % TSD | **0.150** |

**It orders all four exactly as he does.** Four points is not a threshold, but it
is the first measurement that responds to the thing he has been describing since
the mosaic discussion began, and it is row-wise-within-column — neither purely
per-column nor purely per-copy, which is why every earlier attempt missed it.

Note `pom_SINE_1` scores highest of all, above the acknowledged mosaic. That is
consistent with its other weaknesses (A-box 3 mismatches, TSD 35 %, no tail) and
suggests it deserves the same scrutiny.

**Needs more of his judgements to become a threshold.** The obvious source is
the 22 Hydra candidates, which are unjudged and from a fourth phylum.

---

## 2026-09-02 — pom_SINE_1 structure, and both candidates vs RepBase

Alignments rebuilt at 100 bp flanks, justified and degapped, published as
`NEW__pom100_*__*.degap.aln.fa`. Flank gaps:

| | raw | justified | degapped |
|---|---|---|---|
| SINE_0 top100 | 0.69 / 0.68 | 0.07 / 0.04 | **0.07 / 0.03** |
| SINE_1 top100 | 0.74 / 0.79 | 0.19 / 0.45 | **0.19 / 0.22** |
| SINE_1 rand100 | 0.78 / 0.87 | 0.38 / 0.70 | **0.23 / 0.30** |

### The two candidates side by side

| | pom_SINE_0 | pom_SINE_1 |
|---|---|---|
| length | 175 bp | 262 bp |
| AnnoSINE source | tRNA | 7SL RNA |
| genomic copies | 537 | 339 |
| **A-box** | pos 2-12, **2 mismatches** | pos 0-10, **3 mismatches** |
| **B-box** | pos 54-62 `GTTCAACTC`, **0 mm** | pos 50-58 `GTTCAAGAC`, **0 mm** |
| **TSD** | **86 % of copies**, median 13 bp | **35 % of copies**, median 8 bp |
| Pol III terminator | 24 % strong, 8 % moderate | 16 % strong, 9 % moderate |
| 3' tail | poly-C, 8 units | none detected |

**Both have a perfect B-box** — an intact RNA-pol-III internal promoter, which is
the hardest part to acquire by chance.

**SINE_1 is the weaker of the two**: its A-box has 3 mismatches rather than 2,
only 35 % of copies carry a TSD against 86 %, and no 3' tail repeat is detected.
Real but less convincing, and worth his eye on the degapped alignment.

### Neither has a relative or a LINE partner

Searched at `word_size 7`, no filtering, against **49,007 RepBase entries**,
of which **1,784 L1** and **7,117 other LINE families** (RTE, CR1, Jockey, L2,
Rex, R1/R2, Nimb, Penelope):

| query | whole RepBase | best LINE hit | where on the LINE |
|---|---|---|---|
| pom_SINE_0 | 3 hits, e 0.25-0.89 | 28 bp, e 0.14 | **1,246 bp** from the 3' end |
| pom_SINE_1 | 1 hit, 19 bp, e 0.39 | 28 bp, e 0.22 | **3,000-4,600 bp** from the 3' end |

**The "similarity to L1 tails" does not survive a targeted search.** Every hit is
14-28 bp with e >= 0.14, and none is near a LINE 3' terminus. A genuine SINE-LINE
partnership puts the SINE's tail on the LINE's last ~50 bp; these sit thousands
of bases inside. It is background.

Combined with the SINEBase survey above — 23-55 % of that database is isolated —
being unrelated to known elements is a normal property of a real SINE, not
evidence against one.

**Open question, unresolved by sequence alone:** if neither has a LINE partner in
RepBase, what mobilises them? Molluscan LINEs are poorly represented there, so
the partner may simply be undescribed — findable by searching the snail genome
itself for a LINE whose 3' end matches these SINEs' tails.

---

## 2026-09-02 — pom_SINE_0 structure, and the SINEBase isolation survey

### pom_SINE_0 has a complete SINE structure

He asked whether it has TSDs or boxes. It has all of them:

| feature | detail |
|---|---|
| **A-box** | positions 2-12, `TGGCGCAACAC`, 2 mismatches |
| **B-box** | positions 54-62, `GTTCAACTC`, **0 mismatches** |
| **TSD** | present in **86 % of copies**, median **13 bp** |
| also detected | tRNA region, conserved core, terminator, tail repeat |

A perfect B-box and TSDs on 86 % of copies means an intact RNA-pol-III internal
promoter and direct repeats flanking each insertion. That is strong internal
evidence, and it is consistent with his reading: a real SINE with no relationship
to known families. 537 genomic copies.

Rebuilt with 100 bp flanks (was 400) on his instruction — published as
`NEW__pom100_SINE_0__*`.

### SINEBase all-vs-all: how much of it is isolated?

The only local SINEBase was 958 x 50 bp fragments already reduced to 85 %
non-redundant — useless for this question. Downloaded the real bank from
sines.eimb.ru: **231 entries, 230 unique ids** (one duplicate header).

blastn all-vs-all, `word_size 7`, self-hits removed:

| significance cut | isolated | share |
|---|---|---|
| e <= 1e-20 | **127** | 55 % |
| e <= 1e-10 | 94 | 41 % |
| **e <= 1e-5** | **53** | **23 %** |
| e <= 1e-3 | 16 | 7 % |
| e <= 0.01 | 5 | 2 % |

**The threshold decides the answer**, which is itself the finding. At a strict
cut more than half of SINEBase has no relative in it; at a loose cut almost
everything is connected to something.

Connectivity of the rest: **median 2 neighbours, maximum 15**. So SINEBase is
not one connected family — it is mostly small clusters plus a large isolated
tail.

The 53 isolated at e<=1e-5 are in `sinebase_isolated.txt`: `ACar-1`, `AFC-3`,
`Asu1`, `Au`, `Bm1`, `BmSE`, `BraS-I`, `CQ-2`, `Cre-3`, `Cry`, `CucuS-II`,
`DR-2`, `ERI-2`, `EuphS-I`, `FabaS-I/III/V/VI/VII/IX`, `Feilai`, `IscNinDC`,
`Lj-2`, `LF`, `Lm1`, `Mad-1`, `Mad-2`, `NV-2`, `NymS-I`, `OK`, `p-SINE1`,
`PinS-I`, `PoaS-I`, `PoaS-II`, `Ruka`, `SaliS-III`, `SB1`, `SB3`, `SB4`,
`SB12`, and 15 more.

**Bearing on pom_SINE_0:** an unrelated SINE is not unusual. Somewhere between a
quarter and a half of the established database has no significant similarity to
anything else in it, so "no similarity to known SINEs" is a normal property of a
real SINE family, not a reason for doubt.

---

## 2026-09-02 — First prospective results: candidates with no answer key

`candidate_to_aln.py` is new and is the step that never existed: candidate
consensus + genome -> ~100-locus alignments with flanks, in the three views he
uses (top hits / random sample / everything). Plan section 6.1, working.

### AnnoSINE_v2 on three genomes with no SINE library

| genome | phylum | candidates found |
|---|---|---|
| *Pomacea canaliculata* (snail) | Mollusca | **2** |
| *Acanthaster planci* (starfish) | Echinodermata | **2** |
| *Hydra vulgaris* | Cnidaria | **22** (20 tRNA-derived, 1 5S, 1 unknown) |

Against Timema's 55. He read the low counts correctly: snail and starfish are
genuinely SINE-poor, not a failed run — all seven AnnoSINE steps completed and
the candidates found are substantial, not junk.

### Verdicts on the four snail/starfish candidates

| candidate | derived from | genomic copies | top100 | rand100 |
|---|---|---|---|---|
| `pom_SINE_0` | tRNA | **537** | **100.0** clean | 93.7 TRUNCATED_COPIES |
| `pom_SINE_1` | 7SL RNA | **339** | **100.0** clean | 86.2 TRUNCATED_COPIES |
| `aca_SINE_0` | tRNA | 176 | **100.0** | 95.9 SUBFAMILY_NOTE |
| `aca_SINE_1` | 5S rRNA | 11 | **26.3** | 26.3 |

Three accepted, one rejected. **Nothing can check this except his eye** — there
is no curated library, no prior annotation, no RepeatMasker track for these
genomes. Published at `newspecies.html`.

Note the pattern in both snail candidates: **top100 clean at 100.0, rand100
flagged TRUNCATED_COPIES**. The best-scoring hits are full-length; a random
sample includes fragments. That is expected for a real family and is itself
weak evidence these are genuine.

### Two more genomes running

Chosen on his instruction — *"find other genome presumably containing sines but
make sure its not properly annotated/studied"*:

- *Acipenser ruthenus* (sterlet sturgeon) 1.9 Gb
- *Amblyraja radiata* (thorny skate) 2.6 Gb

Both are groups where SINEs are documented in the literature while these
assemblies carry no repeat library.

---

## 2026-09-01 — R11. Flanks in the aligner: a hypothesis of mine, DISPROVED

He asked: *"why you insist on aligned flanks in your examples?"*

The two corpora were built differently. `kit/extend.py` (400 bp) aligns the
element with the consensus and attaches raw flanks afterwards — its own comment
says *"no flank in the aligner at all"*. But `aln_v2`, which every review
alignment came from, put element and flanks through MAFFT together and de-gapped
the flanks afterwards.

**My hypothesis:** letting unalignable flanks into the aligner distorts the
ELEMENT alignment, because MAFFT optimises one score across the whole row.

**Tested and disproved.** Re-aligning the element alone and attaching raw flanks
changes the element alignment barely at all:

| alignment | element gaps, flanks in aligner | element only |
|---|---|---|
| `NEGCHIM__ccr__g1_180seqs` | 0.752 | 0.789 |
| `NEGCHIM__ccr__g3_71seqs` | 0.765 | 0.769 |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | 0.451 | 0.434 |
| `NEGTRUNC5__saq__s2_38seqs` | 0.663 | 0.654 |
| `NEGCHIM__saq__s8_225seqs` | 0.738 | 0.747 |

Two go slightly worse, three slightly better - noise. MAFFT was already
effectively ignoring the flanks. **So this is not the cause of what he is
seeing**, and the question remains open: after justification each copy's flank
is raw sequence butted against the element with padding on the outer side only,
because a FASTA alignment must be rectangular. Whether that padding is what he
reads as "aligned" is unresolved and was asked back to him rather than guessed
at a third time.

`realign.py` is kept on DRAGEN in case element-only rebuilding is wanted for
another reason, but it does not fix this.

---

## 2026-09-01 — R10. The mixture test rebuilt: 5/9 → 9/9, but it now fires 3x more often

### The old measure was invalid, not miscalibrated

Group **length** difference, tested against nine of his judgements:

| | rel_len |
|---|---|
| his mixtures | 0.122, 0.133, 0.138, 0.309 |
| his non-mixtures | 0.089, 0.119, **0.228** |

`s8_225seqs` is "not mixture" at 0.228 — above three of the four real mixtures.
The ranges overlap completely, so **no threshold exists**. Re-fitting was the
wrong response; the measure had to be replaced.

### What replaced it

Group **identity** difference — which is what he actually describes: *"top proper
about half longer similarity sequences ... bottom more discordant shorter ones"*.

| | d_ident |
|---|---|
| his mixtures | 0.128 – 0.215 |
| his non-mixtures | 0.012 – 0.124 |

No overlap. Threshold 0.126.

### The grouping mattered as much as the measure

A first attempt applied that threshold to `subfamily_split`'s grouping and got
7/9. Under that grouping CAS and AluYh9 — both of which he calls mixtures —
score **0.013 and 0.012**, the lowest of all nine. `subfamily_split` is driven by
sequence structure and finds the wrong halves.

Splitting instead on each copy's own **identity, length and coverage** (the SVD
split in `row_groups`) recovers 0.209 and 0.128. **9 of 9.**

### Regression, and one number that needs his eye

| | n | mean | extreme | crossing 50 |
|---|---|---|---|---|
| true negatives | 30 | 1.59 | max 47.0 | 0 |
| true positives | 432 | 99.96 | min 90.1 | 0 |
| curated `tim/` | 42 | 98.29 | min 71.3 | 0 |

Unchanged. **But the flag now fires on 120 of 834 sets (14%), up from 39 (4.7%).**

| class | fires | of | |
|---|---|---|---|
| MIXED10 | 28 | 28 | 100% |
| MIXED30 | 28 | 28 | 100% |
| NEGTRUNC5 | 28 | 28 | 100% |
| HUM | 14 | 64 | 22% |
| TIM | 13 | 42 | 31% |

MIXED10/MIXED30 are *deliberately contaminated*, so firing there is arguably
right — contamination is two things in one set. NEGTRUNC5 at 100% is more
doubtful: he called those *"more like a SINE"*, not mixtures.

**This flag withholds the verdict**, so 120 deferred sets is a large behavioural
change resting on nine judgements. Needs his eye on a NEGTRUNC5 and a MIXED30
before it can be trusted at this rate.

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
