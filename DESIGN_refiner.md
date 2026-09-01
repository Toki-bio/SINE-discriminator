# From scorer to refiner: what it would take

Sergei's question, restated precisely: not adjusting *flank width*, but adjusting
**the element's own edges — inward or outward — to where copy-to-copy similarity
actually stops**, for the normal case where the supplied consensus does not
describe the real element. And: is this ready to become a tool that sits
downstream of SINEderella's de novo scan, AnnoSINE2, RepeatModeler, RepeatMasker
or HiTE, all of which emit *inaccurate preliminary candidates*?

This document is the considered answer, with the measurements behind it.

---

## 1. What exists today, honestly

The tool is a **scorer**, not a refiner. It answers "is this set a SINE family,
and if not, what went wrong". It has been validated hard:

- 43-point empty gap between 30 true negatives (max 47.0) and 432 true
  positives (min 90.1)
- 42/42 on the manually curated `tim/` ground truth whose boundaries were
  independently confirmed on all 24 sides
- 95.6 % of 5,306 human benchmark copies corroborated by RepeatMasker's
  independent hg38 annotation

But every one of those runs was given an alignment somebody else built, with a
consensus somebody else decided. **It never moves a boundary and never writes a
corrected consensus.** That is the gap between what exists and what is being
asked for.

The only boundary operation present is `core_window()`: an **inward-only**
search for a better-supported sub-span, which fires as a *flag* and adjusts the
support measurement, but does not rewrite the consensus or re-align.

---

## 2. The measurement needed already exists — and the naive version is wrong

`profiles.profile()` computes pairwise copy-to-copy identity across three zones
in one consistent metric, with flanks indexed outward from each copy's **own
ungapped edge**. That is the correct coordinate system for edge finding: no
aligner is involved outside the element, so the background is the true one.

I probed it: estimate the element plateau and the genomic background, then find
where the curve crosses the midpoint, scanning outward on each side.

**Left edge: works.** On sets whose consensus is known good, the detected left
edge sits within 1–3 bp of the consensus start (`AluY` +1, `AluSx1` +3,
`SINE_17` +2, curated `t1_1` +1).

**Right edge: the naive version fails, and fails the same way every time.** It
reported shifts of −68 to −84 bp on *every* good set. The cause is not noise:

| set | element body | last 90 bp |
|---|---|---|
| `AluY` | pair identity **0.932**, A+T 0.263 | pair identity **0.266**, A+T 0.633 |
| curated `t1_1` | pair identity **0.960** | pair identity **0.305**, A+T 0.700 |

**Copy-to-copy identity collapses in the 3′ tail because poly-A run lengths vary
between copies.** An identity-based edge finder therefore amputates the poly-A
tail from every SINE it touches — deleting a defining feature of the class.

This is the single most important constraint on the design, and it is the reason
this cannot be built as "find where similarity stops".

### What the right edge actually needs

The 3′ boundary must be found with a **composition-aware rule**, not an identity
rule:

- identity falls **and** A+T rises together → still inside the element, in its
  tail; keep going
- identity falls **and** A+T stays at genomic composition → the element has
  ended

`annotate.py` already detects terminators (`TTTT` strong; `TCTTT`/`TGTTT`/
`TATTT` moderate), simple repeats and A-rich runs, and `profile()` already
returns an `at` track. The ingredients are present; the rule is not written.

---

## 3. The other half: is the consensus wrong, or is there no element?

`NO_ELEMENT` currently conflates two situations that need opposite responses:

- **(a) no element** — the copies do not agree with each other either
- **(b) the consensus is wrong** — the copies agree with each other, but not
  with what they were given

(b) is the *normal* state of an upstream de novo candidate and is recoverable.
(a) is a rejection.

Separating them needs two numbers that are both already computed but never
compared: the copy-to-copy plateau (`pair_id`) and identity to the supplied
consensus (`cons_id`), both against measured background.

Measured, with a coverage requirement on the background estimate:

| set | pair | cons | bg | pair − bg | current verdict |
|---|---|---|---|---|---|
| `POS__saq__s5_5seqs` | 0.964 | 0.973 | 0.264 | **0.699** | 100.0 |
| curated `t1_1` | 0.973 | 0.982 | 0.300 | **0.674** | 100.0 |
| `AluY` | 0.905 | 0.914 | 0.249 | **0.656** | 100.0 |
| `MIR1_Amn` | 0.528 | 0.510 | 0.242 | 0.286 | 100.0 |
| `SINE_43` | 0.501 | 0.502 | 0.285 | **0.216** | **0.0 NO_ELEMENT** |
| `SINE_21` | 0.463 | 0.494 | 0.381 | 0.082 | 78.5 |
| `LmeSINE1c` | 0.367 | 0.429 | 0.216 | 0.151 | 69.2 |
| `NEGRAND__dmo__r00` | 0.355 | 0.425 | 0.261 | 0.094 | 0.0 |
| `NEGRAND__saq__r00` | 0.335 | 0.402 | 0.274 | 0.061 | 0.0 |

Two things follow, and one of them corrects an assumption I made out loud:

1. **`pair − cons` is ~0 everywhere.** I expected sets where copies agree with
   each other far more than with their consensus, and there are none in this
   corpus. So case (b) — a consensus that is simply *wrong* rather than
   mis-bounded — does not occur here. It will occur with upstream tools, whose
   consensuses are far rougher than curated ones, but it is not yet demonstrated
   and must not be designed for as if it were.
2. **`pair − bg` is a real gradient the tool currently ignores.** `SINE_43`
   sits at 0.216 — over twice the `NEGRAND` level (0.06–0.09) — and is being
   rejected outright as `NO_ELEMENT`. There is weak but genuine copy-to-copy
   structure there. Whether that is a degraded family or an artifact is exactly
   the judgement a refiner should surface rather than silently reject.

A methodological note: my first pass at this reported `MIR1_Amn` background as
0.573 — higher than its element — which would have been a serious finding. It
was wrong: its right flanks are all under 15 bp and its left flanks median
41 bp, so the far-flank window had almost no copies in it. **The same
small-sample inflation the scorer was already fixed for.** With a coverage
requirement the background is 0.242, matching the scorer. Any edge-finding code
must carry that guard from the start.

---

## 4. What a refiner would actually look like

```
candidate consensus            (AnnoSINE2 / RepeatModeler / HiTE / RepeatMasker /
                                SINEderella de novo)
  │
  ├─ search the genome, take N copies with 400 bp flanks
  ├─ align, consensus as row 1
  │
  ├─ MEASURE   cliff, decay, per-position identity + composition
  │
  ├─ ADJUST    left edge  : identity crossing        (works today)
  │            right edge : identity + A+T rule      (NOT written)
  │            inward trim: core_window()            (exists, flag only)
  │            outward    : decay says it continues  (detected, not acted on)
  │
  ├─ CURATE    prune contamination                   (exists)
  │            split mixtures and recurse            (detected, not acted on)
  │
  ├─ RE-DERIVE consensus from the corrected span
  └─ REPEAT until the boundary stops moving
       │
       └─ emit: corrected consensus + verdict + evidence + copy list
```

Of those eight operations the tool has **three** working (measure, prune,
inward-trim-as-flag), **two** detected but not acted on (outward extension,
mixture splitting), and **three** entirely missing (right-edge rule, consensus
re-derivation, the iteration loop).

---

## 5. Readiness, stated separately

**As a judgement engine: ~85 %.** The hard, novel part — deciding whether a
locus set is a family, with named sub-cases instead of a binary — is built and
tested against synthetic truth, two independent annotations, an expert eye
review, and 1.9 M RepeatMasker calls.

**As the tool Sergei described: ~40 %.** The plumbing that makes it a pipeline
stage rather than a study does not exist:

- no input contract for a bare consensus + genome (upstream tools emit
  consensuses, not alignments)
- no consensus re-derivation, so nothing improves
- no iteration, so nothing converges
- no right-edge rule, and the naive one is actively harmful

None of that is research risk. It is engineering on top of a validated core,
with one genuine unknown — the composition-aware 3′ rule — which is a day of
work to prototype and can be validated immediately against 64 human Dfam
families whose true 3′ ends are known.

---

## 6. Build order, riskiest-assumption-first

1. **The 3′ edge rule.** The only real unknown. Validate against Dfam:
   the corrected right edge must land within a few bp of the curated consensus
   end on all 64 human families, without eating the poly-A tail. If this does
   not work, the refiner cannot be built as conceived, so it goes first.
2. **The wrapper**: consensus + genome → loci → 400 bp extraction → alignment.
   Mostly assembly of existing pieces (`bench_extend.py`, blastn, bedtools,
   MAFFT); SINEderella's step1/step2 already do a version of this.
3. **Consensus re-derivation + the loop**, with an explicit convergence rule and
   a hard iteration cap. Success criterion: starting from a deliberately
   mis-bounded consensus, converge to the curated one.
4. **Mixture splitting**, recursing on each part. `HETEROGENEOUS_SELECTION`
   already identifies the sets and the groups.
5. **Outward extension** for `FRAGMENT_OF_LONGER` / `ELEMENT_CONTINUES`.

Step 1 is the go/no-go. Steps 2–5 are only worth starting once it holds.

---

## 7. What would make this fail

Stated so it can be checked rather than discovered late:

- **The 3′ rule may not generalise.** It will be tuned on Alu and MIR tails. A
  SINE family with no poly-A tail, or with a tail that is not A-rich, breaks it.
  The three-phylum prospective run (Mollusca, Echinodermata, Cnidaria) is the
  test, because none of those families is known.
- **Iteration can drift.** Re-deriving a consensus from copies selected by that
  consensus is a feedback loop and can converge on a subfamily, or on nothing.
  It needs a fixed reference point — the original candidate — and a cap.
- **Two thresholds in the current scorer are overfit** (`HETEROGENEOUS_SELECTION`
  at 0.12 from three examples; the overextension tail bound from one). A refiner
  that acts on them will amplify that, where a scorer that only reports them
  does not.
