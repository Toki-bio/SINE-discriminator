# SINE discriminator — the plan

---

## 1. The problem

Across 72 alignments you have described the **properties of an alignment that
matter for telling a SINE from a non-SINE**. Your words are in `calls.tsv`.

The tool collapses all of them into one number from 0 to 100, and that number
discards most of what you described. Proof: four alignments where you described
clearly different properties all come out as exactly **100.0** — one whose
consensus is too long, one needing internal subgrouping, one with an unstable
pre-tail region, and one that is simply clean.

And four where the number is wrong outright: ERI e1-4 → 0.2, `AluSx_short_` →
0.1, `LmeSINE1c` → 69.2, and the lowest alignment you called clean → 68.3.

---

## 2. The properties you described

Collected from your own words. **This list is the specification.**

| property | what you said |
|---|---|
| flanks alignable or not | "its flanks are good — not alignable" |
| left edge sharpness | "right end is fine but left one is really bad" |
| right edge sharpness | "right end is a bit wobbly unstable" |
| where the right edge actually is | "if its right edge is defined at gggagat… it becomes quite clear" |
| copies agree with each other | throughout — the element is visible across copies |
| divergence level | "pure clear SINE with high divergence" — divergence is not a fault |
| tail behaviour | "long insertions in the tail region but if gaps are removed it becomes consolidated at right edge" |
| a subset of copies differs | "11 lower sequences need manual attention"; "top proper about half … then bottom more discordant shorter ones" |
| consensus length vs the real element | "consensus needs refinement and shortening" |
| enough copies to judge at all | "very weak signal over too few copies"; one copy → "wtf???" |
| presentation quality | "not presented properly with aligned gappy flanks"; "needs … degapped left flank" |
| flank uniqueness genome-wide | "requires post-processing with proving uniqness of at least some flanks on whole-genome level" |
| internal instability in a region | "pre-tail region containing an island of instability" |

Two of these are not statements about the alignment but instructions about what
to do next: **split it and re-run each part**, and **get more copies before
deciding**.

---

## 3. What the tool should report

Each property, measured, with its value — not a single number.

Which properties are already measurable, and which are not, is what §4 finds
out.

---

## 4. Testing

Everything is checked against your 72 judgements, never against sequences I
generated.

### Test 1 — which measurements capture which property *(DONE)*

19 measurements computed on all 72 alignments, in `test1_vars.tsv`. Three of
them already carry information the score throws away:

- **Edge sharpness** (`edge_drop` — how far copy-to-copy similarity falls at the
  element boundary) separates alignments you called clean (0.52) from ones you
  called grey (0.15), badly presented (0.18) or wrong-consensus (0.13) — even
  though all of those score 89–100 and are identical by score.
- **Flank similarity** (`pair_bg`) explains ERI e1-4: its flanks are 0.57
  similar against a normal 0.25. The copies share flanking sequence. The score
  says "no element"; the truth is "the flanks are not independent".
- **Coverage minimum** (`cov_min` = 0.000) marks every alignment you called
  badly presented or unusable.

### Test 2 — what is still not measured

For each property in §2 with no measurement behind it, look at the alignments
where you described it and work out what to compute. The obvious gaps:

- **a subset of copies differing** — everything current is per-column; you were
  describing groups of rows
- **where the right edge actually is** — the naive rule fails, because
  copy-to-copy identity collapses in the poly-A tail (AluY: 0.93 in the body,
  0.27 in the last 90 bp, while A+T rises 0.26 → 0.63). Needs a rule that
  watches base composition, not identity alone.
- **internal instability in a region** — a local measure, not a global one

### Test 3 — build them and re-measure

**Gate:** a new measurement earns its place only if it captures a property
nothing captured before, without spoiling one that already worked.

### Test 4 — more judgements, only where needed

Nine properties rest on a single alignment each. Where a threshold cannot be set
from one example, show you the alignments that would settle it — and only those.

### Test 5 — check on alignments it was not tuned on

**Gate:** if it only works where it was tuned, it is not a tool.

---

## 5. Then: the rules

For each property — the measurement, its range, how much to trust it given how
many examples support it, and what it means for what to do next. Plus, stated
openly, the properties still not measurable.

---

## 6. Then: the tool

1. **One command**: candidate consensus + genome → hits → several ~100-locus
   alignments with flanks (random sample, best hits) → the properties, measured,
   per alignment. This automates your steps 1–4.
2. **The actions that can be automated**: split a mixture and re-run each part;
   shorten an over-long consensus; re-extract with longer flanks; prune
   contamination; de-gap flanks for display.
3. **The subfam peel loop** — written, never run. Group consensus built from the
   original member sequences, not from the consensi-of-50s.
4. **Repeat until stable**: refine the consensus, re-search, re-measure.
5. **Results page**: per run — what flowed through, the families found with
   their measured properties, a picture per family of identity and A+T along the
   sequence with supplied and corrected edges marked, what was left over, and
   provenance.

---

## 7. Order

| # | step | gate |
|---|---|---|
| 1 | Test 1 | **done** |
| 2 | Test 2 — name what is missing | — |
| 3 | Test 3 — build and re-measure | must capture something previously uncaptured |
| 4 | Test 4 — more judgements where thin | enough examples to set each threshold |
| 5 | Test 5 — unseen alignments | matches your judgement where it was not tuned |
| 6 | the rules | — |
| 7 | the tool | rules validated |

Running in parallel, no gate: AnnoSINE2 on snail, starfish and hydra — three
genomes with no SINE library, the first test on families nobody has described.

---

## 8. Rules of work

- Your judgements are the standard. The curated `tim/` set is the standard.
  **Human Dfam is not** — it contains fragments, an ancestral gene, mixtures and
  single-copy entries.
- Record your words exactly as written, before interpreting them.
- Compute on the servers.
- **Do not invent category names for what you described.** Measure the
  properties.
