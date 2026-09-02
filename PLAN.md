# SINE discriminator — the plan

*Updated 2026-09-02. Status against every step is marked; nothing below is
aspirational unless it says so.*

---

## 1. The problem

Across your judgements — now **113 recorded in `calls.tsv`**, 91 on the Tal and
curated sets plus 22 on hydra — you have described the **properties of an
alignment that matter for telling a SINE from a non-SINE**. Your words are the
specification.

The original complaint was that the tool collapsed all of them into one number
and threw the rest away: four alignments with clearly different properties all
came out at exactly **100.0**.

**That is fixed.** The tool now reports every reason that fires, with its
measurement and what to do about it, and three of those reasons *cap* the score
rather than decorate it. On your 22 hydra judgements it agrees on 20, is
undecided on 1 (`hyd_SINE_2`, too few copies to assess), and disagrees on 1
(`hyd_SINE_5`).

Two faults found later were worse than the scoring problem, because they were
upstream of everything:

- **The consensus was never a consensus.** AnnoSINE hands over a single genomic
  locus per family. Every search, alignment and score was made against one copy.
  `hyd_SINE_0`'s covered 79 bp of a 208 bp element, and the missing 130 bp was
  measured as flank — which is what made the tool report no element at all on a
  family you called a good SINE.
- **Nothing ever chose an orientation.** About half of every genome's families
  were reverse-complemented, because blastn reports both strands and MAFFT
  `--adjustdirection` only makes copies agree with each other.

Both are now fixed for every genome, not one alignment.

---

## 2. The properties you described

Collected from your own words. **This list is the specification.** The right
column is what measures it today.

| property | what you said | measured by |
|---|---|---|
| flanks alignable or not | "its flanks are good — not alignable" | `flank_bg`, flank decay |
| left edge sharpness | "right end is fine but left one is really bad" | **not measured** |
| right edge sharpness | "right end is a bit wobbly unstable" | **not measured** |
| where the right edge actually is | "if its right edge is defined at gggagat… it becomes quite clear" | core window, composition-aware 3′ rule |
| copies agree with each other | the element visible across copies | `id_all`, `n_core` |
| divergence level | "pure clear SINE with high divergence" | reported, never penalised |
| tail behaviour | "long insertions in the tail region…" | composition-aware 3′ rule |
| a subset of copies differs | "top proper about half … then bottom more discordant" | `SUBFAMILY_NOTE` + region consistency |
| consensus length vs the real element | "consensus needs refinement and shortening" | `CONSENSUS_OVEREXTENDED` / `UNDEREXTENDED`, and the two-stage rebuild |
| enough copies to judge at all | "very weak signal over too few copies" | `assessable: false` |
| presentation quality | "needs … degapped left flank" | justify + trim |
| flank uniqueness genome-wide | "proving uniqness of at least some flanks on whole-genome level" | `FLANK_ISLANDS` |
| internal instability in a region | "pre-tail region containing an island of instability" | `patch2d` |
| **microsatellite content** | "combination of microsatellites to me, not sine — needs new criteria on microsatellite content?" | `MICROSATELLITE_ELEMENT` / `_FLANK` |
| **mosaic vs subfamilies** | "inconsistent mixture of conserved and discordant columns/spots" | region consistency (reported, cannot reject) |
| **orientation** | "many of those you gave me need reverse complement" | box B position + A-rich tail |
| **consensus fullness** | "make sure that consensuses match real data" | two-stage build |

Two entries are instructions rather than statements: **split it and re-run each
part**, and **get more copies before deciding**. Both exist as actions now.

---

## 3. What the tool should report — *done*

Each property, measured, with its value, not a single number. Fifteen reason
codes, each with the number behind it and what to do next. Three of them
(`MICROSATELLITE_ELEMENT`, `SMALL_CORE`, and a small core fraction) cap the
score at 45 so a number can never read as acceptance while the text beneath says
the evidence is absent.

---

## 4. Testing

Everything checked against your judgements and the curated sets, never against
sequences I generated.

### Test 1 — which measurements capture which property *(done)*

19 measurements over the judged alignments. Edge sharpness, flank similarity and
coverage minimum each carried information the score discarded.

### Test 2 — what is still not measured *(done, and it named the gaps)*

Three gaps identified and since closed: groups of rows rather than columns
(`patch2d`, region consistency), the real 3′ edge under a poly-A tail
(composition-aware rule), and local instability. One gap named and **still
open**: left/right edge sharpness as separate properties.

### Test 3 — build them and re-measure *(done)*

**Gate: a new measurement earns its place only if it captures something nothing
captured before, without spoiling what worked.** Passed by `patch2d`,
`FLANK_ISLANDS`, microsatellite content and the consensus-length rules. Failed
by four attempts at edge sharpness — recorded as failures in `FINDINGS.md`
rather than shipped.

### Test 4 — more judgements, only where needed *(done for hydra)*

22 hydra judgements settled the microsatellite thresholds, which no labelled
class could have set: every class in the 673-set corpus has a median element
microsatellite content of **0.000**.

### Test 5 — check on alignments it was not tuned on *(done)*

**Gate: if it only works where it was tuned, it is not a tool.**

- 673-set labelled corpus: true positives min 81.5, true negatives max 0.2
- human Alu: 55 of 64 accepted, 9 not assessable, **none rejected**
- Timema: 48 accepted, 5 rejected with reasons named
- 15 genomes with no SINE library

---

## 5. The rules — *done, in `ALGORITHM_NOTES.md`*

Twenty-one numbered rules, each with the evidence behind it. The recurring one,
broken four separate times in four places: **an absent measurement must never be
scored as guilt.**

---

## 6. The tool

1. **One command** — candidate consensus + genome → hits → ~100-locus alignments
   with flanks → measured properties. `run_species.py`. **Done.**
2. **Automated actions** — prune contamination, de-gap and trim flanks for
   display, shorten an over-long consensus, extend a short one. **Done.**
3. **The subfam peel loop** — group consensus built from the original member
   sequences, not from consensi-of-50s. **Still not run.**
4. **Repeat until stable** — refine, re-search, re-measure. **Done**, and the
   shape of it matters: stage A bootstrap once, then `conse` over the top 100 by
   bitscore once. Iterating stage B creeps outward, 203 → 252 → 269 bp.
5. **Results page** — per genome, per family, with the measured properties and
   the alignments. **Done**, and now carries a build timestamp.

---

## 7. Order

| # | step | state |
|---|---|---|
| 1 | Test 1 | done |
| 2 | Test 2 — name what is missing | done |
| 3 | Test 3 — build and re-measure | done; edge sharpness still open |
| 4 | Test 4 — judgements where thin | done for hydra |
| 5 | Test 5 — unseen alignments | done |
| 6 | the rules | done |
| 7 | the tool | done except the subfam peel loop |

**Now running:** every genome rebuilt with a full-length oriented consensus, and
AnnoSINE on the remaining scorpion species.

---

## 8. What is still wrong

- **`hyd_SINE_5`** — "weak right end, mosaic left end". Five measurements tried,
  none separates it from families you accept. The property is real and unmeasured.
- **Left and right edge sharpness as separate properties** — never measured.
- **The subfam peel loop** — written, never run.
- **Flank islands cannot be interpreted** — the measure works, but whether an
  island means a shared larger repeat, an insertion preference, or an assembly
  artefact needs the island sequences searched against the genome.
- **`ERI__eri__e2-3`** — the tool now scores it 45.0 where you lean SINE with
  caution. Defensible as "cannot assert", but it crosses the line on a set you
  accept.

---

## 9. Rules of work

- Your judgements are the standard. The curated `tim/` set is the standard.
  **Human Dfam is not.**
- Record your words exactly as written, before interpreting them.
- Compute on the servers.
- **Do not invent category names for what you described.** Measure the properties.
- Save findings after every exchange: current work in `FINDINGS.md`, general
  rules in `ALGORITHM_NOTES.md`.
- Use your tools rather than reinventing them — `conse`, `SINE_consensus`,
  `SubFam`, `sine_scan.sh`, `extract_alignments.sh`. Every time I have written my
  own version of one of these it has been worse.
