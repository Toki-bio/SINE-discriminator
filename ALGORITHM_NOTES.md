# Algorithm notes — general findings

What holds across cases, as opposed to what happened on a particular day. The
day-by-day record is `FINDINGS.md`; his method is `METHOD.md`; this file is what
the algorithm should be built on.

---

## 1. A rejection has several independent reasons, and they must be reported separately

`aca_SINE_1` (starfish) is the clearest case so far. Sergei's reading: *"bad and
it includes mosaicism ... not only mosaic, but also weak identity, too long
identity with bad flanks"*. All confirmed by measurement:

| reason | measured | comparison |
|---|---|---|
| weak identity | id_all **0.552** | 0.72–0.83 in the accepted ones |
| **flanks not independent** | flank background **0.450** | 0.27–0.31 in the accepted ones |
| almost no boundary | cliff **0.102** | 0.45–0.53 in the accepted ones |
| mosaic | patch2d **0.143** | 0.061–0.097 in the accepted ones |

The tool scored it **26.3 with no flags at all** — the right number for the
wrong reason, and no explanation. **A verdict without its reasons is not the
answer he needs**, and this is the same defect as four of his properties all
scoring exactly 100.

**Rule: every reason that fires should be reported, not just the one that
happened to move the score most.**

## 2. Flank background is a first-class signal, not a normalisation constant

`aca_SINE_1`'s flank background of 0.450 against a genomic 0.25 says its copies
sit in related contexts — they are not independent insertions. That is one of
his decisive criteria (*"its flanks are good — not alignable"*), and it happened
to be visible only because the cliff computation needs it.

It should be reported in its own right whenever it exceeds background by a clear
margin, with the interpretation attached: satellite, segmental duplication, or
copies inside one larger repeat.

## 3. The mosaic is row-wise WITHIN column, which is why every earlier measure missed it

Neither per-column nor per-copy statistics detect it:

- column jumpiness: mosaic 0.261, legitimate candidate 0.239 — no separation
- block structure (runs vs chance): the *legitimate* set is blockier — backwards
- mean identity: says "less conserved", which is a different statement

What responds is **per-copy variation across regions** (`patch2d`): for each
copy, its concordance with the consensus in each 20 bp window, then how much
that varies *within* that copy. A clean family is uniformly concordant; a mosaic
has copies that are concordant in some regions and discordant in others, with
the regions differing between copies.

Ordered all four new candidates exactly as he judged them (0.061 / 0.097 /
0.143 / 0.150). Four points, so not yet a threshold.

## 4. AnnoSINE copy counts are a floor, not an estimate

On the two snail candidates, AnnoSINE reported 537 and 339 genomic copies. A
full SINEderella run on the same consensuses found **6,016 and 2,848** — 11x and
8x more.

So a low AnnoSINE copy count is not evidence against a candidate, and any
threshold on copy number must be applied after a proper search, never to the
discovery tool's own count.

## 5. Assignment leak transfers across species as a real-family signal

In Timema, the candidates Sergei judged real had leak 0.00–0.18 % while the ones
he judged noisy ran 65–98 %. Both snail candidates come back at **0.00 %**.

Leak does not answer sine/not-sine on its own — he was explicit that subfamily
ambiguity is a different question — but zero leak across thousands of copies is
strong evidence the locus set is coherent.

## 6. Being unrelated to known elements is normal, not suspicious

All-vs-all blastn over SINEBase (230 sequences): **23 % isolated at e ≤ 1e-5,
55 % at e ≤ 1e-20**, median 2 neighbours, maximum 15. SINEBase is small clusters
plus a large isolated tail, not one connected family.

So "no similarity to known SINEs" is a property a real SINE family routinely
has. It must never count against a candidate.

## 7. A SINE–LINE partnership shows at the LINE's 3′ terminus or not at all

Searching a candidate against 1,784 L1 and 7,117 other LINE entries produced
hits of 14–28 bp at e ≥ 0.14, all sitting **1,200–4,600 bp inside** the LINEs.
A genuine partnership puts the SINE's tail on the LINE's last ~50 bp.

**Test the position of the hit, not its existence.** Short internal hits to long
LINEs are background and will always be found.

## 8. Structural evidence is per-copy, and its strength varies

For the two snail candidates:

| | A-box | B-box | TSD |
|---|---|---|---|
| `pom_SINE_0` | 2 mismatches | **0 mismatches** | **86 % of copies** |
| `pom_SINE_1` | 3 mismatches | **0 mismatches** | **35 % of copies** |

A perfect B-box is the hardest part to acquire by chance, so both have an intact
pol-III promoter. But TSD fraction separates them sharply, and it agrees with
every other measure of which is the stronger candidate. **TSD fraction is worth
more than TSD presence.**

## 9. Display faults are not data faults, and the two must not be confused

Flanks that look bad have had three different causes so far, only one of which
was real:

1. **Column width set by one outlier copy** — one copy with a 125 bp flank pads
   seventy others with gaps. Fixed by capping the width to what the copies
   actually use.
2. **A percentile cap is not enough on skewed lengths**, and a constant floor
   cannot help when copies carry ~10 bases. Fixed by choosing the width that
   keeps the panel under 25 % gaps, with the floor following the copies' median.
3. **The copies genuinely have no flanks** — no presentation can fix this, and
   the tool should say so instead of showing a wall of dashes.

A hypothesis that flanks in the aligner distort the element alignment was
**tested and disproved**: re-aligning the element alone changes element gaps by
at most 0.04, in both directions.

## 10. Flank similarity ISLANDS — a signal the score is blind to

Sergei on `aca_SINE_0`: *"aligned left flank has very faint but non-random
islands of similarity"*, and then *"these islands can be far off on the flanks
(both directions)"*.

He is right, and my first reading was wrong: I dismissed the far-upstream ones
as low-coverage artifacts. With a null that controls for how many copies are
present at each column, they are strongly significant — on `aca_SINE_0`,
**49 patches covering 426 columns, up to z = 69**, on both sides and hundreds
of bases out from the element. The coverage control is the whole difference:
without it a column where 12 copies of 100 happen to agree looks like noise
alongside one where 98 do.

Scanned uniformly over all 50 measurable alignments of the prospective corpus
(z > 8 sustained over at least 6 columns, against a null that controls for
coverage):

| alignment | island columns | flank_bg | score |
|---|---|---|---|
| `hyd_SINE_9` | 739 | 0.914 | 0.0 |
| **`aca_SINE_0`** | **426** | **0.271** | **100.0** |
| `aca_SINE_1` (rejected) | 424 | 0.450 | 26.3 |
| `hyd_SINE_17` | 303 | 0.543 | 100.0 |
| `hyd_SINE_0` | 139 | 0.614 | 0.0 |
| everything else | <= 75 | ~0.30 | |

**`aca_SINE_0` scores a perfect 100 and has as much localised flank similarity
as the set he rejected.** The score cannot see it.

### Islands must be measured on the LONG flank, never the display alignment

Same shape as the flank-decay rule, and it bites hard here:

| `aca_SINE_0` measured on | island columns |
|---|---|
| the full 400 bp flank | **426** |
| the trimmed display alignment | **6** |

Trimming cuts the flank panel to the width the copies actually fill. That is
right for judging a boundary and it is what fixed the walls of dashes, but it
throws away exactly the far-out columns Sergei said the islands sit in. Two
flank geometries, two purposes, and they must never be taken from one file.

### Only one of the two high-island sets is invisible to the score

I first wrote that `hyd_SINE_9` was a second case. It is not, and the
distinction is the whole point:

| | island columns | flank_bg | score | seen? |
|---|---|---|---|---|
| `hyd_SINE_9` | 739 | **0.914** | 0.0, SHARED_FLANKS | **yes** |
| `aca_SINE_0` | 426 | **0.271** | 100.0, no flag | **no** |

Taking "invisible" to mean over 150 island columns with a flank average below
0.40, `aca_SINE_0` is the **only** set in the whole 50-alignment corpus.

`hyd_SINE_9`'s flanks are similar *everywhere*, so an average catches it and
the existing rule already rejects it. `aca_SINE_0`'s flanks are ordinary on
average and similar only in patches. **The islands measure earns its place
precisely because it separates these two, which the average cannot.**

### Why the existing flank measure misses it

`flank_bg` is a single average over the whole flank. `aca_SINE_0`'s is 0.271 -
indistinguishable from background - because a few hundred island columns are
diluted by a thousand ordinary ones. **An average cannot detect a localised
signal.** The islands have to be found as islands.

### Validated against his own calls, and it is a CAUTION signal

Run over the whole labelled corpus (673 sets), measuring the **fraction** of
flank inside an island rather than the raw count. Of the 55 sets he has judged:

| his call | n | island fraction |
|---|---|---|
| plain `SINE` | 31 | **<= 0.067** |
| the three flank-caution calls | 3 | **0.19 - 0.37** |
| `NOT_SINE` / `UNUSABLE` / `MOSAIC` | 7 | **>= 0.50** |

Three bands with wide empty gaps. The middle band is exactly the three sets he
asked for flank checks on and nothing else:

- `ERI__eri__e2-4` - *"requires post-processing with proving uniqness of at
  least some flanks on whole-genome level"* (0.195, flank_bg 0.281)
- `ERI__eri__e2-3` - *"with caution and need manual reinspection"* (0.194, 0.371)
- `NEGSEGDUP__eri__r00` - *"unclear situation with left end"* (0.368, 0.326)

**All three have an ordinary flank average.** `NEGSEGDUP__eri__r00` is a
labelled segmental duplication whose flank_bg is 0.326 - the average cannot
tell it from a clean family, and only the islands can.

By class the measure lands on segmental duplications (0.871), LINE ORFs (0.819)
and satellites (0.523), against 0.025 for POS. Every one of the 11 sets with
flank_bg >= 0.40 also has island fraction >= 0.15, so nothing is lost by using
it. It fires on **14 of 273** sets, all in those classes plus hedgehog ERI, and
on **none** of POS, MIXED10, MIXED30, NEGJITTER, NEGSPLICE, NEGTRUNC5, NEGRAND
or NEGCHIM.

**It must be reported, never subtracted.** He still leans SINE on two of the
three calibration sets; penalising would turn `e2-4` (97.4, which he accepts
with a caveat) into a rejection. `verdict.py` emits `FLANK_ISLANDS` at fraction
>= 0.10 and changes no score - verified, 0 of 273 scores moved.

### Where the prospective candidates land

`aca_SINE_0` reads **0.171** at flank_bg 0.271: between his plain-SINE ceiling
and his caution floor, 2.5x the highest plain SINE he judged. It is the only
one of the 52 new-species alignments above 0.10 with an ordinary flank average.

And the signal is in the **top hits only** - `top100` 0.171, `rand100` under
0.03, while `hyd_SINE_9` is 0.331 on both views. So it is not a family sitting
in a shared context; it is a family whose most similar copies do. That points
at a subset inside a larger duplication, and the check is to search those
particular flanks against the starfish genome.

### What it might mean, unresolved

**Still to do:** pull `aca_SINE_0`'s top-hit flanks and search them against the
starfish genome, which is the only thing that will say which of the four
explanations above is right.

## 11. A background measured next to the element is not a background

The single most damaging class of bug in this tool, hit three separate times:
comparing the element against something that contains the element.

`hyd_SINE_0` is the clearest case. Its consensus covers 79 bp of a 208 bp
element, so ~150 bp of SINE sat in what the code called flank. Flank identity
came out **0.614 against an element identity of 0.560** — the background above
the signal — and the tool concluded there was no element, on a family Sergei
reads as a good SINE and RepBase calls SINE2-2B_HM.

**Rule: take the background from far enough out that the element cannot reach
it.** The flank-decay profile already walks outward in 25 bp steps, so its far
end is unrelated DNA wherever the consensus happens to stop. Both the cliff and
the support threshold now use that value, and nothing else in the corpus moved.

Same shape as the earlier failures: the cliff measured over core copies only
(circular), and the flank background measured on element-only alignments
(1.000). Whenever a denominator can contain the numerator, it will.

## 12. Every measurement has its own geometry, and they must not share a file

Three measurements now need three different views of the same loci:

| measurement | needs | breaks if given the other |
|---|---|---|
| boundary by eye | flank trimmed to what copies fill | a wall of dashes |
| flank islands | full 400 bp | 426 island columns become 6 |
| flank decay, microsatellite | full 400 bp | far field never reached |

The trimmed alignment is a **display artefact**. Nothing should ever be measured
on it. This has now caused two separate wrong answers, so it is a rule rather
than an observation.

## 13. Microsatellite content is an axis nothing else in the corpus tests

Sergei asked for it by name on `hyd_SINE_7` (*"combination of microsatellites to
me, not sine - not caught by your filters!"*) and `hyd_SINE_17` (*"surrounded by
long 2-nt microsatellites"*).

Measured as the fraction of ungapped sequence inside a tandem repeat of period
1–6, median over copies, separately for element and flank:

| | element | flank |
|---|---|---|
| `hyd_SINE_7` | **0.517** | 0.000 |
| `hyd_SINE_17` | **0.229** | **0.253** |
| hydra SINEs he accepted | <= 0.092 | 0.000 |
| **every class of the 673-set corpus** | **0.000** | 0.000 |

The last row is the point: POS, NEGSAT, NEGSEGDUP, NEGLINEORF and all the
synthetic classes sit at zero. **No labelled class exercises this property**, so
no threshold could have been fitted from the corpus — the bands exist only
because he judged the hydra sets. A curated negative set only covers the failure
modes someone thought to build.

Two numbers, not one: repeat *in the element* means the consensus collects loci
by repeat content rather than descent; repeat *in the flank* means the copies
sit in tracts where flanks cannot show independence and mapping is unreliable.

## 14. A reason that cannot change the verdict is decoration

Note 1 said every reason that fires should be reported. The inverse turned out to
be just as wrong: three reasons fired and the score ignored them.

- `hyd_SINE_7`: 52 % simple repeat, scored **77**
- `hyd_SINE_6`: 25 of 100 copies in the core, scored **93**
- `hyd_SINE_2`: 8 copies, scored **95**, its own flag text reading *"this score
  is an impression rather than a measurement"*

He rejected all three. A number that reads as acceptance while the text beneath
it says the evidence is not there is worse than no number at all.

**Rule: some reasons cap the score rather than nudging it.** Not a subtraction —
there is no principled amount to subtract — a cap at 45, just under the
acceptance line, saying the question cannot be answered yes on this evidence.
Currently `MICROSATELLITE_ELEMENT`, `SMALL_CORE`, `INSUFFICIENT_COPIES`.

## 15. When the tool can fix the fault, it should fix it

`CONSENSUS_UNDEREXTENDED` told the user to extend the consensus and re-run. It
now does it:

| stage | length | score |
|---|---|---|
| AnnoSINE seed | 110 bp | 0.0 NO_ELEMENT |
| extended by the decay distance | 310 bp | 0.2, flagged **too long** |
| trimmed to the supported window | **211 bp** | **100.0 clean** |

RepBase's SINE2-2B_HM is 208 bp. The loop converges because the two consensus
rules are mirrors: extension is bounded by where similarity reaches background,
and over-extension is caught by the window the copies actually support. Neither
rule needed changing to make the loop work.

This is Sergei's second starting situation — "a library needing adjustment" —
and it is the one place where the tool can do the adjustment rather than hand
back advice.

## 16. Evidence against and evidence missing are different answers

The score answers "how much evidence that this is a SINE family". Three
situations look alike in a single number and must not:

| | example | right answer |
|---|---|---|
| evidence against | 52 % of the element is simple repeat | reject |
| evidence missing | 8 copies, no flanks | **cannot say** |
| evidence for | 93 of 100 copies in a clean core | accept |

Collapsing the middle into a low score put `HUM__hum__AluYb9` - a real Alu
subfamily - at 45.0 because it has few copies. The scorer now carries
`assessable: false` for `INSUFFICIENT_COPIES` and `NO_FLANKS_PRESENT`, and the
page shows it as a third state rather than a rejection. Human rejections went
from 1 to 0.

This is the same rule as "an absent measurement must never be scored as guilt",
which has now been broken four separate times in four different places. It is
worth checking every new reason against it before adding it.

## 17. Unmeasured is not zero, and it is not the maximum either

Three bugs in flank decay, found together when the profile was regenerated over
the whole corpus for the first time, all one mistake:

- a satellite was called on **distance alone**, so a set at 0.255 flank identity
  - background - was called a satellite because its noisy profile never
  confidently dropped
- an unmeasurable distance **defaulted to the 400 bp maximum**, which collapsed
  the uniqueness term to zero
- the distance was **claimed beyond where the profile was readable**: copies
  that carry only 50 bp of flank were credited with 400 bp of similarity

Together these rejected 10 sets including a positive and a set he calls a good
sine, dropping the minimum POS score from 90.1 to 0.2.

**Rule: a measurement may only claim the range it actually covered.** Where the
data runs out, the answer is "not measured", which feeds nothing. And a
classification that depends on a distance must also require the thing being
measured to be present - distance means nothing until there is elevated
similarity to decay from.

## 18. A consensus threshold is a fraction of the sequences, never a count

Building family consensuses from his de novo scan, I reused SubFam's own call,
`cons -plurality 18`. Six of the eleven consensuses then got **zero blastn hits
against the very genome their loci came from**, and none scored above 50.

`-plurality 18` asks for eighteen sequences to agree. Inside SubFam that runs on
chunks of about fifty, where eighteen is a little over a third - a sensible
majority. Applied to a family of twenty it demands near-unanimity, so most
columns come out N; stripping the N's then splices non-adjacent columns into a
chimera that matches nothing.

Rebuilt as a plain column majority - the most common base wherever at least half
the copies have one - the same loci give 47 usable consensuses, all SINE-sized.

**Any parameter borrowed from another stage has to be re-read as a fraction of
whatever it is now being applied to.** The same trap as measuring a background
next to the element: a number that was right in one context, carried unchanged
into another.

## 19. Grouping a de novo locus set: use the provenance, not the clustering

SubFam on the 3,408 loci from the scorpion de novo scan gives 69 chunk
consensuses whose median pairwise identity is **0.205** - near random. That is
not SubFam failing. It chunks by input order after a rough reorder, so it works
when the input is already one family, which is exactly how SINEderella uses it.
A de novo sweep is the opposite case: those loci came from 792 different query
families, so every 50-locus chunk spans many of them.

For a de novo set the grouping is already in the data - **the query that found
each locus**. Grouping by query name gives 47 families of 20+ loci, all
SINE-sized once their consensuses are built properly.

## 20. A discovery tool's candidate list is redundant, and it must be collapsed first

AnnoSINE proposed **160** candidates for zebrafish. All-vs-all blastn gives 424
redundant pairs, and single-link clustering at 80 % identity over 60 % of the
shorter sequence collapses them to **87 families** — one cluster alone holding
**55 near-identical seeds**, with others of 12 and 8, and 83 singletons.

So **46 % of that candidate list is the same families proposed over and over.**
Sergei suspected it on sight: *"maybe zebrafish 160 candidates can be
deduplicated? too much redundancy suspected."*

This matters beyond compute. Every per-class count, every "n of N accepted", and
any threshold fitted on a candidate list is distorted by whichever family
happened to be proposed 55 times. Deduplicate before scoring, and pick the
representative by genomic copy count rather than by order.

It also cuts real time: each candidate costs a blastn against the whole genome,
so on a 1.7 Gb assembly the redundant 73 were about an hour of wasted search.

Related to note 4 — a discovery tool's copy counts are a floor, not an estimate —
and the same lesson from the other side: its candidate *list* is not a family
list either, in both directions.

## 21. A process-matching pattern that can match itself will kill the wrong thing

`pgrep -f "script.py"` matches its own command line, so it reports a dead job as
alive. Already recorded, already fixed with `pgrep -f "script[.]py"`.

The same trap in the other direction is worse: **`pkill -f auto_pipeline.sh`
matched the ssh command that contained that string and killed the session before
the rest of the command ran.** It happened three times tonight - once killing a
zebrafish run mid-flight, twice silently dropping the `sed` that was meant to
follow.

Two rules:

- bracket a character in every `pgrep`/`pkill` pattern: `auto[_]pipeline.sh`
- never put `pkill` and the work it is clearing the way for in the same remote
  command; the kill can take the shell with it


---

## 2026-09-02 — peel on diagnostic columns, not identity

MANUAL §6.1.6 defines a subfamily as a **shared diagnostic pattern** — a
synapomorphy — and says the accumulation of private per-copy mutations is
*noise on top of the subfamily signal*. Average pairwise identity measures
exactly that noise, so the identity-based peel loop was answering the wrong
question. Rebuilt on features (`peel_features.py`).

**A group is a block of co-occurring (position, character) features**, and a
sequence joins if it carries ≥60 % of the block. That tolerance is the point:
one noisy column no longer expels a member, which is what makes the exact
co-occurrence sets in `SINEClusterer` fragment on 600 noisy chunks.

**Seeded growth, not single linkage, over the features.** Linking features
transitively chains them the same way single linkage chained t3/t4/t5 into one
143-chunk cluster — measured: union-find gave 5 blocks and placed 111 of 596
sequences. Taking the widest unused feature as a seed and attaching only what
agrees *with that seed* gave 26 blocks and placed 503.

### Result on Timema, against his t1–t8

| | identity-based | feature-based |
|---|---|---|
| placed | 111 / 596 | **503 / 596** |
| weighted purity | — | **0.958** |
| t3 | never recovered | 105 @ 96 %, 42 @ 98 % |
| groups | 5 of 8 | all of t1, t2, t3, t6, t7, t8 |
| re-align | never fired | fired twice |

Per group: t1 81 @100 %, t2 29 @100 %, t3 105 @96 %, t7 32 @100 %, t8 56 @100 %,
t6 7 @100 %. Residue 93, mostly t1.

**t1 comes out as four separate pure clusters** — 81, 72, 20, 16, all 100 % —
which is his own v4 curation splitting it into `t1_1 t1-2 t1-3 t1-4`. The
substructure is real and the method finds it without being told.

**His re-alignment rule is load-bearing, and this is the first run that showed
it.** Round 2, after 330 sequences left and the remainder was re-aligned, found
a t3 block of **169 co-occurring features** at 98 % purity. In round 1 that block
did not exist. Everything peeled in round 1 under the identity method, so the
effect could never be observed.

### Correction to the ground truth itself

The chunk labels were being built by matching member coordinates against
`assignment_full.tsv` with only `(+)` and `(-)` suffixes — silently dropping the
23.9 % of loci whose strand is `(+,-)`. Adding that suffix raised mean chunk
purity **0.953 → 0.975** and moved chunks between t3 and t4. Earlier t3/t4
conclusions were drawn on labels distorted by that gap.

---

## 2026-09-02 — re-running the page, and two silent-failure bugs

Re-scored the whole corpus with the current code to update the site. Result
first, because it is the boring one: **not a single set changed score.** 273
sets, published vs fresh run, zero differences. The lowest true positive is
81.5 of 132, exactly what the page already claims. The graph *data* was not
obsolete — the prose around it was.

### Bug 1 — the scorer silently degrades to "no evidence"

`verdict.py` loads four side inputs, each in a bare `try/except → {}`:

```python
try:    _DECAY   = json.load(open("flankdecay.json"))
except: _DECAY   = {}
```
…and the same for `flank_islands.json`, `microsat.json`, `region_part.json`.

A missing or *wrong* file is therefore indistinguishable from "this set has no
flank signal", and no-evidence scores **high**. Running the corpus without them
put every negative at 100 and gave a separation of −55 — with no error, no
warning, nothing in the log.

Worse, the files exist under two naming schemes. The corpus versions are
`flankdecay_corpus.json`, `islands_corpus.json`, `microsat_corpus.json`,
`corpus_regionpart.json` (673/673 coverage of `aln_c`). The generic names that
`verdict.py` actually opens hold the **new-species** data — 119 keys, **0 of
which match any corpus set**. Anyone who runs the scorer in that directory
scores the corpus with the flank evidence switched off and is told nothing.

This is the same failure shape as the zero-byte `genome.sizes` that made all 150
zebrafish alignments contain only the consensus row: an input that is absent or
wrong reads as an input that is legitimately empty. **Load these with an
explicit key-overlap check against the corpus and refuse to score at 0 %
overlap.**

### Bug 2 — step7's background is deflated by soft-masking

`step7_boundary_refine.sh` samples 300 random 70 bp genomic regions and compares
them **positionally and case-sensitively**:

```python
ident = (1 - sum(1 for i in range(min_len) if a[i] != b[i]) / min_len) * 100
```

The Timema genome is **48.4 % soft-masked**, so lowercase-vs-uppercase counts as
a mismatch. Measured on `genome.clean.fa`:

| | mean identity | elevated_frac(>45 %) |
|---|---|---|
| case-sensitive (as written) | **13.640 %** | 0.000 % |
| case-insensitive | **27.176 %** | 0.100 % |

13.640 reproduces the `background_identity_pct` in his Timema
`boundary_refinement.tsv` exactly, and scorpion's 14.586 is the same effect.
The true background is ~27 %, near the 0.25 the discriminator assumes.

**Scope of the damage: the reported number, not the calls.** The decision
variable is `elevated_frac`, which is 0.000 case-sensitive against 0.100
case-insensitive, both far under the 5 % threshold — so boundaries are
unaffected. But `background_identity_pct` is wrong in every
`boundary_refinement.tsv` written so far. One `.upper()` fixes it.

### What was actually obsolete on the page

- **"Subfamily structure — unreliable — confounded by age and copy number."**
  Superseded: 503 of 596 Timema chunks at 0.958 purity on diagnostic columns.
- **"Natural negatives — missing — no LINEs, satellites or processed
  pseudogenes yet."** Contradicted by the page's own data, which already
  contains **44** of them: 14 LINE, 16 chimera, 8 mosaic, 3 satellite,
  3 segmental duplication.
- **No subfamily section at all**, despite that being the strongest result.

Added a Subfamilies section, fixed both rows, added the nav entry. Verified by
rendering the page in headless Chrome: no console errors, table fits, nav has
the new entry.

---

## 2026-09-02 — the alignments were never built with `sear`

`grep -rn "sear " *.py *.sh` over the whole repo returns **nothing**. Every
alignment in this project was built by `candidate_to_aln.py`, which is my
reimplementation of his search:

```python
FLANK = 400
blastn -query q -subject genome -evalue 1e-10 -word_size 11 -dust no
       -outfmt '6 sseqid sstart send pident length bitscore sstrand'
```

then a hand-written single-pass overlap dedupe and `bedtools slop`.

**`sear -b100` already does this, better, and was sitting there the whole time:**
top N by bitscore → ±50 bp flanks → `getfasta -s` → query prepended as row 1 →
`mafft --localpair --maxiterate 1000 --ep 0.123`. Correctly oriented by
construction. What differs:

| | his `sear` | my `candidate_to_aln.py` |
|---|---|---|
| engine | `ssearch36`, Smith-Waterman | `blastn -word_size 11` |
| hit filter | span ≥80 % of query length, ≥65 % identity | `evalue 1e-10`, my own `hi-lo >= 40` |
| merge | `bedtools merge -s` | hand-written single-pass overlap loop |
| flank | 50 | **400** |
| top-100 alignment | built in, `-b100` | rebuilt by hand |
| orientation | stranded extraction, consensus row 1 | `--adjustdirection`, lost |

That single `FLANK = 400` is upstream of everything in this conversation:
`conse` on 400 bp-flanked copies calls **541 bases for a 253 bp element**, and
the over-extension then propagates through
`measure_c.py:91` / `justify_all.py:24`, which take the element window as the
consensus's first-to-last non-gap column with no validation.

**Sixth reimplementation.** SubFam, conse, sine_consensus.sh,
step7_boundary_refine.sh, step8a, and now sear.

### His calibration call on NEGTRUNC5__saq__s5_5seqs

He read the alignment by eye and said the consensus should start at **AGTTCGA**.
Measured, copy occupancy through the consensus:

| consensus base | occupancy |
|---|---|
| 1 | 0.34 |
| 25 | 0.51 |
| 43 | 0.56 |
| 61 | 0.70 |
| **64 = AGTTCGA** | **0.77** |
| 70–88 | 0.82 / 0.74 / 0.77 / 0.79 |

Occupancy rises monotonically and **plateaus at AGTTCGA**. My `OCC >= 0.50`
threshold fires at base 25 — 39 bases too early. **The boundary he sees is where
the occupancy curve flattens, not where it crosses a fixed fraction.**
`GTTCGA` is the tRNA B-box core, so these copies are 5′-truncated back to the
B box, which is exactly what the NEGTRUNC5 class is constructed to be.

His second observation on the same file — "flanks are incompletely sampled in
many copies and the left end of the flank is hardly identifiable" — is the same
fact from the other side: below 0.5 occupancy half the copies have no base
there, so there is no flank edge to see.

### Flank sizes

**50 L / 70 R**, his instruction. (`extract_alignments.sh` uses 30/70;
`step8a_extract_alignments.sh` uses `BASE_UP=50`, `BASE_DOWN=70`. The step8a
values are the ones to use.)

### What this means for the rebuild

Not a patch to the alignments. The product has to be rebuilt from the search up:

1. `sear -b100 <consensus> <genome>` for the top-100 view, `sear` + sampling for
   rand100 — his tool, his filters, his orientation handling.
2. Flanks 50 L / 70 R.
3. Keep a separate 400 bp extraction **only** for decay/island/microsatellite
   measurement, which is what `run_species.py` always said it was for — and
   never publish those files to the viewer.
4. Boundary from the occupancy plateau, calibrated against his AGTTCGA call.
5. Block on consensus length vs median copy length, so a 541-base consensus for
   a 253 bp element cannot ship.

---

## 2026-09-02 — peel tested against the v4 curated partition (13 subfamilies)

Harder test than the 8-group one, and closer to the saq case his own assist tool
fails on. Timema v4 `run_20260823_103133` has its own 596 chunk consensuses and
its own assignment, so the labels are self-consistent — no cross-run coordinate
mismatch (labelling the 2026-08-21 chunks with v4 copies matched **zero**, the
same failure that blocked saq).

Chunk labelling: 597/597 at mean purity **0.882** — that is the ceiling for any
chunk-level method here.

### Result: 493 of 595 placed, weighted purity 0.826 (94 % of the ceiling)

| v4 group | recovered as | purity |
|---|---|---|
| `t1_1` | 48 + 22 + 10 | 100 / 100 / 90 % |
| `t1-2` | 45 + 12 + 16 | 100 % each |
| `t1-3` | 25 | 96 % |
| `t1-4` | 41 + 11 | 100 / 91 % |
| `t2` | 22 | 100 % |
| `t8-1` | 16 | 100 % |
| `t8-2` | 32 | 100 % |
| `t6-1` | 10 | 100 % |

**All four t1 subfamilies separate.** The method was never told t1 has
substructure; on the 8-group labels it produced four pure t1 clusters, and
against v4 those correspond to `t1_1`, `t1-2`, `t1-3`, `t1-4`.

For scale, `cluster_assist.js` collapses **9** saq subfamilies into two blobs
plus a 6-sequence cluster.

### Two honest failures

- **`t3-1`/`t3-2` do not separate**: one 119-chunk cluster at 58 %. The t3 pair
  resists at every level — it is also the dominant multi-label CONFLICT set
  (`t3,t4,t5`, 69 585 loci) in the original run. Two independent methods and his
  own step1 all say these overlap.
- **`t1-2` fragments**: three pure pieces (45, 12, 16), one impure blob of 53 at
  49 %, and 89 of its 193 chunks left unplaced. It is the largest group and the
  peel shaves pieces off it instead of taking it whole — the 50 % size cap
  (`MAX_FRAC`) forbids a cluster larger than half the pool, and 193 of 595 is
  close enough to that limit to matter. Worth testing with the relaxed-upper-bound
  retry firing earlier.

Residue 102 chunks, 89 of them `t1-2`.

---

## 2026-09-03 — t3-1/t3-2: the peel was blind to indels

He looked at the subfam alignment and saw it immediately:

```
t3-1   CTCTGACAAGACAT-----------------
t3-2   CTCTGACAAGACGTGTATTAATCACAA
```

A ~15 bp insertion, present in t3-2, absent in t3-1. Measured on the v4
alignment: **14 columns sit at 90-97 % gap in t3-1 against 2 % in t3-2** - a
cleaner discriminator than any base-level column (21 columns differ in majority
base, and those alone were not enough).

`peel_features.py` skipped gaps entirely, copied verbatim from `SINEClusterer`'s
"gaps are not diagnostic". **That makes an indel-defined subfamily
unrepresentable**: t3-1's defining feature IS the gap. MANUAL 6.1.6 names
"small indels" *first* in its definition of a subfamily, so the rule I inherited
contradicts the definition it was meant to implement.

Fix: gap is a feature state alongside A/C/G/T. The global-dominance filter still
counts bases only - a column mostly gap across the alignment but carrying a
distinctive base subset must survive, or a minority insertion is filtered before
it can be used.

### Result on the v4 13-group partition

| | gaps skipped | gap as a state |
|---|---|---|
| chunks placed | 493 / 595 | **584 / 595** |
| weighted purity | 0.826 | **0.926** |
| t3-1 | one 119-chunk cluster at 58 % | **63 @ 100 %** |
| t3-2 | (same cluster) | **53 @ 94 %** |
| t1-2 | 45+12+16 pure, 53 @ 49 %, 89 residue | **88 @ 100 %** |
| residue | 102 | **11** |

Also 100 % pure: t1_1 (46), t1-3 (25), t2 (22), t8-1 (16), t8-2 (31).

**Still failing:** t1-4 fragments into five pieces (24 @ 67 %, 14 @ 79 %,
14 @ 86 %, 11 @ 100 %, 10 @ 90 %) and accounts for the entire 11-chunk residue;
t6-1 sits at 69 %.

**What this says about the method, not just the bug:** every earlier conclusion
about t3 was drawn with the discriminating signal switched off. The identity
peel called t3 "internally incoherent at 0.685"; the feature peel called
t3-1/t3-2 "the one group that fails". Both were measuring a set whose actual
defining character had been discarded before measurement. His `t345.al` grouping
is still evidence that t3/t4/t5 are close - but "cannot be separated" was my
artefact, not their property.

---

## 2026-09-03 — three of the viewer's ideas, applied to the subfam peel

Taken from MSA-viewer, not invented:

1. **Trim before clustering.** `getSeqsForClustering()` applies the soft trim
   boundaries first, so ragged ends cannot manufacture features. Ported
   `getTrimBoundaries` (sliding window, gap fraction) and run it per round on the
   feature-finding view only — members are still peeled from untrimmed rows. His
   saved preset's 0.6/0.6 with a 20-column window, not the library's 0.50/0.80/15.
2. **Sweep instead of a single run** (`analyzeClusterability`), with an
   "everything loosest" run, because "a single run cannot distinguish settings
   that are too strict from data with no structure to find".
3. Parameters made overridable from the environment so the sweep can drive them.

### What trimming alone bought

| | purity | t1-4 |
|---|---|---|
| no trim | 0.926 | 24 @ 67 % |
| trimmed | **0.937** | **34 @ 91 %** |

### What the sweep says about the three open groups

26 runs over MIN_SET, MIN_BLOCK, FEAT_JACCARD, CARRY, EXCL_MIN, MAX_FRAC, trim
on/off, and everything-loosest:

- **t3-1 / t3-2 — solved and robust.** 63 @ 100 % and 53 @ 94 % in 22 of 26
  settings. They only re-merge at `MIN_BLOCK 8` (120 @ 57 %) and under
  everything-loosest. This is no longer a borderline result.
- **t6-1 — not unseparable, under-resolved.** 16 @ 69 % at default, but a pure
  **8–9 @ 100 %** core appears whenever stringency rises (`MIN_BLOCK 2`,
  `FEAT_JACCARD 0.60` or `0.75`, `CARRY 0.75`). The default cluster of 16 is
  contaminated; the real group has a clean core of about 8 with a soft edge.
- **t1-4 — the one genuine failure.** 88 chunks in truth, and it fragments at
  *every* setting: 34 @ 91 %, 41 @ 88 %, 29 @ 93 %, 43 @ 81 %, 98 @ 78 %,
  19 @ 53 %, 16 @ 56 %, absent entirely at `FEAT_JACCARD 0.75`. No setting
  recovers it whole. This one is worth his eye.

### Two findings about my own parameters

- **`EXCL_MIN` is inert.** 0.15 / 0.25 / 0.35 / 0.50 give byte-identical output
  on every column. A knob I invented that does nothing — the other filters
  dominate it. Remove or re-think it, do not tune it.
- **Best settings are not the defaults I chose.** `MIN_SET 8` gives 591 placed at
  **0.941** with a residue of 4, against 583 / 0.937 / 12 at `MIN_SET 5`.
  `FEAT_JACCARD 0.75` reaches 0.977 purity but places only 306 of 597 — purity
  bought by placing less, which the sweep makes visible and a single run would
  have hidden.

Full table: `sweep_table.txt` in this repo.

---

## Calibration table — his eye against measured column separation

The core idea of the project: he shows an alignment, says what he sees, and I
record the call and compute variables against it. This is the table. Two entries
so far, both on Timema v4, measured at **copy level** (one representative copy
per chunk, taken from the 50-member `.bnk` banks) with an **occupancy guard** —
a column only counts as a base difference when both groups are at least 30 %
present, otherwise "all of A has T, none of B does" is just where B's copies end.

| pair | his call | cols \|diff\| ≥0.8 | best gap col | best base col | cols ≥0.5 |
|---|---|---|---|---|---|
| **t3-1 / t3-2** | saw it immediately, named the insertion | **14** | **0.95** (0.97 vs 0.02) | — | 14 |
| **t1-4 / t1-2** | "faint difference in 3' part by small deletion and few other places. still very borderline" | **0** | 0.66 (5' end, cols 4-16) | 0.55 (col 345) | 18 gap + 4 base |

**Working threshold from these two points:** a separation he names instantly
carries features at **\|diff\| ≥ 0.9**; one he calls borderline tops out around
**0.5-0.66 with nothing above 0.8**. Needs more of his calls to firm up, but it
is the first quantitative handle on "is this a subfamily".

His "small deletion in the 3' part" is real and measurable: in the 3' half the
best base column is 0.55 (col 345) and the best gap column 0.51 (col 441) - the
deletion is present in about half of t1-4 and absent from about 80 % of t1-2.

### What this says about the peel, and it is not a bug

`CARRY = 0.60` asks a sequence to carry 60 % of a block's features. When no
feature exceeds 0.66, block membership is inherently unstable, so t1-4 shattering
into eight pieces is **the correct response to a signal that weak**, not
something to tune away. The sweep already showed no setting recovers it whole;
this says why.

**Consequence for the method:** the peel needs to report *confidence*, not just a
partition. A group whose best feature is 0.5 should come out labelled grey-zone,
the way he labelled it by eye - not silently fragmented into eight clusters that
look like eight findings.

### Two method lessons from the same exchange

- **Use the general subfam alignment as the source.** He read t1-4's weak
  separation off the full SINEderella-hit subfam alignment, where the competing
  groups are present. Extracting one group and re-aligning it alone destroys the
  comparison that makes a separation visible - which is what I did first, and he
  could not read it.
- **Drop to the banks for a grey-zone call.** Chunk consensuses are
  `cons -plurality 18` over 50 copies: a majority rule, so a diagnostic carried
  by a minority of copies within each chunk cannot survive it. Any weak or
  partially penetrant difference is invisible at chunk level regardless of
  parameters. The `.bnk` member banks are one level down and keep it.
- **The nearest competitor falls out of the peel for free.** I guessed t1-3 from
  the screenshot layout; he corrected it to t1-2. My own output already said so -
  the only t1-4 chunks that leaked into another group's clusters went to t1-2
  clusters (2+2+1), none to t1-3. Leak destination identifies the nearest
  competitor at no extra cost.

---

## MAFFT's ordering — read from the C source

Why: he pointed at MAFFT's sequence ordering and said "separation is faint and
real at the same time — mafft has its good intuition. Can you learn how mafft
uses its magic?"

Source: `GSLBiotech/mafft`, `core/disttbfast.c` and `core/mltaln9.c`. The
`disttbfast.js` in MSA-viewer is minified Emscripten output with no readable
algorithm — go to the C.

### The distance, exactly

```
tuplesize = 6                                        disttbfast.c:214
shared    = commonsextet_p(A, B)                     mltaln9.c:14871
bunbo     = MIN(selfscore_i, selfscore_j)            disttbfast.c:4053
lenfac    = 1 / ( shorter/longer * 0.1
                  + 2500/(longer + 2500)
                  + 0.01 )                           disttbfast.c:4049
dist      = (1 - shared/bunbo) * lenfac * 2.0        disttbfast.c:4057
```

`D6LENFACA 0.01, D6LENFACB 2500, D6LENFACC 2500, D6LENFACD 0.1`
(`disttbfast.c:49`). Amino acid and 10-mer modes have their own constants.

**`commonsextet_p` is a multiset intersection.** Walking B's k-mer occurrences,
an occurrence counts only while B has not yet used up A's copies of that k-mer:

```c
tmp = memo[point]++;
if( tmp < table[point] ) value++;
```

`seq_grp_nuc` drops any non-ACGT before k-merising, and marks a sequence
unusable if it is shorter than `tuplesize`.

Ported and verified in `mafft_dist.py`: self-distance 0.0; on a 24 bp toy, one
substitution costs 0.574 and a 5 bp deletion 0.662.

### Three properties worth taking

1. **Normalised by MIN(self), not by a union.** A short sequence wholly contained
   in a longer one scores similarity 1.0. This is containment, not Jaccard, and
   it is tolerant of truncation — which is what SINE copies of varying
   completeness need. Nothing in my peel has that tolerance.
2. **It aggregates weak evidence.** No single column has to be diagnostic; a
   difference that is faint everywhere but consistent still accumulates. That is
   the mechanism behind "faint and real at the same time".
3. **An indel is weighted by its k-mer footprint, not its column count.** A
   deletion of length L destroys L+k−1 overlapping k-mers; a substitution
   destroys at most k. A 5 bp indel therefore counts ~10 against a SNP's ~6 —
   indels outweigh substitutions automatically, which is the weighting MANUAL
   §6.1.6 asks for ("small indels/SNPs", indels named first). **My peel counts a
   5 bp indel as 5 features and a SNP as 1** — under-weighting indels by roughly
   a factor of two relative to MAFFT, in the wrong direction.

**His edge rule falls out of the same arithmetic.** K-mers spanning a ragged edge
are mostly unique, so they enter `bunbo` (the denominator) but never `shared` —
ragged edges inflate distance systematically. That is a mechanical reason to trim
or downweight edges *before* computing a distance, not after, and it matches his
instruction: edges are used in comparison only when they carry firm evidence well
above the noise level characteristic of simple repeats.

### Tested, and it does NOT work as a group statistic

Mean pairwise MAFFT distance, on representative copies, against his two calls:

| pair | within A | within B | between | separation |
|---|---|---|---|---|
| t3-1 / t3-2 — he named it instantly | 0.562 | **0.994** | 0.888 | **−0.105** |
| t1-4 / t1-2 — he called it faint but real | 1.697 | 1.775 | 1.781 | +0.006 |

**It fails on the easy pair.** t3-2's internal spread (0.994) exceeds its distance
to t3-1 (0.888), so a group-mean test says t3 is not separable — the opposite of
what he saw in one glance, and of what the 14 gap columns at |diff| 0.95 say.

**Conclusion: the aggregate distance is not the magic, and I measured the wrong
thing.** He said *sorting*. MAFFT's ordering comes from the guide tree, which
joins nearest neighbours — local structure. A group mean averages exactly that
away. A heterogeneous group with a tight sub-lineage looks bad by mean and fine
by tree.

### Next step, not yet done

Test the **tree**, not the mean: build the UPGMA/nearest-neighbour tree on this
distance and ask whether t3-1/t3-2 fall into separate clades and t1-4/t1-2 into
partially separate ones. The viewer already exposes this as "guide-tree groups",
cutting the 6-mer tree into N groups and **reporting the height of the next
merge** — a ready-made isolation measure, better than the one I invented.

The likely synthesis, to be tested rather than assumed: use diagnostic columns to
*propose* a split, and the guide-tree cut height to say how isolated it is —
identity for isolation, synapomorphy for identity-of-group.

---

## The measure that reproduces his eye: adjacency in MAFFT's leaf order

Three attempts on the same two calibration pairs. Only the third works.

| test | t3-1/t3-2 (he named it instantly) | t1-4/t1-2 ("faint and real") |
|---|---|---|
| **mean 6-mer distance**, between minus within | **−0.105** — says not separable | +0.006 |
| **guide-tree cut** at k=2 | splits 119 vs 2 — peels an outlier | splits 279 vs 1 |
| **adjacency in leaf order** | **0.858** vs chance 0.512 | **0.785** vs chance 0.569 |

Normalised to a 0-1 scale, where 0 is a random ordering and 1 is perfectly
contiguous — `(observed − chance) / (max − chance)`, with `max = (n−2)/(n−1)`
for two groups:

| pair | his call | **score** |
|---|---|---|
| t3-1 / t3-2 | named it at a glance | **0.72** |
| t1-4 / t1-2 | "faint and real at the same time" | **0.51** |

**Both well above chance, t3 clearly stronger.** That is his reading, as a number.

### Why the first two failed and this one does not

- A **group mean** asks "are these two clouds apart". It fails when a group is
  internally heterogeneous, which t3-2 is: its within-group distance (0.994) is
  larger than its distance to t3-1 (0.888), so the mean says unseparable while
  the eye says obvious.
- A **top-level tree cut** asks "does the highest split fall here". On diverged
  copies the highest splits are driven by outliers — at k=2 it peels off two
  sequences and leaves 119 together.
- **Adjacency** asks "do members of a group sit next to each other in the
  ordering". It survives internal heterogeneity, because it only ever looks at
  neighbours. That is exactly what he is reading off the screen when he says the
  sorting shows it.

This is the same statistic that measured SubFam's chunk ordering — adjacent
chunks share a subfamily 0.916 of the time against 0.313 for a random pair. The
same tool answers both questions and I did not connect them until now.

### Standing use

Report an **ordering-adjacency score** alongside any proposed split, as the
confidence figure the peel currently lacks. Working bands from these two points,
to be firmed up with more of his calls:

- **≥ 0.70** — a separation he names without hesitating
- **~0.50** — real but faint; the grey zone, report as such rather than splitting
- **≈ 0.0** — no structure in the ordering

Caveats worth keeping: two calibration points is not a calibration; the score
depends on the k-mer size and on edge trimming (ragged edges inflate distance
systematically, since edge k-mers enter the denominator but never the shared
count); and it measures *ordering*, so it says a difference is consistent, not
what the difference is. Pair it with the diagnostic columns, which say what.

---

## Variable-sites-only, borrowed from the viewer — it sharpens the score

His observation: the viewer has a variable-sites-only mode, and "at even 15 %
difference it shows a much more readable picture of difference between upper and
lower sequences".

### The exact rule (script.js:5555-5570)

```
covered   = sequences whose OWN first..last non-gap span covers this column
counts    = character tally over those sequences, gaps counted as '-'
diffCount = covered - max(counts)          # size of the minority
keep column if diffCount >= ceil(pct/100 * nSeq)
```

Two details already in it that I had to rediscover the hard way:

- **The occupancy guard.** A sequence that does not reach a column is not
  counted, so ragged ends cannot manufacture a difference. This is exactly the
  guard I had to add by hand for the TACAT test after my first pass reported 31
  "perfect" discriminators that were only where t1-2 copies end.
- **Internal gaps are a character.** An indel inside a sequence's own span makes
  the column variable — the same fix that made t3-1/t3-2 separable in the peel.

Note the threshold is a fraction of **all** sequences, not of the covered ones.

### Measured: it improves both calibration pairs and widens the gap

Ordering-adjacency score, restricted to variable sites:

| threshold | columns (t3) | **t3** (he named it instantly) | columns (t1-4) | **t1-4/t1-2** (faint but real) | gap |
|---|---|---|---|---|---|
| all sites | 591 | 0.652 | 557 | 0.497 | 0.16 |
| 5 % | 94 | 0.826 | 265 | 0.514 | 0.31 |
| 10 % | 34 | **0.896** | 249 | 0.514 | 0.38 |
| **15 %** | 25 | **0.896** | 214 | 0.547 | **0.35** |
| 20 % | 22 | 0.878 | 173 | 0.581 | 0.30 |
| 30 % | 19 | too few | 106 | 0.564 | — |

**t3 rises from 0.65 to 0.90 on 25 of 591 columns** — a 24-fold reduction that
makes the signal sharper rather than noisier. t1-4 rises more modestly, 0.50 to
0.55 at 15 %, peaking at 0.58 at 20 %.

His 15 % is a good default: near the top for both, and the gap between his
"obvious" and his "borderline" call is widest in the 10-15 % band.

### Revised calibration bands (variable sites at 15 %)

| his call | score |
|---|---|
| named at a glance (t3-1/t3-2) | **0.90** |
| "faint and real at the same time" (t1-4/t1-2) | **0.55** |

Better separated than the 0.72 / 0.51 obtained on all sites.

### The part worth keeping beyond this measurement

At 15 % the filter leaves ~25 columns for t3 — the same scale as the 14 gap
columns found by comparing the labelled groups directly. **So the variable-sites
filter is an unsupervised diagnostic-position finder.** `step4_diagnostic.py`'s
MI / RandomForest / KL weighting needs a partition to compute against; this needs
nothing. That makes it usable at the point in the peel where no partition exists
yet — propose candidate positions with the variable-site filter, cluster on them,
then hand the resulting partition to step4_diagnostic for proper weighting.

---

## The separator between variable sites is evidence, not bookkeeping

His point: "one difficult to account but need to be considered is the separator
inserted in place of ignored similarity columns". The viewer records these as
`_varSiteHiddenRanges` and draws them as breakpoints with a count.

He is right, and it changes two things.

### 1. It corrects my own arithmetic

I reported t3-1/t3-2 as having "14 gap columns at |diff| 0.95". Those columns
were **207-218, contiguous** — that is **one indel**, not fourteen independent
diagnostics. Collapsing runs into events is the honest count.

It also means `MIN_BLOCK = 3` in the peel is trivially satisfied by any indel of
3 bp or more while believing it has three separate pieces of evidence. A block
built from one indel and a block built from three dispersed substitutions are
scored identically and should not be.

### 2. The separator distribution is itself a discriminator

Variable sites at the 15 % threshold, collapsed into events (consecutive sites
merged):

| | sites | events | median separator | single-column events | runs >=3 |
|---|---|---|---|---|---|
| **t3-1/t3-2** (he named it instantly) | 25 | **14** | **14 cols** | 12 of 14 | 2 |
| **t1-4/t1-2** (faint but real) | 214 | **78** | **3 cols** | 24 of 78 | 41 |

**t3's evidence is dispersed and independent**: twelve isolated single columns a
median of 14 columns apart, plus one 8-column indel. That is what a set of
distributed synapomorphies looks like.

**t1-4's is packed**: a median of 3 columns between events and 41 runs of three
or more. Its 214 variable sites are really about 78 events, and those events sit
in local clusters — which is the signature of alignment-uncertainty patches and
hypervariable regions rather than evidence spread across the element.

The sites-to-events ratio says the same thing more compactly: **1.8 for t3, 2.7
for t1-4**.

### What to do with it

- **Count events, not columns**, everywhere a block or a feature set is scored.
- **Report the separator distribution** alongside a proposed split. Dispersed
  evidence (large median separator, mostly single-column events) is stronger than
  the same number of variable sites packed into a few patches.
- This is a second axis independent of the adjacency score: adjacency says *how
  consistently* the groups sort apart, the separator distribution says *whether
  the evidence is independent*. t1-4 scores moderately on the first (0.55) and
  badly on the second (median separator 3) — which is a fuller description of
  "faint" than either number alone.

Still to check: whether the packed t1-4 events coincide with regions where MAFFT
itself is uncertain, which would confirm the alignment-artefact reading. The
viewer's realign-and-compare would answer it.

---

## Event stability tested — my artefact hypothesis was wrong, and the answer is better

I predicted t1-4's packed events were alignment artefacts and would move under
re-alignment. **Tested and rejected.** Both pairs realigned five ways — localpair
with and without `--ep 0.123`, `--auto`, `--retree 2`, `--globalpair` — with
events mapped into a fixed reference sequence's coordinates so alignments of
different width can be compared:

| | events | present in all 5 alignments |
|---|---|---|
| t3-1/t3-2 | 13 | **13 (100 %)** |
| t1-4/t1-2 | 31 | **31 (100 %)** |

Neither is alignment noise. Column counts drift a little (t1-4: 212-228) and
event boundaries shift by a base or two, but with 2 bp slop every event survives
every aligner.

### A measurement error this exposed

I earlier reported t1-4 as having **78 events**. That was counted in
**alignment-column** space, where gap columns split one event into several. In
**reference coordinates** it is **31**. Count events in a reference frame, not in
alignment columns.

### The real contrast, in reference coordinates

| | reference | variable cols | events | sites/event | median separator | largest event |
|---|---|---|---|---|---|---|
| **t3-1/t3-2** | 282 bp | 25 (**9 %**) | 13 | 1.9 | **19 bp** | 1 bp |
| **t1-4/t1-2** | 263 bp | 213 (**81 %**) | 31 | 6.9 | **1 bp** | **51 bp** |

**t1-4 and t1-2 differ over 81 % of the element, and still nothing reaches
|diff| 0.66.** That is not a faint difference between two tight groups — it is
two pools that are each internally diverse and mutually overlapping, with no
consistent signal to separate them. Every t1-4 event is stable and real; none is
diagnostic.

t3 is the mirror image: difference at only 9 % of positions, but consistent (up
to 0.95) and dispersed a median of 19 bp apart.

### The rule this gives

- **few, strong, dispersed, consistent** → a real subfamily split
- **many, weak, packed, inconsistent** → one heterogeneous group, do not split

This is stronger support for his "merge t1-2 and t1-4" than anything offered
before, and it arrives from the opposite direction to my guess. The packedness is
not artefact; it is what pervasive low-level divergence looks like when two
groups have no synapomorphy between them.

**Add to the score:** fraction of the reference that is variable. t3 9 %,
t1-4 81 %. A genuine split should be *sparse* and consistent; pervasive variation
with no strong column is the signature of one pool.

---

## Fraction-variable tested across pairs — it does NOT work. `best` does.

He asked whether the 9 %-vs-81 % figure holds up as a discriminator. Tested on
every pair, and on three pairs at copy level. **It does not.**

### At copy level, the three pairs whose status he has judged

| pair | his call | `frac_var` | `best` |
|---|---|---|---|
| t3-1 / t3-2 | named the insertion at a glance — **split** | 0.08 | **0.96** |
| t1-4 / t1-2 | "more plausible to merge" — **reject** | 0.39 | **0.66** |
| t2 / t8-2 | distant, obviously distinct | **0.70** | **1.00** |

`frac_var` runs 0.08 → 0.39 → 0.70 with no relation to the calls, and the most
clearly distinct pair has the **highest** value. It measures evolutionary
**distance**, not separability: two distant subfamilies differ nearly everywhere
*and* carry consistent diagnostics. Pervasive variation is not evidence against a
split.

`best` — the strongest single consistent column, max over gap-fraction difference
and occupancy-guarded majority-base difference — tracks his calls exactly.

**Threshold, from his instruction to reject anything below the t1-4 case:**

| `best` | reading |
|---|---|
| **>= 0.85** | split — a consistent diagnostic exists |
| 0.70 – 0.85 | untested band |
| **<= 0.70** | merge — nothing consistent, whatever the distance |

One pass over columns. Cheaper than the adjacency score, the tree, and the
k-mer distance combined.

### Where my earlier framing went wrong

"9 % vs 81 %" was a real measurement attributed to the wrong cause. Two further
errors in it, both stated for the record: the 81 % used one reference sequence's
full length as denominator, where the shared element span gives 0.39; and the
contrast I read as sparse-vs-pervasive was really consistent-vs-inconsistent.

### The chunk consensuses lie — measured

Surveying all 45 pairs of the 10 well-populated v4 groups **in the general subfam
alignment of chunk consensuses**:

| pair | frac_var | best | reading at chunk level | truth |
|---|---|---|---|---|
| t3-1 / t3-2 | 0.03 | 0.97 | split | correct |
| **t1-2 / t1-4** | **0.13** | **0.86** | **split** | **wrong — he merges it** |

At chunk level t1-2/t1-4 looks *sparse and consistent*. At copy level the same
pair is `best 0.66`. `cons -plurality 18` over 50 copies smooths each chunk, so
an overlapping pair acquires a clean apparent difference that its copies do not
support.

**This is not a refinement of his instruction to use the banks — it is the reason
for it.** A grey-zone call made on chunk consensuses gets the wrong answer.
Chunk level is for proposing candidates; copy level is for judging them.

All 45 chunk-level pairs also show `best` at 0.97-1.00 with one exception, which
is another way of saying the same thing: at chunk level almost everything looks
separable.

---

## saq recreated from scratch — 9 groups, three consensuses recovered exactly

He was asked for the s1-s9 chunk labels and answered "i dont know, recreate from
scratch?". So: the peel run on saq's 600 chunk consensuses
(`run_20260425_110449`), unsupervised, no labels.

### The size cap was the binding constraint

Default settings give 4 groups and 157 unplaced. The sweep found why:

| setting | groups | placed | residue |
|---|---|---|---|
| default (`MAX_FRAC 0.50`) | 4 | 443 | 157 |
| `MIN_BLOCK 2` | 8 | 595 | 5 |
| **`MAX_FRAC 0.90`** | **9** | **598** | **2** |

`MAX_FRAC` caps a group at a fraction of the **remaining** pool, and after each
peel the pool shrinks — so a legitimately large group is refused once the pool is
small enough. saq's s8 is 225 of 600 (37 %) and s7g 172 (29 %), both under 50 %
of the *original* pool but not of what remains. This is inherited from
`SINEClusterer`'s 50 % cap, which has a relaxed-upper-bound retry for exactly this
case; mine gated it too tightly.

Also: saq yields only **79 features** in round 1 against Timema's **616** from a
similar number of chunks. Its subfamilies differ far more subtly, which is why it
is the hard case.

### Result: 9 groups, matched to his by consensus identity

| mine | chunks | best match | identity | runner-up |
|---|---|---|---|---|
| grp05 | **43** | **s3_43seqs** | **1.000** | s7g 0.969 |
| grp06 | 40 | **s2_38seqs** | **1.000** | s9 0.947 |
| grp02 | 109 | **s8_225seqs** | **1.000** | s4 0.980 |
| grp04 | 70 | s1_30seqs | 0.996 | s9 0.988 |
| grp03 | 99 | s7g_172seqs | 0.996 | s9 0.981 |
| grp09 | 8 | s7g_172seqs | 0.985 | s9 0.984 |
| grp07 | 39 | s7g_172seqs | 0.981 | s3 0.977 |
| grp01 | 155 | s4_60seqs | 0.968 | **s5 0.968 (tied)** |
| grp08 | 35 | s8_225seqs | 0.964 | s3 0.945 |

Sizes: mine `155 109 99 70 43 40 39 35 8` (598 placed, 2 residue) against his
`225 172 60 43 38 30 20 7 5`.

**Three consensuses recovered exactly** (s2, s3, s8), two more at 0.996 (s1,
s7g). `grp05` matches `s3_43seqs` at 43 chunks *and* an identical consensus.

For scale: `cluster_assist.js` on this same data gives 2 blobs of 293 and 307
plus one 6-sequence cluster.

### Two failure modes, both structural

- **Over-splits the two largest.** s7g goes to three pieces (99, 39, 8) and s8 to
  two (109, 35). The peel removes a group as soon as a block is exclusive enough,
  so a large group with internal structure is taken in pieces across rounds and
  never reassembled. There is no merge step.
- **Loses the three smallest** — s5 (5 chunks), s6 (7), s9 (20) never appear.
  `MIN_GROUP = 5` and `MIN_SET = 5` put a floor under what can be found, and s5
  at 5 chunks sits exactly on it.

Both are fixable in principle: a merge pass that reconsiders whether two peeled
groups belong together (using the `best` test at copy level), and a lower floor
for late rounds when little is left.

**Not yet judged by him.** The matching above is against his *consensuses*, which
is indirect — his actual chunk-to-group assignment is still not on disk.

---

## Correction: the `subfam/*.al` files are NOT all the same 200 chunks

I reported that saq's nine `.al` files "all contain the same 200 chunks", and
concluded they were per-group views of one shared sample rather than a partition.
**Wrong, and the error was mine.**

`run_subfam_per_sf.sh` runs SubFam **separately for each subfamily**, on that
subfamily's own sampled copies. Each run numbers its chunks `input_001…` from
scratch. So nine files carrying `input_001.bnk` through `input_200.bnk` hold nine
**different** sets of chunk consensuses that happen to share names.

The set-overlap test I ran compared names, so it reported 100 % overlap and I
believed it. The test was measuring a naming collision.

**Rule: never compare SubFam chunk identities by name across runs.** `input_042`
means "the 42nd chunk of whichever run produced this file". Compare by sequence,
or by the genomic coordinates of the members in the `.bnk`.

What this means for saq: each `.al` *is* per-group data — the chunk consensuses of
that group's own copies, plus that group's consensus. It is not the assignment of
the original 600 chunks to nine groups, because those chunk numbers belong to a
different SubFam run. That assignment is still being looked for.

---

## Where the saq chunk assignment is — searched, and what the search established

He said "i never delete important data, try harder to find it or tell me how to
do from start", and later "dragen?".

**I could not find it.** What the search established positively:

- **KIT `.bash_eternal_history` is the decisive evidence.** Between
  `SINEderella GCA_004024925.1… tal.bnk` at **11:04** on 25 Apr 2026 (which
  created `run_20260425_110449`) and `tal.bnk` appearing at **18:21** holding
  exactly the nine `s*` consensuses, there are no relevant commands on KIT at
  all. The curation was not done on the server.
- Every file under the saq tree in that window (12:15-12:30) is pipeline output
  from run_110449.
- **DRAGEN**: the only saq files are my own corpus derivatives from 2026-08-31.

So the grouping was made visually in MSA-viewer and what came back to the server
was the consensus bank, not the assignment. That points at the browser's
localStorage or a local download — neither reachable from here. **If it is on
another machine, or exported under a name I did not try, that changes the
answer.**

### `asSINEment` cannot reconstruct it — measured

Running his own assignment engine, the nine consensuses as references against the
600 chunk consensuses:

| subfamily | chunks assigned |
|---|---|
| s2_38seqs | 32 |
| s3_43seqs | 21 |
| s6_7seqs | 7 |
| s5_5seqs | 4 |
| s7g_172seqs | 1 |
| **total** | **65 of 600** |

**The 10/10 unanimity rule is why.** His nine consensuses are 96-99 % identical to
each other, so across ten stochastic `ssearch36` cycles the winner flips and
unanimity almost never holds. The rule is designed for genomic copies competing
between *distinct families*; it does not transfer to chunk consensuses competing
between near-identical *subfamilies*.

This is a real limit worth remembering: **asSINEment is not a tool for assigning
chunks to subfamilies.** The peel is.

### How to redo it from start, in order of cost

1. **Re-export from MSA-viewer** if the session survives in that browser profile —
   it keeps cluster state in localStorage and can download it.
2. **Re-curate visually** — load `run_20260425_110449`'s
   `genome.clean_step1/subfam_input/input.clw.al` (600 chunks plus the
   g-consensuses), group by eye, export each group. This is the original route.
3. **Correct my reconstruction** — 9 groups, 598 of 600 placed, three consensuses
   recovered at identity 1.000, over-splitting s7g and s8 and missing the three
   smallest. A starting partition to fix, not an answer.

---

## His manual peel, walked step by step — round 1 on saq (2026-09-03)

He walked me through what he actually does, instead of me inferring it. Recorded
verbatim as a procedure, because this is the loop the whole discriminator is meant
to replace.

**Step 1 — throw away the garbage rows.** "last3 sequences are garbage and need to
be thrown away". In the 609-row curation alignment those were `input_587.bnk`
(368 bp), `input_599.bnk` (155 bp), `input_600.bnk` (92 bp) — against a 254 bp
median and a 262 bp 95th percentile. MAFFT's `--reorder` puts them last.
`input_587` alone was generating **119 columns of pure gap in 608 of 609 rows**,
29 % of the alignment width.

*The operational rule: length outliers at the bottom of the reorder are removed
before anything else is judged.* They are `split -l 50` remainders and fused
chunks, not subfamily material.

**Step 2 — peel the clearly-defined block.** "a group of sequences below s8 are
clearly defined, so i copy them into new alignment, and remove from original".
His boundary: from `CONS__s8_225seqs` (row 462) down to `input_592.bnk` (row 605)
= **143 chunks**.

**Step 3 — consensus of the peeled block.** "i copy sequences from s8 to 592 to new
viewalign window and use its internal consensus creation tool - copy consensus
button - with default threshold 50% or lower it if some positions need it".

So the consensus tool is **MSA-viewer's Copy Consensus**, not `conse`.

**Step 4 — continue with the rest.** Remainder is realigned and the loop repeats.

### The falsification this produced

Before he gave the boundary I computed adjacent-row MAFFT 6-mer distance jumps
below s8 and proposed the two largest as candidate boundaries: row 529 (a block of
67) and row 566 (a block of 104), jumps 0.418 and 0.446. **He drew the boundary at
605 — all 143.** So those jumps are *internal structure inside a single subfamily*,
and their magnitude is indistinguishable from a real boundary's.

**Adjacent-row distance jumps do not locate subfamily boundaries.** Rejected, on
his call, the same way mean pairwise distance was rejected on t3.

### What did validate: three consensus routes agree exactly

On the 143 peeled chunks:

| route | ungapped len | identity to his `s8_225seqs` |
|---|---|---|
| `conse` (EMBOSS `cons`) at plurality 25 / 35 / 50 % | 254 | **1.0000** |
| MSA-viewer Copy Consensus at threshold 50 / 40 / 35 / 30 % | 254 | **1.0000** |

254 of 254 columns, zero differences, zero gap differences, at every threshold.

Two things follow. First, **the peel found exactly his s8** — this is the first
step of the reconstruction that reproduces his curation at the sequence level, not
just approximately. Second, **threshold-insensitivity is itself a signal**: a
well-defined group gives the same consensus from 30 % to 50 %, so needing to lower
the threshold (which he named as a possibility) marks a group that is *not* clean.
That is a free, cheap statistic to record per block.

Note the count: 143 chunks here against the 225 encoded in his bank name, and the
consensus is still byte-identical. Consensus identity therefore does **not**
constrain group membership tightly — worth remembering before using consensus
agreement as evidence that a grouping is right.

### MSA-viewer's Copy Consensus, ported exactly

`script.js:2469`, `_computeConsensusCharForColumn`. Defaults from `DEFAULTS`
(`script.js:339`): threshold 50, minCoverage 30, consType `normal`, fallbackMode
`gap`.

```
col      = every row at pos, uppercased, missing -> '-'
nonGap   = col without '-' and '.'
if nonGap empty                     -> '-'
if len(nonGap)/len(col) < 0.30      -> '-'          # coverage floor
counts over nonGap only
freq = maxCount / len(col)                          # DENOMINATOR INCLUDES GAPS
if freq >= threshold -> top base, ACGT preferred, ties alphabetical
else                 -> '-'                         # fallbackMode 'gap'
```

**The denominator including gaps is the part that differs from EMBOSS `cons`**, and
the source comments it as deliberate: "for consistency with display". A column that
is mostly gap can never clear the threshold. Ported in `viewer_cons.py`; it agrees
with `conse` on this block, but the two are not the same rule and will diverge on
gappy columns.

### Round 1 output

| file | rows | cols |
|---|---|---|
| `CURATE__saq__s8_peeled_143.aln.fa` | 143 chunks + his consensus | 271 |
| `CURATE__saq__round2_after_s8.aln.fa` | 454 chunks + 8 landmark consensuses | 289 |

The round-2 alignment is clean: core columns 1..285 of 289, left overhang 1, right
overhang 3. Dropping `input_587` removed the entire 119-column artefact.

## Round 2 on saq: s6 peeled, and the exclusion rule made concrete

His call, verbatim: "Next comes s6, few sequence with 2 point differences - T at
about 200-201 positon, and T deleted after position 225 (approximate) - find it
yourself? BUT 3 lower ones are not in this group, i just discard them (possible
very minor subfamily) - 437,438 (similar) and 342 (probably some transitional one
- subfam has this problem in 50s of loci where it is composed of sequences from
adjacent subfamilies)".

**His coordinates are ungapped element positions, not alignment columns.** Both
landed exactly:

| his call | ungapped pos | aln col | the 7 | all 444 others |
|---|---|---|---|---|
| "T at about 200-201" | 200 | 223 | `T` | `C` 444/444 |
| "T deleted after position 225" | 225 | 251 | `-` | `T` 444/444 |

Two further columns are equally unanimous and he did not name them: ungapped 188
(`G`, 0/444 outside) and ungapped 55 (`A`, 1.6 % outside). So the group carries
**four** diagnostics, not two. He called it on the two he happened to see; the
others are there.

### The exclusion is the valuable part

`input_437`, `input_438`, `input_342` carry the *outside* state at all four
diagnostics. They share **zero** derived columns with the seven. But by identity:

| | mean identity to the 7 | to 120 outside rows |
|---|---|---|
| input_437 | 0.9152 | 0.7900 |
| input_438 | 0.8959 | 0.7720 |
| input_342 | 0.9180 | 0.7833 |
| the 7, internally | **0.9939** | — |

**0.92 identity to a group that is 0.994 internally is still an exclusion.** Every
similarity threshold sweeps these three in; only the diagnostic columns keep them
out. This is the cleanest demonstration so far that MANUAL 6.1.6 is not a stylistic
preference but the operative criterion, and it is a calibration point with a number
attached: the gap that matters is within-group 0.994 against candidate 0.92, and
identity ranks the candidates correctly while still giving the wrong answer.

His note on why they exist is worth keeping: a chunk of 50 loci can be "composed of
sequences from adjacent subfamilies", so a chunk consensus can be genuinely
transitional rather than merely noisy. That is a property of the chunking, not of
the genome, and it means some chunks have no correct assignment.

They are kept in `CURATE__saq__s6_discarded_3.aln.fa` rather than deleted - he
called 437/438 a possible very minor subfamily.

### Second exact reproduction

Viewer-rule consensus of the 7, at thresholds 50/40/35/30 %: identity **1.0000**
to his `CONS__s6_7seqs` over all 250 columns, and the count matches the label
exactly (7 chunks, `s6_7seqs`). Threshold-insensitive again.

Running score: two peels, two consensuses reproduced byte-exactly.

| round | peeled | chunks | label says | consensus identity |
|---|---|---|---|---|
| 1 | s8 | 143 | 225 | 1.0000 over 254 |
| 2 | s6 | 7 | 7 | 1.0000 over 250 |

## Round 3 on saq: s5, the first case where the threshold mattered

His call: "s5 is amidst s4 sequences but they form a group, from s5 consensus down
to 296." In r3 that is rows 391-395: `input_386`, `input_387`, `input_421`,
`input_385`, `input_296` — five chunks, matching `s5_5seqs`.

### The boundary is confirmed by the only decisive test available

Consensus of the block against his `CONS__s5_5seqs`, by the viewer rule:

| block | thr 50 % | thr 40 % | thr 35 % | thr 30 % |
|---|---|---|---|---|
| all 5, his boundary | 0.9960 | **1.0000** | **1.0000** | **1.0000** |
| 4, dropping 296 | 0.9960 | 0.9960 | 0.9960 | 0.9960 |

The single difference is ungapped position 247, where his consensus has `T`.
`input_421` and `input_296` carry that `T` and the other three carry a gap. At
threshold 50 % two of five cannot call it, so it renders as a gap; at 40 % it
calls `T` and the match is exact. **Drop 296 and it is one of four, which no
threshold down to 30 % can call, so the block can never reproduce his consensus.**

Two conclusions. His boundary at 296 is right, and he must have lowered the
threshold on this block — which is exactly the option he described and the first
block where it was load-bearing. *Recording the threshold a block needed is
therefore not a diagnostic of messiness alone; it is part of the block's
definition.*

### Identity in the 0.92-0.93 band carries no information about membership

`input_296` is the identity outlier of its own group and sits closer to the
neighbouring s4 rows than to its own members:

| | to own group | to the 88 s4-block rows |
|---|---|---|
| input_386 | 0.9774 | 0.9412 |
| input_387 | 0.9774 | 0.9412 |
| input_421 | 0.9714 | 0.9410 |
| input_385 | 0.9734 | 0.9364 |
| **input_296** | **0.9300** | **0.9599** |

Set that beside the previous round:

| chunk | mean identity to the group | call | diagnostics carried |
|---|---|---|---|
| input_342 | 0.918 | **excluded** | 0 of 4 |
| input_296 | 0.930 | **included** | 3 of 3 unanimous |

Nearly the same identity, opposite calls, and identity puts 296 *further* from its
own group than 342 was from the one it was thrown out of. **No similarity threshold
can reproduce both calls.** What separates them is diagnostic carriage: 296 carries
every unanimous diagnostic of its group and merely dissents at two columns that are
consequently 4-of-5 (ungapped 100, where the outside is C 439/439, and 235); 342
carried the outside state at all four of s6's.

The operative membership rule, stated computably:

    member      <- carries the group's unanimous diagnostics
    non-member  <- carries the outside state at all of them
    identity is a ranking aid and nothing more

### Why he sees some diagnostics and not others — the colouring rule

He said his two positions "were directly indicated by colouring in current
alignment view". `_computeConservationForColumn` (script.js:4508) and
`applyConservationShadeClass` (6387): conservation = maxCount / nonGapCount, and a
base is shaded black/dark/light only if it is in the column's majority set AND
conservation clears the threshold. A minority base in a near-invariant column falls
through to `other` — pale against black.

Measured on the four s6 diagnostics:

| ungapped | group base | column conservation | visible? |
|---|---|---|---|
| 200 | T | 0.985 | yes — he named it |
| 225 | - | 1.000 | yes — he named it |
| 188 | G | 0.652 | no — mid-shaded, no contrast |
| 55 | A | 0.657 | no — mid-shaded, no contrast |

**The colouring reveals a diagnostic only when the column is otherwise
near-invariant.** A perfectly unanimous group diagnostic in a column that is
already biallelic across the family is invisible, however real it is. That is the
systematic blind spot of curation by eye, it is measurable, and it is the clearest
statement yet of what an automated discriminator adds: not better judgement on the
columns he sees, but the columns the shading cannot show him.

His own note on the gap after ungapped 55 — "also part informative considering
existence of edge cases or 50s forced to consensus" — is a caution on the same
columns: a gap there can be a chunk artefact rather than an indel, because a chunk
of 50 loci can be forced to a consensus it does not really have.

### Running score: three peels, three consensuses reproduced byte-exactly

| round | peeled | chunks | label says | identity | threshold needed |
|---|---|---|---|---|---|
| 1 | s8 | 143 | 225 | 1.0000 over 254 | any, 30-50 |
| 2 | s6 | 7 | 7 | 1.0000 over 250 | any, 30-50 |
| 3 | s5 | 5 | 5 | 1.0000 over 247 | **<= 40** |
