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
