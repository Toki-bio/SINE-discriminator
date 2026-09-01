# Automating the SubFam group-picking loop

## The manual process being replaced

Sergei's current workflow, in his words:

> initial search merged into multitude of all genomic hits then requires manual
> selection of sine candidates in subfam alignment... I do this by looking at
> alignment, pick up groups of similar sequences, copy-pasting them to separate
> alignment and if it looks good then after realignment I create their consensus
> and use it for future re-scan and refinement in next sinederella run. After
> selecting more obvious groups in subfam, I continue with removing them from
> subfam alignment and continue with remaining ones, until no groups remain
> (only noisy ungroupable sequence).

So: **cluster → inspect → accept → derive consensus → peel → repeat**, with the
"inspect" step being expert eyeballing. That eyeballing is precisely what the
discriminator was built to replace.

## Sergei's refinement, which is the important part

> after subfam inspection if the group of sequences is selected, the algorithm
> should not create consensus from these "consensi of 50s" used in subfam but
> look back into those original 50sequence chunks used for subfam consensuses.

A consensus built from chunk-consensuses is a consensus of consensuses: it
inherits each chunk's smoothing and loses the real variation. The members are
still on disk, so there is no reason to accept that loss.

## Verified data structure

`subfam_work/` (and per-run `genome.clean_step1/subfam_input/`) contains:

```
input_007.bnk        50 original sequences, headers carrying genomic coordinates
                     e.g. Lm1_5p_1-50::CM115237.1:112020019-112020166(+)
input_007.bnk.cons   the chunk consensus, header ">input_007.bnk"
```

- 162 chunks, 161 chunk-consensuses, **exactly 50 sequences per chunk**,
  16,042 member sequences total (Timema, `run_20260821_132226`).
- **chunk-consensus → its 50 members is a filename lookup.** Nothing needs to be
  reconstructed.
- Members carry genomic coordinates, so they can be **re-extracted from the
  genome with real flanks** — which the discriminator needs and the raw
  hit-regions (~147 bp, element only) do not have.
- Per-run copies persist: 604–607 `.bnk` files in each of the two 2026-08-21
  runs.

## Why this can be validated rather than guessed

`tim/` is ground truth for this exact task. Sergei clustered these **same 161
chunk-consensuses by hand into 8 groups** (`t1`–`t8`), later refined to 6
(`t1/t2/t345/t6/t7/t8`) and then to the current 14-subfamily v4. So automated
clustering has a known target: it should recover something close to that
partition, from the same input.

That makes step 1 a measurable experiment, not a design preference.

## The algorithm

```
INPUT: chunk consensuses (161) + their member chunks + genome

repeat:
  1. CLUSTER the remaining chunk-consensuses by pairwise identity
  2. rank clusters by tightness x size  ("more obvious groups" first)
  3. for the top cluster:
       a. gather the union of member sequences from its chunks
       b. re-extract those loci from the genome WITH 400bp flanks
       c. align (MAFFT), derive a consensus from the MEMBERS
       d. SCORE with the discriminator, using that consensus as row 1
       e. if accepted -> emit consensus + copy list; remove those chunks
          if rejected -> mark the cluster unusable; remove it from consideration
  4. stop when no remaining cluster is accepted
OUTPUT: a set of consensuses ready for the next SINEderella rescan,
        plus the residue that formed no acceptable group
```

Two things make this more than a clustering script:

- **The discriminator is the accept oracle.** "If it looks good" becomes a
  score plus named flags, on an alignment with real flanks — which is strictly
  more information than the flankless subfam view Sergei is eyeballing now.
- **Consensus comes from members, not chunk-consensuses**, per his point above.

## Build order

1. **Cluster the 161 and compare to `t1`–`t8`.** Pure measurement, no genome
   needed. If automated clustering cannot approximately recover his partition,
   the rest is not worth building — so this is the go/no-go.
2. **Member gathering + flank re-extraction** for one accepted cluster; confirm
   the discriminator scores it as a real family.
3. **The peel loop**, with a stopping rule and a cap.
4. **Compare end to end**: run the whole loop on the Timema chunks and compare
   the emitted consensuses against `t1`–`t8` and against v4's 14.

## Open questions for Sergei

These change the design and I would rather ask than assume:

1. **Accept threshold.** Should a cluster be accepted on discriminator score
   alone (and at what value), or do you want specific flags to be
   disqualifying regardless of score — e.g. `HETEROGENEOUS_SELECTION`, which
   would suggest the cluster should be split rather than accepted?
2. **Granularity.** Your first pass gave 8 groups, v4 gives 14. Should the loop
   aim for the coarse split (and leave subfamily separation to a later pass), or
   go straight to the finer one?
3. **Residue.** When no cluster is acceptable, do you want the leftovers dumped
   for inspection, or is "noise" a sufficient answer?

None of these block step 1, which is why step 1 starts now.
