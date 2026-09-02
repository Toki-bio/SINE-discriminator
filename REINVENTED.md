# Tools of his I rewrote instead of calling

Asked directly on 2026-09-02: "maybe there are more reinvented scripts that you
are missing?" Yes. This is the audit, from reading his scripts on KIT
(`/data/V/toki/bin/`) and DRAGEN (`/staging/tmp/SINEderella/`), not from
guessing by filename.

**Check this file before writing any new script.**

## Confirmed duplicates

| his tool | what it does (his own header) | what I wrote instead |
|---|---|---|
| **`sear`** | top-N by bitscore, ±50 bp flanks, `getfasta -s`, query as row 1, `mafft --localpair --maxiterate 1000 --ep 0.123`; ssearch36, hit ≥80 % of query length at ≥65 % identity | `candidate_to_aln.py` — `blastn -word_size 11 -evalue 1e-10`, hand-written overlap dedupe, `FLANK = 400` |
| **`CutByCons`** | "Trims multiple fasta sequences from both ends according to a query" — 23 lines | `boundary.py`, `fix_alignments.py`, `trim_flanks.py`, `copy_boundary.py` |
| **`Trim.sh`, `Trim53.sh`** | trims a locus at a position *relative to the consensus* in pairwise alignment | the same four scripts again |
| **`100sim.sh`** | "takes consensus and a bank of 100 randomly selected SINEs… sample, align with MAFFT, trim to consensus coords, calculate avg identity" | `measure_c.py` identity block + the rand100 path in `candidate_to_aln.py` |
| **`SimilarCopies`** | query → most similar copies at 0.99 length / 99 % identity, caps at 300, returns an MSA | the top100/rand100 views in `candidate_to_aln.py` |
| **`PairClust.sh`** | "defines clusters of multiple or pairwise similarity in a list of pairs" — single-linkage from a pair list | `single_link()` in `peel_step1.py`; `cluster_chunks.py`, `subfam_cluster.py` |
| **`SubFam`, `SubSort`** | chunked consensus clustering | `subfam_group.py`, `subfam_peel.py` |
| **`cluster_assist.js` / `SINEClusterer`** | peel loop on shared (position, character) features | `peel_loop.py`, `peel_features.py` |
| **`conse`** | EMBOSS cons at pct % plurality, N→`-` | consensus calls inside `measure_c.py`, `consensus_check.py` |
| **`sine_consensus.sh`, `_smart`, `_pairwise`** | bootstrapped / variance-flagged / pairwise-projected consensus | `build_consensuses.sh`, `rebuild_consensus.py`, `consensus_two_stage.py`, `extend_consensus.py`, `refine_iterate.py`, `refine_pass2.py` |
| **`ReadExtend`** | "iteratively extend short reads coverage given a starting sequence" | `extend_consensus.py` |
| **`step7_boundary_refine.sh`** | empirical background from random genomic regions, threshold bg + 5 %, walked per subfamily × side × population | `edges.py`, `edge_probe.py`, `edge_quality.py`, `boundary.py` |
| **`step8a_extract_alignments.sh`** | top100 / rand100 / subfam with step7 flanks, consensus forced to row 1 | `build_viewer_data.py`, `newsp_all.py` |
| **`postprocess_flanks()`** in `extract_alignments.sh` | flank gaps deleted, lowercased, butted to the element edge | `justify_all.py`, `fix_alignments.py` |
| **`tsd1.sh` … `tsd33.sh`** | TSD detection via ssearch36 | the `tsd_frac` computation in `verdict.py` |
| **`trf`, `trf2bed`, `UnitFinder.sh`, `SatComp.sh`** | tandem repeat finding; "represent locus sequence as a graphical scheme of its parts" | `microsat.py`, `segmap.py`, `mosaic_kmer.py` |

## Not duplicates — his tools with no equivalent of mine

Worth knowing they exist, because they cover things the discriminator claims are
missing or partial:

- `Sindel`, `Sindel2022`, `SindelAll` — cross-genome dimorphic indel comparison
- `SINE_orth_loc.bash`, `ComPair.sh`, `loccomp.sh` — orthologous locus comparison, PM/MP/SINE
- `SINE_walker`, `sine_walker.py`, `LINE_walker`, `line_walker.py` — walking elements
- `sine_scan.sh` + `SINE-de-novo-genome-scan/` — library-free de novo search
- `coordinates_by_consensus*.sh` — coordinate mapping through a consensus
- `dotplot`, `FaBest`, `FaSort`, `sample`, `liner`, `refa`, `msf`, `tribes`
- `maft`, `mafta`, `maftm`, `maftr` — his mafft wrappers
- `ssearch36tofasta*`, `ssearch36_gaps.awk`, `ssearch36_rename*.awk` — ssearch36 parsers

Also: many `sear` variants exist (`sear100`, `sear100blast`, `sear1step`,
`sear2k`, `sear2kloop*`, `searFULLblast`, `seart*`, `sea`) — check which variant
fits before assuming plain `sear`.

## The pattern

Every one of these was written because I did not look first. The cost is not
duplicated effort — it is that my version is worse in ways I do not notice:
`candidate_to_aln.py`'s `FLANK = 400` is why `conse` calls a 541-base consensus
for a 253 bp element, which is why the element window is wrong, which is why the
flanks render ragged and the mosaic statistic saturates. One unchecked
reimplementation produced every symptom in that chain.
