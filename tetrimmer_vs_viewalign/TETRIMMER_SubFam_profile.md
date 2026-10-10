# TEtrimmer, read at the source: what it does, and how it relates to SubFam

Source-level profile of TEtrimmer v1.7.4 (https://github.com/qjiangzhao/TEtrimmer, commit 8954274,
2026-07-24; paper: Qian et al. 2025 Nat Commun, doi:10.1038/s41467-025-63889-y). File paths are
under `tetrimmer/`. Everything is read from the code; the lines that decide the comparison were
checked twice (`MSAcluster.py:366`, `TEtrimmer.py:100`, `functions.py:445-446`, `MSAcluster.py:584`,
no `random.seed` anywhere). Items marked *inferred* are deductions.

## 1. What TEtrimmer is

A library-curation tool: for each consensus of a de novo library (RepeatModeler/EDTA) it finds the
copies in a genome, extends them to the true element boundaries, cleans the alignment, builds a
consensus, classifies it and grades it (Perfect / Good / Reco_check / Need_check). Its clustering
step exists to keep a consensus from being a chimera of variants, not to enumerate subfamilies.

Entry: `TEtrimmer.py:416 main`. `--genome_file` is `required=True` (`TEtrimmer.py:100`). `--input_file`
is optional only because RepeatModeler2 is run de novo without it (`:774-790`). **There is no entry
point for a FASTA of copies or a pre-made alignment**: `clean_and_cluster_MSA` takes an `input_msa`
but still needs a BED of genome coordinates (`MSAcluster.py:426-437`), and every later step
re-extracts from the genome (`bedtools slop/getfasta`, `functions.py:568-576`).

## 2. Pipeline per input consensus (`analyze.py:1080-1689`)

1. Type-specific settings (`:1155-1177`): SINE → `ex_step=200, max_ext=1400, min_blast_len=50`.
2. `blastn -max_target_seqs 10000 -evalue 1e-40 -qcov_hsp_perc 15`, HSP length `> min_blast_len`
   (`functions.py:253-300`); no identity filter; one HSP = one BED line.
3. Fewer than `min_seq_num` (10) hits → low-copy path and return.
4. **Cap**: if more than `--max_msa_lines` (100) HSPs, keep the `--top_msa_lines` (100) **longest**
   (`functions.py:445-446`), then fill randomly — unseeded (`random.shuffle`, no `random.seed` in the
   repo). With defaults the MSA is "the 100 longest HSPs".
5. MAFFT `--quiet --nuc --retree 1` (`functions.py:679-683`), no `--thread`; MUSCLE is tried first
   but is not in the conda env, so MAFFT runs (*inferred*).
6. **Clustering** (section 3) → per-cluster BED files.
7. Per cluster: iterative extension (`boundarycrop.py:109-273`), re-alignment, column cleaning,
   boundary definition by MSA conservation + genome BLAST coverage + poly-A / terminal repeats
   (`boundarycrop.py:406-593`), a second clustering pass that keeps only the largest cluster
   (`boundarycrop.py:1019-1049`), consensus, grading, PDFs.

## 3. The clustering step (`MSAcluster.py:483-627`)

- Column cleaning: any base below 2 % within its column becomes a gap (`:81-100`).
- Only if ≥ 20 sequences (`:552`): keep **divergent columns** — major allele ≤ 0.80 (≤ 0.85 when ≥ 90
  sequences) (`:112-122`) — plus conserved gap blocks when ≥ 30 such columns exist (`:138-146`).
  Clustering is attempted only when more than max(50, 5 % of the length) divergent columns survive
  (`:158-162`; 100 when the alignment is over 2,000 columns).
- IQ-TREE on that column subset, model `K2P+I` (`:239-248`), single thread; pairwise **patristic
  distance** from the tree (`:193-226`).
- `DBSCAN(eps=0.1, min_samples=3, metric='precomputed')` (`:366`, scikit-learn). eps and
  min_samples are not on the command line.
- Acceptance (`:283-338`): a cluster needs ≥ `min_seq_num` (10) members; clusters are ranked by size
  and the top `--max_cluster_num` (5) kept; DBSCAN noise is kept as one more cluster only if ≥ 12
  sequences and ≥ 60 % of all, else dropped. If the largest cluster has < 18 members the whole
  clustering is discarded (`:584`).
- Fallback (`:599-612`): one cluster on the full gap-filtered alignment with outliers removed.

So for a ~190 bp SINE the gate "more than 50 divergent columns" is rarely met and the fallback
(single cluster) is the likely path (*inferred*).

## 4. Consensus and outputs

- Consensus: per column, count A/C/G/T only, winner if its frequency ≥ threshold else `N`
  (`functions.py:54-152`); final threshold `cons_thr` 0.7 (0.65 for divergent inputs)
  (`boundarycrop.py:1877-1878`). Gap columns were removed earlier, so consensuses are gap-free.
- `--gap_thr` and `--gap_nul_thr` are accepted but never used (`boundarycrop.py:283-284, 948-949`).
- Files: `TEtrimmer_consensus.fasta` (header `>{name}#{type}\tTEtrimmer\t{grade}`; cluster k of
  input X is `X_0k`, `boundarycrop.py:800-805`; an input whose clusters all fail is written back
  unchanged, `analyze.py:1637-1641`); `TEtrimmer_consensus_merged.fasta` after cd-hit-est 0.9 then
  0.85 (`analyze.py:552-723`); `summary.txt`; per-consensus final MSA `{name}#{type}.fa` with headers
  `chr:start-end(strand)` (the copy→cluster map), `.cluster.fa` (the first-round 100-copy MSA),
  with `--debug` the raw `*_g_{i}.bed` DBSCAN clusters and the `_pat_MSA.fa.treefile`.
- Dependencies (`TEtrimmer_env.yml`): bedtools, blast, cd-hit, emboss, hmmer, iqtree, mafft,
  pfam_scan + Pfam-A (mandatory, auto-download), RepeatMasker, RepeatModeler, samtools, trf, R,
  scikit-learn, biopython. One process per input consensus; every external tool single-threaded.
  README runtime: 48 cores, D. rerio 3,504 families → 2.3 h.

## 5. Hard-coded thresholds that matter for a comparison

| parameter | default | where | effect |
|---|---|---|---|
| evalue / qcov_hsp | 1e-40 / 15 % | functions.py:261-262 | copy discovery |
| min_blast_len | 150 (SINE 50) | analyze.py:1165-1175 | HSP length |
| min_seq_num | 10 | TEtrimmer.py:582 | min hits, min cluster |
| max_msa_lines / top_msa_lines | 100 / 100 | functions.py:454-495 | **copies per MSA, longest first** |
| sequences to attempt clustering | 20 | MSAcluster.py:552 | |
| divergent column | major allele ≤ 0.8 (0.85 if ≥ 90 seqs) | MSAcluster.py:112-122 | feature columns |
| min divergent columns | max(50, 5 % L); 100 if L > 2000 | MSAcluster.py:158-162 | gate |
| tree model | K2P+I | MSAcluster.py:242 | distance |
| DBSCAN | eps 0.1, min_samples 3 | MSAcluster.py:366 | **not exposed** |
| max_cluster_num | 5 | TEtrimmer.py:234-241 | **cluster cap** |
| largest cluster must be ≥ | 18 | MSAcluster.py:584 | else no clustering |
| noise kept as cluster | ≥ 12 and ≥ 60 % | MSAcluster.py:309 | |
| cons_thr | 0.7 / 0.65 | boundarycrop.py:1878 | consensus |
| cd-hit-est merge | 0.9 then 0.85 | analyze.py:567-672 | merged library |
| random seed | none | — | fill is irreproducible |

## 6. Relation to SubFam

| | TEtrimmer | SubFam |
|---|---|---|
| purpose | one curated, boundary-correct consensus per library entry | the variety of a family on one screen, and a proxy alignment |
| input | genome (required) + seed consensus | any FASTA of related sequences |
| copies used | ≤ 100 HSPs, the longest | all |
| grouping | DBSCAN on patristic distance of divergent columns; eps 0.1 fixed; ≤ 5 clusters; floors 10 / 18 / 20 | consecutive chunks of N along a k-mer UPGMA order; every copy in a chunk |
| copies without a group | dropped (noise) | none |
| consensus | A/C/G/T majority ≥ 0.7, gap-free | EMBOSS plurality 0.36, gap columns dropped |
| order / variety | clusters ranked by size, unordered | the order is the output |
| determinism | unseeded fill; IQ-TREE seed time-based (*inferred*) | deterministic |
| what it has that SubFam lacks | boundary extension, TSD / poly-A / terminal-repeat logic, grading, PDFs | — |

The overlap is one auxiliary step: splitting a 100-copy alignment into at most five groups to
protect the consensus. SubFam starts from the whole family and does not pick groups. For the paper,
one paragraph (RELATED_WORK §A, PLAN §5.4) states this; a run is optional.

## 7. Recipe for a fair comparison on mouse B2 (if wanted)

Inputs: a mouse genome FASTA (GRCm39, a single chromosome is acceptable), a seed library with
headers `>B2_Mm1a#SINE/B2` etc. (parsed at `analyze.py:372-383`; type SINE sets the SINE
parameters), Pfam-A available.

```
conda env create -f TEtrimmer_env.yml && conda activate TEtrimmer
python tetrimmer/TEtrimmer.py -i B2_seeds.fa -g GRCm39.fa -o out_B2 -t 8 \
    --max_msa_lines 100 --top_msa_lines 100 --max_cluster_num 5 --debug --export_coverage
```

Read: `out_B2/TEtrimmer_consensus.fasta` (≤ 5 entries per seed), `summary.txt`
(`output_MSA_seq_n`), `TEtrimmer_for_proof_curation/Annotations_*/<name>.fa` (final cluster MSA;
headers give the copy→cluster map by genomic interval — compare by overlap, coordinates are
post-extension), with `--debug` `Multiple_sequence_alignment/*_g_{i}.bed` (raw DBSCAN clusters).
Then align its consensuses into the SubFam B2 example as `PUBLISHED_TEtrimmer_<k>` rows.

To use more copies or clusters: `--max_msa_lines N --top_msa_lines N --max_cluster_num K` have no
upper bound in code (`functions.py:468`, `MSAcluster.py:321`); the floors (10, 18, 20), the divergent-
column gate and eps 0.1 / min_samples 3 / K2P+I cannot be changed without editing `MSAcluster.py`.
A "no genome" run by concatenating copies with N-spacers into a fake genome is possible, but BLAST
selection, longest-first capping and extension would still reshape the copy set (*inferred*).
