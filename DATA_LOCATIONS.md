# Where everything is

*Verified 2026-09-02 17:21 by `verify_locations.py`. Re-run it rather than trusting age.*

A location is **server plus absolute path**. A bare run ID is not a location -
recording run IDs without paths is exactly how the Timema runs were declared
lost while sitting on DRAGEN.

| what | server | path | state | note |
|---|---|---|---|---|
| human SINEderella run | `DRAGEN` | `/staging/tmp/sinederella_benchmark/runs/human/run_20260822_231705` | ok | benchmark run, Alu subfamilies |
| zebrafish SINEderella run | `DRAGEN` | `/staging/tmp/sinederella_benchmark/runs/zebrafish/run_20260822_231705` | ok | benchmark run |
| Timema SINEderella run (timb benchmark) | `DRAGEN` | `/staging/tmp/sinederella_benchmark/runs/timema/run_20260822_231705` | ok | the AnnoSINE-55 benchmark set |
| Timema curated v4 (tim, live) | `DRAGEN` | `/staging/tmp/timema_sines/v4_curated/run_20260823_103133` | ok | 597 chunks, his refined 14-subfamily partition - ground truth |
| Timema curated v3 | `DRAGEN` | `/staging/tmp/timema_sines/v3_curated/run_20260823_063444` | ok |  |
| Timema run_20260821_132226 | `DRAGEN` | `/staging/tmp/timema_sines/run_20260821_132226` | ok | 598 chunks, 59,636 members, his 8-group partition |
| Timema run_20260821_145119 | `DRAGEN` | `/staging/tmp/timema_sines/run_20260821_145119` | ok |  |
| hedgehog (eri) SINEderella run | `KIT` | `/data/W/toki/Genomes/Mammalia/Eulipotyphla/Erniacidae/run_20260820_221537` | ok | on KIT, not DRAGEN - this is why searches on DRAGEN found nothing |
| snail (pom) SINEderella run | `DRAGEN` | `/staging/tmp/newsp/pom/run_20260901_131644` | ok | both snail families |
| Tal subfam chunks (toc) | `KIT` | `/data/W/toki/Tal/genomes/toc/sf` | ok | 600 chunk consensuses + 600 member chunks of 50 = 30,000 copies |
| Tal consensus bank | `KIT` | `/data/W/toki/Tal/consi.bnk` | ok | TAL + the 6 groups he peeled by hand - ground truth for the peel loop |
| Tal starting query | `KIT` | `/data/W/toki/Tal/Tal.q` | ok | one consensus |
| Tal manual subfamilies | `KIT` | `/data/W/toki/SINE_disc/aln_c` | ok | POS__saq__s1-s9, ccr g1-g7+a_ccr, teu t1-t6, dmo d1-d5 |
| labelled corpus (674 alignments) | `KIT` | `/data/W/toki/SINE_disc/aln_c` | ok | raw, 400 bp flanks |
| labelled corpus, justified | `KIT` | `/data/W/toki/SINE_disc/aln_v2` | ok |  |
| human + Timema benchmark alignments | `KIT` | `/data/W/toki/SINE_disc/bench_in` | ok | 64 HUM + 55 TIMB |
| SINEderella toolkit | `DRAGEN` | `/staging/tmp/SINEderella` | ok | MANUAL.md, README.md, all step scripts, SubFam |
| SINE_consensus repo | `DRAGEN` | `/staging/tmp/SINE_consensus_repo` | ok | sine_consensus.sh, sine_pairwise_consensus.sh |
| conse | `DRAGEN` | `/usr/bin/conse` | ok | EMBOSS cons wrapper, 35 % plurality |
| SubFam | `KIT` | `/data/V/toki/bin/SubFam` | ok | also /usr/local/bin/SubFam on DRAGEN |
| his de novo scan | `KIT` | `/data/W/toki/scorpio/denovo_scan/sine_scan.sh` | ok |  |
| AnnoSINE env | `KIT` | `/data/V/toki/envs/annosine` | ok | built 2026-09-02 |
| AnnoSINE env | `DRAGEN` | `/staging/tmp/kmer_experiment/envs/annosine` | ok |  |
| discriminator working dir | `DRAGEN` | `/staging/tmp/sinedisc` | ok |  |
| discriminator working dir | `KIT` | `/data/W/toki/SINE_disc` | ok |  |
| scorpion genomes (9 species) | `KIT` | `/data/V/toki/Genomes/Scorpions` | ok | one best assembly per species, 14.6 GB |
| scorpion genomes (working copies) | `DRAGEN` | `/staging/tmp/scorpions` | ok |  |
| Timema genome | `KIT` | `/data/W/toki/Genomes/lower/Arthropoda/Timema/timema.fna` | ok | 1.24 GB, bwa-indexed |
| prospective genomes | `DRAGEN` | `/staging/tmp/newsp` | ok | hyd pom aca stu ska zeb + rebuilds |
| RepBase blast db | `KIT` | `/data/W/toki/SINE_disc/repbase.nhr` | ok | blastn -db repbase |
| SINEBase | `KIT` | `/data/W/toki/SINEbase` | ok |  |

## Connecting

- **DRAGEN** — `plink -ssh copilot@100.104.25.22 -batch -noagent -i <ppk>`;
  strip the 9-line banner with `tail -n +10`; never write to `/home`.
- **KIT** — one tunnel per session, then everything through it:
  `plink -batch -i <ppk> -N -L 2223:localhost:22 toki@85.89.102.78`, then
  `plink -batch -i <ppk> -P 2223 toki@127.0.0.1 '<cmd>'`. Never open a direct
  connection per command (fail2ban ban, 2026-08-10).
- `/data/W` on KIT is ~98 % full. Put new data on `/data/V`.

## Searching for something not listed

`find A -o B` without an explicit `-print` can return nothing while both
match. Always write `find ... \( -a -o -b \) -print`, or run two finds.
A negative result from a malformed search is not a negative result.
