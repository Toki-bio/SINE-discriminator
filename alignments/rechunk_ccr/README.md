# ccr locus-level re-chunk (his plan, 2026-10-05)

Source: the 50-locus member files `input_NNN.bnk` of the ccr SubFam run on KIT (`/data/W/toki/Genomes/Mammalia/Eulipotyphla/ccr/run_20260517_205955`),
for two of his groups: new13 (13 chunks, 650 loci) and g5 (5 chunks, 250 loci). All loci of a group re-aligned together
(`mafft --localpair --maxiterate 1000 --ep 0.123 --nuc`), then re-chunked by `MSA-viewer/tests/kmer/rechunk.js` (no outliers found; loci ordered by
similarity; 50 per chunk). Files: `<g>.chunks.aln.fa` = the new chunk consensi, `<g>.loci.tsv` = which new chunk each locus went to
(locus names are `<old chunk>|<genome coordinates>`), `new13.loci.aln.fa` = the re-aligned loci.
