# SINE discriminator

Formalising the manual "is this a real SINE family?" judgement into measurable
statistics. Report: **https://toki-bio.github.io/SINE-discriminator/**

- `HANDOFF.md` — project state, read this first
- `FINDINGS.md` — chronological record (§3–§5 partly retracted by §11)
- `SINE_discriminator_spec.md` — the original design
- `measure_c.py` — MEASURE(), consensus-anchored
- `profiles.py` — per-nucleotide positional tracks
- `alignments/` — 200 consensus-anchored MAFFT alignments, viewable in
  [MSA-viewer](https://toki-bio.github.io/MSA-viewer/)
- `kit/` — the extraction and alignment scripts that run on KIT

Corpus: 28 curated Talpidae subfamilies across four genomes, plus synthetic
negatives and a 400-set resampling side quest.
