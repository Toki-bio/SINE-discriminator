# Calibration cases for the peel rule

Each case: the peel loop split one of your groups; **X** is the piece it split off, **Y** the rest of your group.
Please mark each case **real** (X is its own group), **together** (X belongs with Y) or **unsure**.

In each alignment the X rows come first (`X|name`), then Y (`Y|name`). The two top rows are markers, not sequences:
`EVIDENCE_clean` has a letter only at columns where every X sequence has that character and no Y sequence does;
`EVIDENCE_1or2_exceptions` the same with one or two exceptions (named in the table below each case). Column numbers are the
viewer's ruler numbers for that file.

| case | your group | X | Y | clean substitutions | clean indels | with 1 exception | with 2 | distance in X / in Y / X-Y | your call |
|---|---|---|---|---|---|---|---|---|---|
| [01](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case01_ccr_g3.aln.fa&title=case%2001) | ccr g3 | 29 | 108 | 0 | 0 | 0 | 0 | 0.020 / 0.037 / 0.053 | |
| [02](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case02_ccr_g3.aln.fa&title=case%2002) | ccr g3 | 33 | 104 | 0 | 0 | 1 | 3 | 0.017 / 0.034 / 0.056 | |
| [03](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case03_ccr_g2.aln.fa&title=case%2003) | ccr g2 | 2 | 67 | 0 | 0 | 1 | 0 | 0.031 / 0.055 / 0.073 | |
| [04](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case04_oma_SINE21.aln.fa&title=case%2004) | oma SINE21 | 3 | 9 | 0 | 0 | 4 | 3 | 0.181 / 0.097 / 0.141 | |
| [05](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case05_ccr_g2.aln.fa&title=case%2005) | ccr g2 | 3 | 66 | 1 | 0 | 2 | 0 | 0.014 / 0.056 / 0.058 | |
| [06](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case06_oma_SINE24.aln.fa&title=case%2006) | oma SINE24 | 8 | 16 | 0 | 1 | 1 | 5 | 0.074 / 0.025 / 0.076 | |
| [07](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case07_oma_SINE22.aln.fa&title=case%2007) | oma SINE22 | 2 | 65 | 1 | 2 | 3 | 6 | 0.579 / 0.193 / 0.459 | |
| [08](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case08_oma_SINE21.aln.fa&title=case%2008) | oma SINE21 | 3 | 9 | 1 | 1 | 10 | 17 | 0.036 / 0.117 / 0.130 | |
| [09](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case09_oma_SINE5.aln.fa&title=case%2009) | oma SINE5 | 5 | 12 | 3 | 1 | 8 | 5 | 0.042 / 0.032 / 0.083 | |
| [10](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case10_oma_SINE25.aln.fa&title=case%2010) | oma SINE25 | 2 | 29 | 7 | 2 | 1 | 1 | 0.000 / 0.050 / 0.180 | |
| [11](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case11_oma_SINE18.aln.fa&title=case%2011) | oma SINE18 | 3 | 5 | 8 | 1 | 3 | 3 | 0.039 / 0.130 / 0.146 | |
| [12](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case12_oma_long12.aln.fa&title=case%2012) | oma long12 | 3 | 9 | 3 | 1 | 24 | 1 | 0.161 / 0.036 / 0.136 | |
| [13](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/peel_calib/case13_oma_SINE19.aln.fa&title=case%2013) | oma SINE19 | 4 | 6 | 9 | 2 | 29 | 56 | 0.358 / 0.379 / 0.506 | |

`~` in an EVIDENCE row marks a gap in X where Y has bases (a deletion in X).

## Columns per case

### Case 01: ccr g3, X = 29, Y = 108

- no column at all with 2 or fewer exceptions

### Case 02: ccr g3, X = 33, Y = 104

- column 211: gap (deletion in X) - exceptions: input_052.bnk (in X, has T); input_064.bnk (in X, has C)
- column 215: T - exceptions: input_052.bnk (in X, has C)
- column 217: C - exceptions: input_052.bnk (in X, has T); input_064.bnk (in X, has T)
- column 249: A - exceptions: input_115.bnk (in X, has C); input_116.bnk (outside X, has A)

### Case 03: ccr g2, X = 2, Y = 67

- column 73: A - exceptions: input_227.bnk (outside X, has A)

### Case 04: oma SINE21, X = 3, Y = 9

- column 141: A - exceptions: input_390.bnk (in X, has G)
- column 142: T - exceptions: input_390.bnk (in X, has C)
- column 220: A - exceptions: input_390.bnk (in X, has T); input_387.bnk (outside X, has A)
- column 222: T - exceptions: input_391.bnk (in X, has -)
- column 276: T - exceptions: input_390.bnk (in X, has A); input_389.bnk (outside X, has T)
- column 300: gap (deletion in X) - exceptions: input_386.bnk (outside X, has -); input_383.bnk (outside X, has -)
- column 301: gap (deletion in X) - exceptions: input_386.bnk (outside X, has -)

### Case 05: ccr g2, X = 3, Y = 66

- column 206: G - exceptions: input_325.bnk (outside X, has G)
- column 233: T - exceptions: input_317.bnk (outside X, has T)
- column 253: A - clean

### Case 06: oma SINE24, X = 8, Y = 16

- column 22: gap (deletion in X) - clean
- column 23: G - exceptions: input_553.bnk (in X, has -)
- column 25: A - exceptions: input_553.bnk (in X, has -); input_539.bnk (outside X, has A)
- column 37: A (Y mostly gapped: insertion in X) - exceptions: input_550.bnk (in X, has -); input_553.bnk (in X, has -)
- column 38: A (Y mostly gapped: insertion in X) - exceptions: input_550.bnk (in X, has -); input_553.bnk (in X, has -)
- column 39: A (Y mostly gapped: insertion in X) - exceptions: input_550.bnk (in X, has -); input_553.bnk (in X, has -)
- column 77: C - exceptions: input_550.bnk (in X, has T); input_553.bnk (in X, has -)

### Case 07: oma SINE22, X = 2, Y = 65

- column 2: gap (deletion in X) - exceptions: input_154.bnk (outside X, has -); input_378.bnk (outside X, has -)
- column 3: gap (deletion in X) - exceptions: input_154.bnk (outside X, has -); input_378.bnk (outside X, has -)
- column 4: gap (deletion in X) - exceptions: input_154.bnk (outside X, has -); input_378.bnk (outside X, has -)
- column 5: gap (deletion in X) - exceptions: input_154.bnk (outside X, has -); input_378.bnk (outside X, has -)
- column 117: A - clean
- column 137: T - exceptions: input_197.bnk (outside X, has T)
- column 156: gap (deletion in X) - clean
- column 157: gap (deletion in X) - exceptions: input_230.bnk (outside X, has -)
- column 225: A - exceptions: input_208.bnk (outside X, has A)
- column 241: gap (deletion in X) - clean
- column 242: gap (deletion in X) - clean
- column 243: gap (deletion in X) - exceptions: input_180.bnk (outside X, has -); input_226.bnk (outside X, has -)
- column 247: gap (deletion in X) - exceptions: input_226.bnk (outside X, has -); input_378.bnk (outside X, has -)

### Case 08: oma SINE21, X = 3, Y = 9

- column 226: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 232: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 233: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 234: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 235: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 236: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 237: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 238: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 239: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 240: gap (deletion in X) - clean
- column 241: gap (deletion in X) - clean
- column 257: A - exceptions: input_380.bnk (outside X, has A); input_391.bnk (outside X, has A)
- column 259: G - exceptions: input_380.bnk (outside X, has G)
- column 264: T - exceptions: input_380.bnk (outside X, has T)
- column 267: C - exceptions: input_386.bnk (in X, has G)
- column 268: G - exceptions: input_386.bnk (in X, has T); input_380.bnk (outside X, has G)
- column 273: A - exceptions: input_380.bnk (outside X, has A); input_391.bnk (outside X, has A)
- column 276: C - clean
- column 282: T - exceptions: input_380.bnk (outside X, has T); input_391.bnk (outside X, has T)
- column 284: G - exceptions: input_380.bnk (outside X, has G)
- column 287: G - exceptions: input_380.bnk (outside X, has G)
- column 289: T - exceptions: input_380.bnk (outside X, has T)
- column 290: T - exceptions: input_380.bnk (outside X, has T); input_383.bnk (outside X, has T)
- column 291: A - exceptions: input_380.bnk (outside X, has A)
- column 296: G - exceptions: input_386.bnk (in X, has -)
- column 297: A - exceptions: input_386.bnk (in X, has -); input_389.bnk (outside X, has A)
- column 299: T - exceptions: input_386.bnk (in X, has -)
- column 300: A - exceptions: input_386.bnk (in X, has -)
- column 303: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)
- column 304: gap (deletion in X) - exceptions: input_380.bnk (outside X, has -); input_391.bnk (outside X, has -)

### Case 09: oma SINE5, X = 5, Y = 12

- column 17: T - exceptions: input_139.bnk (in X, has C); input_141.bnk (in X, has C)
- column 38: T - exceptions: input_140.bnk (in X, has -); input_141.bnk (in X, has G)
- column 42: A - exceptions: input_140.bnk (in X, has G)
- column 73: T - clean
- column 156: A - clean
- column 157: A - clean
- column 167: gap (deletion in X) - exceptions: input_141.bnk (in X, has T)
- column 168: gap (deletion in X) - exceptions: input_141.bnk (in X, has T)
- column 169: gap (deletion in X) - exceptions: input_141.bnk (in X, has G)
- column 171: C - exceptions: input_141.bnk (in X, has T)
- column 190: gap (deletion in X) - exceptions: input_140.bnk (in X, has C); input_141.bnk (in X, has T)
- column 191: gap (deletion in X) - exceptions: input_140.bnk (in X, has A)
- column 192: gap (deletion in X) - clean
- column 193: gap (deletion in X) - clean
- column 194: gap (deletion in X) - clean
- column 195: gap (deletion in X) - clean
- column 197: gap (deletion in X) - exceptions: input_152.bnk (outside X, has -)
- column 198: gap (deletion in X) - exceptions: input_152.bnk (outside X, has -)
- column 199: gap (deletion in X) - exceptions: input_153.bnk (outside X, has -); input_152.bnk (outside X, has -)
- column 200: gap (deletion in X) - exceptions: input_153.bnk (outside X, has -); input_152.bnk (outside X, has -)

### Case 10: oma SINE25, X = 2, Y = 29

- column 33: T - exceptions: input_411.bnk (outside X, has T)
- column 38: gap (deletion in X) - clean
- column 39: gap (deletion in X) - clean
- column 40: gap (deletion in X) - clean
- column 41: T - exceptions: input_415.bnk (outside X, has T); input_395.bnk (outside X, has T)
- column 44: T - clean
- column 45: T - clean
- column 93: A - clean
- column 94: A - clean
- column 96: A - clean
- column 100: A - clean
- column 104: A - clean

### Case 11: oma SINE18, X = 3, Y = 5

- column 17: T - clean
- column 35: G - clean
- column 37: A - clean
- column 41: gap (deletion in X) - clean
- column 47: T - clean
- column 50: C - exceptions: input_286.bnk (outside X, has C); input_287.bnk (outside X, has C)
- column 53: T - exceptions: input_291.bnk (outside X, has T); input_293.bnk (outside X, has T)
- column 59: A - clean
- column 64: T - exceptions: input_289.bnk (in X, has A)
- column 66: T - exceptions: input_291.bnk (outside X, has T); input_293.bnk (outside X, has T)
- column 81: A - exceptions: input_293.bnk (outside X, has A)
- column 92: T - clean
- column 97: T - exceptions: input_290.bnk (in X, has C)
- column 131: A - clean
- column 155: C - clean

### Case 12: oma long12, X = 3, Y = 9

- column 41: C - exceptions: input_223.bnk (in X, has T)
- column 60: A - exceptions: input_224.bnk (in X, has G)
- column 79: gap (deletion in X) - exceptions: input_223.bnk (in X, has C)
- column 90: T (Y mostly gapped: insertion in X) - exceptions: input_214.bnk (outside X, has T)
- column 132: A - exceptions: input_214.bnk (outside X, has A)
- column 135: C - exceptions: input_214.bnk (outside X, has C)
- column 155: C - clean
- column 167: T - exceptions: input_214.bnk (outside X, has T)
- column 189: G - exceptions: input_213.bnk (in X, has -); input_214.bnk (outside X, has G)
- column 243: C - exceptions: input_214.bnk (outside X, has C)
- column 398: A - exceptions: input_213.bnk (in X, has G)
- column 402: A - exceptions: input_223.bnk (in X, has G)
- column 420: T - exceptions: input_223.bnk (in X, has A)
- column 421: A - clean
- column 422: T - clean
- column 433: gap (deletion in X) - exceptions: input_223.bnk (in X, has T)
- column 434: gap (deletion in X) - exceptions: input_223.bnk (in X, has T)
- column 435: gap (deletion in X) - exceptions: input_223.bnk (in X, has T)
- column 436: gap (deletion in X) - exceptions: input_223.bnk (in X, has T)
- column 437: gap (deletion in X) - exceptions: input_223.bnk (in X, has A)
- column 438: gap (deletion in X) - exceptions: input_223.bnk (in X, has T)
- column 439: gap (deletion in X) - clean
- column 440: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 441: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 442: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 443: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 444: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 445: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)
- column 446: gap (deletion in X) - exceptions: input_214.bnk (outside X, has -)

### Case 13: oma SINE19, X = 4, Y = 6

- column 8: gap (deletion in X) - exceptions: input_302.bnk (outside X, has -); input_306.bnk (outside X, has -)
- column 9: gap (deletion in X) - clean
- column 10: gap (deletion in X) - clean
- column 11: gap (deletion in X) - clean
- column 12: gap (deletion in X) - clean
- column 13: gap (deletion in X) - clean
- column 14: gap (deletion in X) - clean
- column 15: gap (deletion in X) - clean
- column 16: gap (deletion in X) - clean
- column 17: gap (deletion in X) - clean
- column 18: gap (deletion in X) - clean
- column 19: gap (deletion in X) - clean
- column 20: gap (deletion in X) - clean
- column 21: gap (deletion in X) - exceptions: input_272.bnk (in X, has A)
- column 22: gap (deletion in X) - exceptions: input_272.bnk (in X, has T)
- column 36: T - exceptions: input_272.bnk (in X, has -)
- column 38: C (Y mostly gapped: insertion in X) - clean
- column 39: A (Y mostly gapped: insertion in X) - exceptions: input_302.bnk (outside X, has A); input_303.bnk (outside X, has A)
- column 40: G (Y mostly gapped: insertion in X) - exceptions: input_302.bnk (outside X, has G); input_303.bnk (outside X, has G)
- column 44: C - exceptions: input_306.bnk (outside X, has C)
- column 45: A - exceptions: input_302.bnk (outside X, has A)
- column 46: A (Y mostly gapped: insertion in X) - exceptions: input_306.bnk (outside X, has A); input_307.bnk (outside X, has A)
- column 47: T (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has -)
- column 48: T (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has -)
- column 49: gap (deletion in X) - exceptions: input_315.bnk (in X, has T); input_316.bnk (in X, has A)
- column 52: A - exceptions: input_272.bnk (in X, has -)
- column 68: T (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has A); input_317.bnk (in X, has -)
- column 69: A (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has T); input_317.bnk (in X, has -)
- column 70: A (Y mostly gapped: insertion in X) - exceptions: input_315.bnk (in X, has T); input_317.bnk (in X, has -)
- column 71: A (Y mostly gapped: insertion in X) - exceptions: input_317.bnk (in X, has -); input_302.bnk (outside X, has A)
- column 72: T (Y mostly gapped: insertion in X) - exceptions: input_317.bnk (in X, has A); input_302.bnk (outside X, has T)
- column 75: T - exceptions: input_272.bnk (in X, has A); input_316.bnk (in X, has -)
- column 77: A - exceptions: input_272.bnk (in X, has T); input_316.bnk (in X, has -)
- column 81: A - exceptions: input_272.bnk (in X, has T); input_317.bnk (in X, has G)
- column 82: A - exceptions: input_272.bnk (in X, has T)
- column 85: T - exceptions: input_272.bnk (in X, has C)
- column 90: A - clean
- column 97: A - exceptions: input_272.bnk (in X, has G); input_306.bnk (outside X, has A)
- column 98: G - clean
- column 102: T - exceptions: input_272.bnk (in X, has A)
- column 109: C - clean
- column 120: A (Y mostly gapped: insertion in X) - exceptions: input_303.bnk (outside X, has A)
- column 121: A (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has G); input_315.bnk (in X, has C)
- column 122: A (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has T)
- column 123: T (Y mostly gapped: insertion in X) - exceptions: input_315.bnk (in X, has A); input_316.bnk (in X, has A)
- column 124: T (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has G); input_315.bnk (in X, has C)
- column 125: A (Y mostly gapped: insertion in X) - exceptions: input_315.bnk (in X, has G)
- column 126: T (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has -)
- column 131: gap (deletion in X) - exceptions: input_272.bnk (in X, has T)
- column 132: gap (deletion in X) - exceptions: input_272.bnk (in X, has A)
- column 133: gap (deletion in X) - exceptions: input_272.bnk (in X, has A)
- column 134: gap (deletion in X) - exceptions: input_272.bnk (in X, has A)
- column 135: gap (deletion in X) - exceptions: input_272.bnk (in X, has A)
- column 137: gap (deletion in X) - exceptions: input_272.bnk (in X, has A); input_315.bnk (in X, has T)
- column 146: C - exceptions: input_272.bnk (in X, has T); input_317.bnk (in X, has G)
- column 148: T - exceptions: input_317.bnk (in X, has A); input_302.bnk (outside X, has T)
- column 149: T - clean
- column 158: gap (deletion in X) - exceptions: input_303.bnk (outside X, has -); input_307.bnk (outside X, has -)
- column 159: gap (deletion in X) - exceptions: input_303.bnk (outside X, has -); input_307.bnk (outside X, has -)
- column 160: A - exceptions: input_315.bnk (in X, has G); input_316.bnk (in X, has G)
- column 161: A - exceptions: input_272.bnk (in X, has G); input_307.bnk (outside X, has A)
- column 169: T - exceptions: input_272.bnk (in X, has A)
- column 179: T - exceptions: input_307.bnk (outside X, has T)
- column 182: gap (deletion in X) - exceptions: input_272.bnk (in X, has T); input_307.bnk (outside X, has -)
- column 183: gap (deletion in X) - exceptions: input_272.bnk (in X, has C); input_307.bnk (outside X, has -)
- column 184: gap (deletion in X) - exceptions: input_272.bnk (in X, has T); input_307.bnk (outside X, has -)
- column 185: gap (deletion in X) - exceptions: input_272.bnk (in X, has T); input_307.bnk (outside X, has -)
- column 197: A (Y mostly gapped: insertion in X) - exceptions: input_315.bnk (in X, has -); input_316.bnk (in X, has -)
- column 202: T - exceptions: input_272.bnk (in X, has A)
- column 204: T - exceptions: input_272.bnk (in X, has C)
- column 210: A - exceptions: input_272.bnk (in X, has T); input_317.bnk (in X, has G)
- column 215: G - clean
- column 220: A - exceptions: input_302.bnk (outside X, has A); input_305.bnk (outside X, has A)
- column 222: A - exceptions: input_272.bnk (in X, has T); input_305.bnk (outside X, has A)
- column 226: C - clean
- column 234: T - clean
- column 242: T (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has -)
- column 243: T (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has -); input_317.bnk (in X, has A)
- column 251: A - exceptions: input_316.bnk (in X, has G); input_304.bnk (outside X, has A)
- column 259: gap (deletion in X) - exceptions: input_272.bnk (in X, has T); input_315.bnk (in X, has G)
- column 272: T (Y mostly gapped: insertion in X) - exceptions: input_304.bnk (outside X, has T); input_306.bnk (outside X, has T)
- column 273: T (Y mostly gapped: insertion in X) - exceptions: input_272.bnk (in X, has C)
- column 274: A (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has -)
- column 280: T - clean
- column 281: G - exceptions: input_306.bnk (outside X, has G); input_305.bnk (outside X, has G)
- column 286: G - clean
- column 297: T - exceptions: input_315.bnk (in X, has A); input_317.bnk (in X, has A)
- column 298: A - exceptions: input_317.bnk (in X, has G); input_305.bnk (outside X, has A)
- column 304: T - exceptions: input_302.bnk (outside X, has T)
- column 308: gap (deletion in X) - exceptions: input_317.bnk (in X, has G); input_307.bnk (outside X, has -)
- column 309: A - exceptions: input_303.bnk (outside X, has A)
- column 310: T - exceptions: input_302.bnk (outside X, has T); input_305.bnk (outside X, has T)
- column 311: T - exceptions: input_302.bnk (outside X, has T); input_303.bnk (outside X, has T)
- column 314: A - exceptions: input_302.bnk (outside X, has A)
- column 316: gap (deletion in X) - exceptions: input_317.bnk (in X, has A); input_307.bnk (outside X, has -)
- column 318: A - exceptions: input_304.bnk (outside X, has A)
- column 322: G - exceptions: input_272.bnk (in X, has A); input_304.bnk (outside X, has G)
- column 325: G - exceptions: input_304.bnk (outside X, has G); input_306.bnk (outside X, has G)
- column 327: A - exceptions: input_304.bnk (outside X, has A); input_306.bnk (outside X, has A)
- column 329: A - exceptions: input_315.bnk (in X, has G); input_316.bnk (in X, has G)
- column 330: A - exceptions: input_272.bnk (in X, has G); input_302.bnk (outside X, has A)
- column 332: A - exceptions: input_316.bnk (in X, has T)
- column 339: A (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has T); input_303.bnk (outside X, has A)
- column 340: T (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has A)
- column 343: T (Y mostly gapped: insertion in X) - exceptions: input_315.bnk (in X, has A); input_316.bnk (in X, has -)
- column 344: T (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -); input_317.bnk (in X, has A)
- column 345: A (Y mostly gapped: insertion in X) - exceptions: input_316.bnk (in X, has -)
