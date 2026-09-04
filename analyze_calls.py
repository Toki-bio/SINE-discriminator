import io, csv, re
from collections import Counter, defaultdict

R = list(csv.DictReader(io.open('C:/work/SINE_discriminator/calls.tsv', encoding='utf-8'), delimiter='\t'))

# ---- 1. did the synthetic negatives work as negatives? -----------------------
NEGP = ('NEGSPLICE','NEGCHIM','NEGTRUNC5','NEGMOSAIC','NEGLINE','NEGLINEORF',
        'NEGRAND','NEGLOWBIT','NEGSEGDUP','NEGSAT')
print("=" * 78)
print("1. SYNTHETIC NEGATIVES: did construction predict his verdict?")
print("=" * 78)
tot_neg = 0
agree = 0
per = defaultdict(Counter)
for r in R:
    p = r['set'].split('__')[0]
    if p in NEGP:
        tot_neg += 1
        per[p][r['call']] += 1
        if r['call'].startswith('NOT_SINE') or r['call'] in ('UNUSABLE_UNALIGNED',):
            agree += 1
print("  sets built to be negative: %d" % tot_neg)
print("  of those he actually called negative/unusable: %d (%.0f%%)" % (agree, 100.0*agree/tot_neg))
print()
for p in sorted(per, key=lambda z: -sum(per[z].values())):
    tot = sum(per[p].values())
    neg = sum(v for k, v in per[p].items() if k.startswith('NOT_SINE') or k == 'UNUSABLE_UNALIGNED')
    print("  %-12s n=%-3d  he called negative: %d   %s" % (p, tot, neg, dict(per[p])))

# ---- 2. where do the REAL negatives come from? ------------------------------
print()
print("=" * 78)
print("2. WHERE THE USABLE NEGATIVES ACTUALLY COME FROM")
print("=" * 78)
neg = [r for r in R if r['call'].startswith('NOT_SINE')]
print("  %d NOT_SINE* verdicts:" % len(neg))
for r in neg:
    print("    %-34s %-22s %s" % (r['set'][:34], r['call'], r['verbatim'][:70]))

# ---- 3. what does he actually talk about? -----------------------------------
print()
print("=" * 78)
print("3. WHAT HE TALKS ABOUT (term frequency across all 113 verbatims)")
print("=" * 78)
TERMS = {
 'edge/end definition': r'\bedge|\bend\b|\bends\b',
 'flank handling (degap/length)': r'flank|degap|gappy',
 'consensus wrong/too long': r'consensus',
 'copy number too few': r'too few|not enough|few copies|underrepresented',
 'mosaic': r'mosaic',
 'mixture of subfamilies': r'mixtur|subfamil|subgroup|separat',
 'length (too long/short)': r'too long|too short|longer|shorter',
 'microsatellite / tandem': r'microsat|tandem|satellite',
 'uniqueness / genome check': r'uniq|genome level|whole-genome',
 'alignment quality / presentation': r'align|presented|prepared|realign',
 'divergence / age': r'diverg|ancient|old\b|young',
 'nested / other repeats': r'nested|LINE|CR1|Ginger|repeat',
}
low = [(r['verbatim'] or '').lower() for r in R]
for k, pat in sorted(TERMS.items(), key=lambda z: -sum(1 for t in low if re.search(z[1], t))):
    n = sum(1 for t in low if re.search(pat, t))
    print("  %-36s %3d / 113  (%2.0f%%)" % (k, n, 100.0*n/113))

# ---- 4. verdicts that are complaints about INPUT, not biology ---------------
print()
print("=" * 78)
print("4. IS IT ABOUT THE ELEMENT, OR ABOUT HOW IT WAS SHOWN TO HIM?")
print("=" * 78)
INPUT_FAULT = ('ANCIENT_BADLY_PRESENTED','BADLY_PRESENTED','UNUSABLE_UNALIGNED',
               'UNUSABLE','NO_FLANKS','SINE_CONSENSUS_TOO_LONG','SINE_CONSENSUS_WRONG',
               'SINE_CONSENSUS_TOO_LONG_NOT_MIXTURE','GREY_OR_BADLY_PREPARED',
               'TOO_FEW_AND_BAD_LEFT_FLANK','MIXTURE_SPLIT_FIRST','SINE_NEEDS_SUBGROUPING',
               'MIXTURE_TWO_SUBFAMILIES','SINE_TOO_FEW')
c = Counter(r['call'] for r in R)
n_input = sum(v for k, v in c.items() if k in INPUT_FAULT)
print("  verdicts that are wholly or partly a complaint about the INPUT: %d of 113 (%.0f%%)"
      % (n_input, 100.0*n_input/113))
for k in sorted(INPUT_FAULT):
    if c.get(k):
        print("    %-38s %d" % (k, c[k]))
print()
print("  pure biology verdicts (SINE / NOT_SINE / GREY only): %d"
      % sum(v for k, v in c.items() if k in ('SINE','NOT_SINE','GREY')))
