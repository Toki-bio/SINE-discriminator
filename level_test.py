"""Can any statistic tell a SUBFAMILY from a CLADE above it or a FRAGMENT below it?

Test bed: teu_subfam_input.aln.fa, six blocks of 200 chunks whose membership is
known by construction. Build three classes of candidate block and ask whether any
of the four proposed level-criteria separates them.

  SUBFAMILY  each of the 6 true groups                       (the right level)
  CLADE      unions of 2 and 3 true groups                   (one level too coarse)
  FRAGMENT   contiguous halves and quarters of a true group  (one level too fine)

A usable level criterion must score SUBFAMILY differently from BOTH others.
"""
import io, random, statistics, itertools
from collections import Counter

S = 'C:/Users/T/AppData/Local/Temp/claude/c--work-Raks-COI/e2cef7cb-0c8d-4c4f-af8b-cb95bed931f5/scratchpad/'
GAPS = '-.'

def rd(p):
    n, s, c, b = [], [], None, []
    for l in io.open(p, encoding='utf-8'):
        l = l.rstrip()
        if l.startswith('>'):
            if c is not None: s.append(''.join(b))
            c = l[1:].split()[0]; n.append(c); b = []
        else: b.append(l.strip())
    s.append(''.join(b)); return n, s

n, s = rd(S + 'teu_subfam_input.aln.fa')
L = len(s[0])
BOUNDS = [(1,201,'t1'),(202,402,'t2'),(403,603,'t3'),(604,804,'t4'),(805,1005,'t5'),(1006,1206,'t6')]
groups = {}
for a, b, lab in BOUNDS:
    groups[lab] = [i for i in range(a-1, b) if 'CONSENSUS' not in n[i]]
ALL = sorted(sum(groups.values(), []))

def col_counts(idx):
    out = []
    for j in range(L):
        c = Counter(s[i][j].upper() for i in idx)
        out.append(c)
    return out

TOTAL = col_counts(ALL)

def diagnostics(mem, f_in=0.90, f_out=0.02):
    ms = set(mem); nout = len(ALL) - len(mem)
    got = []
    for j in range(L):
        c = Counter(s[i][j].upper() for i in mem)
        ch, cnt = c.most_common(1)[0]
        if cnt / float(len(mem)) < f_in: continue
        fo = (TOTAL[j].get(ch, 0) - cnt) / float(nout)
        if fo <= f_out: got.append((j, ch))
    return got

def pid(a, b):
    m = t = 0
    for x, y in zip(a, b):
        if x in GAPS and y in GAPS: continue
        t += 1
        if x.upper() == y.upper(): m += 1
    return m / float(t) if t else 0.0

def within_pairs(mem, k=300, seed=0):
    rnd = random.Random(seed)
    return [pid(s[rnd.choice(mem)], s[rnd.choice(mem)]) for _ in range(k)]

def bimodality(v):
    """Sarle's bimodality coefficient. >0.555 suggests bimodal/multimodal."""
    nn = len(v)
    m = statistics.mean(v); sd = statistics.pstdev(v)
    if sd == 0: return 0.0
    g1 = sum((x-m)**3 for x in v)/nn/sd**3
    g2 = sum((x-m)**4 for x in v)/nn/sd**4 - 3.0
    return (g1*g1 + 1.0) / (g2 + 3.0*((nn-1)**2)/float((nn-2)*(nn-3)))

def stability(mem, reps=8, drop=0.10, seed=1):
    """Fraction of the full block's diagnostics that survive dropping 10% of rows."""
    base = set(diagnostics(mem))
    if not base: return 0.0
    rnd = random.Random(seed); keep_frac = []
    for _ in range(reps):
        sub = rnd.sample(mem, int(round(len(mem)*(1-drop))))
        d = set(diagnostics(sub))
        keep_frac.append(len(base & d) / float(len(base)))
    return statistics.mean(keep_frac)

def consensus(mem, thr, cov=0.30):
    out = []
    N = len(mem)
    for j in range(L):
        col = [s[i][j].upper() for i in mem]
        nz = [x for x in col if x not in GAPS]
        if not nz or len(nz)/float(N) < cov: out.append('-'); continue
        c = Counter(nz); ch, cnt = c.most_common(1)[0]
        out.append(ch if cnt/float(N) >= thr else '-')
    return ''.join(out)

def thr_sensitivity(mem):
    """How much does the consensus move between threshold 50% and 30%?
    His signal: a clean block gives the same consensus across that range."""
    a = consensus(mem, 0.50); b = consensus(mem, 0.30)
    diff = sum(1 for x, y in zip(a, b) if x != y)
    span = sum(1 for x in a if x not in GAPS) or 1
    return diff / float(span)

# ---- build the three classes ------------------------------------------------
cands = []
for lab, mem in groups.items():
    cands.append(('SUBFAMILY', lab, mem))
for a, b in itertools.combinations(sorted(groups), 2):
    cands.append(('CLADE2', a+'+'+b, groups[a]+groups[b]))
for combo in itertools.combinations(sorted(groups), 3):
    cands.append(('CLADE3', '+'.join(combo), sum((groups[c] for c in combo), [])))
rnd = random.Random(7)
for lab, mem in groups.items():
    cands.append(('FRAG_half', lab+'/2a', mem[:len(mem)//2]))
    cands.append(('FRAG_half', lab+'/2b', mem[len(mem)//2:]))
    q = len(mem)//4
    cands.append(('FRAG_quarter', lab+'/4a', mem[:q]))
    cands.append(('FRAG_quarter', lab+'/4b', mem[q:2*q]))

print("%-13s %-12s %5s %6s %8s %8s %8s %8s %8s" %
      ("class","block","n","ndiag","dens/100","cohgap","bimod","stabil","thrsens"))
rows = []
for cls, name, mem in cands:
    d = diagnostics(mem)
    w = within_pairs(mem)
    others = [i for i in ALL if i not in set(mem)]
    rnd2 = random.Random(3)
    bet = [pid(s[rnd2.choice(mem)], s[rnd2.choice(others)]) for _ in range(300)]
    gap = statistics.mean(w) - statistics.mean(bet)
    dens = 100.0*len(d)/len(mem)
    bim = bimodality(w)
    stab = stability(mem)
    ts = thr_sensitivity(mem)
    rows.append((cls, name, len(mem), len(d), dens, gap, bim, stab, ts))
    print("%-13s %-12s %5d %6d %8.2f %8.4f %8.3f %8.3f %8.4f" %
          (cls, name, len(mem), len(d), dens, gap, bim, stab, ts))

print()
print("=" * 92)
print("SEPARATION: does any statistic distinguish SUBFAMILY from CLADE and FRAGMENT?")
print("=" * 92)
def grp(c):
    return [r for r in rows if r[0] == c]
stats = [("ndiag",3),("density/100",4),("cohesion gap",5),("bimodality",6),
         ("stability",7),("thr-sensitivity",8)]
classes = ["SUBFAMILY","CLADE2","CLADE3","FRAG_half","FRAG_quarter"]
for label, idx in stats:
    print("\n  %s" % label)
    for c in classes:
        v = [r[idx] for r in grp(c)]
        if not v: continue
        print("    %-13s n=%-3d min %8.3f  median %8.3f  max %8.3f" %
              (c, len(v), min(v), statistics.median(v), max(v)))
    sub = [r[idx] for r in grp("SUBFAMILY")]
    oth = [r[idx] for r in rows if r[0] != "SUBFAMILY"]
    lo, hi = min(sub), max(sub)
    inside = sum(1 for x in oth if lo <= x <= hi)
    print("    -> non-subfamily blocks falling inside the subfamily range: %d of %d (%.0f%%)"
          % (inside, len(oth), 100.0*inside/len(oth)))
