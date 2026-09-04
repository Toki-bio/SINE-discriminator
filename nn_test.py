"""Derived from the failure of the four candidates.

Every statistic computed ON a block fails to tell a subfamily from a fragment of
one, because a fragment of a homogeneous group is itself homogeneous. The
information has to be in the block's relationship to its NEAREST outside neighbour,
not its mean.

  fragment   -> each member's nearest neighbour is usually its own sibling, which
                lies OUTSIDE the block. Containment low.
  subfamily  -> nearest neighbour is a fellow member. Containment high.
  clade      -> containment also high, but threshold-sensitivity already catches
                clades.

So the hypothesis is that the PAIR (nearest-neighbour containment, threshold
sensitivity) covers both directions of error.
"""
import io, random, statistics, itertools
import numpy as np
from collections import Counter

S = 'C:/Users/T/AppData/Local/Temp/claude/c--work-Raks-COI/e2cef7cb-0c8d-4c4f-af8b-cb95bed931f5/scratchpad/'

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
pos = {g: k for k, g in enumerate(ALL)}

# encode: 0-3 = ACGT, 4 = gap
M = np.full((len(ALL), L), 4, dtype=np.int8)
for k, i in enumerate(ALL):
    a = np.frombuffer(s[i].upper().encode(), dtype=np.uint8)
    for q, ch in enumerate("ACGT"):
        M[k][a == ord(ch)] = q
print("matrix", M.shape)

N = M.shape[0]
ID = np.zeros((N, N), dtype=np.float32)
for k in range(N):
    same = (M[k] == M).sum(axis=1)
    bothgap = ((M[k] == 4) & (M == 4)).sum(axis=1)
    ID[k] = (same - bothgap) / np.maximum(L - bothgap, 1)
np.fill_diagonal(ID, -1.0)
NN = ID.argmax(axis=1)                      # index of each row's nearest neighbour
print("nearest-neighbour map built")

truth = {}
for lab, mem in groups.items():
    for i in mem: truth[pos[i]] = lab
nn_same = sum(1 for k in range(N) if truth[k] == truth[NN[k]])
print("sanity: nearest neighbour is in the SAME true group for %d of %d rows (%.1f%%)"
      % (nn_same, N, 100.0*nn_same/N))

GAPS = '-.'
def consensus(memk, thr, cov=0.30):
    sub = M[memk]
    out = []
    for j in range(L):
        col = sub[:, j]
        nz = col[col != 4]
        if nz.size == 0 or nz.size / float(len(memk)) < cov: out.append('-'); continue
        vals, cnts = np.unique(nz, return_counts=True)
        mx = cnts.max()
        out.append("ACGT"[vals[cnts.argmax()]] if mx / float(len(memk)) >= thr else '-')
    return ''.join(out)

def thr_sens(memk):
    a = consensus(memk, 0.50); b = consensus(memk, 0.30)
    diff = sum(1 for x, y in zip(a, b) if x != y)
    span = sum(1 for x in a if x not in GAPS) or 1
    return diff / float(span)

def nn_containment(memk):
    ms = set(memk)
    return sum(1 for k in memk if NN[k] in ms) / float(len(memk))

cands = []
for lab, mem in groups.items():
    cands.append(('SUBFAMILY', lab, [pos[i] for i in mem]))
for a, b in itertools.combinations(sorted(groups), 2):
    cands.append(('CLADE2', a+'+'+b, [pos[i] for i in groups[a]+groups[b]]))
for combo in itertools.combinations(sorted(groups), 3):
    cands.append(('CLADE3', '+'.join(combo), [pos[i] for i in sum((groups[c] for c in combo), [])]))
for lab, mem in groups.items():
    mk = [pos[i] for i in mem]
    cands.append(('FRAG_half', lab+'/2a', mk[:len(mk)//2]))
    cands.append(('FRAG_half', lab+'/2b', mk[len(mk)//2:]))
    q = len(mk)//4
    cands.append(('FRAG_quarter', lab+'/4a', mk[:q]))
    cands.append(('FRAG_quarter', lab+'/4b', mk[q:2*q]))

rows = []
print()
print("%-13s %-12s %5s %10s %10s" % ("class","block","n","NN-contain","thr-sens"))
for cls, name, memk in cands:
    c = nn_containment(memk); t = thr_sens(memk)
    rows.append((cls, name, len(memk), c, t))
    print("%-13s %-12s %5d %10.3f %10.4f" % (cls, name, len(memk), c, t))

print()
print("=" * 70)
print("SEPARATION by class")
print("=" * 70)
for label, idx in (("NN-containment", 3), ("thr-sensitivity", 4)):
    print("\n  %s" % label)
    for c in ["SUBFAMILY","CLADE2","CLADE3","FRAG_half","FRAG_quarter"]:
        v = [r[idx] for r in rows if r[0] == c]
        if v: print("    %-13s n=%-3d min %7.3f median %7.3f max %7.3f"
                    % (c, len(v), min(v), statistics.median(v), max(v)))

print()
print("=" * 70)
print("COMBINED RULE: subfamily iff NN-containment >= A and thr-sensitivity <= B")
print("=" * 70)
best = None
for A in [x/100.0 for x in range(50, 100)]:
    for B in [x/1000.0 for x in range(5, 120, 1)]:
        tp = sum(1 for r in rows if r[0]=='SUBFAMILY' and r[3] >= A and r[4] <= B)
        fp = sum(1 for r in rows if r[0]!='SUBFAMILY' and r[3] >= A and r[4] <= B)
        fn = 6 - tp
        acc = (tp + (59 - fp)) / 65.0
        if best is None or (tp, -fp) > (best[0], -best[1]):
            best = (tp, fp, A, B, acc)
print("  best: %d/6 subfamilies accepted, %d/59 non-subfamilies wrongly accepted"
      % (best[0], best[1]))
print("  at NN-containment >= %.2f and thr-sensitivity <= %.3f  (accuracy %.3f)"
      % (best[2], best[3], best[4]))
print("  wrongly accepted:", [r[1] for r in rows if r[0]!='SUBFAMILY' and r[3] >= best[2] and r[4] <= best[3]])
print("  missed:", [r[1] for r in rows if r[0]=='SUBFAMILY' and not (r[3] >= best[2] and r[4] <= best[3])])
