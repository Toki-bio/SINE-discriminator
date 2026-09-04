"""The sister-exclusion fix, in a form computable WITHOUT labels.

Excluding "the nearest sister" needs the answer. But the same effect is available
blind: measure f_out against the nearest OUTSIDE ROWS only, since a block's nearest
competitors are its sister group whether or not we can name it.

  f_out_local = fraction of the nearest K% of non-members carrying the state

A subfamily should still have diagnostics under f_out_local, because its private
derived states are absent even from its nearest competitor.
A fragment should NOT, because its nearest competitors are its own siblings, which
carry every state it has.

This is the crucial asymmetry the global f_out throws away.
"""
import io, itertools, statistics
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
groups = {lab: [i for i in range(a-1, b) if 'CONSENSUS' not in n[i]] for a, b, lab in BOUNDS}
labs = sorted(groups)
ALL = sorted(sum(groups.values(), []))
pos = {g: k for k, g in enumerate(ALL)}
N = len(ALL)

M = np.full((N, L), 4, dtype=np.int8)
for k, i in enumerate(ALL):
    a = np.frombuffer(s[i].upper().encode(), dtype=np.uint8)
    for q, ch in enumerate("ACGT"):
        M[k][a == ord(ch)] = q
ID = np.zeros((N, N), dtype=np.float32)
for k in range(N):
    same = (M[k] == M).sum(axis=1)
    bg = ((M[k] == 4) & (M == 4)).sum(axis=1)
    ID[k] = (same - bg) / np.maximum(L - bg, 1)
np.fill_diagonal(ID, -1.0)

def diag_local(memk, k_frac=0.20, f_in=0.90, f_out=0.02):
    ms = np.zeros(N, dtype=bool); ms[list(memk)] = True
    out = np.where(~ms)[0]
    if out.size == 0: return []
    # nearest outside rows: rank non-members by max identity to any member
    prox = ID[np.ix_(out, list(memk))].max(axis=1)
    k = max(1, int(round(k_frac * out.size)))
    near = out[np.argsort(-prox)[:k]]
    sub = M[list(memk)]; nb = M[near]
    got = []
    for j in range(L):
        col = sub[:, j]
        vals, cnts = np.unique(col, return_counts=True)
        mx = cnts.max(); ch = vals[cnts.argmax()]
        if mx / float(len(memk)) < f_in: continue
        fo = (nb[:, j] == ch).sum() / float(k)
        if fo <= f_out: got.append((j, int(ch)))
    return got

def diag_global(memk, f_in=0.90, f_out=0.02):
    ms = np.zeros(N, dtype=bool); ms[list(memk)] = True
    out = np.where(~ms)[0]
    sub = M[list(memk)]; nb = M[out]
    got = []
    for j in range(L):
        col = sub[:, j]
        vals, cnts = np.unique(col, return_counts=True)
        mx = cnts.max(); ch = vals[cnts.argmax()]
        if mx / float(len(memk)) < f_in: continue
        fo = (nb[:, j] == ch).sum() / float(out.size)
        if fo <= f_out: got.append((j, int(ch)))
    return got

cands = []
for lab in labs:
    cands.append(('SUBFAMILY', lab, [pos[i] for i in groups[lab]]))
for a, b in itertools.combinations(labs, 2):
    cands.append(('CLADE2', a+'+'+b, [pos[i] for i in groups[a]+groups[b]]))
for combo in itertools.combinations(labs, 3):
    cands.append(('CLADE3', '+'.join(combo), [pos[i] for i in sum((groups[c] for c in combo), [])]))
for lab in labs:
    mk = [pos[i] for i in groups[lab]]
    cands.append(('FRAG_half', lab+'/2a', mk[:len(mk)//2]))
    cands.append(('FRAG_half', lab+'/2b', mk[len(mk)//2:]))
    q = len(mk)//4
    cands.append(('FRAG_quarter', lab+'/4a', mk[:q]))
    cands.append(('FRAG_quarter', lab+'/4b', mk[q:2*q]))

rows = []
print("%-13s %-12s %5s %8s %10s" % ("class","block","n","global","local20%"))
for cls, name, memk in cands:
    g = len(diag_global(memk)); l = len(diag_local(memk))
    rows.append((cls, name, len(memk), g, l))
    print("%-13s %-12s %5d %8d %10d" % (cls, name, len(memk), g, l))

print()
print("=" * 72)
for label, idx in (("diagnostics vs ALL (current criterion)", 3),
                   ("diagnostics vs NEAREST 20% (proposed)", 4)):
    print("\n  %s" % label)
    for c in ["SUBFAMILY","CLADE2","CLADE3","FRAG_half","FRAG_quarter"]:
        v = [r[idx] for r in rows if r[0] == c]
        if v: print("    %-13s n=%-3d min %4d median %6.1f max %4d  |  with >=1: %d/%d"
                    % (c, len(v), min(v), statistics.median(v), max(v),
                       sum(1 for x in v if x >= 1), len(v)))
    sub = [r[idx] for r in rows if r[0] == 'SUBFAMILY']
    oth = [r for r in rows if r[0] != 'SUBFAMILY']
    for T in (1, 3, 5, 8):
        tp = sum(1 for x in sub if x >= T)
        fp = sum(1 for r in oth if r[idx] >= T)
        print("      threshold >=%d diagnostics: %d/6 subfamilies, %d/%d others wrongly accepted"
              % (T, tp, fp, len(oth)))
