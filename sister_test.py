"""Why do t1 and t2 have ZERO diagnostics each, while t1+t2 has 16?

Hypothesis: the criterion measures f_out against EVERYTHING outside the block. A
subfamily's private derived states are shared with its nearest sister, so the
sister's presence in the "outside" set destroys them. Two sisters therefore each
score zero while their union scores high -- which looks like evidence that the
union is the real unit, but is actually an artefact of the denominator.

If true, computing diagnostics against everything EXCEPT the nearest sister should
resurrect them.
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

def diagnostics(mem, out, f_in=0.90, f_out=0.02):
    got = []
    for j in range(L):
        c = Counter(s[i][j].upper() for i in mem)
        ch, cnt = c.most_common(1)[0]
        if cnt / float(len(mem)) < f_in: continue
        co = Counter(s[i][j].upper() for i in out)
        if co.get(ch, 0) / float(len(out)) <= f_out: got.append((j, ch))
    return got

def pid(a, b):
    m = t = 0
    for x, y in zip(a, b):
        if x in '-.' and y in '-.': continue
        t += 1
        if x.upper() == y.upper(): m += 1
    return m / float(t) if t else 0.0

# consensus of each group, to find nearest sister
import random
def cons(mem, thr=0.35):
    out = []
    for j in range(L):
        col = [s[i][j].upper() for i in mem]
        nz = [x for x in col if x not in '-.']
        if not nz: out.append('-'); continue
        c = Counter(nz); ch, cnt = c.most_common(1)[0]
        out.append(ch if cnt/float(len(mem)) >= thr else '-')
    return ''.join(out)

C = {lab: cons(groups[lab]) for lab in labs}
print("consensus identity between true groups:")
print("      " + "  ".join("%-6s" % l for l in labs))
sister = {}
for a in labs:
    row = []
    best = (-1, None)
    for b in labs:
        v = pid(C[a], C[b]) if a != b else 1.0
        row.append(v)
        if a != b and v > best[0]: best = (v, b)
    sister[a] = best[1]
    print("%-5s " % a + "  ".join("%-6.4f" % v for v in row) + "   nearest sister: %s (%.4f)" % (best[1], best[0]))

print()
print("%-6s %-8s %14s %26s" % ("group", "sister", "diag vs ALL", "diag vs ALL-minus-sister"))
for a in labs:
    out_all = [i for i in ALL if i not in set(groups[a])]
    out_ns = [i for i in ALL if i not in set(groups[a]) and i not in set(groups[sister[a]])]
    d1 = diagnostics(groups[a], out_all)
    d2 = diagnostics(groups[a], out_ns)
    print("%-6s %-8s %14d %26d" % (a, sister[a], len(d1), len(d2)))

print()
print("and the unions, for contrast:")
for a, b in itertools.combinations(labs, 2):
    mem = groups[a] + groups[b]
    out = [i for i in ALL if i not in set(mem)]
    d = diagnostics(mem, out)
    if d: print("   %-8s %d diagnostics vs everything else" % (a+'+'+b, len(d)))
