"""Island fraction against his own judgements - the only ground truth that counts.

Question: is a high island fraction a REJECT signal or a CAUTION signal? If it
lines up with the calls where he asked for extra checking rather than the ones
he rejected outright, it belongs in the report as a caution, not in the score
as a penalty.
"""
import json
from collections import defaultdict
import numpy as np

isl = json.load(open("islands_corpus.json"))
ver = json.load(open("verdicts.json"))

rows = []
for line in open("calls.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 4 or f[0] == "set":
        continue
    s = f[0]
    if s not in isl:
        continue
    v = ver.get(s, {})
    rows.append((s, f[2], isl[s]["frac"], isl[s]["cols"],
                 v.get("flank_bg"), v.get("score"), f[3]))

print("judged sets with an island measurement: %d" % len(rows))
by = defaultdict(list)
for r in rows:
    by[r[1]].append(r)

print()
print("%-26s %4s %9s %9s %8s" % ("his call", "n", "med frac", "max frac", "med bg"))
print("-" * 62)
for c in sorted(by, key=lambda c: -np.median([x[2] for x in by[c]])):
    v = by[c]
    b = [x[4] for x in v if x[4] is not None]
    print("%-26s %4d %9.3f %9.3f %8s"
          % (c, len(v), np.median([x[2] for x in v]), max(x[2] for x in v),
             "-" if not b else "%.3f" % np.median(b)))

print()
print("every judged set with island fraction >= 0.10, in order:")
print("%-30s %-26s %6s %7s %6s" % ("set", "his call", "frac", "bg", "score"))
for r in sorted([x for x in rows if x[2] >= 0.10], key=lambda x: -x[2]):
    print("%-30s %-26s %6.3f %7s %6s"
          % (r[0][:30], r[1], r[2],
             "-" if r[4] is None else "%.3f" % r[4],
             "-" if r[5] is None else "%.1f" % r[5]))
    print("      %s" % r[6][:96])
