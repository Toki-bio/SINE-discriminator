#!/usr/bin/env python3
"""Collapse a discovery tool's candidate list before anything is scored.

AnnoSINE proposed 160 candidates for zebrafish; all-vs-all blastn collapses them
to 87 families, with one cluster holding 55 near-identical seeds. Sergei spotted
it immediately: "maybe zebrafish 160 candidates can be deduplicated? too much
redundancy suspected."

It costs more than compute. Every per-genome count, and any threshold ever
fitted on a candidate list, is skewed by whichever family happened to be
proposed fifty-five times.

Two candidates are the same family when they align at >= 80 % identity over
>= 60 % of the shorter one. Groups are single-link, and the representative is the
one with the most genomic copies (AnnoSINE's own blast_count), longest first on
a tie - not whichever came first in the file.
"""
import os
import re
import subprocess
import sys
from collections import defaultdict

MIN_ID = 80.0
MIN_COV = 0.60


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    for line in open(path):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    tmp = dst + ".db"

    names, seqs = read_fa(src)
    ids = [n.split()[0] for n in names]
    copies, length = {}, {}
    for n, s in zip(names, seqs):
        cid = n.split()[0]
        m = re.search(r"blast_count:(\d+)", n)
        copies[cid] = int(m.group(1)) if m else 0
        length[cid] = len(s)

    if len(ids) < 2:
        open(dst, "w").write(open(src).read())
        print("%d candidate, nothing to collapse" % len(ids))
        return

    subprocess.run("makeblastdb -in %s -dbtype nucl -out %s" % (src, tmp),
                   shell=True, capture_output=True)
    r = subprocess.run(
        "blastn -query %s -db %s -evalue 1e-5 -word_size 11 -dust no "
        "-outfmt '6 qseqid sseqid pident length qlen slen' -num_threads 16"
        % (src, tmp), shell=True, capture_output=True, text=True)

    lab = {i: i for i in ids}

    def find(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x

    merges = 0
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6 or f[0] == f[1]:
            continue
        q, s = f[0], f[1]
        if q not in lab or s not in lab:
            continue
        pid, aln, ql, sl = float(f[2]), int(f[3]), int(f[4]), int(f[5])
        if pid >= MIN_ID and aln >= MIN_COV * min(ql, sl):
            a, b = find(q), find(s)
            if a != b:
                lab[a] = b
                merges += 1

    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)
    reps = set()
    for g in groups.values():
        reps.add(sorted(g, key=lambda k: (-copies[k], -length[k]))[0])

    keep = False
    with open(dst, "w") as out:
        for line in open(src):
            if line.startswith(">"):
                keep = line[1:].split()[0] in reps
            if keep:
                out.write(line)

    sizes = sorted((len(g) for g in groups.values()), reverse=True)
    print("%d candidates -> %d families (%d merges); largest clusters %s"
          % (len(ids), len(reps), merges, sizes[:5]))
    for ext in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto", ".njs"):
        try:
            os.remove(tmp + ext)
        except OSError:
            pass


if __name__ == "__main__":
    main()
