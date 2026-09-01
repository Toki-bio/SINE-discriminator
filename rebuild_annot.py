#!/usr/bin/env python3
"""Regenerate annotations for every set the viewer actually shows.

annotate.main() only looked in aln_c/, so the Timema and human picks - which
live in bench/ and tim_bench/ - got no feature annotations at all and would
render with an empty track. This walks the viewer payload instead and searches
every corpus directory.
"""
import io
import json
import os
import sys

sys.path.insert(0, ".")
import annotate as A

DIRS = ["aln_c", "aln_v2", "bench", "tim_bench"]


def find(setname):
    for d in DIRS:
        p = os.path.join(d, setname + ".aln.fa")
        if os.path.exists(p):
            return p
    return None


def main():
    viewer = json.load(open("viewer_data.json"))
    prof = {}
    if os.path.exists("profiles.json"):
        prof = json.load(open("profiles.json"))
    else:
        # profiles live inside _embed.js now; pull them back out
        for line in io.open("_embed.js", encoding="utf-8").read().split("\n"):
            if line.startswith("const PROFILES="):
                prof = json.loads(line[len("const PROFILES="):].rstrip().rstrip(";"))
                break

    out = {}
    missing, failed = [], []
    for d in viewer:
        s = d["set"]
        p = find(s)
        if p is None:
            missing.append(s)
            continue
        try:
            a = A.annotate(p, prof.get(s))
        except Exception:
            failed.append(s)
            continue
        if a:
            out[s] = a
    json.dump(out, open("annotations.json", "w"), separators=(",", ":"))
    print("annotated %d of %d sets" % (len(out), len(viewer)))
    if missing:
        print("no alignment found: %d %s" % (len(missing), missing[:4]))
    if failed:
        print("annotate failed:   %d %s" % (len(failed), failed[:4]))
    th = [s for s in out if s.startswith(("TIM", "HUM"))]
    print("timema/human annotated: %d" % len(th))


if __name__ == "__main__":
    main()
