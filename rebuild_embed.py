#!/usr/bin/env python3
"""Rebuild _embed.js from the current data files.

_embed.js is six single-line `const NAME=...;` declarations. It was previously
assembled by hand, which is how it drifted: the page was serving verdicts and
captions from before several scorer fixes, and the alignment/profile arrays
carried only the 24 talpid picks with no Timema or human sets at all.

This script makes that reproducible. Print only counts - the arrays contain
whole alignments and dumping one floods the terminal.
"""
import io
import json
import os
import sys

sys.path.insert(0, ".")

OUT = "_embed.js"
NAMES = ["ALIGNMENTS", "PROFILES", "ANNOTATIONS", "COMMENTARY", "VERDICTS", "GLOSSARY"]


def load_json(path, default):
    if not os.path.exists(path):
        return default
    try:
        return json.load(open(path))
    except Exception:
        return default


def existing():
    """Current values, so anything we do not regenerate is preserved."""
    cur = {}
    if not os.path.exists(OUT):
        return cur
    for line in io.open(OUT, encoding="utf-8").read().split("\n"):
        for n in NAMES:
            p = "const %s=" % n
            if line.startswith(p):
                try:
                    cur[n] = json.loads(line[len(p):].rstrip().rstrip(";"))
                except Exception:
                    pass
    return cur


def main():
    cur = existing()
    viewer = load_json("viewer_data.json", [])
    verdicts = load_json("verdicts.json", {})

    # ALIGNMENTS is the viewer payload itself: one entry per picked set.
    alignments = viewer

    # PROFILES is keyed by set name. Recompute for every picked set so the
    # Timema and human picks get graphs rather than silently rendering blank.
    profiles = dict(cur.get("PROFILES", {}))
    try:
        import profiles as P
        maker = getattr(P, "profile", None) or getattr(P, "build", None) \
            or getattr(P, "compute", None) or getattr(P, "make", None)
    except Exception:
        maker = None

    made, kept, failed = 0, 0, 0
    for d in alignments:
        st = d["set"]
        src = None
        for cand in ("aln_v2", "bench", "tim_bench"):
            f = os.path.join(cand, st + ".aln.fa")
            if os.path.exists(f):
                src = f
                break
        if src is None:
            failed += 1
            continue
        if maker is None:
            if st in profiles:
                kept += 1
            else:
                failed += 1
            continue
        try:
            profiles[st] = maker(src)
            made += 1
        except Exception:
            if st in profiles:
                kept += 1
            else:
                failed += 1

    out = {
        "ALIGNMENTS": alignments,
        "PROFILES": profiles,
        "ANNOTATIONS": load_json("annotations.json", cur.get("ANNOTATIONS", {})),
        "COMMENTARY": load_json("commentary.json", cur.get("COMMENTARY", {})),
        "VERDICTS": verdicts,
        "GLOSSARY": cur.get("GLOSSARY", []),
    }

    with io.open(OUT, "w", encoding="utf-8") as fh:
        for n in NAMES:
            fh.write("const %s=%s;\n" % (n, json.dumps(out[n], separators=(",", ":"))))

    print("alignments      %d" % len(alignments))
    print("profiles        %d  (recomputed %d, kept %d, failed %d)"
          % (len(profiles), made, kept, failed))
    print("verdicts        %d" % len(verdicts))
    print("annotations     %d" % len(out["ANNOTATIONS"]))
    print("commentary      %d" % len(out["COMMENTARY"]))
    print("glossary        %d" % len(out["GLOSSARY"]))
    tim = [d["set"] for d in alignments if d["set"].startswith(("TIM", "HUM"))]
    print("timema/human picks: %d" % len(tim))
    print("with profiles:      %d" % sum(1 for s in tim if s in profiles))


if __name__ == "__main__":
    main()
