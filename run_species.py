#!/usr/bin/env python3
"""One command per genome: AnnoSINE candidates in, judged alignments out.

Everything that was a sequence of hand-typed steps, in the order that matters:

  1  search the genome with each candidate consensus, extract with 400 bp
     flanks, align                                        (candidate_to_aln)
  2  measure on those LONG flanks - decay, islands, microsatellite
  3  justify and trim the flanks for viewing              (justify + trim)
  4  score
  5  where the consensus is too short, EXTEND it, re-search, and where the
     extension overshoots, trim to the window the copies support - then score
     again. hyd_SINE_0 goes 0.0 NO_ELEMENT -> 100.0 clean through this loop,
     ending at 211 bp against RepBase's 208 bp for SINE2-2B_HM.

Steps 2 and 3 use different geometry on purpose: measurements need the full
400 bp, the viewer needs the flank cut to what the copies fill. Taking both
from one file is what hid the flank islands.
"""
import glob
import json
import os
import subprocess
import sys

BASE = "/staging/tmp/sinedisc"
NEWSP = "/staging/tmp/newsp"
sys.path.insert(0, BASE)


def sh(cmd, quiet=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=BASE)
    if not quiet:
        print(r.stdout[-2000:] or r.stderr[-800:])
    return r


def side_tables(aln_dir, prefix):
    """decay / islands / microsatellite, measured on the raw long flanks."""
    sh("python3 flankdecay.py %s > /dev/null 2>&1 && mv flankdecay.json %s_decay.json"
       % (aln_dir, prefix))
    sh("python3 islands_corpus.py %s %s_islands.json > /dev/null 2>&1" % (aln_dir, prefix))
    sh("python3 microsat.py %s %s_msat.json > /dev/null 2>&1" % (aln_dir, prefix))
    for src, dst in (("%s_decay.json" % prefix, "flankdecay.json"),
                     ("%s_islands.json" % prefix, "flank_islands.json"),
                     ("%s_msat.json" % prefix, "microsat.json")):
        p = os.path.join(BASE, src)
        if not os.path.exists(p):
            json.dump({}, open(os.path.join(BASE, dst), "w"))
            continue
        d = json.load(open(p))
        out = {}
        for k, v in d.items():
            out[k] = v
            out["NEW__%s.clean" % k] = v
        json.dump(out, open(os.path.join(BASE, dst), "w"))


def clean(aln_dir, clean_dir):
    import justify_all as J
    import trim_flanks as T
    os.makedirs(clean_dir, exist_ok=True)
    for f in glob.glob(os.path.join(clean_dir, "*")):
        os.remove(f)
    n = 0
    for f in sorted(glob.glob(os.path.join(aln_dir, "*.aln.fa"))):
        b = os.path.basename(f).replace(".aln.fa", "")
        tmp = os.path.join(clean_dir, "_j_" + b + ".aln.fa")
        out = os.path.join(clean_dir, "NEW__" + b + ".clean.aln.fa")
        try:
            if not J.justify(f, tmp):
                continue
            if not T.trim(tmp, out):
                os.replace(tmp, out)
            elif os.path.exists(tmp):
                os.remove(tmp)
            n += 1
        except Exception as exc:
            print("   clean failed on %s: %s" % (b, exc))
    return n


def score(clean_dir, out_json):
    sh("python3 verdict.py %s %s > /dev/null 2>&1" % (clean_dir, out_json))
    p = os.path.join(BASE, out_json)
    return json.load(open(p)) if os.path.exists(p) else {}


def refine(sp, verdicts, genome, aln_dir):
    """Repair the consensuses the scorer says are the wrong length.

    Under-extended: extend by the decay distance and re-search.
    Then over-extended: cut to the window the copies actually support.
    Both readings come from the scorer itself, so nothing here is a guess.
    """
    todo = []
    for k, d in verdicts.items():
        if "flags" not in d:
            continue
        if "top100" not in k:
            continue
        tag = k.replace("NEW__", "").replace(".clean", "").replace("__top100", "")
        codes = [f["code"] for f in d["flags"]]
        if "CONSENSUS_UNDEREXTENDED" in codes:
            todo.append((tag, "short", None))
        elif "CONSENSUS_OVEREXTENDED" in codes and d["score"] < 90:
            cw = [f["n"] for f in d["flags"] if f["code"] == "CONSENSUS_OVEREXTENDED"][0]
            todo.append((tag, "long", cw))
    if not todo:
        return {}, []
    print("   refining %d consensus%s: %s"
          % (len(todo), "" if len(todo) == 1 else "es",
             ", ".join("%s(%s)" % (t, w) for t, w, _ in todo)))

    ext_dir = os.path.join(NEWSP, sp, "ext")
    os.makedirs(ext_dir, exist_ok=True)
    made = []
    for tag, why, cw in todo:
        out = os.path.join(ext_dir, "%s.new.fa" % tag)
        if why == "short":
            r = sh("python3 extend_consensus.py %s %s_decay.json %s %s"
                   % (aln_dir, sp, ext_dir, tag))
            src = os.path.join(ext_dir, "%s.ext.fa" % tag)
            if os.path.exists(src):
                os.replace(src, out)
                made.append(out)
        else:
            # cut the existing consensus to the supported window
            seq, keep = [], False
            p = os.path.join(aln_dir, "%s__top100.aln.fa" % tag)
            for line in open(p):
                if line.startswith(">"):
                    if keep:
                        break
                    keep = "CONSENSUS_" in line.upper()
                    continue
                if keep:
                    seq.append(line.strip())
            s = "".join(seq).replace("-", "")
            lo, hi = cw
            if hi - lo >= 60:
                open(out, "w").write(">%s_ref\n%s\n" % (tag, s[lo:hi + 1]))
                made.append(out)
    if not made:
        return {}, []

    ref_aln = os.path.join(NEWSP, sp, "refaln")
    os.makedirs(ref_aln, exist_ok=True)
    for f in glob.glob(os.path.join(ref_aln, "*")):
        os.remove(f)
    for m in made:
        sh("python3 candidate_to_aln.py %s %s %s %sR" % (m, genome, ref_aln, sp))
    if not glob.glob(os.path.join(ref_aln, "*.aln.fa")):
        return {}, []
    side_tables(ref_aln, sp + "_ref")
    ref_clean = os.path.join(NEWSP, sp, "refclean")
    clean(ref_aln, ref_clean)
    return score(ref_clean, "%s_ref_verdicts.json" % sp), [t for t, _, _ in todo]


def main():
    sp = sys.argv[1]
    genome = os.path.join(NEWSP, sp, "genome.fna")
    seeds = os.path.join(NEWSP, sp, "anno", "Seed_SINE.fa")
    aln = os.path.join(NEWSP, sp, "aln")
    cln = os.path.join(NEWSP, sp, "clean")

    if not glob.glob(os.path.join(aln, "*.aln.fa")):
        os.makedirs(aln, exist_ok=True)
        print("1. searching the genome with each candidate")
        sh("python3 candidate_to_aln.py %s %s %s %s" % (seeds, genome, aln, sp), quiet=False)
    n_aln = len(glob.glob(os.path.join(aln, "*.aln.fa")))
    print("   %d alignments" % n_aln)

    print("2. decay / islands / microsatellite on the 400 bp flanks")
    side_tables(aln, sp)
    print("3. justify and trim for viewing")
    print("   %d cleaned" % clean(aln, cln))
    print("4. scoring")
    v = score(cln, "%s_verdicts.json" % sp)
    print("   %d scored, %d at 50 or above"
          % (len(v), sum(1 for d in v.values() if d.get("score", 0) >= 50 and "top100" in d.get("set", ""))))
    print("5. repairing consensus lengths")
    rv, tags = refine(sp, v, genome, aln)
    # Extending stops where similarity reaches background, which overshoots the
    # element edge; the scorer then names the window the copies support. Cutting
    # to it is what finishes the repair - without this second pass hyd_SINE_0
    # stalls at its 310 bp overshoot instead of reaching 211 bp and 100.0.
    if rv:
        sh("python3 refine_pass2.py %s" % sp, quiet=False)
    if rv:
        for k, d in sorted(rv.items()):
            if "top100" not in k:
                continue
            print("   %-44s %5.1f  %s" % (k.replace("NEW__", "")[:44], d["score"],
                                          ",".join(f["code"] for f in d["flags"]) or "clean"))
    json.dump({"verdicts": v, "refined": rv, "refined_tags": tags},
              open(os.path.join(BASE, "%s_all.json" % sp), "w"))
    print("done: %s_all.json" % sp)


if __name__ == "__main__":
    main()
