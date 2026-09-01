"""Second pass of the consensus repair.

Extending an under-length consensus overshoots - it stops where similarity
reaches background, which is past the element edge. The scorer then says
CONSENSUS_OVEREXTENDED and names the window the copies actually support. Cutting
to that window is the step that finishes the job: hyd_SINE_0 goes 110 bp (0.0)
-> 310 bp (over-extended) -> 211 bp scoring 100.0, against RepBase 208 bp for
SINE2-2B_HM.

run_species.py did only the first pass, so every repaired consensus was left
sitting at the overshoot.
"""
import glob, json, os, subprocess, sys

BASE = "/staging/tmp/sinedisc"
NEWSP = "/staging/tmp/newsp"
sys.path.insert(0, BASE)
import run_species as R

sp = sys.argv[1]
genome = os.path.join(NEWSP, sp, "genome.fna")
rv_path = os.path.join(BASE, "%s_ref_verdicts.json" % sp)
if not os.path.exists(rv_path):
    print("  %s: nothing was refined" % sp); sys.exit(0)
rv = json.load(open(rv_path))

ref_aln = os.path.join(NEWSP, sp, "refaln")
todo = []
for k, d in rv.items():
    if "top100" not in k or "flags" not in d:
        continue
    cw = [f["n"] for f in d["flags"] if f["code"] == "CONSENSUS_OVEREXTENDED"]
    if not cw:
        continue
    tag = k.replace("NEW__", "").replace(".clean", "").replace("__top100", "")
    todo.append((tag, cw[0]))
if not todo:
    print("  %s: no overshoot to trim" % sp); sys.exit(0)

out_dir = os.path.join(NEWSP, sp, "ext2")
os.makedirs(out_dir, exist_ok=True)
made = []
for tag, (lo, hi) in todo:
    p = os.path.join(ref_aln, "%s__top100.aln.fa" % tag)
    if not os.path.exists(p):
        print("  no alignment for %s" % tag); continue
    seq, keep = [], False
    for line in open(p):
        if line.startswith(">"):
            if keep: break
            keep = "CONSENSUS_" in line.upper(); continue
        if keep: seq.append(line.strip())
    s = "".join(seq).replace("-", "")
    if hi - lo < 60:
        continue
    f = os.path.join(out_dir, "%s.fa" % tag)
    open(f, "w").write(">%s_p2\n%s\n" % (tag, s[lo:hi + 1]))
    made.append(f)
    print("  %-30s %d bp -> %d bp (positions %d-%d)" % (tag, len(s), hi - lo + 1, lo, hi))
if not made:
    sys.exit(0)

aln2 = os.path.join(NEWSP, sp, "ref2aln")
os.makedirs(aln2, exist_ok=True)
for f in glob.glob(os.path.join(aln2, "*")): os.remove(f)
for m in made:
    R.sh("python3 candidate_to_aln.py %s %s %s %sP" % (m, genome, aln2, sp))
if not glob.glob(os.path.join(aln2, "*.aln.fa")):
    print("  no alignments built"); sys.exit(0)
R.side_tables(aln2, sp + "_ref2")
cl2 = os.path.join(NEWSP, sp, "ref2clean")
R.clean(aln2, cl2)
v2 = R.score(cl2, "%s_ref2_verdicts.json" % sp)
for k, d in sorted(v2.items()):
    if "top100" not in k: continue
    print("  -> %-44s %5.1f  n_core %3d/%3d  %s"
          % (k.replace("NEW__", "")[:44], d["score"], d["n_core"], d["n"],
             ",".join(f["code"] for f in d["flags"]) or "clean"))
