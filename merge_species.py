"""Combine every species run into the two files the page builder reads.

newsp_verdicts_all.json  every view of every candidate, with the REFINED verdict
                         substituted where the consensus was repaired, and the
                         original kept alongside so the repair is visible
seeds.json               seed RNA and AnnoSINE copy count per candidate, read
                         from each genome Seed_SINE.fa header rather than typed
"""
import glob, json, os, re

SP = ["hyd", "pom", "aca", "stu", "ska"]
BASE = "/staging/tmp/sinedisc"
NEWSP = "/staging/tmp/newsp"

out, refined = {}, {}
for sp in SP:
    p = os.path.join(BASE, "%s_verdicts.json" % sp)
    if os.path.exists(p):
        out.update(json.load(open(p)))
    for suffix in ("_ref_verdicts.json", "_ref2_verdicts.json"):
        r = os.path.join(BASE, sp + suffix)
        if os.path.exists(r):
            for k, v in json.load(open(r)).items():
                refined[k] = v

# map a refined set back to the candidate it came from:
#   NEW__<sp>R_<tag>_ref__top100.clean  or  ..._ext__top100.clean
for k, v in refined.items():
    base = k.replace("NEW__", "").replace(".clean", "")
    m = re.match(r"([a-z]+[RP]?)_(.+?)__(top100|rand100|all)$", base)
    if not m:
        continue
    tag, view = m.group(2), m.group(3)
    # strip every repair-pass wrapper back to the original candidate name
    for _ in range(4):
        tag = re.sub(r"^[a-z]{3}[RP]_", "", tag)
        tag = re.sub(r"_(ref|ext|p2)$", "", tag)
    orig = "NEW__%s__%s.clean" % (tag, view)
    if orig in out:
        out[orig]["refined"] = {kk: v.get(kk) for kk in
                                ("score", "n", "n_core", "core_frac", "all_identity",
                                 "flank_bg", "flags")}
    out[k] = v

json.dump(out, open(os.path.join(BASE, "newsp_verdicts_all.json"), "w"),
          separators=(",", ":"))

seeds = {}
for sp in SP:
    f = os.path.join(NEWSP, sp, "anno", "Seed_SINE.fa")
    if not os.path.exists(f):
        continue
    for line in open(f):
        if not line.startswith(">"):
            continue
        cid = line[1:].split()[0]
        cnt = re.search(r"blast_count:(\d+)", line)
        rna = line.rstrip().split("|")[-1]
        seeds["%s_%s" % (sp, cid)] = [rna.replace("_", " "),
                                      int(cnt.group(1)) if cnt else None]
json.dump(seeds, open(os.path.join(BASE, "seeds.json"), "w"), indent=0, sort_keys=True)

print("merged %d verdict entries (%d refined), %d seeds"
      % (len(out), len(refined), len(seeds)))
for sp in SP:
    n = sum(1 for k in out if k.replace("NEW__", "").startswith(sp + "_"))
    print("  %-4s %d" % (sp, n))
