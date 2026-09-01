"""Islands must be measured on the FULL 400 bp flank, not the trimmed display.

Trimming cuts the flank to the width the copies actually fill, which is right
for looking at a boundary and wrong for this: it throws away exactly the far-off
columns Sergei said the islands sit in. aca_SINE_0 reads 839 island columns at
full width and 6 after trimming.
"""
import glob, json, os, sys
sys.path.insert(0, "/staging/tmp/sinedisc")
import newsp_all as N

out = {}
for sp in ("pom", "aca", "hyd"):
    for f in sorted(glob.glob("/staging/tmp/newsp/%s/aln/*.aln.fa" % sp)):
        b = os.path.basename(f).replace(".aln.fa", "")
        v = N.island_scan(f)
        if v:
            out[b] = v
json.dump(out, open("/staging/tmp/sinedisc/flank_islands_raw.json", "w"),
          indent=1, sort_keys=True)
print("raw-flank island scan: %d sets" % len(out))
