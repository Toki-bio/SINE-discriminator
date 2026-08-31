#!/usr/bin/env python3
"""Build and score alternative flank-handling strategies for an alignment.

The problem: MAFFT is asked to align 70 bp of genuinely non-homologous flank per
copy, so it scatters those bases over 3x as many columns as there are bases.
That is not just ugly - any statistic computed on aligned flank columns is
inflated, because the aligner has co-located whatever happened to match.

Strategies produced (element alignment preserved except where stated):
  v1_current      as MAFFT left it
  v2_justify      flanks de-gapped and pushed against the element, no gaps at all
  v3_op3          flanks de-gapped then re-aligned, moderate gap penalty
  v4_op10         flanks de-gapped then re-aligned, high gap penalty
  v5_op20         flanks de-gapped then re-aligned, very high gap penalty
  v6_whole_op10   everything re-aligned from scratch at high gap penalty

Scoring. The honest reference is the ungapped, edge-anchored flank identity - no
aligner involved. A strategy is good when identity measured ON ITS ALIGNMENT
matches that reference instead of exceeding it, while the element stays compact
and well aligned.
"""
import os, sys, json, subprocess, shutil
import numpy as np
import measure_c as M

SETS = ["POS__ccr__g5_7seqs", "POS__saq__s5_5seqs", "NEGRAND__dmo__r00"]
OUT = "variants"


def split(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n][0]
    nz = np.where(A[ci] != M.GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    keep = [i for i in range(len(names)) if i != ci]
    rows = []
    for i in keep:
        r = A[i]
        lf = r[:lo]
        lf = lf[lf != M.GAP]
        rf = r[hi + 1:]
        rf = rf[rf != M.GAP]
        rows.append({"name": names[i], "lf": lf, "el": r[lo:hi + 1], "rf": rf})
    cons = {"name": names[ci], "lf": np.array([], dtype=np.int8),
            "el": A[ci][lo:hi + 1], "rf": np.array([], dtype=np.int8)}
    return cons, rows


def s(arr):
    return "".join("ACGT-"[c] for c in arr)


def write_aln(path, cons, rows, lfblocks, rfblocks):
    wL = max(len(x) for x in lfblocks.values()) if lfblocks else 0
    wR = max(len(x) for x in rfblocks.values()) if rfblocks else 0
    with open(path, "w") as fh:
        fh.write(">%s\n%s%s%s\n" % (cons["name"], "-" * wL, s(cons["el"]), "-" * wR))
        for i, r in enumerate(rows):
            L = lfblocks.get(i, "")
            R = rfblocks.get(i, "")
            fh.write(">%s\n%s%s%s\n" % (r["name"], L.rjust(wL, "-"),
                                        s(r["el"]), R.ljust(wR, "-")))


def mafft(seqs, op, ep, tag):
    """Align a list of (name, string) with the given gap penalties."""
    if not any(len(x[1]) > 0 for x in seqs):
        return {i: "" for i in range(len(seqs))}
    tmp_in = os.path.join(OUT, "_%s.in.fa" % tag)
    tmp_out = os.path.join(OUT, "_%s.out.fa" % tag)
    with open(tmp_in, "w") as fh:
        for i, (nm, sq) in enumerate(seqs):
            fh.write(">%d\n%s\n" % (i, sq if sq else "N"))
    cmd = ("mafft --retree 2 --maxiterate 0 --quiet --thread 2 --op %s --ep %s %s > %s"
           % (op, ep, tmp_in, tmp_out))
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write(r.stderr[-500:])
        return None
    out, nm = {}, None
    for line in open(tmp_out):
        line = line.rstrip()
        if line.startswith(">"):
            nm = int(line[1:].split()[0].replace("_R_", ""))
            out[nm] = ""
        elif nm is not None:
            out[nm] += line.upper()
    return out


def score(path):
    """Objective comparison, plus the honest ungapped reference."""
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n][0]
    nz = np.where(A[ci] != M.GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, ci, axis=0)

    def colpid(block):
        vals = []
        for c in range(block.shape[1]):
            b = block[:, c]
            b = b[b != M.GAP]
            if len(b) >= 6:
                cnt = np.bincount(b, minlength=4)[:4].astype(float)
                vals.append((cnt * (cnt - 1)).sum() / (len(b) * (len(b) - 1)))
        return float(np.mean(vals)) if vals else float("nan")

    Lf, El, Rf = C[:, :lo], C[:, lo:hi + 1], C[:, hi + 1:]
    fl_cols = Lf.shape[1] + Rf.shape[1]
    fl_bases = np.median((Lf != M.GAP).sum(axis=1) + (Rf != M.GAP).sum(axis=1))
    el_bases = np.median((El != M.GAP).sum(axis=1))
    aligned_flank = np.nanmean([colpid(Lf), colpid(Rf)])

    # honest reference: ungapped, walking outward from each copy's own edge
    lefts, rights = [], []
    for i in range(C.shape[0]):
        l = C[i, :lo]
        lefts.append(l[l != M.GAP][::-1])
        r = C[i, hi + 1:]
        rights.append(r[r != M.GAP])

    def ung(seqs):
        vals = []
        for off in range(1, 71):
            col = np.array([x[off - 1] if len(x) >= off else M.GAP for x in seqs],
                           dtype=np.int8)
            b = col[col != M.GAP]
            if len(b) >= 6:
                cnt = np.bincount(b, minlength=4)[:4].astype(float)
                vals.append((cnt * (cnt - 1)).sum() / (len(b) * (len(b) - 1)))
        return float(np.mean(vals)) if vals else float("nan")

    return {"width": int(A.shape[1]),
            "flank_cols_per_base": round(fl_cols / max(1, fl_bases), 2),
            "elem_cols_per_base": round(El.shape[1] / max(1, el_bases), 2),
            "flank_id_aligned": round(aligned_flank, 3),
            "flank_id_ungapped": round(np.nanmean([ung(lefts), ung(rights)]), 3),
            "elem_col_id": round(colpid(El), 3)}


def main():
    os.makedirs(OUT, exist_ok=True)
    report = {}
    for st in SETS:
        src = "aln_c/%s.aln.fa" % st
        if not os.path.exists(src):
            continue
        cons, rows = split(src)
        made = {}

        shutil.copy(src, "%s/%s__v1_current.aln.fa" % (OUT, st))
        made["v1_current"] = "%s/%s__v1_current.aln.fa" % (OUT, st)

        # v2: de-gapped, pushed against the element
        lfb = {i: s(r["lf"][::-1])[::-1] for i, r in enumerate(rows)}
        rfb = {i: s(r["rf"]) for i, r in enumerate(rows)}
        p = "%s/%s__v2_justify.aln.fa" % (OUT, st)
        write_aln(p, cons, rows, lfb, rfb)
        made["v2_justify"] = p

        # v3-v5: de-gapped then re-aligned at increasing gap penalties.
        # Left flanks are reversed before aligning so the edge adjacent to the
        # element is the anchored end, then reversed back.
        for tag, op, ep in (("v3_op3", 3.0, 0.5), ("v4_op10", 10.0, 1.0),
                            ("v5_op20", 20.0, 2.0)):
            la = mafft([(str(i), s(r["lf"])) for i, r in enumerate(rows)], op, ep, "l" + tag)
            ra = mafft([(str(i), s(r["rf"])) for i, r in enumerate(rows)], op, ep, "r" + tag)
            if la is None or ra is None:
                continue
            p = "%s/%s__%s.aln.fa" % (OUT, st, tag)
            write_aln(p, cons, rows,
                      {i: la.get(i, "") for i in range(len(rows))},
                      {i: ra.get(i, "") for i in range(len(rows))})
            made[tag] = p

        # v6: everything re-aligned from scratch at a high gap penalty
        whole = [(cons["name"], s(cons["el"]).replace("-", ""))] + \
                [(r["name"], s(r["lf"]) + s(r["el"]).replace("-", "") + s(r["rf"]))
                 for r in rows]
        wa = mafft(whole, 10.0, 1.0, "w" + st[:12])
        if wa:
            p = "%s/%s__v6_whole_op10.aln.fa" % (OUT, st)
            with open(p, "w") as fh:
                for i, (nm, _) in enumerate(whole):
                    fh.write(">%s\n%s\n" % (nm, wa.get(i, "")))
            made["v6_whole_op10"] = p

        report[st] = {k: score(v) for k, v in made.items()}

    json.dump(report, open("variant_report.json", "w"), indent=1)
    for st, r in report.items():
        print("\n=== %s" % st)
        print("%-16s %7s %10s %10s %11s %11s %9s"
              % ("variant", "width", "flank c/b", "elem c/b", "flank id", "ungapped", "elem id"))
        for k in ["v1_current", "v2_justify", "v3_op3", "v4_op10", "v5_op20", "v6_whole_op10"]:
            if k not in r:
                continue
            d = r[k]
            print("%-16s %7d %10.2f %10.2f %11.3f %11.3f %9.3f"
                  % (k, d["width"], d["flank_cols_per_base"], d["elem_cols_per_base"],
                     d["flank_id_aligned"], d["flank_id_ungapped"], d["elem_col_id"]))


if __name__ == "__main__":
    main()
