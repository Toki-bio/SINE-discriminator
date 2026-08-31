#!/usr/bin/env python3
"""Ground-truth test: damage the consensus edge by a known amount, then ask the
search to recover it. If the reported d_best does not track the damage, the
method does not work and nothing downstream of it should be believed.

The damage is applied by moving the anchor inward, which is exactly what a
consensus trimmed too short would look like.
"""
import sys
import numpy as np
import measure_c as M
import edges as E


def damaged_profiles(path, trimL, trimR):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n][0]
    cons = A[ci]
    nzc = np.where(cons != M.GAP)[0]
    # pretend the consensus started trimL later and ended trimR earlier
    nz = nzc[trimL:len(nzc) - trimR] if trimR else nzc[trimL:]
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, ci, axis=0)
    n = C.shape[0]
    inside = np.array([E.pid(C[:, c]) for c in nz])
    lefts, rights = [], []
    for i in range(n):
        l = C[i, :lo]
        lefts.append(l[l != M.GAP][::-1])
        r = C[i, hi + 1:]
        rights.append(r[r != M.GAP])

    def fp(seqs, kmax):
        return np.array([E.pid(np.array([s[o - 1] if len(s) >= o else M.GAP
                                         for s in seqs], dtype=np.int8))
                         for o in range(1, kmax + 1)])
    return inside, fp(lefts, E.SEARCH + E.W + 5), fp(rights, E.SEARCH + E.W + 5)


def main():
    sets = ["POS__saq__s5_5seqs", "POS__saq__s1_30seqs", "POS__teu__t2_75seqs",
            "POS__dmo__d4_266seqs", "POS__ccr__g3_71seqs"]
    trims = [0, 10, 20, 40]
    print("RECOVERY TEST - consensus deliberately trimmed, can the search find it back?")
    print("d_best should equal -trim (negative = extend outward by that much)\n")
    print("%-24s %6s | %-22s %-22s" % ("set", "trim", "LEFT d_best (want)", "RIGHT d_best (want)"))
    print("-" * 80)
    errL, errR = [], []
    for s in sets:
        for t in trims:
            ins, lf, rf = damaged_profiles("aln_c/%s.aln.fa" % s, t, t)
            L = E.scan_edge(ins, lf)
            R = E.scan_edge(ins[::-1], rf)
            if not L or not R:
                continue
            errL.append(L["d_best"] + t)
            errR.append(R["d_best"] + t)
            print("%-24s %6d | %6d (%-4d) %-8s %6d (%-4d)"
                  % (s[:24], t, L["d_best"], -t, "", R["d_best"], -t))
    print("-" * 80)
    print("left  error: mean %+.1f bp, median %+.1f, |err|<=5 in %d/%d"
          % (np.mean(errL), np.median(errL), sum(abs(np.array(errL)) <= 5), len(errL)))
    print("right error: mean %+.1f bp, median %+.1f, |err|<=5 in %d/%d"
          % (np.mean(errR), np.median(errR), sum(abs(np.array(errR)) <= 5), len(errR)))


if __name__ == "__main__":
    main()
