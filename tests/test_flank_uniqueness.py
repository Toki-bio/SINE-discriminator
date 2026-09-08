#!/usr/bin/env python3
"""Tests for flank_uniqueness.scan — synthetic and SINE10 smoke."""
import hashlib
import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from flank_uniqueness import scan, side_report, cluster_labels, identity_matrix
import numpy as np

GT = os.path.join(
    REPO, "alignments", "oma", "oma_SINE10_top100.aln.fa")


def write_aln(path, left_w, elem, right_w, copies):
    """Build a gapped MSA: [left flank][element upper][right flank]."""
    cons_seq = "-" * left_w + elem + "-" * right_w
    with open(path, "w") as fh:
        fh.write(">%s\n%s\n" % ("CONSENSUS", cons_seq))
        for h, left, right in copies:
            fh.write(">%s\n%s\n" % (h, left + elem + right))


def _uniq_flank(i, width, salt=0):
    raw = hashlib.sha256(("%d:%d" % (i, salt)).encode()).digest()
    out = []
    while len(out) < width:
        for b in raw:
            out.append("ACGT"[b % 4])
            if len(out) >= width:
                break
        raw = hashlib.sha256(raw).digest()
    return "".join(out)


class SyntheticFlankClusters(unittest.TestCase):
    def test_detects_shared_block(self):
        shared = ("ATCG" * 20)[:80]
        elem = "G" * 50
        right = "a" * 30
        copies = []
        for i in range(20):
            copies.append(("shared_%d" % i, shared.lower(), right))
        for i in range(20):
            uniq = ("ACGTN" * 20 + str(i) * 5)[:80]
            copies.append(("uniq_%d" % i, uniq.lower(), right))
        td = tempfile.mkdtemp()
        p = os.path.join(td, "oma_test_rand100.aln.fa")
        write_aln(p, 80, elem, 30, copies)
        r = scan(p)
        self.assertTrue(r["left"]["measured"])
        self.assertGreaterEqual(r["left"]["largest_cluster"], 10)
        self.assertGreater(r["left"]["shared_copy_frac"], 0.3)
        self.assertEqual(r["tier"], "rand100")
        self.assertIn(r["worst_flag"], ("medium", "high"))

    def test_all_unique_low_flag(self):
        elem = "G" * 50
        copies = []
        for i in range(30):
            copies.append(("u%d" % i, _uniq_flank(i, 80, 0), _uniq_flank(i, 60, 1)))
        td = tempfile.mkdtemp()
        p = os.path.join(td, "oma_test_top100.aln.fa")
        write_aln(p, 80, elem, 60, copies)
        r = scan(p)
        self.assertTrue(r["left"]["measured"])
        self.assertTrue(r["right"]["measured"])
        self.assertIsNone(r["worst_flag"])
        self.assertGreater(r["left"]["n_clusters"], 20)
        self.assertGreater(r["right"]["n_clusters"], 20)

    def test_identical_right_flanks_flags(self):
        """All copies share one right flank — must not pass as 'unique'."""
        elem = "G" * 50
        shared_r = ("ATCG" * 15)[:60]
        copies = []
        for i in range(25):
            copies.append(("u%d" % i, _uniq_flank(i, 80, 0), shared_r))
        td = tempfile.mkdtemp()
        p = os.path.join(td, "oma_test_top100.aln.fa")
        write_aln(p, 80, elem, 60, copies)
        r = scan(p)
        self.assertTrue(r["left"]["measured"])
        self.assertIsNone(next((f for f in r["flags"] if f["side"] == "left"), None))
        rfl = next(f for f in r["flags"] if f["side"] == "right")
        self.assertEqual(rfl["severity"], "high")
        self.assertEqual(r["worst_flag"], "high")


class ClusterLabels(unittest.TestCase):
    def test_single_linkage(self):
        M = np.array([
            [1.0, 0.9, 0.1],
            [0.9, 1.0, 0.1],
            [0.1, 0.1, 1.0],
        ])
        lb = cluster_labels(M, thr=0.55)
        self.assertEqual(lb[0], lb[1])
        self.assertNotEqual(lb[0], lb[2])


class SINE10GroundTruth(unittest.TestCase):
    @unittest.skipUnless(os.path.isfile(GT), "missing SINE10 alignment")
    def test_independent_flanks_no_flag(self):
        r = scan(GT)
        self.assertEqual(r["set"], "oma_SINE10_top100")
        self.assertEqual(r["tier"], "top100")
        self.assertTrue(r["left"]["measured"])
        self.assertTrue(r["right"]["measured"])
        self.assertIsNone(r["worst_flag"])
        for side in ("left", "right"):
            s = r[side]
            self.assertGreaterEqual(s["unique_frac"], 0.85)
            self.assertLessEqual(s["largest_cluster"], 3)
            self.assertLess(s["largest_cluster_frac"], 0.10)


if __name__ == "__main__":
    unittest.main()
