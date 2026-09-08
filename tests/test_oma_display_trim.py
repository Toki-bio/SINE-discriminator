#!/usr/bin/env python3
"""Regression: oma_SINE10_top100 display trim ground truth."""
import os
import shutil
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from fix_alignments import read_fa
from trim_display_flanks import compute_window, trim_inplace

GT_PRETRIM = os.path.join(
    REPO, "alignments", "oma", "fixtures",
    "oma_SINE10_top100_pretrim.aln.fa")
GT = os.path.join(
    REPO, "alignments", "oma", "oma_SINE10_top100.aln.fa")
WANT_LEFT, WANT_RIGHT = 430, 969
WANT_ELEM_LEN = 243


class SINE10DisplayTrim(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.isfile(GT_PRETRIM):
            raise unittest.SkipTest("missing %s" % GT_PRETRIM)
        if not os.path.isfile(GT):
            raise unittest.SkipTest("missing %s" % GT)

    def _copy(self):
        td = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(td, ignore_errors=True))
        path = os.path.join(td, "oma_SINE10_top100.aln.fa")
        shutil.copy(GT_PRETRIM, path)
        return path

    def test_window_before_trim(self):
        names, seqs = read_fa(GT_PRETRIM)
        left, right, diag = compute_window(seqs, 0, mode="occupancy")
        self.assertIsNotNone(left)
        self.assertEqual(left, WANT_LEFT)
        self.assertEqual(right, WANT_RIGHT)
        self.assertGreaterEqual(diag["median_left_bp"], 25)
        self.assertGreaterEqual(diag["median_right_bp"], 25)

    def test_trim_preserves_element(self):
        path = self._copy()
        names, seqs = read_fa(path)
        elem_before = "".join(c for c in seqs[0] if c.isupper())
        d = trim_inplace(path, mode="occupancy")
        self.assertTrue(d.get("elem_unchanged"))
        self.assertEqual(d.get("elem_len_before"), d.get("elem_len_after"))
        names2, seqs2 = read_fa(path)
        self.assertEqual(len(seqs2[0]), WANT_RIGHT - WANT_LEFT + 1)

    def test_published_trimmed_width(self):
        """Published file is already display-trimmed."""
        names, seqs = read_fa(GT)
        self.assertEqual(len(seqs[0]), WANT_RIGHT - WANT_LEFT + 1)
        elem = "".join(c for c in seqs[0] if c.isupper())
        self.assertEqual(len(elem), WANT_ELEM_LEN)

    def test_hybrid_same_on_sine10(self):
        names, seqs = read_fa(GT_PRETRIM)
        lo, ro, _ = compute_window(seqs, 0, mode="occupancy")
        lh, rh, _ = compute_window(seqs, 0, mode="hybrid")
        self.assertEqual((lo, ro), (lh, rh))


class TrimGuards(unittest.TestCase):
    def test_refuse_short_flank(self):
        elem = "A" * 40
        rows = []
        for i in range(12):
            rows.append(">c%d" % i)
            rows.append("-" * 50 + elem + "-" * 50)
        td = tempfile.mkdtemp()
        path = os.path.join(td, "t.aln.fa")
        with open(path, "w") as fh:
            fh.write("\n".join(rows) + "\n")
        names, seqs = read_fa(path)
        left, right, diag = compute_window(seqs, 0, mode="occupancy", min_flank_bp=25)
        self.assertIsNone(left)
        self.assertEqual(diag["reason"], "flank_too_short")
        import shutil
        shutil.rmtree(td)


if __name__ == "__main__":
    unittest.main()
