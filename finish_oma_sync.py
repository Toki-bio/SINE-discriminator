#!/usr/bin/env python3
"""Extract downloaded tarball and patch oma_report (local only)."""
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

SITE = Path(__file__).resolve().parent
TAR = SITE / "_oma_site_sync2.tgz"
if not TAR.is_file():
    TAR = SITE / "_oma_site_sync.tgz"
ALN = SITE / "alignments" / "oma"
PLOTS = ALN / "plots"
DROP = {"oma_big76", "oma_group34", "oma_sub515"}


def extract() -> Path:
    ext = SITE / "_oma_sync_extract"
    if ext.exists():
        shutil.rmtree(ext)
    for p in ALN.glob("oma_*.aln.fa"):
        p.unlink()
    for p in PLOTS.glob("*"):
        if p.is_file():
            p.unlink()
    ALN.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    with tarfile.open(TAR, "r:gz") as tf:
        tf.extractall(ext)
    src = ext / "alignments"
    for p in src.glob("oma_oma_*.aln.fa"):
        out = ALN / p.name.replace("oma_oma_", "oma_", 1)
        if any(d in out.name for d in DROP):
            continue
        shutil.copy2(p, out)
    for p in src.glob("oma_*.aln.fa"):
        if p.name.startswith("oma_oma_"):
            continue
        if any(d in p.name for d in DROP):
            continue
        shutil.copy2(p, ALN / p.name)
    plots_src = src / "oma" / "plots"
    if plots_src.is_dir():
        for p in plots_src.iterdir():
            if p.is_file():
                shutil.copy2(p, PLOTS / p.name)
    summary = SITE / "_summary.by_subfam.tsv"
    shutil.copy2(
        ext / "step2/step2_output_rerun_20260908_run1/summary.by_subfam.tsv",
        summary,
    )
    shutil.rmtree(ext)
    n = len(list(ALN.glob("oma_*_top100.aln.fa")))
    print(
        f"extracted {n} top100 alignments, "
        f"{len(list(PLOTS.glob('*_pctid.tsv')))} pctid"
    )
    return summary


def patch_counts(report: Path, summary: Path) -> None:
    counts = {}
    for i, line in enumerate(summary.read_text(encoding="utf-8").splitlines()):
        if i == 0 or not line.strip():
            continue
        p = line.split("\t")
        if p[0] in DROP:
            continue
        try:
            counts[p[0]] = int(p[3])
        except (IndexError, ValueError):
            pass
    text = report.read_text(encoding="utf-8")

    def repl(m: re.Match) -> str:
        sf = m.group(1)
        if sf in DROP:
            return ""
        n = counts.get(sf)
        if n is None:
            return m.group(0)
        return (
            f"<tr><td>{sf}</td><td>{m.group(2)}</td><td>{m.group(3)}</td>"
            f"<td class=\"hl\">{n:,}</td></tr>"
        )

    text = re.sub(
        r"<tr><td>(oma_[^<]+)</td><td>(\d+)</td><td>(\d+)</td>"
        r"<td class=\"hl\">\d+</td></tr>",
        repl,
        text,
    )
    text = text.replace("(26 families)", "(23 families)")
    text = text.replace("26 subfamilies", "23 subfamilies")
    text = text.replace("26 families", "23 families")
    report.write_text(text, encoding="utf-8")


def main() -> None:
    if not TAR.is_file():
        sys.exit(f"missing {TAR}")
    summary = extract()
    report = SITE / "oma_report.html"
    patch_counts(report, summary)
    py = sys.executable
    subprocess.run(
        [py, str(SITE / "patch_oma_report.py"), str(report), str(PLOTS)],
        check=True,
    )
    subprocess.run([py, str(SITE / "apply_divergence_simple.py")], check=True)
    subprocess.run(
        [py, str(SITE / "inject_oma_aln_section.py"), str(report), str(ALN)],
        check=True,
    )
    n = len(list(ALN.glob("oma_*_top100.aln.fa")))
    print(f"patched {report} ({n} subfamilies in alignments dir)")


if __name__ == "__main__":
    main()
