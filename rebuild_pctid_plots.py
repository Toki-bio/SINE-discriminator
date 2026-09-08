#!/usr/bin/env python3
"""Regenerate pctid plots on DRAGEN, fetch, patch oma_report.html."""
import subprocess
import sys
from pathlib import Path

SITE = Path(__file__).parent
PLINK = r"C:\Program Files\PuTTY\plink.exe"
KEY = r"C:\Users\T\.ssh\id_ed25519.ppk"
HOST = "copilot@100.104.25.22"
RUN = "/staging/tmp/scorpions/oma/run_oma"
S2 = f"{RUN}/step2/step2_output"
FRAG = f"{RUN}/report_fragments"
PATCH = "/staging/tmp/sinedisc/patch_oma_report.py"


def plink(cmd: str) -> None:
    subprocess.run(
        [PLINK, "-ssh", HOST, "-batch", "-noagent", "-i", KEY, cmd],
        check=True,
    )


def main() -> None:
    subprocess.run(
        [sys.executable, str(SITE / "_upload_dragen_scripts.py"),
         str(SITE / "patch_oma_report.py"), PATCH],
        check=True,
    )
    plink(f"python3 {PATCH} x {S2} {S2}/plots {FRAG}")
    subprocess.run(
        [sys.executable, str(SITE / "fetch_dragen_file.py"),
         f"{FRAG}/violins.js", str(SITE / "report_fragments" / "violins.js")],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(SITE / "fetch_dragen_file.py"),
         f"{FRAG}/divergence_extra.js",
         str(SITE / "report_fragments" / "divergence_extra.js")],
        check=True,
    )
    js = (SITE / "report_fragments" / "divergence_extra.js").read_text(encoding="utf-8")
    marker = "Plotly.newPlot('plot_pctid_kde',"
    if marker not in js:
        raise SystemExit("plot_pctid_kde not in divergence_extra.js")
    kde_line = js[js.find(marker):].split("\n", 1)[0].strip()
    if not kde_line.endswith(";"):
        kde_line += ";"
    kde_path = SITE / "report_fragments" / "divergence_kde.js"
    kde_path.write_text(kde_line + "\n", encoding="utf-8")
    subprocess.run(
        [sys.executable, str(SITE / "patch_plot_block.py"),
         "plot_pctid_kde", str(kde_path)],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(SITE / "patch_violins_in_report.py")],
        check=True,
    )
    print("rebuilt oma_report.html")


if __name__ == "__main__":
    main()
