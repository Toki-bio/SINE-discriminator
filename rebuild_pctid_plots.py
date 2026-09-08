#!/usr/bin/env python3
"""Regenerate divergence picker on DRAGEN and patch oma_report.html."""
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
         f"{FRAG}/oma_divergence.js", str(SITE / "oma_divergence.js")],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(SITE / "apply_divergence_simple.py")],
        check=True,
    )
    print("done")


if __name__ == "__main__":
    main()
