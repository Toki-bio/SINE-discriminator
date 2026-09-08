#!/usr/bin/env python3
"""Download a file from DRAGEN via base64 (avoids banner in text cat)."""
import base64
import re
import subprocess
import sys

PLINK = r"C:\Program Files\PuTTY\plink.exe"
KEY = r"C:\Users\T\.ssh\id_ed25519.ppk"
HOST = "copilot@100.104.25.22"


def download(remote: str, local: str) -> None:
    r = subprocess.run(
        [PLINK, "-ssh", HOST, "-batch", "-noagent", "-i", KEY,
         "base64 -w0 %s" % remote],
        capture_output=True, check=True,
    )
    raw = r.stdout.decode("ascii", errors="ignore")
    # strip DRAGEN login banner junk before base64 payload
    m = re.search(r"([A-Za-z0-9+/=]{40,})$", raw.replace("\n", "").replace("\r", ""))
    if not m:
        raise SystemExit("no base64 payload for %s" % remote)
    data = base64.b64decode(m.group(1))
    open(local, "wb").write(data)
    print("wrote", local, len(data), "bytes")


if __name__ == "__main__":
    download(sys.argv[1], sys.argv[2])
