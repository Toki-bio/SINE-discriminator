#!/usr/bin/env python3
"""Replace a Plotly.newPlot(...) block in oma_report.html."""
import sys
from pathlib import Path

REPORT = Path(__file__).parent / "oma_report.html"


def replace_block(plot_id: str, js_path: Path) -> None:
    js = js_path.read_text(encoding="utf-8").strip()
    if not js.endswith(";"):
        js += ";"
    js += "\n"
    text = REPORT.read_text(encoding="utf-8")
    marker = "Plotly.newPlot('%s'," % plot_id
    start = text.find(marker)
    if start < 0:
        raise SystemExit("block not found: %s" % plot_id)
    end = text.find(";\n", start)
    if end < 0:
        end = text.find(";", start)
    if end < 0:
        raise SystemExit("block end not found: %s" % plot_id)
    REPORT.write_text(text[:start] + js + text[end + 2 :], encoding="utf-8")
    print("replaced", plot_id)


if __name__ == "__main__":
    replace_block(sys.argv[1], Path(sys.argv[2]))
