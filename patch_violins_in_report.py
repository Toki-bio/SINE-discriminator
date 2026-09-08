#!/usr/bin/env python3
"""Replace plot_sim_violins Plotly block in oma_report.html."""
import re
from pathlib import Path

REPORT = Path(__file__).parent / "oma_report.html"
VIOLINS = Path(__file__).parent / "report_fragments" / "violins.js"


def main():
    js = VIOLINS.read_text(encoding="utf-8").strip()
    if not js.endswith(";"):
        js += ";"
    js += "\n"
    text = REPORT.read_text(encoding="utf-8")
    marker = "Plotly.newPlot('plot_sim_violins',"
    start = text.find(marker)
    if start < 0:
        raise SystemExit("plot_sim_violins block not found")
    end = text.find(";\n", start)
    if end < 0:
        end = text.find(";", start)
    if end < 0:
        raise SystemExit("plot_sim_violins block end not found")
    text = text[:start] + js + text[end + 2 :]
    REPORT.write_text(text, encoding="utf-8")
    print("replaced violins in", REPORT)


if __name__ == "__main__":
    main()
