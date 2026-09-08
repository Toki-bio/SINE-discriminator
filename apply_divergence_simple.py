#!/usr/bin/env python3
"""Apply simple divergence section + oma_divergence.js to oma_report.html."""
import re
import sys
from pathlib import Path

from patch_oma_report import divergence_section_html

REPORT = Path(__file__).parent / "oma_report.html"
DIV_JS = Path(__file__).parent / "oma_divergence.js"


def strip_old_plots(text: str) -> str:
    for plot_id in ("plot_pctid_kde", "plot_sim_violins"):
        marker = "Plotly.newPlot('%s'," % plot_id
        while marker in text:
            start = text.index(marker)
            end = text.find(";\n", start)
            if end < 0:
                end = text.find(";", start)
            if end < 0:
                break
            text = text[:start] + text[end + 2:]
    return text


def main() -> None:
    js_src = DIV_JS
    if len(sys.argv) > 1:
        js_src = Path(sys.argv[1])
    js = js_src.read_text(encoding="utf-8")
    DIV_JS.write_text(js, encoding="utf-8")

    text = REPORT.read_text(encoding="utf-8")
    text = re.sub(
        r'<section class="card" id="divergence">.*?</section>\s*',
        divergence_section_html(),
        text,
        count=1,
        flags=re.DOTALL,
    )
    text = strip_old_plots(text)
    if "oma_divergence.js" not in text:
        text = text.replace(
            '<script src="oma_report_plots.js"></script>',
            '<script src="oma_divergence.js"></script>\n'
            '<script src="oma_report_plots.js"></script>',
        )
    REPORT.write_text(text, encoding="utf-8")
    print("patched", REPORT)


if __name__ == "__main__":
    main()
