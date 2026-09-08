#!/usr/bin/env python3
"""Apply spline divergence plot (variant 3) to oma_report.html."""
import re
import sys
from pathlib import Path

from patch_oma_report import divergence_section_html

REPORT = Path(__file__).parent / "oma_report.html"
PLOT_JS = Path(__file__).parent / "oma_divergence.js"

PLOT_IDS = (
    "plot_divergence",
    "plot_div_kde",
    "plot_pctid_kde",
    "plot_div_kde_filtered",
    "plot_sim_violins",
)


def strip_plotly_calls(text: str) -> str:
    for plot_id in PLOT_IDS:
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


def strip_nav_gallery_link(text: str) -> str:
    text = re.sub(
        r'\s*<span style="opacity:\.45">&bull;</span>\s*\n'
        r'\s*<a href="oma_divergence_gallery\.html">Chart variants \(pick one\)</a>\n',
        "\n",
        text,
        count=1,
    )
    return text


def inject_plot_script(text: str, plot_js: str) -> str:
    plot_js = plot_js.strip()
    if not plot_js.endswith(";"):
        plot_js += ";"
    block = "<script>\n%s\n</script>\n" % plot_js
    text = strip_plotly_calls(text)
    text = text.replace('<script src="oma_divergence.js"></script>\n', "")
    anchor = '<script src="oma_report_plots.js"></script>'
    if anchor in text:
        if "Plotly.newPlot('plot_div_kde'," not in text:
            text = text.replace(anchor, block + anchor)
    else:
        text = text.replace("</body></html>", block + "</body></html>")
    return text


def main() -> None:
    js_src = PLOT_JS
    if len(sys.argv) > 1:
        js_src = Path(sys.argv[1])
    plot_js = js_src.read_text(encoding="utf-8")

    text = REPORT.read_text(encoding="utf-8")
    text = re.sub(
        r'<section class="card" id="divergence">.*?</section>\s*',
        divergence_section_html(),
        text,
        count=1,
        flags=re.DOTALL,
    )
    text = strip_nav_gallery_link(text)
    text = inject_plot_script(text, plot_js)
    REPORT.write_text(text, encoding="utf-8")
    print("patched", REPORT)


if __name__ == "__main__":
    main()
