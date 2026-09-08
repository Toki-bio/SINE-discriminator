#!/usr/bin/env python3
"""Keep only ssearch36 %identity KDE in oma_report divergence section."""
import re
from pathlib import Path

REPORT = Path(__file__).parent / "oma_report.html"

NEW_INTRO = (
    '<p class="intro"><b>Metric:</b> divergence = 100 &minus; ssearch36 '
    '%identity to the subfamily consensus (same as the Gallery histograms). '
    'One KDE curve per subfamily (up to 3,000 assigned copies sampled); '
    'click legend entries to toggle.</p>\n    '
    '<h3>ssearch36 %identity divergence (step4)</h3>\n    '
    '<div class="plot" id="plot_pctid_kde"></div>'
)

NEW_SOURCE = (
    '<p class="small muted">Source: step4 <code>*_pctid.tsv</code> on assigned copies.</p>'
)


def main():
    text = REPORT.read_text(encoding="utf-8")

    text, n = re.subn(
        r'<section class="card" id="divergence">.*?<div class="plot" id="plot_div_kde"></div>.*'
        r'<div class="plot" id="plot_pctid_kde"></div>',
        '<section class="card" id="divergence">\n'
        '    <h2>Divergence from consensus &mdash; per copy</h2>\n    '
        + NEW_INTRO,
        text,
        count=1,
        flags=re.DOTALL,
    )
    if n != 1:
        raise SystemExit("divergence section pattern not found (%d)" % n)

    text = re.sub(
        r'<p class="small muted">Source: <code>sim_scores\.tsv</code>.*?</p>',
        NEW_SOURCE,
        text,
        count=1,
        flags=re.DOTALL,
    )

    for plot_id in ("plot_div_kde_filtered", "plot_div_kde"):
        pat = re.compile(
            r"Plotly\.newPlot\('%s',.*?\);\n" % re.escape(plot_id),
            re.DOTALL,
        )
        text, c = pat.subn("", text, count=1)
        if c != 1:
            raise SystemExit("missing Plotly block for %s" % plot_id)

    REPORT.write_text(text, encoding="utf-8")
    print("patched", REPORT)


if __name__ == "__main__":
    main()
