#!/usr/bin/env python3
"""Build oma_divergence_gallery.html — 10 binned divergence chart variants."""
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from patch_oma_report import SF_PALETTE, bin_divergence, load_pctid_by_sf

BIN = 1.0
HOVER = (
    "%{fullData.name}<br>divergence ~%{x:.0f}%<br>copies %{y:,d}<extra></extra>"
)


def bin_series(by_sf: Dict[str, List[float]]) -> Dict[str, Tuple[List[float], List[float]]]:
    out = {}
    for sf, vals in by_sf.items():
        counts = bin_divergence(vals, BIN)
        if not counts:
            continue
        bins = sorted(counts)
        x = [round(b + BIN / 2.0, 1) for b in bins]
        y = [counts[b] for b in bins]
        out[sf] = (x, y)
    return out


def x_max(series: Dict[str, Tuple[List[float], List[float]]]) -> float:
    hi = 0.0
    for x, _ in series.values():
        if x:
            hi = max(hi, max(x))
    return min(100.0, max(5.0, math.ceil((hi + BIN) / 5.0) * 5.0))


def base_layout(title: str, xmax: float, height: int = 380) -> dict:
    return {
        "title": {"text": title, "font": {"size": 14}},
        "xaxis": {
            "title": "Divergence (100 − %identity to consensus)",
            "range": [0, xmax],
            "dtick": 5,
        },
        "yaxis": {"title": "Copies", "rangemode": "tozero"},
        "height": height,
        "margin": {"t": 48, "r": 16, "b": 48, "l": 52},
        "showlegend": True,
        "legend": {"font": {"size": 9}},
    }


def line_traces(series, sf_sorted, line_kw=None, marker=False, opacity=1.0):
    traces = []
    for i, sf in enumerate(sf_sorted):
        if sf not in series:
            continue
        x, y = series[sf]
        kw = {"width": 2, "color": SF_PALETTE[i % len(SF_PALETTE)]}
        if line_kw:
            kw.update(line_kw)
        if opacity < 1:
            c = kw["color"]
            if c.startswith("#") and len(c) == 7:
                kw["color"] = c
        tr = {
            "type": "scatter",
            "mode": "lines+markers" if marker else "lines",
            "x": x,
            "y": y,
            "name": sf,
            "line": kw,
            "opacity": opacity,
            "hovertemplate": HOVER,
        }
        if marker:
            tr["marker"] = {"size": 4}
        traces.append(tr)
    return traces


def build_variants(by_sf: Dict[str, List[float]]) -> List[Tuple[str, str, dict]]:
    series = bin_series(by_sf)
    sf_sorted = sorted(series.keys())
    xmax = x_max(series)
    top8 = sorted(by_sf.keys(), key=lambda s: len(by_sf[s]), reverse=True)[:8]
    top10 = sorted(by_sf.keys(), key=lambda s: len(by_sf[s]), reverse=True)[:10]

    variants = []

    # 1 — Tal ccr reference
    variants.append((
        "1 — Tal ccr (frequency polygon)",
        "Straight lines through 1% bin centres; one line per subfamily; "
        "legend toggles. This is what Tal ccr/toc reports use.",
        {
            "data": line_traces(series, sf_sorted),
            "layout": base_layout("Tal ccr style — straight lines", xmax),
        },
    ))

    # 2 — Markers at each bin
    variants.append((
        "2 — Lines + markers",
        "Same bins; a visible dot at each (divergence, copy-count) point.",
        {
            "data": line_traces(series, sf_sorted, marker=True),
            "layout": base_layout("Lines with markers at bin centres", xmax),
        },
    ))

    # 3 — Spline smooth through bin points
    variants.append((
        "3 — Smooth spline through bin points",
        "Plotly spline through the same 1% bin counts — smooth curve, "
        "not KDE (still Y = copy count).",
        {
            "data": line_traces(series, sf_sorted, line_kw={"shape": "spline"}),
            "layout": base_layout("Spline-smoothed frequency polygon", xmax),
        },
    ))

    # 4 — Step (horizontal-vertical)
    variants.append((
        "4 — Step line (hv)",
        "Classic step chart for binned counts: flat over each 1% bin, "
        "then step to next count.",
        {
            "data": line_traces(series, sf_sorted, line_kw={"shape": "hv"}),
            "layout": base_layout("Step line (hv) — bin-native", xmax),
        },
    ))

    # 5 — Semi-transparent overlay
    variants.append((
        "5 — Semi-transparent overlay",
        "All subfamilies; 35% opacity so overlapping lines remain readable.",
        {
            "data": line_traces(series, sf_sorted, opacity=0.35),
            "layout": base_layout("Transparent overlay — all subfamilies", xmax),
        },
    ))

    # 6 — Stacked bars (step6 barmode)
    bar_traces = []
    for i, sf in enumerate(sf_sorted):
        if sf not in series:
            continue
        x, y = series[sf]
        bar_traces.append({
            "type": "bar",
            "x": x,
            "y": y,
            "name": sf,
            "width": BIN * 0.92,
            "marker": {"color": SF_PALETTE[i % len(SF_PALETTE)]},
            "hovertemplate": HOVER,
        })
    variants.append((
        "6 — Stacked bars (1% bins)",
        "Histogram bars stacked by subfamily — step6_report barmode stack.",
        {
            "data": bar_traces,
            "layout": {**base_layout("Stacked 1% bins", xmax),
                       "barmode": "stack", "bargap": 0.02},
        },
    ))

    # 7 — Top 8 only (less clutter)
    variants.append((
        "7 — Top 8 subfamilies only",
        "Largest 8 by copy count; spline smooth.",
        {
            "data": line_traces(
                series, [s for s in sf_sorted if s in top8],
                line_kw={"shape": "spline"}),
            "layout": base_layout("Top 8 — spline", xmax),
        },
    ))

    # 8 — Relative % within subfamily
    rel_traces = []
    for i, sf in enumerate(sf_sorted):
        if sf not in series:
            continue
        x, y = series[sf]
        n = sum(y) or 1
        rel_traces.append({
            "type": "scatter",
            "mode": "lines",
            "x": x,
            "y": [round(100.0 * v / n, 2) for v in y],
            "name": sf,
            "line": {"width": 2, "color": SF_PALETTE[i % len(SF_PALETTE)],
                     "shape": "spline"},
            "hovertemplate": (
                "%{fullData.name}<br>divergence ~%{x:.0f}%"
                "<br>% of copies %{y:.1f}<extra></extra>"),
        })
    variants.append((
        "8 — Relative frequency (% within subfamily)",
        "Y = percent of that subfamily's copies in each bin — "
        "compare shape when copy totals differ wildly.",
        {
            "data": rel_traces,
            "layout": {
                **base_layout("Relative % — spline", xmax),
                "yaxis": {"title": "% of subfamily copies", "rangemode": "tozero"},
            },
        },
    ))

    # 9 — Filled area (to zero)
    area_traces = []
    for i, sf in enumerate(sf_sorted):
        if sf not in series:
            continue
        x, y = series[sf]
        col = SF_PALETTE[i % len(SF_PALETTE)]
        area_traces.append({
            "type": "scatter",
            "mode": "lines",
            "x": x,
            "y": y,
            "name": sf,
            "fill": "tozeroy",
            "fillcolor": col,
            "line": {"width": 1, "color": col, "shape": "spline"},
            "opacity": 0.25,
            "hovertemplate": HOVER,
        })
    variants.append((
        "9 — Filled area (spline)",
        "Area under each subfamily's curve; low opacity.",
        {
            "data": area_traces,
            "layout": base_layout("Filled area — spline", xmax),
        },
    ))

    # 10 — Small multiples: top 10 subfamilies
    mini = []
    for j, sf in enumerate(top10):
        if sf not in series:
            continue
        x, y = series[sf]
        mini.append({
            "id": "mini_%d" % j,
            "sf": sf,
            "n": len(by_sf[sf]),
            "fig": {
                "data": [{
                    "type": "scatter",
                    "mode": "lines",
                    "x": x,
                    "y": y,
                    "line": {"width": 2, "color": "#4C72B0", "shape": "spline"},
                    "fill": "tozeroy",
                    "fillcolor": "rgba(76,114,176,0.15)",
                    "hovertemplate": HOVER.replace("%{fullData.name}<br>", ""),
                }],
                "layout": {
                    "title": {"text": "%s (%d)" % (sf, len(by_sf[sf])),
                              "font": {"size": 11}},
                    "xaxis": {"range": [0, xmax], "dtick": 10, "title": ""},
                    "yaxis": {"rangemode": "tozero", "title": ""},
                    "height": 220,
                    "margin": {"t": 28, "r": 8, "b": 28, "l": 36},
                    "showlegend": False,
                },
            },
        })
    variants.append((
        "10 — Small multiples (top 10)",
        "One spline line graph per large subfamily — no overlap clutter.",
        {"minis": mini, "grid": True},
    ))

    return variants


def render_html(variants: List[Tuple[str, str, dict]]) -> str:
    blocks = []
    scripts = []
    for i, (title, blurb, fig) in enumerate(variants):
        if fig.get("grid"):
            cells = []
            for m in fig["minis"]:
                cells.append(
                    '<div class="mini-cell"><div class="plot" id="%s"></div></div>'
                    % m["id"])
                scripts.append(
                    "Plotly.newPlot(%r,%s,%s,{displayModeBar:false});\n"
                    % (m["id"], json.dumps(m["fig"]["data"]),
                       json.dumps(m["fig"]["layout"])))
            plot_html = '<div class="mini-grid">' + "".join(cells) + "</div>"
        else:
            pid = "vplot_%d" % i
            plot_html = '<div class="plot" id="%s"></div>' % pid
            scripts.append(
                "Plotly.newPlot(%r,%s,%s,{displayModeBar:true});\n"
                % (pid, json.dumps(fig["data"]), json.dumps(fig["layout"])))
        blocks.append(
            '<section class="card variant">'
            "<h3>%s</h3><p class=\"intro\">%s</p>%s</section>"
            % (title, blurb, plot_html))

    return """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Oma divergence — 10 chart variants</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root { --bg:#f6f7f9; --card:#fff; --text:#1a1a1a; --muted:#666; }
body { font-family: system-ui, sans-serif; background: var(--bg); color: var(--text);
  max-width: 960px; margin: 0 auto; padding: 16px; line-height: 1.45; }
header { margin-bottom: 20px; }
.card { background: var(--card); border-radius: 8px; padding: 14px 16px;
  margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,.08); }
.intro { color: var(--muted); font-size: .9rem; }
.plot { width: 100%%; min-height: 200px; }
.mini-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; }
@media (min-width: 720px) { .mini-grid { grid-template-columns: repeat(3, 1fr); } }
.mini-cell .plot { height: 220px; }
a { color: #1a5fb4; }
</style>
</head><body>
<header>
<h1>Divergence chart variants — pick one</h1>
<p>Data: oma assigned copies, divergence = 100 − ssearch36 %%identity, <b>1%% bins</b>,
Y = copy count (except variant 8). Built for comparison — tell me the number you want
on <a href="oma_report.html">oma_report.html</a>.</p>
<p class="intro">Best-practice refs: frequency polygon (connect bin midpoints);
Tal ccr/toc reports; avoid KDE when Y should be counts (DataRekha, OnlineStatBook).</p>
</header>
%s
<script>
%s
</script>
</body></html>
""" % ("".join(blocks), "".join(scripts))


def main():
    plots_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "oma_divergence_gallery.html"
    if not plots_dir or not plots_dir.is_dir():
        raise SystemExit("usage: build_divergence_gallery.py <plots_dir> [out.html]")
    by_sf = load_pctid_by_sf(plots_dir)
    if not by_sf:
        raise SystemExit("no pctid data in %s" % plots_dir)
    html = render_html(build_variants(by_sf))
    out.write_text(html, encoding="utf-8")
    print("wrote", out, "subfamilies", len(by_sf))


if __name__ == "__main__":
    main()
