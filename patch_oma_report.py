#!/usr/bin/env python3
"""Patch oma_report.html: header, dual divergence plots, gallery, plot JS."""
import base64
import json
import math
import random
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

SF_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
]

SIM_FLOOR_PCT = 45.0  # similarity floor for filtered bitscore KDE


def kde_curve(vals: List[float], n_pts: int = 300) -> Tuple[List[float], List[float]]:
    vals = [float(v) for v in vals if v == v]
    n = len(vals)
    if n < 2:
        return [], []
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / (n - 1)
    std = math.sqrt(max(var, 1e-12))
    h = std * n ** (-0.2)
    lo = min(vals) - 2.5 * h
    hi = max(vals) + 2.5 * h
    step = (hi - lo) / (n_pts - 1)
    x_arr = [lo + i * step for i in range(n_pts)]
    c = 1.0 / (n * h * math.sqrt(2 * math.pi))
    y_arr = [
        c * sum(math.exp(-0.5 * ((xi - v) / h) ** 2) for v in vals)
        for xi in x_arr
    ]
    return x_arr, y_arr


def stratified_sample_sim(sim_path: Path, assign_path: Path,
                          per_group: int = 3000,
                          min_sim: float = None) -> Dict[str, List[float]]:
    seq_to_sf: Dict[str, str] = {}
    with assign_path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[4] == "assigned":
                seq_to_sf[parts[0]] = parts[1]
    rng = random.Random(42)
    by_sf: Dict[str, List[float]] = {}
    counts: Dict[str, int] = {}
    with sim_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            sid = parts[0]
            sf = seq_to_sf.get(sid)
            if not sf:
                continue
            try:
                v = float(parts[3]) * 100.0
            except ValueError:
                continue
            if min_sim is not None and v < min_sim:
                continue
            counts[sf] = counts.get(sf, 0) + 1
            buf = by_sf.setdefault(sf, [])
            if len(buf) < per_group:
                buf.append(v)
            else:
                j = rng.randrange(counts[sf])
                if j < per_group:
                    buf[j] = v
    return by_sf


def load_pctid_by_sf(plots_dir: Path) -> Dict[str, List[float]]:
    by_sf: Dict[str, List[float]] = {}
    for p in sorted(plots_dir.glob("*_pctid.tsv")):
        sf = p.name.replace("_pctid.tsv", "")
        vals = []
        with p.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        vals.append(float(parts[1]))
                    except ValueError:
                        pass
        if vals:
            by_sf[sf] = vals
    return by_sf


def fig_divergence_kde(by_sf: Dict[str, List[float]], title: str,
                       xaxis_title: str = None) -> dict:
    sf_sorted = sorted(by_sf.keys())
    traces = []
    for i, sf in enumerate(sf_sorted):
        vals = [max(0.0, 100.0 - v) for v in by_sf[sf]]
        if not vals:
            continue
        x, y = kde_curve(vals)
        if not x:
            continue
        traces.append({
            "type": "scatter", "mode": "lines",
            "x": [round(max(0.0, v), 3) for v in x],
            "y": [round(v, 6) for v in y],
            "name": sf,
            "line": {"color": SF_PALETTE[i % len(SF_PALETTE)], "width": 2},
            "hovertemplate": (
                "%{fullData.name}<br>divergence %{x:.1f}%"
                "<br>density %{y:.5f}<extra></extra>"),
        })
    layout = {
        "title": title,
        "xaxis": {"title": xaxis_title or (
            "Divergence (100 − bitscore/self-bits × 100%)"),
                  "rangemode": "nonnegative"},
        "yaxis": {"title": "Density"},
        "legend": {"title": {"text": "Subfamily (click to toggle)"}},
        "height": 460,
        "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
    }
    return {"data": traces, "layout": layout}


def bin_divergence(vals: List[float], bin_width: float = 1.0) -> Dict[float, int]:
    """Divergence bins: 100 − pctid, 1% bins (Tal step4 / ccr report style)."""
    counts: Dict[float, int] = {}
    for v in vals:
        d = max(0.0, 100.0 - float(v))
        bin_start = math.floor((d + 1e-9) / bin_width) * bin_width
        bin_start = round(bin_start, 6)
        counts[bin_start] = counts.get(bin_start, 0) + 1
    return counts


def fig_pctid_divergence(by_sf: Dict[str, List[float]]) -> dict:
    """Binned copy counts — same metric as step4 gallery PNGs, not KDE."""
    sf_sorted = sorted(by_sf.keys())
    traces = []
    max_bin_end = 0.0
    bin_width = 1.0
    for i, sf in enumerate(sf_sorted):
        counts = bin_divergence(by_sf[sf], bin_width)
        if not counts:
            continue
        bins = sorted(counts)
        max_bin_end = max(max_bin_end, max(bins) + bin_width)
        traces.append({
            "type": "scatter",
            "mode": "lines",
            "x": [round(b + bin_width / 2.0, 3) for b in bins],
            "y": [counts[b] for b in bins],
            "name": sf,
            "line": {"color": SF_PALETTE[i % len(SF_PALETTE)], "width": 2},
            "hovertemplate": (
                "%{fullData.name}<br>divergence ~%{x:.0f}%"
                "<br>copies %{y:,d}<extra></extra>"),
        })
    x_range_max = min(100.0, max(5.0, math.ceil(max_bin_end / 5.0) * 5.0))
    return {
        "data": traces,
        "layout": {
            "title": "ssearch36 %identity divergence (step4 — Tal gallery metric)",
            "xaxis": {"title": "Divergence (100 − %identity to consensus)",
                      "rangemode": "nonnegative", "range": [0, x_range_max]},
            "yaxis": {"title": "Copies"},
            "legend": {"title": {"text": "Subfamily (click to toggle)"}},
            "height": 460,
            "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
            "uirevision": "oma-pctid-hist",
        },
    }


def fig_pctid_violins(by_sf: Dict[str, List[float]]) -> dict:
    """Per-subfamily violin — same layout as step6_report.py."""
    sf_sorted = sorted(by_sf.keys())
    traces = []
    for sf in sf_sorted:
        vals = [round(max(0.0, 100.0 - v), 2) for v in by_sf[sf]]
        if not vals:
            continue
        hi = max(vals)
        traces.append({
            "type": "violin",
            "y": vals,
            "name": sf,
            "box": {"visible": True},
            "meanline": {"visible": True},
            "points": False,
            "spanmode": "hard",
            "span": [0, hi + 1],
        })
    return {
        "data": traces,
        "layout": {
            "title": "ssearch36 %identity divergence per subfamily (step4)",
            "yaxis": {
                "title": "Divergence (100 − %identity to consensus)",
                "rangemode": "nonnegative",
                "range": [0, None],
            },
            "xaxis": {
                "title": "Subfamily",
                "tickangle": -45,
            },
            "height": 520,
            "showlegend": False,
            "margin": {"t": 60, "r": 20, "b": 140, "l": 70},
            "uirevision": "oma-pctid-violins",
        },
    }


def plotly_js(div_id: str, fig: dict) -> str:
    return "Plotly.newPlot(%r, %s, %s);\n" % (
        div_id,
        json.dumps(fig["data"]),
        json.dumps(fig["layout"]),
    )


def patch_header(text: str) -> str:
    old = re.search(
        r'(<header>\s*<h1>.*?</h1>\s*)<div class="sub">.*?</div>'
        r'(?:\s*<p class="small".*?</p>)?(\s*</header>)',
        text, re.DOTALL)
    if not old:
        return text
    new_sub = (
        '<div class="sub">Genome: <i>Olivierus martensii</i> (Chinese scorpion) '
        '&middot; assembly <code>GCA_000484575.1</code> '
        '(<i>M_martensii_Version_1</i>, ~901&nbsp;Mb WGS, AYEL010 contigs) '
        '&middot; Query consensuses: <code>oma_seeds.fa</code> '
        '(23 AnnoSINE seeds, 26 curated families) &middot; '
        'Generated 2026-09-07 12:17 UTC</div>'
    )
    return text[:old.start(1)] + old.group(1) + new_sub + old.group(2) + text[old.end(2):]


def patch_overview_genome(text: str) -> str:
    return text.replace(
        "searched <b>this genome</b>",
        "searched <b><i>Olivierus martensii</i> "
        "(GCA_000484575.1, ~901&nbsp;Mb)</b>",
        1,
    )


GALLERY_IMG_BASE = (
    "https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/"
    "alignments/oma/plots"
)


def build_gallery_section(plots_dir: Path, copy_counts: Dict[str, int],
                          min_copies: int = 50,
                          img_base_url: str = None) -> str:
    if not plots_dir.is_dir() and not img_base_url:
        return ""
    blocks = []
    png_names = sorted(p.name for p in plots_dir.glob("*_divergence.png")) if plots_dir.is_dir() else []
    if img_base_url and not png_names:
        # allow building from known sf list when PNGs live on GitHub only
        png_names = [f"{sf}_divergence.png" for sf in sorted(copy_counts)
                     if copy_counts.get(sf, 0) >= min_copies]
    for png_name in png_names:
        sf = png_name.replace("_divergence.png", "")
        n = copy_counts.get(sf, 0)
        if n and n < min_copies:
            continue
        nuc_name = f"{sf}_nucfreq.png"
        div_block = f"<div class='subfam-block'><h3>{sf}</h3>"
        if img_base_url:
            div_src = f"{img_base_url}/{png_name}"
            nuc_src = f"{img_base_url}/{nuc_name}"
            div_block += (
                f"<img src='{div_src}' alt='{sf} divergence histogram' "
                f"onclick=\"document.getElementById('lightbox-img').src=this.src;"
                f"document.getElementById('lightbox').classList.add('active')\">"
            )
            if plots_dir.is_dir() and (plots_dir / nuc_name).is_file():
                div_block += (
                    f"<img src='{nuc_src}' alt='{sf} nucfreq' "
                    f"onclick=\"document.getElementById('lightbox-img').src=this.src;"
                    f"document.getElementById('lightbox').classList.add('active')\">"
                )
        else:
            png = plots_dir / png_name
            if not png.is_file():
                continue
            nuc = plots_dir / nuc_name
            data = base64.b64encode(png.read_bytes()).decode("ascii")
            div_block += (
                f"<img src='data:image/png;base64,{data}' "
                f"alt='{sf} divergence histogram' "
                f"onclick=\"document.getElementById('lightbox-img').src=this.src;"
                f"document.getElementById('lightbox').classList.add('active')\">"
            )
            if nuc.is_file():
                data2 = base64.b64encode(nuc.read_bytes()).decode("ascii")
                div_block += (
                    f"<img src='data:image/png;base64,{data2}' "
                    f"alt='{sf} nucfreq' "
                    f"onclick=\"document.getElementById('lightbox-img').src=this.src;"
                    f"document.getElementById('lightbox').classList.add('active')\">"
                )
        div_block += "</div>"
        blocks.append(div_block)
    if not blocks:
        return ""
    return (
        "<section class='card' id='gallery'>"
        "<h2>Per-subfamily diagnostic plots (step4)</h2>"
        "<p class='intro'>Binned histogram: divergence = 100 &minus; ssearch36 "
        "%identity to the subfamily consensus. Stacked bars: per-position "
        "nucleotide frequencies in the copy alignment. "
        f"Families with &lt;{min_copies} assigned copies omitted.</p>"
        "<div class='subfam-grid'>" + "".join(blocks) + "</div></section>"
    )


def load_copy_counts_from_summary(summary: Path) -> Dict[str, int]:
    counts = {}
    if not summary.is_file():
        return counts
    with summary.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            p = line.split("\t")
            if len(p) >= 4:
                try:
                    counts[p[0]] = int(p[3])
                except ValueError:
                    pass
    return counts


def export_fragments(data_dir: Path, plots_dir: Path,
                     min_copies: int = 50) -> dict:
    html = ""
    js = ""
    gallery = ""
    pctid: Dict[str, List[float]] = {}
    if data_dir and data_dir.is_dir():
        pctid_dir = plots_dir if plots_dir and plots_dir.is_dir() else data_dir / "plots"
        pctid = load_pctid_by_sf(pctid_dir)
        if pctid:
            js = plotly_js("plot_pctid_kde", fig_pctid_divergence(pctid))
    if plots_dir and plots_dir.is_dir():
        counts = load_copy_counts_from_summary(
            data_dir / "summary.by_subfam.tsv" if data_dir else Path())
        gallery = build_gallery_section(
            plots_dir, counts, min_copies=min_copies,
            img_base_url=GALLERY_IMG_BASE)
    return {"html": html, "js": js, "gallery": gallery,
            "violins_js": plotly_js("plot_sim_violins", fig_pctid_violins(pctid))
            if pctid else ""}


def patch_divergence_section(text: str, data_dir: Path) -> Tuple[str, bool]:
    sim = data_dir / "sim_scores.tsv"
    assign = data_dir / "assignment_full.tsv"
    plots = data_dir / "plots"
    if not sim.is_file() or not assign.is_file():
        return text, False

    all_sim = stratified_sample_sim(sim, assign)
    filt_sim = stratified_sample_sim(sim, assign, min_sim=SIM_FLOOR_PCT)
    pctid = load_pctid_by_sf(plots) if plots.is_dir() else {}

    extra_html = (
        '<h3>Bitscore divergence — copies above similarity floor '
        f'({SIM_FLOOR_PCT:.0f}% of self-bitscore)</h3>'
        '<p class="intro">Same bitscore-based metric as above, but copies with '
        'similarity below the step2 floor are excluded. Values below ~45% often '
        'reflect marginal assignments, not real subfamily divergence.</p>'
        '<div class="plot" id="plot_div_kde_filtered"></div>'
    )
    extra_js = plotly_js("plot_div_kde_filtered", fig_divergence_kde(
        filt_sim,
        "Bitscore divergence — filtered (sim ≥ %.0f%%)" % SIM_FLOOR_PCT,
    ))

    if pctid:
        extra_html += (
            '<h3>ssearch36 %identity divergence (step4)</h3>'
            '<p class="intro">Nucleotide %identity from ssearch36 against the '
            'subfamily consensus — the same metric as Tal&rsquo;s per-subfamily '
            'binned histograms in the Gallery. Compare with the bitscore plots above '
            'and pick which you trust.</p>'
            '<div class="plot" id="plot_pctid_kde"></div>'
        )
        extra_js += plotly_js("plot_pctid_kde", fig_pctid_divergence(pctid))

    # Insert after first KDE plot div
    text, n = re.subn(
        r'(<div class="plot" id="plot_div_kde"></div>)',
        r'\1' + extra_html,
        text, count=1)
    if n != 1:
        return text, False

    # Append Plotly calls before closing script tag (after plot_div_kde line)
    text, n2 = re.subn(
        r"(Plotly\.newPlot\('plot_div_kde',",
        extra_js + r"\1",
        text, count=1)
    return text, n2 == 1


def embed_gallery(text: str, plots_dir: Path, copy_counts: Dict[str, int],
                  min_copies: int = 50) -> Tuple[str, int]:
    gallery = build_gallery_section(
        plots_dir, copy_counts, min_copies, img_base_url=GALLERY_IMG_BASE)
    if not gallery:
        return text, 0
    n_blocks = gallery.count("subfam-block")
    if 'id="gallery"' in text or "id='gallery'" in text:
        return text, n_blocks
    text, n = re.subn(
        r'(<section class="card" id="pca">)',
        gallery + r"\1",
        text, count=1)
    if 'href="#gallery"' not in text:
        text = text.replace(
            '<a href="#pca">PCA</a>',
            '<a href="#gallery">Gallery</a>\n    <a href="#pca">PCA</a>',
            1)
    return text, n_blocks if n else 0


def parse_copy_counts(text: str) -> Dict[str, int]:
    counts = {}
    for m in re.finditer(
            r"<tr><td>(oma_[^<]+)</td><td>\d+</td><td>\d+</td>"
            r"<td class=\"hl\">(\d+)</td>", text):
        counts[m.group(1)] = int(m.group(2))
    return counts


def ensure_plot_js(text: str) -> str:
    if "oma_report_plots.js" in text:
        return text
    text = text.replace(
        "</body></html>",
        '<script src="oma_report_plots.js"></script>\n</body></html>',
    )
    return text


def merge_fragments(report: Path, frag_dir: Path) -> None:
    text = report.read_text(encoding="utf-8")
    div_html = (frag_dir / "divergence_extra.html").read_text(encoding="utf-8")
    div_js = (frag_dir / "divergence_extra.js").read_text(encoding="utf-8")
    gallery = (frag_dir / "gallery.html").read_text(encoding="utf-8")
    if div_html and 'plot_div_kde_filtered' not in text:
        text, n = re.subn(
            r'(<div class="plot" id="plot_div_kde"></div>)',
            r'\1' + div_html, text, count=1)
        marker = "Plotly.newPlot('plot_div_kde',"
        if n and div_js and marker in text:
            text = text.replace(marker, div_js + marker, 1)
    elif div_html and 'plot_pctid_kde' not in text and 'plot_pctid_kde' in div_html:
        m = re.search(
            r'(<h3>ssearch36 %identity divergence \(step4\)</h3>.*'
            r'<div class="plot" id="plot_pctid_kde"></div>)',
            div_html, re.DOTALL)
        if m:
            text = text.replace(
                '<div class="plot" id="plot_div_kde_filtered"></div>',
                '<div class="plot" id="plot_div_kde_filtered"></div>' + m.group(1),
                1)
        pct_marker = "Plotly.newPlot('plot_pctid_kde',"
        if pct_marker in div_js and pct_marker not in text:
            end = div_js.find(";\n", div_js.index(pct_marker))
            if end == -1:
                end = div_js.find(";", div_js.index(pct_marker))
            pct_js = div_js[div_js.index(pct_marker):end + 1] + "\n"
            kde_marker = "Plotly.newPlot('plot_div_kde',"
            if kde_marker in text:
                text = text.replace(kde_marker, pct_js + kde_marker, 1)
    if gallery and "id='gallery'" not in text and 'id="gallery"' not in text:
        text, _ = re.subn(
            r'(<section class="card" id="pca">)',
            gallery + r"\1", text, count=1)
        if 'href="#gallery"' not in text:
            text = text.replace(
                '<a href="#pca">PCA</a>',
                '<a href="#gallery">Gallery</a>\n    <a href="#pca">PCA</a>',
                1)
    text = ensure_plot_js(text)
    report.write_text(text, encoding="utf-8")


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--merge":
        merge_fragments(Path(sys.argv[2]), Path(sys.argv[3]))
        print("merged fragments into", sys.argv[2])
        return

    report = Path(sys.argv[1])
    data_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    plots_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else None
    out_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else None

    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        fragments = export_fragments(data_dir, plots_dir)
        (out_dir / "divergence_extra.html").write_text(
            fragments["html"], encoding="utf-8")
        (out_dir / "divergence_extra.js").write_text(
            fragments["js"], encoding="utf-8")
        (out_dir / "gallery.html").write_text(
            fragments.get("gallery", ""), encoding="utf-8")
        if fragments.get("violins_js"):
            (out_dir / "violins.js").write_text(
                fragments["violins_js"], encoding="utf-8")
        print("wrote fragments to", out_dir)
        return

    text = report.read_text(encoding="utf-8")
    text = patch_header(text)
    text = patch_overview_genome(text)
    text = ensure_plot_js(text)
    if data_dir and data_dir.is_dir():
        text, ok = patch_divergence_section(text, data_dir)
        print("divergence patch:", "ok" if ok else "skipped")
    counts = parse_copy_counts(text)
    if plots_dir and plots_dir.is_dir():
        text, n = embed_gallery(text, plots_dir, counts, min_copies=50)
        print("gallery blocks:", n)
    report.write_text(text, encoding="utf-8")
    print("patched", report)


if __name__ == "__main__":
    main()
