#!/usr/bin/env python3
"""Fix embedded pctid KDE + violins in oma_report.html without DRAGEN."""
import json
import re
from pathlib import Path

REPORT = Path(__file__).parent / "oma_report.html"

PLOTLY_RE = re.compile(
    r"Plotly\.newPlot\('(?P<id>[^']+)',\s*(?P<data>\[.*?\]),\s*(?P<layout>\{.*?\})\);",
    re.DOTALL,
)


def parse_plotly_call(text: str, plot_id: str):
    marker = "Plotly.newPlot('%s'," % plot_id
    start = text.find(marker)
    if start < 0:
        raise SystemExit("missing %s" % plot_id)
    i = start + len(marker)
    depth = 0
    data_start = i
    in_str = False
    esc = False
    for j in range(i, len(text)):
        c = text[j]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
            continue
        if c == '"':
            in_str = True
            continue
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                data_end = j + 1
                break
    else:
        raise SystemExit("data end not found for %s" % plot_id)
    k = data_end
    while k < len(text) and text[k] in " \t\n\r,":
        k += 1
    if text[k] != "{":
        raise SystemExit("layout start not found for %s" % plot_id)
    depth = 0
    in_str = False
    esc = False
    layout_start = k
    for j in range(k, len(text)):
        c = text[j]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
            continue
        if c == '"':
            in_str = True
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                layout_end = j + 1
                break
    else:
        raise SystemExit("layout end not found for %s" % plot_id)
    end = text.find(";", layout_end)
    data = json.loads(text[data_start:data_end])
    layout = json.loads(text[layout_start:layout_end])
    return start, end + 1, data, layout


def fix_kde(data, layout):
    for tr in data:
        if tr.get("type") == "scatter" and tr.get("x"):
            tr["x"] = [round(max(0.0, float(v)), 3) for v in tr["x"]]
    xa = layout.setdefault("xaxis", {})
    xa["rangemode"] = "nonnegative"
    hi = 0.0
    for tr in data:
        for v in tr.get("x") or []:
            hi = max(hi, float(v))
    xa["range"] = [0, round(hi * 1.06 + 0.5, 2)]
    layout["uirevision"] = "oma-pctid-kde"
    return data, layout


def fix_violins(data, layout):
    clean = []
    for tr in data:
        if tr.get("type") != "violin":
            continue
        clean.append({
            "type": "violin",
            "y": tr.get("y") or [],
            "name": tr.get("name", ""),
            "box": {"visible": True},
            "meanline": {"visible": True},
            "points": False,
        })
    layout = {
        "title": layout.get(
            "title",
            "ssearch36 %identity divergence per subfamily (step4)"),
        "yaxis": {
            "title": "Divergence (100 − %identity to consensus)",
            "rangemode": "nonnegative",
        },
        "xaxis": {
            "title": "Subfamily",
            "tickangle": -45,
        },
        "height": layout.get("height", 520),
        "showlegend": False,
        "margin": {"t": 60, "r": 20, "b": 140, "l": 70},
        "uirevision": "oma-pctid-violins",
    }
    return clean, layout


def render(plot_id, data, layout):
    return "Plotly.newPlot('%s', %s, %s);" % (
        plot_id,
        json.dumps(data, separators=(",", ":")),
        json.dumps(layout, separators=(",", ":")),
    )


def main():
    text = REPORT.read_text(encoding="utf-8")
    for plot_id, fixer in (
        ("plot_pctid_kde", fix_kde),
        ("plot_sim_violins", fix_violins),
    ):
        start, end, data, layout = parse_plotly_call(text, plot_id)
        data, layout = fixer(data, layout)
        text = text[:start] + render(plot_id, data, layout) + text[end:]
        print("fixed", plot_id, "traces", len(data))
    REPORT.write_text(text, encoding="utf-8")
    print("wrote", REPORT)


if __name__ == "__main__":
    main()
