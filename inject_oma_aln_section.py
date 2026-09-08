#!/usr/bin/env python3
"""Replace alignment section in oma report with Tal-style table + verdict columns."""
import html
import os
import re
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verdict import verdict
from flank_uniqueness import scan as flank_scan
from parse_oma_copy_counts import from_report

MSA = "https://toki-bio.github.io/MSA-viewer/"
RAW_OMA = "https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/oma/"
RAW_ALIGN = "https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/"

VCHIP_CSS = """
.vchip { display: inline-block; font-size: .78rem; font-weight: 600;
         padding: 2px 7px; border-radius: 3px; white-space: nowrap; }
.vchip.v-ok { background: #e2efe9; color: #1f6f5c; }
.vchip.v-edge { background: #f5eed8; color: #7a5c12; }
.vchip.v-warn { background: #f6e8de; color: #a8501d; }
.vchip.v-muted { background: #eceee8; color: #68766f; }
"""


def msa_href(path: str, title: str, base: str = RAW_OMA) -> str:
    url = base + path
    return f"{MSA}?url={quote(url, safe='')}&title={quote(title, safe='')}"


def raw_href(path: str, base: str = RAW_ALIGN) -> str:
    return base + path


def _codes(v):
    return {f["code"] for f in v.get("flags", [])}


def _chip(label, kind, title=""):
    t = f" title='{html.escape(title, quote=True)}'" if title else ""
    return f"<span class='vchip v-{kind}'{t}>{html.escape(label)}</span>"


def status_flanks(v):
    if not v or v.get("error"):
        return _chip("n/a", "muted", v.get("error", "not scored") if v else "")
    codes = _codes(v)
    if "NO_FLANKS_PRESENT" in codes:
        return _chip("No flanks", "muted", "Alignment carries no flank sequence.")
    if "NOT_ISOLATED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "NOT_ISOLATED")
        return _chip("Satellite / dup", "warn", f.get("text", ""))
    if "FRAGMENT_OF_LONGER" in codes or "ELEMENT_CONTINUES" in codes:
        f = next(x for x in v["flags"]
                 if x["code"] in ("FRAGMENT_OF_LONGER", "ELEMENT_CONTINUES"))
        return _chip("Fragment / LINE", "warn", f.get("text", ""))
    if "SHARED_FLANKS" in codes:
        f = next(x for x in v["flags"] if x["code"] == "SHARED_FLANKS")
        return _chip("Shared flanks", "warn", f.get("text", ""))
    if "FLANK_ISLANDS" in codes:
        f = next(x for x in v["flags"] if x["code"] == "FLANK_ISLANDS")
        return _chip("Flank islands", "edge", f.get("text", ""))
    if "FLANKS_UNMEASURED" in codes:
        return _chip("Short-flank view", "muted",
                      "50L/70R publish geometry; no 400 bp decay profile.")
    fb = v.get("flank_bg")
    if fb is None:
        return _chip("Not measured", "muted")
    if fb < 0.32:
        return _chip("Independent", "ok",
                      "Flank pairwise identity %.2f vs ~0.25 background." % fb)
    if fb < 0.42:
        return _chip("Borderline", "edge",
                      "Flank background %.2f — check by eye." % fb)
    return _chip("Raised flanks", "warn",
                  "Flank background %.2f — copies may share context." % fb)


def status_element(v):
    if not v or v.get("error"):
        return _chip("n/a", "muted", v.get("error", "not scored") if v else "")
    codes = _codes(v)
    g = v.get("groups", {}).get("element", 0.0)
    if "NO_ELEMENT" in codes or g < 0.25:
        f = next((x for x in v["flags"] if x["code"] == "NO_ELEMENT"), None)
        return _chip("No element", "warn", (f or {}).get("text", ""))
    if "MICROSATELLITE_ELEMENT" in codes:
        f = next(x for x in v["flags"] if x["code"] == "MICROSATELLITE_ELEMENT")
        return _chip("Microsatellite", "warn", f.get("text", ""))
    if "CONSENSUS_OVEREXTENDED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "CONSENSUS_OVEREXTENDED")
        return _chip("Over-extended", "edge", f.get("text", ""))
    if "CONSENSUS_UNDEREXTENDED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "CONSENSUS_UNDEREXTENDED")
        return _chip("Under-extended", "edge", f.get("text", ""))
    if "SMALL_CORE" in codes:
        f = next(x for x in v["flags"] if x["code"] == "SMALL_CORE")
        return _chip("Small core", "edge", f.get("text", ""))
    if g >= 0.85:
        return _chip("Strong", "ok",
                      "%d/%d copies support the consensus."
                      % (v.get("n_supported", 0), v.get("n", 0)))
    if g >= 0.5:
        return _chip("Supported", "ok",
                      "%d/%d copies support the consensus."
                      % (v.get("n_supported", 0), v.get("n", 0)))
    return _chip("Weak", "edge", "Element group score %.2f." % g)


def status_overall(v):
    if not v or v.get("error"):
        return _chip("n/a", "muted", v.get("error", "not scored") if v else "")
    if v.get("deferred"):
        return _chip("Deferred", "edge", "Mixture — split before a firm call.")
    if not v.get("assessable", True):
        return _chip("Cannot assess", "muted", "Too little evidence for a score.")
    s = float(v.get("score", 0))
    note = "%.0f/100" % s
    if v.get("capped_by"):
        note += " (capped: %s)" % ", ".join(v["capped_by"])
    if s >= 90:
        return _chip("SINE", "ok", note)
    if s >= 75:
        return _chip("SINE (caveats)", "ok", note)
    if s >= 55:
        return _chip("Grey zone", "edge", note)
    return _chip("Not SINE", "warn", note)


def _side_note(side):
    if not side.get("measured"):
        return side.get("reason", "not measured")
    return (
        "%.0f%% unique; largest shared group %.0f%% (%d copies)"
        % (100 * side.get("unique_frac", 0),
           100 * side.get("largest_cluster_frac", 0),
           side.get("largest_cluster", 0))
    )


def status_flank_context(top, rand):
    """Per-side flank clustering: rand100 high, top100 medium."""
    parts = []
    worst = None
    for label, r in (("rand", rand), ("top", top)):
        if not r or r.get("error"):
            continue
        sev = r.get("worst_flag")
        if sev == "high":
            worst = "high"
        elif sev == "medium" and worst != "high":
            worst = "medium"
        for f in r.get("flags", []):
            parts.append("%s %s: %s" % (label, f["side"], f["text"]))
    if not top and not rand:
        return _chip("n/a", "muted", "no alignment")
    if not parts:
        note = ""
        if top and top.get("left", {}).get("measured"):
            note = "Top100 L/R " + _side_note(top["left"]) + "; " + _side_note(top["right"])
        return _chip("Independent", "ok", note or "No shared-flank clusters above threshold.")
    title = " ".join(parts)
    if worst == "high":
        return _chip("Shared context", "warn", title)
    return _chip("Subgroup", "edge", title)


def score_top100(aln_dir: Path, sf: str):
    path = aln_dir / f"{sf}_top100.aln.fa"
    if not path.is_file():
        return None
    try:
        return verdict(str(path))
    except Exception as exc:
        return {"error": str(exc)}


def scan_flanks(aln_dir: Path, sf: str, tier: str):
    suffix = "top100" if tier == "top100" else "rand100"
    path = aln_dir / f"{sf}_{suffix}.aln.fa"
    if not path.is_file():
        return None
    try:
        return flank_scan(str(path))
    except Exception as exc:
        return {"error": str(exc)}


def build(subfams, aln_dir: Path, copy_counts: dict) -> str:
    intro_links = (
        "&nbsp;&nbsp;<a class='aln-link' href='"
        + msa_href("oma_consensuses.aln.fa", "oma all consensi", RAW_ALIGN)
        + "' target='_blank'>All consensi (MSA)</a>"
        "&nbsp;&nbsp;<a class='aln-link' href='"
        + raw_href("oma_consensuses_v2.fa")
        + "' target='_blank'>All consensi (FASTA)</a>"
        "&nbsp;&nbsp;<a class='aln-link green' href='"
        + msa_href("oma_subfam_input.aln.fa",
                   "oma SubFam input (30k loci sample)", RAW_ALIGN)
        + "' target='_blank'>SubFam input (step1, 30k loci)</a>"
    )
    rows = []
    for sf in sorted(subfams):
        t100 = f"{sf}_top100.aln.fa"
        r100 = f"{sf}_rand100.aln.fa"
        sub = f"{sf}_subfam.aln.fa"
        dash = "<span class='muted small'>&mdash;</span>"
        n_copies = copy_counts.get(sf)
        n_cell = (
            f"<td class='hl'>{n_copies:,}</td>"
            if n_copies is not None
            else f"<td class='muted small'>&mdash;</td>"
        )
        c_t = (
            f"<a class='aln-link' href='{msa_href(t100, f'oma {sf} top100')}' "
            f"target='_blank'>top 100</a>"
            if (aln_dir / t100).is_file()
            else dash
        )
        c_r = (
            f"<a class='aln-link orange' href='{msa_href(r100, f'oma {sf} rand100')}' "
            f"target='_blank'>100 random</a>"
            if (aln_dir / r100).is_file()
            else dash
        )
        c_s = (
            f"<a class='aln-link green' href='{msa_href(sub, f'oma {sf} subfam')}' "
            f"target='_blank'>SubFam</a>"
            if (aln_dir / sub).is_file()
            else dash
        )
        v = score_top100(aln_dir, sf)
        fu_t = scan_flanks(aln_dir, sf, "top100")
        fu_r = scan_flanks(aln_dir, sf, "rand100")
        rows.append(
            f"<tr><td><code>{html.escape(sf)}</code></td>"
            f"{n_cell}"
            f"<td>{c_t}</td><td>{c_r}</td><td>{c_s}</td>"
            f"<td>{status_flanks(v)}</td>"
            f"<td>{status_flank_context(fu_t, fu_r)}</td>"
            f"<td>{status_element(v)}</td>"
            f"<td>{status_overall(v)}</td></tr>"
        )
    return (
        "<section class='card' id='alignments'>"
        "<h2>Subfamily Alignments &mdash; open in MSA Viewer</h2>"
        "<p class='intro'>Full publish path: border loop, step7 flank widen, "
        "step8a re-extract + MAFFT, <code>rebuild_consensus_row.py</code> "
        "(copy-column majority + edge trim), then "
        "<code>boundary_justify.py</code>. Verdict columns score each "
        "family&rsquo;s <strong>top&nbsp;100</strong> alignment with "
        "<code>verdict.py</code> on the 50L/70R publish geometry "
        "(no 400&nbsp;bp decay sidecar for oma yet). Consensus is row&nbsp;1."
        + intro_links
        + "</p>"
        "<table class='tbl'><thead><tr><th>Subfamily</th>"
        "<th title='Genome-wide assigned copies (firm + soft) from summary.by_subfam.tsv'>"
        "Copies</th>"
        "<th>Top 100 by bitscore</th><th>100 random copies</th>"
        "<th>SubFam (chunk consensuses)</th>"
        "<th>Flanks</th><th>Flank context</th><th>Element</th><th>Overall</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def ensure_vchip_css(text: str) -> str:
    if ".vchip" in text:
        return text
    return text.replace("</style>", VCHIP_CSS + "</style>", 1)


def main():
    report = Path(sys.argv[1])
    aln_dir = Path(sys.argv[2])
    subfams = sorted(
        p.name.replace("_top100.aln.fa", "")
        for p in aln_dir.glob("oma_*_top100.aln.fa")
        if not p.name.startswith("oma_oma_")
    )
    text = report.read_text(encoding="utf-8")
    text = ensure_vchip_css(text)
    copy_counts = from_report(text)
    new_sec = build(subfams, aln_dir, copy_counts)
    text, n = re.subn(
        r"<section class='card' id='alignments'>.*?</section>",
        new_sec,
        text,
        count=1,
        flags=re.DOTALL,
    )
    if n != 1:
        sys.exit("alignment section not found")
    report.write_text(text, encoding="utf-8")
    print("patched %d subfamilies" % len(subfams))


if __name__ == "__main__":
    main()
