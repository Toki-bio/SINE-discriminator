#!/usr/bin/env python3
"""Inject SINE-discriminator verdict columns into a SINEderella report."""
import argparse
import html
import os
import re
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verdict import verdict  # noqa: E402
from flank_uniqueness import scan as flank_scan  # noqa: E402
from overall import overall_call  # noqa: E402

MSA = "https://toki-bio.github.io/MSA-viewer/"

VCHIP_CSS = """
.vchip { display: inline-block; font-size: .78rem; font-weight: 600;
         padding: 2px 7px; border-radius: 3px; white-space: nowrap; }
.vchip.v-ok { background: #e2efe9; color: #1f6f5c; }
.vchip.v-edge { background: #f5eed8; color: #7a5c12; }
.vchip.v-warn { background: #f6e8de; color: #a8501d; }
.vchip.v-muted { background: #eceee8; color: #68766f; }
"""

# One styled tooltip for the whole page. Native title= tooltips appear only after a delay, vanish
# after a few seconds and never on touch; this moves every title= into data-tip on load and shows it
# at once on hover (the texts themselves are unchanged, so every generator's titles keep working).
TIP_MARK = "sd-tip-js"
TIP_HTML = """<style id="sd-tip-css">
#sd-tip { position: fixed; z-index: 9999; max-width: 460px; background: #1f2933; color: #fff;
          font-size: 12.5px; line-height: 1.45; padding: 7px 10px; border-radius: 5px;
          box-shadow: 0 4px 14px rgba(0,0,0,.25); pointer-events: none; white-space: pre-line;
          display: none; }
[data-tip] { cursor: help; }
a[data-tip] { cursor: pointer; }
</style>
<script id="sd-tip-js">
(function () {
  function init() {
    var tip = document.createElement('div'); tip.id = 'sd-tip'; document.body.appendChild(tip);
    document.querySelectorAll('[title]').forEach(function (el) {
      if (el.closest('svg')) return;
      el.setAttribute('data-tip', el.getAttribute('title')); el.removeAttribute('title');
    });
    function place(e) {
      var x = e.clientX + 14, y = e.clientY + 16, w = tip.offsetWidth, h = tip.offsetHeight;
      if (x + w > window.innerWidth - 8) x = Math.max(8, window.innerWidth - w - 8);
      if (y + h > window.innerHeight - 8) y = Math.max(8, e.clientY - h - 10);
      tip.style.left = x + 'px'; tip.style.top = y + 'px';
    }
    document.addEventListener('mouseover', function (e) {
      var el = e.target.closest && e.target.closest('[data-tip]');
      if (!el) { tip.style.display = 'none'; return; }
      tip.textContent = el.getAttribute('data-tip'); tip.style.display = 'block'; place(e);
    });
    document.addEventListener('mousemove', function (e) { if (tip.style.display === 'block') place(e); });
    document.addEventListener('scroll', function () { tip.style.display = 'none'; }, true);
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
</script>
"""

TOP_RULE = ("Top 100 by bitscore: the firmly assigned copies (10/10 votes, above the bitscore bar) with "
            "the highest bitscore to the consensus. If fewer than 100 are firmly assigned, the plate is "
            "filled up to 100 with soft-assigned copies (found by this query but not unanimous), ranked "
            "by search score and marked [soft] in the row name.")
RAND_RULE = ("100 random copies of the firmly assigned set. If fewer than 100 are firmly assigned, all of "
             "them are used and the plate is filled up to 100 with random soft-assigned copies, marked "
             "[soft] in the row name.")
ROWS_NOTE = ("Row 1 is the consensus rebuilt from these copies; row 2 (<subfamily>_seed_as_searched) is the "
             "consensus the genome was searched with.")


def plate_counts(path: Path):
    """(copies, soft) on a published plate: rows other than the consensus and the seed row."""
    n = soft = 0
    first = True
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        if first:
            first = False
            continue
        if "_seed_as_searched" in line:
            continue
        n += 1
        soft += "[soft]" in line
    return n, soft

# Re-use chip helpers from inject_oma_aln_section (same logic)
from inject_oma_aln_section import (  # noqa: E402
    status_element,
    status_flank_context,
    status_flanks,
    _chip,
)


def msa_href(raw_base: str, filename: str, title: str) -> str:
    url = raw_base.rstrip("/") + "/" + filename
    return (
        f"{MSA}?url={quote(url, safe='')}"
        f"&title={quote(title, safe='')}"
    )


def aln_name(species: str, sf: str, kind: str) -> str:
    """step8a naming: {species}_{subfam}_{kind}.aln.fa"""
    prefix = f"{species}_" if not sf.startswith(f"{species}_") else ""
    return f"{prefix}{sf}_{kind}.aln.fa"


def load_copy_counts(summary_tsv: Path) -> dict:
    counts = {}
    if not summary_tsv.is_file():
        return counts
    for line in summary_tsv.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sf = parts[0].strip()
        try:
            counts[sf] = int(parts[1]) + int(parts[2])
        except ValueError:
            try:
                counts[sf] = int(parts[3])
            except (ValueError, IndexError):
                pass
    return counts


def discover_subfams(aln_dir: Path, species: str) -> list:
    out = set()
    for p in aln_dir.glob("*_top100.aln.fa"):
        name = p.name.replace("_top100.aln.fa", "")
        if name.startswith(f"{species}_"):
            name = name[len(species) + 1:]
        out.add(name)
    return sorted(out)


def score_top100(aln_dir: Path, species: str, sf: str):
    for fn in (aln_name(species, sf, "top100"), f"{sf}_top100.aln.fa"):
        path = aln_dir / fn
        if path.is_file():
            try:
                return verdict(str(path))
            except Exception as exc:
                return {"error": str(exc)}
    return None


def scan_flanks(aln_dir: Path, species: str, sf: str, tier: str):
    kind = "top100" if tier == "top100" else "rand100"
    for fn in (aln_name(species, sf, kind), f"{sf}_{kind}.aln.fa"):
        path = aln_dir / fn
        if path.is_file():
            try:
                return flank_scan(str(path))
            except Exception as exc:
                return {"error": str(exc)}
    return None


def build_section(species: str, subfams, aln_dir: Path,
                  copy_counts: dict, raw_base: str) -> str:
    dash = "<span class='muted small'>&mdash;</span>"
    rows = []
    for sf in sorted(subfams):
        t100 = aln_name(species, sf, "top100")
        r100 = aln_name(species, sf, "rand100")
        sub = aln_name(species, sf, "subfam")
        n = copy_counts.get(sf) or copy_counts.get(f"{species}_{sf}")
        n_cell = (
            f"<td class='hl'>{n:,}</td>" if n is not None
            else "<td class='muted small'>&mdash;</td>"
        )
        def link(fn, label, css="", rule=""):
            path = aln_dir / fn
            if not path.is_file():
                return dash
            tip = ""
            if rule:
                n_c, n_soft = plate_counts(path)
                if n_c < 100:
                    label = label.replace("100", str(n_c))
                if n_soft:
                    label += f" ({n_soft} soft)"
                tip = (f"{rule}\n\nThis plate: {n_c - n_soft} firm + {n_soft} soft = {n_c} copies."
                       + ("" if n_c >= 100 else " The run has no more copies of this subfamily.")
                       + "\n" + ROWS_NOTE)
            href = msa_href(raw_base, fn, f"{species} {sf} {label}") if raw_base else fn
            cls = "aln-link" + (f" {css}" if css else "")
            t = f" title='{html.escape(tip, quote=True)}'" if tip else ""
            return f"<a class='{cls}' href='{href}' target='_blank'{t}>{label}</a>"
        v = score_top100(aln_dir, species, sf)
        top_scan = scan_flanks(aln_dir, species, sf, 'top100')
        rand_scan = scan_flanks(aln_dir, species, sf, 'rand100')
        oc = overall_call(v, top_scan, rand_scan)
        rows.append(
            f"<tr><td><code>{html.escape(sf)}</code></td>{n_cell}"
            f"<td>{link(t100, 'top 100', '', TOP_RULE)}</td>"
            f"<td>{link(r100, '100 random', 'orange', RAND_RULE)}</td>"
            f"<td>{link(sub, 'SubFam', 'green')}</td>"
            f"<td>{status_flanks(v)}</td>"
            f"<td>{status_flank_context(top_scan, rand_scan)}</td>"
            f"<td>{status_element(v)}</td>"
            f"<td>{_chip(oc['label'], oc['kind'], chr(10).join(oc['reasons']))}</td></tr>"
        )
    intro = (
        "<p class='intro'>One row per subfamily. "
        "<strong>Links</strong> open 100-copy MAFFT alignments (50&nbsp;bp left, "
        "70&nbsp;bp right flanks) in the MSA viewer: the <strong>consensus rebuilt from the copies is "
        "row&nbsp;1</strong>, the consensus as searched is row&nbsp;2. A subfamily with fewer than 100 "
        "firmly assigned copies is filled up with soft-assigned ones, marked <code>[soft]</code>; the "
        "button then shows the count (hover for the rule). "
        "Last four columns: automated <strong>SINE-discriminator</strong> checks "
        "on the top-100 view (<strong>Flanks</strong>, <strong>Flank context</strong>, "
        "<strong>Element</strong>, <strong>Overall</strong>).</p>"
        "<details class='legend'><summary>How these alignments were built</summary>"
        "<ul class='small muted'>"
        "<li><code>publish/align_for_publish.sh</code>: step7 → border loop → step8a.</li>"
        "<li><code>rebuild_consensus_row.py</code>, <code>boundary_justify.py</code>, "
        "<code>trim_display_flanks.py</code> (SINE-discriminator).</li>"
        "<li>Verdict chips: <code>verdict.py</code> on 50L/70R publish geometry.</li>"
        "</ul></details>"
    )
    thead = (
        "<th>Subfamily</th><th>Copies</th>"
        f"<th title='{html.escape(TOP_RULE, quote=True)}'>Top 100 by bitscore</th>"
        f"<th title='{html.escape(RAND_RULE, quote=True)}'>100 random copies</th>"
        "<th>SubFam (chunk consensuses)</th>"
        "<th>Flanks</th><th>Flank context</th><th>Element</th>"
        "<th title='Verdict score on the top-100 plate. Below 30 copies no score is given: the call is "
        "Doubtful if the copies show negative evidence, Cannot assess if they show none. Hover a chip "
        "for its reasons.'>Overall</th>"
    )
    return (
        "<section class='card' id='alignments'>"
        "<h2>Subfamily alignments</h2>"
        + intro
        + "<table class='tbl'><thead><tr>" + thead + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def ensure_vchip_css(text: str) -> str:
    if ".vchip" in text:
        return text
    return text.replace("</style>", VCHIP_CSS + "</style>", 1)


def ensure_tooltip(text: str) -> str:
    if TIP_MARK in text:
        return text
    i = text.rfind("</body>")
    return text + TIP_HTML if i == -1 else text[:i] + TIP_HTML + text[i:]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("report", type=Path)
    ap.add_argument("aln_dir", type=Path)
    ap.add_argument("species", help="Species code prefix (ccr, oma, …)")
    ap.add_argument("--raw-base", default="",
                    help="Published raw URL prefix for MSA viewer links")
    ap.add_argument("--summary-tsv", type=Path, default=None,
                    help="summary.by_subfam.tsv for copy counts")
    args = ap.parse_args()

    subfams = discover_subfams(args.aln_dir, args.species)
    if not subfams:
        sys.exit("no *_top100.aln.fa in %s" % args.aln_dir)

    counts = load_copy_counts(args.summary_tsv) if args.summary_tsv else {}
    text = args.report.read_text(encoding="utf-8")
    text = ensure_vchip_css(text)
    text = ensure_tooltip(text)
    new_sec = build_section(args.species, subfams, args.aln_dir,
                            counts, args.raw_base.rstrip("/"))
    pat = r'<section class=[\'"]card[\'"] id=[\'"]alignments[\'"].*?</section>'
    text, n = re.subn(pat, new_sec, text, count=1, flags=re.DOTALL)
    if n != 1:
        # insert before overview if no alignment section yet
        text, n = re.subn(
            r'(<section class="card" id="overview">)',
            new_sec + r"\1",
            text, count=1,
        )
    if n != 1:
        sys.exit("could not insert alignment section")
    args.report.write_text(text, encoding="utf-8")
    print("patched %d subfamilies → %s" % (len(subfams), args.report))


if __name__ == "__main__":
    main()
