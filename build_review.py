#!/usr/bin/env python3
"""Build the Test 4 review page: alignments whose judgement would settle a threshold.

Eleven of the sixteen things Sergei described rest on one or two alignments
each. No threshold can be set from one point. This page shows, for each of
those, the unjudged alignments that already carry the matching signal - so his
answer either confirms the threshold or moves it.

Only alignments where his answer changes something. Not a request to re-judge
anything.
"""
import io
import json

VIEW = "https://toki-bio.github.io/MSA-viewer/?url="
BASE = "https://toki-bio.github.io/SINE-discriminator/alignments/"

WHY = {
 "SINE_CONSENSUS_WRONG":
   ("You said of MIR1_Amn: <em>\"weak old sine, consensus needs refinement and "
    "shortening, otherwise difficult but very real SINE\"</em>. That is the only "
    "example. These four are measured as having the same property - the copies "
    "support a shorter span than the consensus claims.",
    "Is the consensus too long here, or is something else wrong?"),
 "MIXTURE / needs split":
   ("You said of AluYk11: <em>\"truly a mixture ... needs additional work before "
    "verdict\"</em>. Three examples so far. These are measured as splitting into "
    "two groups of copies with different lengths.",
    "Is this one family or two things collected together?"),
 "SINE_TOO_FEW":
   ("You said of AluSx_short_: <em>\"too few copies but looks like SINE ... needs "
    "more evidence from genome\"</em>, and of a one-copy set simply "
    "<em>\"wtf???\"</em>. Where is the line - is 11 copies enough? 15? 22?",
    "Enough copies to judge, or not?"),
 "flanks not independent":
   ("The decisive signal you named is that flanks must NOT align. These are "
    "measured as having flanks that DO resemble each other - satellites and "
    "segmental duplications.",
    "Are these clearly not SINEs, or is the flank similarity misleading?"),
 "GREY (40-70)":
   ("Scores between 40 and 70, where the tool is least certain and where you have "
    "given only two judgements.",
    "SINE, not a SINE, or genuinely undecidable from this alignment?"),
 "contaminated":
   ("Measured as a real family with junk mixed in. You have judged none of these "
    "directly, and there are 56 of them.",
    "Is the contamination call right, and does it change the verdict or just add a note?"),
 "fragment of longer":
   ("Measured as similarity continuing past the annotated boundary. Includes "
    "SINE_25, which you called real - and three LINE fragments, which should look "
    "different.",
    "Does the element end where the consensus says, or continue?"),
}

CSS = """
:root{--ground:#f4f6f2;--surface:#fff;--ink:#131c1a;--muted:#68766f;--rule:#d6dcd4;
--accent:#1f6f5c;--accent-soft:#e2efe9;--warn:#a8501d;--warn-soft:#f6e8de;--mism:#c94f3d;
--sans:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",monospace;
--serif:"IBM Plex Serif",Georgia,serif}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.55}
.wrap{max-width:1000px;margin:0 auto;padding:44px 26px 90px}
.eyebrow{font-family:var(--mono);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);margin:0 0 12px}
h1{font-family:var(--serif);font-weight:700;font-size:38px;line-height:1.12;margin:0 0 16px;max-width:22ch;text-wrap:balance}
h2{font-family:var(--serif);font-weight:700;font-size:22px;margin:34px 0 6px}
p{margin:0 0 12px;max-width:74ch}
a{color:var(--accent)}
.lede{font-size:17px;color:#2c3a36;max-width:68ch}
.rule{height:1px;background:var(--rule);margin:30px 0}
.q{background:var(--surface);border-left:3px solid var(--warn);padding:11px 15px;margin:12px 0;font-weight:600;max-width:74ch}
table{border-collapse:collapse;width:100%;font-size:13.5px;background:var(--surface);border:1px solid var(--rule);margin:10px 0 4px}
th,td{padding:8px 11px;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap}
th{font-family:var(--mono);font-size:10px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);font-weight:500}
td.n{font-family:var(--mono);text-align:right;font-variant-numeric:tabular-nums}
td.name{font-family:var(--mono);font-size:12.5px}
tr:last-child td{border-bottom:0}
.pill{display:inline-block;font-family:var(--mono);font-size:10px;padding:2px 6px;border-radius:2px;margin-right:4px;background:#eceee9;color:var(--muted)}
.lnk{font-family:var(--mono);font-size:11.5px;text-decoration:none;border-bottom:1px solid var(--rule)}
.cap{font-family:var(--mono);font-size:11px;color:var(--muted)}
"""


def main():
    import sys
    src = sys.argv[1] if len(sys.argv) > 1 else "test4_list.json"
    items = json.load(open(src))
    groups = {}
    for label, name, corp, score, flags, n in items:
        groups.setdefault(label, []).append((name, score, flags, n))

    h = []
    h.append("<title>Alignments needing your judgement</title>")
    h.append('<link rel="preconnect" href="https://fonts.googleapis.com">')
    h.append('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
             'family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600'
             '&family=IBM+Plex+Serif:wght@600;700&display=swap">')
    h.append("<style>%s</style>" % CSS)
    h.append('<div class="wrap">')
    h.append('<p class="eyebrow">Test 4 &middot; %d alignments &middot; %d questions</p>'
             % (len(items), len(groups)))
    h.append("<h1>Alignments where your answer would settle something</h1>")
    h.append('<p class="lede">Eleven of the sixteen properties you described rest on '
             'one or two alignments each. A threshold cannot be set from one point. '
             'Each group below is a property with too few examples, and the alignments '
             'the tool measures as having that same property.</p>')
    h.append('<p class="lede">Every alignment opens in the MSA viewer with its '
             'flanks justified <strong>and trimmed to the width the copies actually '
             'use</strong>, so a few long outliers no longer pad the rest with gaps. '
             'On these, flank gaps drop from 0.56 to 0.16, 0.68 to 0.26, 0.67 to 0.07.</p>')
    h.append('<p class="lede">Six at a time. A one-line answer per alignment is '
             'enough.</p>')

    for label, rows in groups.items():
        why, question = WHY.get(label, ("", ""))
        h.append('<div class="rule"></div>')
        h.append("<h2>%s</h2>" % label)
        if why:
            h.append("<p>%s</p>" % why)
        h.append('<div class="q">%s</div>' % question)
        h.append("<table><thead><tr><th>alignment</th><th>copies</th>"
                 "<th>score</th><th>what the tool flags</th><th>open</th></tr></thead><tbody>")
        for name, score, flags, n in rows:
            # link the TRIMMED alignment: flanks justified AND capped to the width
            # the copies actually use, so a few long outliers stop padding the rest
            # with gaps ("flanks are badly degapped and it hinders my estimates")
            full = VIEW + (BASE + name + ".trim.aln.fa").replace(":", "%3A").replace("/", "%2F")
            pills = "".join('<span class="pill">%s</span>' % f for f in flags) or \
                    '<span style="color:#9aa69f">none</span>'
            h.append('<tr><td class="name">%s</td><td class="n">%s</td>'
                     '<td class="n">%s</td><td>%s</td>'
                     '<td><a class="lnk" href="%s" target="_blank" rel="noopener">view</a></td></tr>'
                     % (name, n, "" if score is None else "%.1f" % score, pills, full))
        h.append("</tbody></table>")

    h.append('<div class="rule"></div>')
    h.append('<p class="cap">Built by build_review.py from test4_list.json. '
             'Candidates chosen as unjudged alignments already carrying the signal '
             'for a property that has too few examples.</p>')
    h.append("</div>")
    io.open(sys.argv[2] if len(sys.argv) > 2 else "site/review.html", "w", encoding="utf-8").write("\n".join(h))
    print("wrote site/review.html: %d groups, %d alignments" % (len(groups), len(items)))


if __name__ == "__main__":
    main()
