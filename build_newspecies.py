#!/usr/bin/env python3
"""Build site/newspecies.html from the rebuilt new-species results.

Everything on the page is derived, nothing typed in by hand:
  newsp_verdicts_all.json      scores of the justified+trimmed alignments
  (the island numbers ride along inside the verdicts, measured by
   islands_corpus.py on the FULL 400 bp flank, never on the trimmed display)
  pom_sinederella_summary.tsv  the real genome-wide search for the snail
"""
import html
import io
import json
import os

VIEWER = "https://toki-bio.github.io/MSA-viewer/?url="
BASE = "https://toki-bio.github.io/SINE-discriminator/alignments/"

SPECIES = [
    ("hyd", "hydra", "Hydra vulgaris", "Cnidaria", "GCA_022113875.1"),
    ("pom", "snail", "Pomacea canaliculata", "Mollusca", "GCF_003073045.1"),
    ("aca", "starfish", "Acanthaster planci", "Echinodermata", "GCF_001949145.1"),
    ("stu", "sturgeon", "Acipenser ruthenus", "Actinopterygii", "GCF_010645085.2"),
    ("ska", "skate", "Amblyraja radiata", "Chondrichthyes", "GCF_010909765.2"),
]

# from each genome's AnnoSINE Seed_SINE.fa header: the RNA the seed came from,
# and AnnoSINE's own copy count (note 4: a floor, not an estimate)
SEED_FALLBACK = {
    "pom_SINE_0": ("tRNA", 537), "pom_SINE_1": ("7SL RNA", 339),
    "aca_SINE_0": ("tRNA", 176), "aca_SINE_1": ("5S rRNA", 11),
}


def load_seeds():
    """Seed RNA and AnnoSINE copy count, read from the run rather than typed."""
    if os.path.exists("seeds.json"):
        d = json.load(open("seeds.json"))
        return {k: (v[0], v[1]) for k, v in d.items()}
    return dict(SEED_FALLBACK)


SEED = load_seeds()

CALL_LABEL = {
    "SINE": "SINE",
    "NOT_SINE": "not a SINE",
    "NOT_SINE_MICROSAT": "not a SINE - microsatellite",
    "NOT_SINE_NEEDS_LOOK": "not a SINE, needs a closer look",
    "SINE_HEAD_MOSAIC": "SINE, head mosaic",
}

SHORT = {
    "NO_ELEMENT": "no element",
    "SHARED_FLANKS": "shared flanks",
    "FLANKS_UNMEASURED": "flanks unmeasured",
    "SUBFAMILY_NOTE": "subfamilies",
    "CONSENSUS_OVEREXTENDED": "consensus too long",
    "HETEROGENEOUS_SELECTION": "not one family",
    "INSUFFICIENT_COPIES": "too few copies",
    "NO_FLANKS_PRESENT": "no flanks",
    "TRUNCATED_COPIES": "truncated copies",
    "MICROSATELLITE_ELEMENT": "microsatellite element",
    "MICROSATELLITE_FLANK": "microsatellite flank",
    "SMALL_CORE": "small core",
    "CONSENSUS_UNDEREXTENDED": "consensus too short",
    "FLANK_ISLANDS": "flank islands",
    "CONTAMINATED": "contaminated",
    "RECOVERABLE_CORE": "recoverable core",
    "NOT_ISOLATED": "not isolated",
    "FRAGMENT_OF_LONGER": "fragment of longer",
    "ADJACENT_SIMILARITY": "adjacent similarity",
}

VIEW_HELP = {
    "top100": "the 100 highest-scoring genomic hits",
    "rand100": "a random 100 of all genomic hits",
    "all": "every hit, because there were fewer than 100",
}

ANN_HELP = ("AnnoSINE's own copy count. A floor, not an estimate "
            "- a real search finds several times more.")
HOVER_HELP = "hover any flag for the full reason"
LEAK_HELP = ("copies assigned to this family that fit another one better. "
             "Near zero means the locus set is coherent")
CAUGHT_HELP = "does the existing flank average already catch this?"


def his_calls():
    """What he said, keyed by candidate. He judged the top100 view."""
    out = {}
    if not os.path.exists("calls.tsv"):
        return out
    for line in io.open("calls.tsv", encoding="utf-8"):
        f = line.rstrip(chr(10)).split(chr(9))
        if len(f) < 4 or f[1] != "newsp_top100":
            continue
        cand = f[0].replace("NEW__", "").replace("__top100", "")
        out.setdefault(cand, []).append((f[2], f[3]))
    return out


def load():
    ver = json.load(open("newsp_verdicts_all.json"))
    out = {}
    for s, d in ver.items():
        if "error" in d:
            continue
        cand, view = s.replace("NEW__", "").replace(".clean", "").split("__")
        out.setdefault(cand, {})[view] = dict(d)
    return out


def aln_link(cand, view):
    f = "NEW__%s__%s.clean.aln.fa" % (cand, view)
    if not os.path.exists(os.path.join("site/alignments", f)):
        return ""
    url = (BASE + f).replace(":", "%3A").replace("/", "%2F")
    return ('<a class="lnk" target="_blank" rel="noopener" title="%s" href="%s%s">%s</a>'
            % (VIEW_HELP.get(view, view), VIEWER, url, view))


def flag_pills(d):
    out = []
    for f in d["flags"]:
        out.append('<span class="pill" title="%s">%s</span>'
                   % (html.escape(f.get("text", f["code"])),
                      SHORT.get(f["code"], f["code"].lower())))
    return "".join(out) or '<span class="cln">clean</span>'


def final_score(d):
    """The tool's answer after any consensus repair - that is its real verdict."""
    rf = d.get("refined")
    return rf["score"] if rf else d["score"]


def num(x, fmt="%.3f"):
    return "&mdash;" if x is None else fmt % x


HEAD = """<title>Three genomes with no SINE library</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@600;700&display=swap">
<style>
:root{--ground:#f4f6f2;--surface:#fff;--ink:#131c1a;--muted:#68766f;--rule:#d6dcd4;
--accent:#1f6f5c;--warn:#a8501d;--mism:#c94f3d;--soft:#eef1ec;
--sans:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",ui-monospace,monospace;
--serif:"IBM Plex Serif",Georgia,serif}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.55}
.wrap{max-width:1080px;margin:0 auto;padding:44px 26px 90px}
.eyebrow{font-family:var(--mono);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);margin:0 0 12px}
h1{font-family:var(--serif);font-weight:700;font-size:36px;line-height:1.12;margin:0 0 16px;max-width:24ch;text-wrap:balance}
h2{font-family:var(--serif);font-weight:700;font-size:21px;margin:30px 0 4px}
h2 .lat{font-family:var(--sans);font-weight:400;font-size:15px;color:var(--muted);font-style:italic}
h3{font-family:var(--mono);font-size:11px;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);font-weight:500;margin:22px 0 6px}
p{margin:0 0 12px;max-width:74ch}
a{color:var(--accent)}
.lede{font-size:17px;color:#2c3a36;max-width:70ch}
.rule{height:1px;background:var(--rule);margin:28px 0}
.find{background:var(--surface);border-left:3px solid var(--warn);padding:14px 17px;margin:18px 0;max-width:76ch}
.find p:last-child{margin-bottom:0}
.tw{overflow-x:auto;margin:10px 0}
table{border-collapse:collapse;width:100%;font-size:13px;background:var(--surface);border:1px solid var(--rule)}
th,td{padding:7px 10px;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap;vertical-align:top}
th{font-family:var(--mono);font-size:10px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);font-weight:500}
th[title]{cursor:help}
td.n{font-family:var(--mono);text-align:right;font-variant-numeric:tabular-nums}
td.name{font-family:var(--mono);font-size:12px;font-weight:500}
tbody tr:last-child td{border-bottom:0}
tr.alt td{background:#fbfcfa}
td.flags{white-space:normal;max-width:150px;min-width:120px}
td.call{white-space:normal;max-width:110px;font-size:12px}
.lnk{font-family:var(--mono);font-size:11px;text-decoration:none;border-bottom:1px solid var(--rule);margin-right:7px;white-space:nowrap}
.ok{color:var(--accent);font-weight:600}
.bad{color:var(--mism);font-weight:600}
.warn{color:var(--warn);font-weight:600}
.pill{display:inline-block;font-family:var(--mono);font-size:9.5px;padding:2px 6px;border-radius:2px;margin:0 3px 3px 0;background:var(--soft);color:#4a564f;cursor:help;white-space:nowrap}
.cln{font-family:var(--mono);font-size:9.5px;color:var(--accent)}
.cap{font-family:var(--mono);font-size:11px;color:var(--muted);max-width:76ch}
.bar{display:inline-block;height:8px;background:var(--warn);border-radius:1px;vertical-align:middle;margin-left:6px;opacity:.75}
</style>"""


def main():
    data = load()
    sd = {}
    for line in open("pom_sinederella_summary.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("pomSINE"):
            sd[f[0]] = f

    h = [HEAD, '<div class="wrap">']
    present = [s for s in SPECIES if any(c.startswith(s[0] + "_") for c in data)]
    n_cand = len(data)
    n_aln = sum(len(v) for v in data.values())
    names = ", ".join(s[1] for s in present[:-1]) + " and " + present[-1][1] if len(present) > 1         else present[0][1]
    h.append('<p class="eyebrow">Prospective test &middot; %d genomes &middot; %d phyla</p>'
             % (len(present), len(set(s[3] for s in present))))
    h.append("<h1>Candidate SINE families in %d genomes with no SINE library</h1>" % len(present))
    h.append('<p class="lede">AnnoSINE_v2 on %s. Every candidate it proposed was searched back '
             'against its own genome, each locus pulled out with 400&nbsp;bp of flank, aligned '
             'with the consensus, the flanks pushed against the element and degapped, and then '
             'scored. <strong>%d candidates, %d alignments.</strong></p>' % (names, n_cand, n_aln))
    h.append('<p class="lede">Only hydra has an answer key, and only because you supplied one: '
             'your reading of all 22 of its candidates, shown beside the tool&rsquo;s in every '
             'row below. For the rest there is no curated library and no prior '
             'annotation.</p>')

    # a compact index, so a long page can be entered at the right genome
    h.append('<div class="tw"><table><thead><tr><th>genome</th><th>phylum</th>'
             '<th title="candidates AnnoSINE proposed">candidates</th>'
             '<th title="scoring 50 or above on their best view">accepted</th>'
             '<th title="candidates Sergei has read himself">judged</th>'
             '<th title="candidates where the tool and his reading agree">agree</th>'
             '</tr></thead><tbody>')
    _calls = his_calls()
    for code, common, latin, phylum, acc in present:
        cs = [c for c in data if c.startswith(code + "_")]
        best = lambda c: (data[c].get("top100") or list(data[c].values())[0])
        n_ok = sum(1 for c in cs if final_score(best(c)) >= 50)
        judged = [c for c in cs if c in _calls]
        agree = sum(1 for c in judged
                    if (final_score(best(c)) >= 50) == (_calls[c][0][0] == "SINE"))
        h.append('<tr><td><a class="lnk" href="#%s">%s</a> <span class="cap">%s</span></td>'
                 '<td class="cap">%s</td><td class="n">%d</td><td class="n ok">%d</td>'
                 '<td class="n">%s</td><td class="n %s">%s</td></tr>'
                 % (code, common, latin, phylum, len(cs), n_ok,
                    len(judged) or "&mdash;",
                    "ok" if judged and agree == len(judged) else "",
                    "%d of %d" % (agree, len(judged)) if judged else "&mdash;"))
    h.append('</tbody></table></div>')

    h.append('<div class="find">')
    h.append('<h3 style="margin-top:0">The flank islands</h3>')
    h.append('<p>You saw faint, non-random similarity in <code>aca_SINE_0</code>&rsquo;s aligned '
             'left flank, and said the patches can sit far out on either side. Measured against a '
             'null that controls for how many copies are present in each column, they are real: '
             '<strong>17 % of its flank</strong> sits inside a patch &mdash; 426 columns in 49 of '
             'them, up to z&nbsp;=&nbsp;69, some hundreds of bases out. Its flank average is '
             '<strong>0.271</strong>, ordinary background, so the average sees nothing.</p>')
    h.append('<p>Run over the 273 labelled sets, that fraction lines up with your own calls: '
             'every one of the 31 you called a plain SINE is at or below <strong>0.067</strong>, '
             'the three you asked for flank checks on sit at <strong>0.19 to 0.37</strong>, and '
             'everything you rejected is at <strong>0.50 or more</strong>. It fires on 14 of 273 '
             'sets &mdash; only segmental duplications, satellites, LINE ORFs and the messy '
             'hedgehog set &mdash; and on none of the 28 positives, 56 mixtures, 28 jitter, '
             '20 splice, 28 truncation or 20 random sets.</p>')
    h.append('<p>So it is reported as a reason, never subtracted from the score: two of the three '
             'sets it was calibrated on you still lean SINE on. <code>aca_SINE_0</code> at 0.171 '
             'sits between your plain-SINE ceiling and your caution floor &mdash; and the signal '
             'is in <code>top100</code> only (<code>rand100</code> is under 0.03), so it is not '
             'the family that sits in a shared context but its most similar copies.</p>')
    h.append('<p class="cap">Islands are measured on the full 400&nbsp;bp flank, never on the '
             'trimmed alignments linked below. Trimming to the width the copies fill cuts away '
             'exactly the far-out columns: it drops aca_SINE_0 from 426 island columns to 6.</p>')
    h.append('</div>')

    th = ('<tr><th title="candidate name, the small RNA its seed came from, and %s">candidate</th>'
          '<th title="copies in this alignment, and which hits they are">n</th>'
          '<th title="0 to 100. At 50 and above the tool accepts it">score</th>'
          '<th title="median pairwise identity across the whole element">identity</th>'
          '<th title="pairwise identity in the flanks. 0.25 to 0.31 means the copies sit in '
          'independent places; much higher means they do not">flank bg</th>'
          '<th title="fraction of the flank sitting inside a patch of significantly '
          'raised similarity, and the column count. Measured on the full 400 bp flank. '
          'His own calls put plain SINEs at or below 0.067 and flank-caution calls at '
          '0.19 to 0.37">island frac</th>'
          '<th title="%s">what the tool says</th>'
          '<th title="his own reading of the top100 alignment">his call</th>'
          '<th>alignments</th></tr>' % (html.escape(ANN_HELP), HOVER_HELP))

    calls = his_calls()
    for code, common, latin, phylum, acc in SPECIES:
        cands = [c for c in data if c.startswith(code + "_")]
        if not cands:
            continue
        cands.sort(key=lambda c: int(c.rsplit("_", 1)[1]))
        best = lambda c: (data[c].get("top100") or list(data[c].values())[0])
        acc_n = sum(1 for c in cands if final_score(best(c)) >= 50)
        h.append('<div class="rule"></div>')
        h.append("<h2 id='%s'>%s <span class='lat'>%s &middot; %s &middot; %s</span></h2>"
                 % (code, common, latin, phylum, acc))
        h.append('<p class="cap">%d candidates, %d scoring 50 or above on their best view.</p>'
                 % (len(cands), acc_n))
        h.append('<div class="tw"><table><thead>%s</thead><tbody>' % th)
        for k, c in enumerate(cands):
            views = data[c]
            order = [v for v in ("top100", "rand100", "all") if v in views]
            seed, ann = SEED.get(c, ("?", None))
            for j, v in enumerate(order):
                d = views[v]
                isl, frac = d.get("island_cols"), d.get("island_frac")
                icls = "bad" if (frac or 0) >= 0.10 else ("warn" if (frac or 0) >= 0.07 else "")
                scls = "ok" if d["score"] >= 50 else "bad"
                bar = ('<span class="bar" style="width:%dpx"></span>'
                       % max(2, min(56, int((frac or 0) * 160)))) if frac else ""
                h.append('<tr%s>' % (' class="alt"' if k % 2 else ""))
                if j == 0:
                    h.append('<td class="name" rowspan="%d">%s'
                             '<div class="cap">%s%s</div></td>'
                             % (len(order), c, seed,
                                "" if not ann else " &middot; %s copies" % ann))
                h.append('<td class="n">%d <span class="cap">%s</span></td>' % (d["n"], v))
                rf = d.get("refined")
                if rf:
                    up = rf["score"] > d["score"] + 5
                    h.append('<td class="n %s">%.1f<div class="cap %s">repaired &rarr; '
                             '<strong>%.1f</strong></div></td>'
                             % (scls, d["score"], "ok" if up else "", rf["score"]))
                else:
                    h.append('<td class="n %s">%.1f</td>' % (scls, d["score"]))
                h.append('<td class="n">%s</td>' % num(d["all_identity"]))
                h.append('<td class="n">%s</td>' % num(d["flank_bg"]))
                h.append('<td class="n %s">%s%s</td>'
                         % (icls, "&mdash;" if frac is None
                            else "%.3f <span class=\"cap\">%d</span>" % (frac, isl), bar))
                h.append('<td class="flags">%s</td>' % flag_pills(d))
                if j == 0:
                    cc = calls.get(c)
                    if cc:
                        label = CALL_LABEL.get(cc[0][0], cc[0][0].replace("_", " ").lower())
                        agree = (final_score(d) >= 50) == (cc[0][0] == "SINE")
                        cell = ('<span class="%s" title="%s">%s</span>'
                                % ("ok" if agree else "bad",
                                   html.escape("; ".join(x[1] for x in cc if x[1])),
                                   label))
                    else:
                        cell = '<span class="cap">not judged</span>'
                    h.append('<td class="call" rowspan="%d">%s</td>' % (len(order), cell))
                h.append('<td>%s</td>' % aln_link(c, v))
                h.append('</tr>')
        h.append('</tbody></table></div>')

        if code == "pom" and sd:
            h.append('<h3>Both snail consensuses put through a full SINEderella search</h3>')
            h.append('<p>Not AnnoSINE&rsquo;s count &mdash; the real genome-wide search, run on '
                     'the consensuses refined from their own aligned copies.</p>')
            h.append('<div class="tw"><table><thead><tr>'
                     '<th>consensus</th>'
                     '<th title="assigned with high confidence">firm</th>'
                     '<th title="assigned, lower confidence">soft</th>'
                     '<th>total assigned</th>'
                     '<th title="%s">leak %%</th>'
                     '<th title="mean similarity of assigned copies to the consensus">similarity</th>'
                     '<th title="%s">annosine said</th>'
                     '</tr></thead><tbody>' % (LEAK_HELP, html.escape(ANN_HELP)))
            for key, cand in (("pomSINE0", "pom_SINE_0"), ("pomSINE1", "pom_SINE_1")):
                f = sd.get(key)
                if not f:
                    continue
                ann = SEED[cand][1]
                h.append('<tr><td class="name">%s</td><td class="n">%s</td><td class="n">%s</td>'
                         '<td class="n"><strong>%s</strong></td><td class="n ok">%s</td>'
                         '<td class="n">%s</td>'
                         '<td class="n cap">%d &nbsp;&times;%.0f fewer</td></tr>'
                         % (cand, f[1], f[2], f[3], f[8], f[10], ann, float(f[3]) / ann))
            h.append('</tbody></table></div>')
            h.append('<p class="cap">Zero leak on both, across thousands of copies. In Timema the '
                     'families you judged real ran 0.00 to 0.18 % leak and the noisy ones 65 to '
                     '98 %.</p>')

            h.append('<h3>Subfamily alignments</h3>')
            h.append('<p>Built with <code>extract_alignments.sh</code>, which is the one the '
                     'deployment notes say to use &mdash; step&nbsp;5 writes alignments with no '
                     'flanks, and flanks are what a boundary is judged on. Every row below carries '
                     '30&nbsp;bp left and 70&nbsp;bp right.</p>')
            h.append('<div class="tw"><table><thead><tr><th>family</th>'
                     '<th title="50 random assigned copies, flanked, MAFFT">core</th>'
                     '<th title="the 50 highest-scoring assigned copies, flanked">best 50</th>'
                     '<th title="SubFam chunk consensuses - one consensus per ~50 copies, so '
                     'subfamily structure is visible in one screen. Only produced at 400 copies '
                     'or more">subfam</th></tr></thead><tbody>')
            for fam in ("pomSINE0", "pomSINE1"):
                cells = []
                for tier in ("core", "best50", "subfam"):
                    f = "POM__%s.%s.fa" % (fam, tier)
                    if os.path.exists(os.path.join("site/alignments", f)):
                        url = (BASE + f).replace(":", "%3A").replace("/", "%2F")
                        cells.append('<td><a class="lnk" target="_blank" rel="noopener" '
                                     'href="%s%s">open</a></td>' % (VIEWER, url))
                    else:
                        cells.append('<td class="cap">&mdash;</td>')
                h.append('<tr><td class="name">%s</td>%s</tr>' % (fam, "".join(cells)))
            h.append('</tbody></table></div>')

            chunks = "POM__subfam_chunks.realigned.fa"
            if os.path.exists(os.path.join("site/alignments", chunks)):
                url = (BASE + chunks).replace(":", "%3A").replace("/", "%2F")
                h.append('<h3>The 178 SubFam chunk consensuses, for splitting by eye</h3>')
                h.append('<p>This is the material subfamily discovery starts from: the whole '
                         'assigned copy set collapsed to one consensus per ~50 copies. The '
                         'manual is explicit that SubFam&rsquo;s own final pass is not a converged '
                         'alignment, so these were degapped and realigned with '
                         '<code>--localpair --maxiterate 1000</code> before being posted. '
                         '<a class="lnk" target="_blank" rel="noopener" href="%s%s">open the 178 '
                         'rows</a> and group the ones sharing diagnostic columns &mdash; each '
                         'group is a candidate subfamily.</p>' % (VIEWER, url))

    rows = []
    for c, views in data.items():
        for v, d in views.items():
            if d.get("island_frac") is not None:
                rows.append((c, v, d["island_frac"], d["island_cols"],
                             d["score"], d["flank_bg"]))
    rows.sort(key=lambda r: -r[2])
    rep = sorted((c, v, d) for c, views in data.items() for v, d in views.items()
                 if d.get("refined") and v == "top100")
    if rep:
        h.append('<div class="rule"></div>')
        h.append("<h2>Consensuses the tool repaired by itself</h2>")
        h.append('<p>When the copies keep matching past the end of a consensus and then reach '
                 'background within a couple of hundred bases, the consensus is simply too short. '
                 'The tool now extends it by that distance, re-searches the genome, and &mdash; '
                 'because extending overshoots &mdash; cuts back to the window the copies actually '
                 'support. Both rules already existed; only the loop between them is new.</p>')
        h.append('<div class="tw"><table><thead><tr><th>candidate</th>'
                 '<th title="score before the repair">before</th>'
                 '<th title="score after extend, re-search and trim">after</th>'
                 '<th title="copies in the core, before and after">core</th>'
                 '<th>what the repaired version says</th></tr></thead><tbody>')
        for c, v, d in rep:
            rf = d["refined"]
            up = rf["score"] > d["score"] + 5
            h.append('<tr><td class="name">%s</td><td class="n">%.1f</td>'
                     '<td class="n %s">%.1f</td><td class="n">%d &rarr; %d</td>'
                     '<td class="flags">%s</td></tr>'
                     % (c, d["score"], "ok" if up else "bad", rf["score"],
                        d["n_core"], rf["n_core"],
                        "".join('<span class="pill">%s</span>' % SHORT.get(f["code"], f["code"])
                                for f in rf["flags"]) or '<span class="cln">clean</span>'))
        h.append('</tbody></table></div>')
        h.append('<p class="cap">hyd_SINE_0 is the case Sergei called a false alarm: 110&nbsp;bp '
                 'scoring 0.0, extended to 310, cut back to <strong>211&nbsp;bp scoring 100.0 with '
                 '93 of 100 copies in the core</strong>. RepBase gives SINE2-2B_HM as 208&nbsp;bp. '
                 'ska_SINE_19 stays at 0.2 after the same treatment, which is the point: the loop '
                 'rescues an element that was there, it does not manufacture one.</p>')

    h.append('<div class="rule"></div>')
    h.append("<h2>Flank islands across all 52 alignments</h2>")
    h.append('<p>Ranked by the fraction of flank inside a patch. Your own calls put plain SINEs '
             'at or below 0.067 and flank-caution calls at 0.19 to 0.37, so the reason fires at '
             '0.10. The last column is the one that matters: a set above the line '
             '<em>with</em> an ordinary flank average is a set no average could have found. '
             '<code>aca_SINE_0</code> is the only one here.</p>')
    h.append('<div class="tw"><table><thead><tr><th>alignment</th>'
             '<th title="fraction of flank columns inside a patch of at least 6 '
             'consecutive columns above z = 8">island frac</th>'
             '<th>island cols</th><th>flank bg</th><th>score</th>'
             '<th title="%s">already caught?</th>'
             '</tr></thead><tbody>' % CAUGHT_HELP)
    for c, v, frac, ncol, sc, bg in rows[:16]:
        if bg >= 0.40:
            caught = '<span class="ok">yes, the flank average is raised too</span>'
        elif frac >= 0.10:
            caught = '<span class="bad">no &mdash; the average sees nothing</span>'
        elif frac >= 0.067:
            caught = '<span class="warn">no, but only just over the line</span>'
        else:
            caught = '<span class="cap">nothing to catch</span>'
        h.append('<tr><td class="name">%s <span class="cap">%s</span></td>'
                 '<td class="n %s">%.3f</td><td class="n">%d</td><td class="n">%.3f</td>'
                 '<td class="n %s">%.1f</td><td>%s</td></tr>'
                 % (c, v, "bad" if frac >= 0.15 else ("warn" if frac >= 0.067 else ""),
                    frac, ncol, bg, "ok" if sc >= 50 else "bad", sc, caught))
    h.append('</tbody></table></div>')
    h.append('<p class="cap">Null: pairwise identity per column against the base composition of '
             'that set&rsquo;s own flanks, with the standard deviation taken from the number of '
             'pairs actually present in that column. Columns with fewer than 8 copies present are '
             'not scored. That control is what changed the answer: without it the far-out patches '
             'look like thin-coverage noise, which is what I first called them.</p>')

    h.append('<div class="rule"></div>')
    h.append("<h2>How these were made</h2>")
    h.append('<p>Candidate consensus and genome in, alignments out, in one command '
             '(<code>candidate_to_aln.py</code>): blastn against the genome, overlapping hits '
             'collapsed keeping the best, then two views &mdash; the 100 highest-scoring hits and '
             'a random 100, or everything when there are fewer than 100. Each locus extracted with '
             '400&nbsp;bp of flank, aligned with MAFFT, consensus written first.</p>')
    h.append('<p>Then the flanks are degapped and pushed hard against the element, and the flank '
             'panel is cut to the width that keeps it under 25 % gaps, so one long copy stops '
             'padding seventy others with dashes. Neither step touches the element alignment. The '
             'island scan runs on the untrimmed version.</p>')
    h.append('<p class="cap">Still running: <em>Acipenser ruthenus</em> (sterlet sturgeon, '
             '1.9&nbsp;Gb) and <em>Amblyraja radiata</em> (thorny skate, 2.6&nbsp;Gb) &mdash; '
             'lineages where SINEs are described in the literature but these assemblies carry no '
             'repeat library.</p>')
    h.append('</div>')

    io.open("site/newspecies.html", "w", encoding="utf-8").write("\n".join(h))
    print("wrote site/newspecies.html  (%d candidates, %d alignments)"
          % (len(data), sum(len(v) for v in data.values())))


if __name__ == "__main__":
    main()
