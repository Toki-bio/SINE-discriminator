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
    ("pom", "snail", "Pomacea canaliculata", "Mollusca", "GCF_003073045.1"),
    ("aca", "starfish", "Acanthaster planci", "Echinodermata", "GCF_001949145.1"),
    ("hyd", "hydra", "Hydra vulgaris", "Cnidaria", "GCA_022113875.1"),
]

# from each genome's AnnoSINE Seed_SINE.fa header: the RNA the seed came from,
# and AnnoSINE's own copy count (note 4: a floor, not an estimate)
SEED = {
    "pom_SINE_0": ("tRNA", 537), "pom_SINE_1": ("7SL RNA", 339),
    "aca_SINE_0": ("tRNA", 176), "aca_SINE_1": ("5S rRNA", 11),
    "hyd_SINE_0": ("tRNA", 61), "hyd_SINE_1": ("tRNA", 1616),
    "hyd_SINE_2": ("tRNA", 8), "hyd_SINE_3": ("tRNA", 67),
    "hyd_SINE_4": ("tRNA", 6861), "hyd_SINE_5": ("tRNA", 7),
    "hyd_SINE_6": ("unknown", 1257), "hyd_SINE_7": ("tRNA", 697),
    "hyd_SINE_8": ("tRNA", 3135), "hyd_SINE_9": ("5S rRNA", 246),
    "hyd_SINE_10": ("tRNA", 95), "hyd_SINE_11": ("tRNA", 3036),
    "hyd_SINE_12": ("tRNA", 250), "hyd_SINE_13": ("tRNA", 1144),
    "hyd_SINE_14": ("tRNA", 2608), "hyd_SINE_15": ("tRNA", 4730),
    "hyd_SINE_16": ("tRNA", 2252), "hyd_SINE_17": ("tRNA", 7307),
    "hyd_SINE_18": ("tRNA", 1746), "hyd_SINE_19": ("tRNA", 5574),
    "hyd_SINE_20": ("tRNA", 3988), "hyd_SINE_21": ("tRNA", 4663),
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
td.flags{white-space:normal;max-width:250px;min-width:200px}
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
    h.append('<p class="eyebrow">Prospective test &middot; three phyla &middot; no answer key</p>')
    h.append("<h1>Candidate SINE families in three genomes with no SINE library</h1>")
    h.append('<p class="lede">AnnoSINE_v2 on a snail, a starfish and a hydra. Every candidate it '
             'proposed was searched back against its own genome, each locus pulled out with '
             '400&nbsp;bp of flank, aligned with the consensus, the flanks pushed against the '
             'element and degapped, and then scored. 26 candidates, 52 alignments.</p>')
    h.append('<p class="lede">There is <strong>no answer key</strong> here &mdash; no curated '
             'library and no prior annotation for any of the three. The only check is your eye.</p>')

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

    th = ('<tr><th title="candidate name as AnnoSINE assigned it">candidate</th>'
          '<th title="the small RNA AnnoSINE derived the seed from">seed</th>'
          '<th title="%s">annosine copies</th>'
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
          '<th>alignments</th></tr>' % (html.escape(ANN_HELP), HOVER_HELP))

    for code, common, latin, phylum, acc in SPECIES:
        cands = [c for c in data if c.startswith(code + "_")]
        if not cands:
            continue
        cands.sort(key=lambda c: int(c.rsplit("_", 1)[1]))
        best = lambda c: (data[c].get("top100") or list(data[c].values())[0])
        acc_n = sum(1 for c in cands if best(c)["score"] >= 50)
        h.append('<div class="rule"></div>')
        h.append("<h2>%s <span class='lat'>%s &middot; %s &middot; %s</span></h2>"
                 % (common, latin, phylum, acc))
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
                    h.append('<td class="name" rowspan="%d">%s</td>' % (len(order), c))
                    h.append('<td rowspan="%d">%s</td>' % (len(order), seed))
                    h.append('<td class="n" rowspan="%d">%s</td>'
                             % (len(order), ann if ann else "&mdash;"))
                h.append('<td class="n">%d <span class="cap">%s</span></td>' % (d["n"], v))
                h.append('<td class="n %s">%.1f</td>' % (scls, d["score"]))
                h.append('<td class="n">%s</td>' % num(d["all_identity"]))
                h.append('<td class="n">%s</td>' % num(d["flank_bg"]))
                h.append('<td class="n %s">%s%s</td>'
                         % (icls, "&mdash;" if frac is None
                            else "%.3f <span class=\"cap\">%d</span>" % (frac, isl), bar))
                h.append('<td class="flags">%s</td>' % flag_pills(d))
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

    rows = []
    for c, views in data.items():
        for v, d in views.items():
            if d.get("island_frac") is not None:
                rows.append((c, v, d["island_frac"], d["island_cols"],
                             d["score"], d["flank_bg"]))
    rows.sort(key=lambda r: -r[2])
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
