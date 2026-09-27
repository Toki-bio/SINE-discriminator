#!/usr/bin/env python3
"""The Overall call on a published plate: the verdict score plus what the copies were seen to do.

verdict.py withholds its score below 30 copies (INSUFFICIENT_COPIES -> "Cannot assess"), so that a
small but clean family is not rejected for being small (HUM AluYb9). That rule stays. What it also
did was hide negative evidence the copies DO show: on rsi (2026-09-27) MEG-RL had 20 copies, no
supported element (NO_ELEMENT), a consensus the copies cut ~45 bp short at the 3' end, and left and
right flanks shared by 14-20 % of copies - and the report said only "Cannot assess". The flank
context column (flank_uniqueness.scan) was never part of the Overall call at all.

Rules:
  * assessable (>= 30 copies, flanks present): the verdict score as before.
  * not assessable, with negative evidence in the verdict flags: "Doubtful", reasons listed.
  * not assessable, nothing negative observed: "Cannot assess", as before.
Every call carries its reasons, shown in the tooltip.

The flank-context scan is listed ("also seen") but does not move the call. Tested 2026-09-27 on the
224 published plates scoring >= 75: 30 have rand100 "high" shared flanks, among them hum CAS, FLAM_C,
FRAM and Rhin-1 in six bats - random copies sharing a flank group is common in real families, so it
cannot demote a call on its own.
"""

NEGATIVE = {
    "NO_ELEMENT": "the copies do not support the consensus above background",
    "SMALL_CORE": "too few copies form the core the consensus is built on",
    "MICROSATELLITE_ELEMENT": "the element is largely simple repeat",
    "NOT_ISOLATED": "copies sit in satellites or duplications",
    "SHARED_FLANKS": "copies share flanking sequence",
    "FRAGMENT_OF_LONGER": "the element continues past the consensus (part of a longer repeat)",
    "ELEMENT_CONTINUES": "the element continues past the consensus (part of a longer repeat)",
    "CONSENSUS_OVEREXTENDED": "the consensus is longer than the element the copies support",
}


def _flags(v):
    return {f["code"]: f for f in v.get("flags", [])}


def _context(top, rand):
    """(severity, reasons) from flank_uniqueness scans of top100 and rand100."""
    sev, out = None, []
    for label, r in (("random copies", rand), ("top copies", top)):
        if not r or r.get("error"):
            continue
        for f in r.get("flags", []):
            first = f["text"].split(". ")[0].rstrip(".")
            out.append("%s, %s" % (label, first[0].lower() + first[1:]))
        w = r.get("worst_flag")
        if w == "high" or (w == "medium" and sev is None):
            sev = w
    return sev, out


def overall_call(v, top_scan=None, rand_scan=None):
    """-> {"label", "kind" (ok/edge/warn/muted), "reasons": [str]}"""
    if not v or v.get("error"):
        return {"label": "n/a", "kind": "muted",
                "reasons": [v.get("error", "not scored") if v else "not scored"]}
    if v.get("deferred"):
        return {"label": "Deferred", "kind": "edge",
                "reasons": ["Mixture of two groups - split before a firm call."]}

    fl = _flags(v)
    n = v.get("n", 0)
    neg = [NEGATIVE[c] + (" (%d/%d copies support it)" % (v.get("n_supported", 0), n)
                          if c == "NO_ELEMENT" else "")
           for c in NEGATIVE if c in fl]
    # the two ELEMENT_CONTINUES/FRAGMENT texts are identical; keep one
    neg = list(dict.fromkeys(neg))
    sev, ctx = _context(top_scan, rand_scan)
    also = (["Also seen (does not change the call) - shared flank groups:"] + ["- " + c for c in ctx]
            if ctx else [])

    if not v.get("assessable", True):
        why = []
        if "INSUFFICIENT_COPIES" in fl:
            why.append("only %d copies (a score needs about 30)" % n)
        if "NO_FLANKS_PRESENT" in fl:
            why.append("no flank sequence could be measured")
        head = "Too little evidence for a score: " + "; ".join(why) + "."
        if neg:
            return {"label": "Doubtful", "kind": "warn",
                    "reasons": [head, "What the copies do show is negative:"]
                    + ["- " + r for r in neg] + also}
        return {"label": "Cannot assess", "kind": "muted",
                "reasons": [head, "Nothing negative observed."] + also}

    s = float(v.get("score", 0))
    head = "Score %.0f/100" % s
    if v.get("capped_by"):
        head += " (capped: %s)" % ", ".join(v["capped_by"])
    if s >= 90:
        label, kind = "SINE", "ok"
    elif s >= 75:
        label, kind = "SINE (caveats)", "ok"
    elif s >= 55:
        label, kind = "Grey zone", "edge"
    else:
        label, kind = "Not SINE", "warn"
    reasons = [head + "."]
    reasons += ["- " + r for r in neg] + also
    return {"label": label, "kind": kind, "reasons": reasons}
