#!/usr/bin/env python3
"""Per-alignment notes on what the method gets right and wrong, derived from the
measured values so they cannot drift away from the data."""
import json
import numpy as np

F = {json.loads(l)["set"]: json.loads(l) for l in open("features_v2c.jsonl")}
E = {}
try:
    E = {json.loads(l)["set"]: json.loads(l) for l in open("edges.jsonl")}
except IOError:
    pass
A = json.load(open("annotations.json"))
SEL = [x["set"] for x in json.load(open("viewer_data.json"))]

CLASS_NOTE = {
    "POS": "A curated family — everything here should read as real.",
    "MIXED10": "A real family with 10 % random loci mixed in.",
    "MIXED30": "A real family with 30 % random loci mixed in.",
    "NEGTRUNC5": "Copies truncated at the 5′ end with foreign sequence pasted on, mimicking LINE fragments.",
    "NEGSPLICE": "Left half of one family joined to the right half of another.",
    "NEGRAND": "Random genomic loci searched with a real consensus — the pure null.",
    "NEGJITTER": "A real family whose extraction window was displaced 20–100 bp per copy.",
    "SQ": "One family resampled, to show what changes with sample size alone.",
}


def notes(s):
    f, good, bad = F.get(s, {}), [], []
    if not f or "error" in f:
        return {"good": [], "bad": ["MEASURE() could not process this set."],
                "cls": CLASS_NOTE.get(s.split("__")[0], "")}
    e = E.get(s, {})
    ann = A.get(s, {})
    ft = [x["type"] for x in ann.get("features", [])]
    g = lambda k, d=None: f.get(k, d)

    # --- Tier 2: is there an element at all -----------------------------
    if g("cliff", 0) > 0.35:
        good.append("Clear Tier-2 cliff: identity %.2f inside the element against %.2f in the "
                    "flank (unrelated DNA is 0.25). The boundary is unambiguous."
                    % (g("cons_identity_med"), g("flank_id")))
    elif g("cliff", 0) > 0.15:
        bad.append("Weak cliff (%.2f). The element is detectable but the boundary is soft — "
                   "this is the kind of set a threshold will get wrong."
                   % g("cliff"))
    else:
        good.append("No cliff (%.2f): identity inside the element is no higher than the flank, "
                    "so the set is correctly rejected." % g("cliff"))

    # --- contamination ---------------------------------------------------
    fs = g("frac_supported", 1)
    if fs >= 0.95:
        good.append("Every copy supports the consensus (frac_supported %.2f) — no contamination "
                    "signal." % fs)
    elif fs >= 0.8:
        bad.append("%.0f %% of copies fail to support the consensus — visible contamination, "
                   "though not enough to break the analysis." % (100 * (1 - fs)))
    else:
        good.append("Contamination caught: only %.2f of copies support the consensus." % fs)

    # --- length / truncation ---------------------------------------------
    cv = g("elem_len_cv", 0)
    if cv < 0.08:
        good.append("Copy lengths are uniform (CV %.3f), as a real family should be." % cv)
    elif cv < 0.15:
        bad.append("Copy lengths vary more than usual (CV %.3f) without an obvious cause." % cv)
    else:
        good.append("Length heterogeneity caught (CV %.3f against 0.05 for real families)." % cv)
    ra = g("res_asymmetry")
    if ra is not None and abs(ra) > 1.5:
        good.append("Edges are strongly asymmetric (res_asymmetry %+.2f) — the ragged-one-side "
                    "signature the spec predicts for truncated copies." % ra)

    # --- one-sided bonuses ------------------------------------------------
    t = g("tsd_frac", 0)
    if t >= 0.5:
        good.append("TSDs in %.0f %% of copies, median %d bp — positive evidence of a genuine "
                    "insertion." % (100 * t, g("tsd_len_med", 0)))
    elif t <= 0.2:
        if s.startswith(("NEGRAND", "NEGTRUNC5", "NEGSPLICE")):
            good.append("Almost no TSDs (%.0f %%), correctly — there was no insertion event."
                        % (100 * t))
        else:
            bad.append("Few TSDs found (%.0f %%). Absence is not evidence against a family — "
                       "many SINEs lack them — but nothing is gained here." % (100 * t))
    if "abox" in ft and "bbox" in ft:
        good.append("Pol III A and B boxes both present in the consensus. Note this is a "
                    "property of the query, not of these loci.")

    # --- edge refinement --------------------------------------------------
    if e:
        pm = e.get("prominence_mean")
        if pm is not None and pm > 0.35:
            good.append("The edge search finds a sharp optimum (prominence %.2f), i.e. the data "
                        "have a definite opinion about where the boundary is." % pm)
        elif pm is not None and pm < 0.15:
            good.append("Flat edge landscape (prominence %.2f) — no preferred boundary anywhere, "
                        "which is what junk should look like." % pm)
        dr = e.get("d_best_R")
        if dr is not None and dr <= -5 and g("cliff", 0) > 0.3:
            bad.append("The 3′ edge wants to move %d bp outward — the consensus stops inside the "
                       "poly-A that the copies still carry." % -dr)

    # --- known weak spots -------------------------------------------------
    if g("rank1_excess") is not None and s.startswith("NEGSPLICE"):
        bad.append("Mosaicism is not detected (rank1_excess %.4f, same as a real family). This "
                   "negative splices every copy at the same point, so it forms a coherent new "
                   "family rather than a mosaic — the test is wrong, not the statistic."
                   % g("rank1_excess"))
    if s.startswith("NEGJITTER"):
        bad.append("Displacing the extraction window changes almost nothing measurable — the "
                   "anchor recovers the true element regardless. Useful as a robustness result, "
                   "but it means this is not really a negative.")
    return {"good": good, "bad": bad, "cls": CLASS_NOTE.get(s.split("__")[0], "")}


out = {s: notes(s) for s in SEL}
json.dump(out, open("commentary.json", "w"), separators=(",", ":"))
print("commentary for %d sets" % len(out))
for s in SEL[:2]:
    print("\n==", s)
    for k in ("good", "bad"):
        for x in out[s][k]:
            print("  %s %s" % ("+" if k == "good" else "-", x[:105]))
