#!/usr/bin/env python3
"""Where every dataset in this project actually lives, checked against reality.

This exists because the same failure happened repeatedly: run IDs were recorded
("run_20260821_132226") but never the resolvable path, so a later session -
mine or another - searched the wrong machine, found nothing, and reported the
data lost. It was on DRAGEN the whole time.

Two rules this enforces:

  1. A location is <server>:<absolute path>. A bare run ID is not a location.
  2. Every entry is verified before it is written down. An unverified path is
     marked MISSING, not quietly dropped.

Run it any time: python3 verify_locations.py
It rewrites DATA_LOCATIONS.md from what actually resolves.
"""
import io
import subprocess
import sys
import datetime

KEY = "C:/Users/T/.ssh/id_ed25519.ppk"

# server -> how to run a command there
HOSTS = {
    "DRAGEN": ["plink", "-ssh", "copilot@100.104.25.22", "-batch", "-noagent", "-i", KEY],
    "KIT": ["plink", "-batch", "-noagent", "-i", KEY, "-P", "2223", "toki@127.0.0.1"],
}

# what -> (server, path, what it is)
LOCATIONS = [
    # --- SINEderella runs -------------------------------------------------
    ("human SINEderella run", "DRAGEN",
     "/staging/tmp/sinederella_benchmark/runs/human/run_20260822_231705",
     "benchmark run, Alu subfamilies"),
    ("zebrafish SINEderella run", "DRAGEN",
     "/staging/tmp/sinederella_benchmark/runs/zebrafish/run_20260822_231705",
     "benchmark run"),
    ("Timema SINEderella run (timb benchmark)", "DRAGEN",
     "/staging/tmp/sinederella_benchmark/runs/timema/run_20260822_231705",
     "the AnnoSINE-55 benchmark set"),
    ("Timema curated v4 (tim, live)", "DRAGEN",
     "/staging/tmp/timema_sines/v4_curated/run_20260823_103133",
     "597 chunks, his refined 14-subfamily partition - ground truth"),
    ("Timema curated v3", "DRAGEN",
     "/staging/tmp/timema_sines/v3_curated/run_20260823_063444", ""),
    ("Timema run_20260821_132226", "DRAGEN",
     "/staging/tmp/timema_sines/run_20260821_132226",
     "598 chunks, 59,636 members, his 8-group partition"),
    ("Timema run_20260821_145119", "DRAGEN",
     "/staging/tmp/timema_sines/run_20260821_145119", ""),
    ("hedgehog (eri) SINEderella run", "KIT",
     "/data/W/toki/Genomes/Mammalia/Eulipotyphla/Erniacidae/run_20260820_221537",
     "on KIT, not DRAGEN - this is why searches on DRAGEN found nothing"),
    ("snail (pom) SINEderella run", "DRAGEN",
     "/staging/tmp/newsp/pom/run_20260901_131644", "both snail families"),

    # --- the peel-loop material -------------------------------------------
    ("Tal subfam chunks (toc)", "KIT",
     "/data/W/toki/Tal/genomes/toc/sf",
     "600 chunk consensuses + 600 member chunks of 50 = 30,000 copies"),
    ("Tal consensus bank", "KIT",
     "/data/W/toki/Tal/consi.bnk",
     "TAL + the 6 groups he peeled by hand - ground truth for the peel loop"),
    ("Tal starting query", "KIT", "/data/W/toki/Tal/Tal.q", "one consensus"),

    # --- his manual separations, the standard ------------------------------
    ("Tal manual subfamilies", "KIT",
     "/data/W/toki/SINE_disc/aln_c",
     "POS__saq__s1-s9, ccr g1-g7+a_ccr, teu t1-t6, dmo d1-d5"),
    ("labelled corpus (674 alignments)", "KIT",
     "/data/W/toki/SINE_disc/aln_c", "raw, 400 bp flanks"),
    ("labelled corpus, justified", "KIT",
     "/data/W/toki/SINE_disc/aln_v2", ""),
    ("human + Timema benchmark alignments", "KIT",
     "/data/W/toki/SINE_disc/bench_in", "64 HUM + 55 TIMB"),

    # --- tools -------------------------------------------------------------
    ("SINEderella toolkit", "DRAGEN", "/staging/tmp/SINEderella",
     "MANUAL.md, README.md, all step scripts, SubFam"),
    ("SINE_consensus repo", "DRAGEN", "/staging/tmp/SINE_consensus_repo",
     "sine_consensus.sh, sine_pairwise_consensus.sh"),
    ("conse", "DRAGEN", "/usr/bin/conse", "EMBOSS cons wrapper, 35 % plurality"),
    ("SubFam", "KIT", "/data/V/toki/bin/SubFam", "also /usr/local/bin/SubFam on DRAGEN"),
    ("his de novo scan", "KIT", "/data/W/toki/scorpio/denovo_scan/sine_scan.sh", ""),
    ("AnnoSINE env", "KIT", "/data/V/toki/envs/annosine", "built 2026-09-02"),
    ("AnnoSINE env", "DRAGEN", "/staging/tmp/kmer_experiment/envs/annosine", ""),
    ("discriminator working dir", "DRAGEN", "/staging/tmp/sinedisc", ""),
    ("discriminator working dir", "KIT", "/data/W/toki/SINE_disc", ""),

    # --- genomes -----------------------------------------------------------
    ("scorpion genomes (9 species)", "KIT",
     "/data/V/toki/Genomes/Scorpions", "one best assembly per species, 14.6 GB"),
    ("scorpion genomes (working copies)", "DRAGEN", "/staging/tmp/scorpions", ""),
    ("Timema genome", "KIT",
     "/data/W/toki/Genomes/lower/Arthropoda/Timema/timema.fna", "1.24 GB, bwa-indexed"),
    ("prospective genomes", "DRAGEN", "/staging/tmp/newsp",
     "hyd pom aca stu ska zeb + rebuilds"),
    ("RepBase blast db", "KIT", "/data/W/toki/SINE_disc/repbase.nhr", "blastn -db repbase"),
    ("SINEBase", "KIT", "/data/W/toki/SINEbase", ""),
]


def check(host, path):
    cmd = HOSTS[host] + ["test -e '%s' && echo YES || echo NO" % path]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        return "YES" in r.stdout
    except Exception:
        return None


def main():
    rows = []
    for name, host, path, note in LOCATIONS:
        ok = check(host, path)
        state = "ok" if ok else ("unreachable" if ok is None else "MISSING")
        rows.append((name, host, path, note, state))
        print("%-9s %-40s %s" % (state, name[:40], path))

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    out = ["# Where everything is",
           "",
           "*Verified %s by `verify_locations.py`. Re-run it rather than trusting age.*" % stamp,
           "",
           "A location is **server plus absolute path**. A bare run ID is not a location -",
           "recording run IDs without paths is exactly how the Timema runs were declared",
           "lost while sitting on DRAGEN.",
           "",
           "| what | server | path | state | note |",
           "|---|---|---|---|---|"]
    for name, host, path, note, state in rows:
        mark = {"ok": "ok", "MISSING": "**MISSING**"}.get(state, state)
        out.append("| %s | `%s` | `%s` | %s | %s |" % (name, host, path, mark, note))
    out += ["",
            "## Connecting",
            "",
            "- **DRAGEN** — `plink -ssh copilot@100.104.25.22 -batch -noagent -i <ppk>`;",
            "  strip the 9-line banner with `tail -n +10`; never write to `/home`.",
            "- **KIT** — one tunnel per session, then everything through it:",
            "  `plink -batch -i <ppk> -N -L 2223:localhost:22 toki@85.89.102.78`, then",
            "  `plink -batch -i <ppk> -P 2223 toki@127.0.0.1 '<cmd>'`. Never open a direct",
            "  connection per command (fail2ban ban, 2026-08-10).",
            "- `/data/W` on KIT is ~98 % full. Put new data on `/data/V`.",
            "",
            "## Searching for something not listed",
            "",
            "`find A -o B` without an explicit `-print` can return nothing while both",
            "match. Always write `find ... \\( -a -o -b \\) -print`, or run two finds.",
            "A negative result from a malformed search is not a negative result.",
            ""]
    io.open("DATA_LOCATIONS.md", "w", encoding="utf-8").write("\n".join(out))
    bad = [r for r in rows if r[4] != "ok"]
    print("\n%d locations, %d not ok" % (len(rows), len(bad)))


if __name__ == "__main__":
    main()
