# SINE_discriminator — working rules

These are binding. They exist because each one was violated at least once.

## 1. Compute runs on the servers. Always.

Sergei's instruction from the first message of this project: *"use dragen or any
other server you need"*. It has not changed and does not need restating.

- **DRAGEN** `copilot@100.104.25.22` — 64 cores, 503 GB RAM, 4 TB free on
  `/staging`. Has python3 3.12 + numpy 2.3.2, mafft, blastn, makeblastdb,
  bedtools, samtools, bwa, and AnnoSINE_v2 at
  `/staging/tmp/kmer_experiment/envs/annosine/bin/`.
- **KIT** `toki@85.89.102.78` via the port-2223 tunnel — genomes under
  `/data/W/toki/Genomes/`, hg38 at `/usr/local/genomes/hg38.mfa`.
- The pipeline lives at `/staging/tmp/sinedisc/` on DRAGEN.

**The Windows box is for editing code, building pages, and reporting. Nothing
else.** It has no mafft, no blastn, no bedtools, no genomes.

Before writing any script, decide where it runs. If it shells out to a
bioinformatics tool, or touches a genome, or loops over more than a few dozen
alignments — it runs on a server, and it must be written that way from the
start, not written locally and then moved.

**Re-ask this every time the task changes.** The failure mode is drift: the rule
was honoured at the start of the project and then quietly abandoned when the
work shifted to scorer code, which is legitimately edited locally.

## 2. File transfer to DRAGEN

`pscp` is broken to this host — SFTP fails outright, and `-scp` mode silently
writes nothing. Use:

```bash
gzip -c file | base64 | tr -d '\n' > /tmp/x.b64
split -b 8000 /tmp/x.b64 /tmp/part_
for f in /tmp/part_*; do plink ... "cat >> /remote/x.b64" < "$f"; done
plink ... "tr -d '\n' < /remote/x.b64 | base64 -d | gunzip > /remote/file"
```

**8 KB chunks.** 60 KB chunks truncated silently at 19 KB of 80 KB and the
error surfaced only as a wrong result later.

Better: don't transfer. Keep data and code server-side and pull back only small
summaries.

## 3. Verify, don't assume

- `pgrep -f "script.py"` inside a `bash -c` containing that text **matches its
  own command line** and always returns ≥1. It reported a dead job as alive
  twice. Use `pgrep -f "script[.]py"`.
- DRAGEN's login banner is exactly 9 lines; strip with `tail -n +10`.
- `find /` hangs. Scope searches.
- A "done" claim from any agent, or a progress note saying an edit was made, is
  not evidence. Check the artifact.
- Before reporting a number, check the sample it came from. Small-sample
  inflation has produced two wrong findings in this project (a flank background
  of 1.000 from a few stub bases; 0.573 from a near-empty window).

## 4. Regenerate, don't hand-assemble

`_embed.js` was built by hand and silently drifted: the published page served
verdicts computed with criteria that had since been removed, and captions
asserting conclusions that had been retracted. Use `rebuild_embed.py` and
`rebuild_annot.py`. After any scorer change, regenerate before publishing.

## 5. Read the source of truth first

`hum/LOG.md`, `timb/LOG.md`, `tim/LOG.md` in the Tal repo record what was run,
on which genome, with which flank width, and where the outputs live. Two
sessions were spent re-deriving facts that were written down. `tim/` is the
manually curated ground truth and should have been the first thing scored, not
the last.

## 6. Standing scientific rules

- Element evidence from the short-flank alignment; flank decay from the long
  one. Never both from one geometry.
- An absent measurement is reported as absent, never scored as guilt. Three
  separate bugs were this same error.
- Subfamily ambiguity does not affect the sine/not-sine verdict.

## 7. Save after every exchange — his standing instruction

> "save these intermediate results in your internal documents. do it every chat
> message, every answer not to waste my time. important results save
> separately."

Three files, three purposes. Update the relevant one **before** replying, not at
the end of a session:

| file | holds |
|---|---|
| `FINDINGS.md` | **results that change the tool**, newest first — this is the "important results, separately" file |
| `calls.tsv` | his judgements, verbatim, one row per alignment |
| `HANDOFF.md` | working history, retractions, how things were done |

`PLAN.md` is what to do next; `METHOD.md` is his method in his words.

Do not batch this to the end of a session. He has lost work to that twice.


## 8. Why the rule keeps getting broken, and what removes the failure

He asked: *"why i have to catch you for no reasons decisions?"*

The honest answer: the decision gets made at the moment of action, and a reading
is found that lets the rule be skipped. Packing a tarball "isn't compute".
Running numpy on 72 local files "isn't heavy". Both were rationalisations made
in the second before acting, against a rule already written down.

Vigilance is not the fix. **Removing the decision is.** The data now lives where
the work happens:

**KIT** `/data/W/toki/SINE_disc/` holds `aln_c` (674, the source corpus),
`aln_v2` (673, regenerated from it — never transfer this, it is derivable),
`aln_ext` (64), and the pipeline: `verdict.py`, `measure_c.py`, `profiles.py`,
`flankdecay.py`, `annotate.py`, `prune.py`, `justify_all.py`, `trim_flanks.py`,
`test2_props.py`, `test3_combine.py`, `already_answered.py`, `calls.tsv`.
python3.12 + numpy 1.26.4, mafft, blastn, bedtools, samtools all present.

**DRAGEN** `/staging/tmp/sinedisc/` holds the same pipeline plus the genomes,
AnnoSINE_v2 and the SubFam chunk work.

**Derived data is regenerated on the server, never shipped.** `aln_v2` is
`justify_all.py` applied to `aln_c`; regenerating took one command. Shipping it
would have been 71 MB through a base64 channel.

If a step needs a file that is not on the server, the question is not "how do I
transfer it" but "why is it not there".
