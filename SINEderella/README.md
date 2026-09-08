# SINEderella publish pipeline (mirrored from working tree)

Scripts here wire the **Tal-style report alignment** path for any SINEderella
`run_*` directory. On DRAGEN, set `DISC=/staging/tmp/sinedisc` if you prefer
the server copy; from this repo, `DISC` defaults to the **repository root**
(parent of `SINEderella/`).

## Main entry

```bash
./SINEderella/run_publish_alignments.sh <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]
```

Chain:

1. `step7_boundary_refine.sh` (in `SINEDERELLA_BIN`) — flank widen
2. `needs_border_loop.py` + `border_loop_subfam.py` — element border loop
3. `step8a_extract_alignments.sh` — top100 / rand100 / subfam MSAs
4. `rebuild_consensus_row.py` — row 0 = copy-column majority at CONS_EDGE 0.50
5. `boundary_justify.py` — display justify
6. `trim_display_flanks.py` — drop abandoned end columns (`TRIM_DISPLAY_MODE`)

## Report helpers (repo root)

| script | role |
|---|---|
| `inject_oma_aln_section.py` | Patch alignment table into species report HTML |
| `parse_oma_copy_counts.py` | Copy counts from report composition section |
| `flank_uniqueness.py` | Per-side flank clustering / shared-context flags |
| `scan_oma_flank_uniqueness.py` | Corpus scan → TSV |
| `verdict.py` | Flanks / Element / Overall chips |

## Tests

From repo root:

```bash
python -m unittest discover -s tests -p "test_*.py"
```

## Oma-specific launchers

`launch_oma_*.sh` — one-off DRAGEN jobs for the scorpion peel run; kept as
reference. Override `RUN` and `DISC` when reusing on another species.
