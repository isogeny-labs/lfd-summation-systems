# The sweep (Sections 5 and 6 of the paper)

Exact last fall degrees of 3,005 distinct summation systems from 3,163 seeded cases: prime-field, binary
ordinary and Koblitz curves, with random and symmetric factor bases (seed 20261104).

| File | Content |
|---|---|
| `cases.jsonl`, `cases.csv` | every case: curve, field, factor base, target, depth |
| `systems.jsonl` | the distinct systems (cases that give the same system are merged) |
| `results.jsonl` | one exact result per system: profile, strict falls, LFD, local certificates, planted decomposition when present |
| `summary.json` | LFD histogram, counts by curve kind, the 36 groups of the paper's tables, skipped cases and limitations |
| `base-class-analysis.json` | the base-class analysis (`scripts/analyze_ic_base_classes.py`) |
| `verification.json` | SHA-256 of the sources and of every record |

Regenerate every case and replay every saved closure without writing files (about 40 minutes with 4 workers):

    python scripts/run_ic_curve_sweep.py --check-only --workers 4

The check first compares the hashes in `verification.json` with the current files, then recomputes each system
and compares it with `results.jsonl`. The sweep is a seeded, structured grid: it supports no population
inference, and it measures last fall degrees only, not relation yield or ECDLP cost (see `limitations` in
`summary.json`).
