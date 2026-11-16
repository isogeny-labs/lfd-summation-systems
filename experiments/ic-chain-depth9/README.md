# Two-summand chain systems at depth 9

Ten seeded systems (two controls at depth 7, eight at depth 9) over p = 7681, made with the depth-8 driver.
Every system has exact LFD 5 with strict falls at degrees 4 and 5 (`summary.json`).

| File | Content |
|---|---|
| `cases.jsonl`, `systems.jsonl` | the seeded cases and the distinct systems they give |
| `results.jsonl` | one exact result per system: closure profile, strict falls, LFD, certificate status |
| `summary.json` | counts by depth and by factor-base type, the engine used and its limitations |
| `verification.json` | SHA-256 of the sources and records of the run |
| `replay-20261108.json` | full replay: every recomputed result equals the stored one |

Engine: exact dense closure through degree 4, then the shift-structured closure (once per pivot) at degrees 5
and 6, with rank certificates. An exact value needs the certificate at degree 5; a certificate only at degree 6
gives an interval.

Replay (about 2.5 hours on one core):

    python scripts/run_ic_chain_depth8.py --check-only --depths 9 --primes 7681 \
      --out experiments/ic-chain-depth9 --workers 1 --product-cap 12000 \
      --shift-max-cells 450000000

On macOS this command fails. The driver sets `--product-cap` and `--shift-max-cells` inside `main()`, and worker
processes started with `spawn`, the macOS default, do not see them: they fall back to 16000 and 200,000,000, and
the first depth-9 system stops after about 1.5 hours with "Shift-closure matrix budget exceeded". Linux starts
workers with `fork`, which keeps the values, so there the command runs as written. On macOS, start the driver under `fork` instead (the driver itself is unchanged, since its hash is pinned in
`verification.json`):

    python -c "import multiprocessing, runpy, sys; multiprocessing.set_start_method('fork'); \
      sys.argv[0] = 'scripts/run_ic_chain_depth8.py'; runpy.run_path(sys.argv[0], run_name='__main__')" \
      --check-only --depths 9 --primes 7681 --out experiments/ic-chain-depth9 --workers 1 \
      --product-cap 12000 --shift-max-cells 450000000
