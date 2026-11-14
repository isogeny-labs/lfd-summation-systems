# Replication of the binary deep runs: plan (frozen before the run)

Written 2026-11-12 after the first run (`RESULTS.md`) and before any replication on a rented machine. The owner asked
for one more replication. The only replication records that exist at this point are two local smoke tests of the
wrapper (one k = 11 system on the simulated GPU device, one on the CPU device), both equal to the first run.

## Question

Are the 57 certified records of the first run reproduced by a second execution that shares the engine source but
not the execution path (batch, tile and work sizes; for the smaller systems also the device layer and the linear
algebra library)?

## What is run (one H200 pod, two chains side by side, `pod/repl_bin.py`)

| run | systems | device | settings | first run for comparison |
|---|---|---|---|---|
| r2 | all 57 of `cells.json`, same degree caps | GPU (CuPy), device pool off | `--batch 768 --tile-mb 512 --work-gb 12` | GPU, `--batch 1024 --tile-mb 1024 --work-gb 20` |
| r3 | the 46 systems of `cells_cpu.json` | CPU (numpy, OpenBLAS, 16 threads) | defaults | GPU |

`cells_cpu.json` (`make_repl_cells.py`) holds the systems whose closure at the first run's certificate degree has at
most 60,000 columns. That subset is chosen by cost. The 11 larger systems (the trace-zero cells at k = 12 and 13, one
at k = 11, the subfield cells, and (18, 13)) are repeated on the GPU only: they are the ones with the unexpected
values, and their replication varies the settings but not the device.

## Comparison and pass criterion (`compare_repl.py`)

For every system and every degree D from 0 to the certificate degree: `columns`, `rank`, `codim`, `profile`,
`contains_one`, `closed`, `unit`, `system_hash` and the canonical basis hash `rref_sha256` equal to the first run's
record; on the final record `exact_lfd`, `certificate_degree`, `strict_falls`, `quotient_ranks`, `local_checks`,
`n_points`, `points_sha256`, `status` equal.

- **Pass** for a run: every system of the run equal on every field above.
- **Fail**: any inequality. It is reported per system and field. Nothing is rerun to turn a fail into a pass; a
  failing system is investigated and the first run's record for it is marked as not replicated.
- A system that does not finish (time limit, deadline, pod loss) is reported as not replicated, not as a pass.

The canonical basis depends only on the subspace and the fixed column order, so an equal hash between two executions
with different batch and tile sizes (or devices) means the two computed the same subspace. Where the subspace is the
whole space or zero the hash carries no information; `compare_repl.py` counts the proper subspaces separately.
A negative control of the comparison (one hash altered in a copy) is reported as a problem.

## What this does not test

- A second implementation: both executions run `csd/ic_bounds_bin_big.py`. An error in the algorithm that does not
  depend on batch, tile or device (for example a wrong shift rule, or a wrong certificate criterion) would replicate.
  The guard against that remains the validation against the existing engine at k <= 10 (`validation/`).
- Other instance draws: the systems are the same 57.

## Machine and budget

One H200 pod, lifetime cap 3 hours (at most $15.87 at $5.29/h); the first run took 78 minutes. Results leave the pod
as two archives rewritten every two minutes and copied by the detached loop `mirror.sh`, which terminates the pod when
both chains have ended.
