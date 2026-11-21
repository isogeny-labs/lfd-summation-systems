# Protocol: benchmarking relation finders for index calculus (written before the runs below)

Harness: `csd/bench/`, driver `scripts/bench.py`, tests `tests/test_bench.py`.

## Question

How fast is each point-decomposition solver, how does the cost grow with the factor-base dimension `k`, and which
differences are larger than the noise of the machine? Separately: given a measured per-target cost, how far is the
whole attack from Pollard rho, and what per-target cost would break even (`bench.py model`)?

## Design

| Concern | Control |
|---|---|
| A wrong solver is fast | every decomposition re-added to `R` with the group law; all solvers on one instance compared with the lowest-`oracle_rank` solver that ran (`py_ref`, then `enum`, then `lrf`); a row that fails is `mismatch` / `invalid` and the driver exits non-zero |
| Agreement on the empty set proves nothing | `plant=True` builds `R` as a sum of `m` factor-base points, so at least one decomposition exists |
| Setup counted as search | `t_setup` (tables, descent) and `t_min` / `t_median` of `reps` runs of the search are separate; `t_first` keeps the cold run |
| Hung or crashing kernel | forked child per run, timeout enforced by the parent; status `timeout` / `crash` / `error` is a row, not an abort |
| Drift of the machine favours one solver | solver order within an instance is a seeded random permutation (`pos` recorded) |
| Noise floor unknown | `--aa`: every solver runs twice as arms A and B on the same instance; `summarize` prints the B/A spread, and ratios inside it are not resolved |
| Neighbours, frequency | load average before each run, CPU pinning (`--pin`), governor and CPU model in `<out>.meta.json` |
| Version drift | git commit and dirty flag, source hashes of the kernels and of the harness in the meta file |
| Peak memory | `VmHWM` of the child after resetting the watermark through `/proc/self/clear_refs` (falls back to `ru_maxrss`, recorded) |
| Replication unit | the instance (seed); intervals are percentile bootstraps over seeds, never over repetitions of one instance |

Instances reuse the random stream of `scripts/run_ic_rootfind.py`, so a `Spec` with the same `(n, k, seed, basis)`
rebuilds the instance of a record in `experiments/ic-rootfind/`.

## Run C1: calibration against the recorded sweep

`bench.py run --grid calib --solvers enum,lrf --aa --pin 0` : `m = 3`, `V = <1, w, .., w^(k-1)>`, `n = 3k - 2`,
`k = 7..10`, seeds 0 to 2, `reps = 3`; and the same instances through `--cell` for `k = 11, 12` (enumeration still
finishes). Reference: `experiments/ic-rootfind/scaling-m3.jsonl` (4-core Xeon, minimum of three runs, Python call
included in the small sizes).

Predictions, fixed before the run:

- **P1 (harness reproduces the record).** For `k >= 9`, the median over seeds of `t_min(harness) / t_lrf(record)` is
  in [0.67, 1.5] for LRF and for enumeration. Falsified outside that interval; then the harness (fork, GC off,
  pinning) or the machine differs from the recorded run, and the cause is to be found before any other number is used.
- **P2 (growth).** The fitted slope of `log2 t_min` against `k`, over `k = 9..12`: enumeration 2.0 (4x per step,
  95% interval containing 2), LRF 1.0 (2x per step, interval containing 1). Falsified if either interval excludes
  its value. LRF's per-target cost is flat in `k` (0.2 to 0.6 microseconds in the record) and `|FB|` doubles.
- **P3 (noise floor).** B/A median in [0.95, 1.05]; 5th to 95th percentile inside [0.8, 1.25] for times above 1 ms.
  Falsified otherwise; then differences below the measured band are reported as unresolved.
- **P4 (correctness).** Zero `mismatch`, `invalid`, `crash` or `error` rows in C1 and in C2.
- **P5 (memory).** Peak RSS of `lrf` and `enum` varies by less than 10 MiB across `k = 7..12` (interpreter and numpy
  dominate; the tables are kilobytes). Falsified otherwise.

## Run C2: correctness with decompositions present

`bench.py run --grid smoke --solvers py_ref,enum,lrf,mitm4 --aa` and
`--cell n=13,k=4,m=4,plant=1,seed=0..2`: every row must be `ok`, every planted instance must have `n_decs >= 1`,
and all four solvers must return the same set (checked by the harness, recorded as `decs_sha`).

## Run M: the cost model

`bench.py model --preset ecc2k130 --m 3 --from-records experiments/ic-rootfind/n131.jsonl` and `--m 4`, `--m 5`.
Model, not measurement. Assumptions are parameters: `|FB| = 2^k`, `mu = |FB|^m / (m! 2^n)`, `t_la_ns = 10`
(assumed cost of one sparse-matrix entry update), `t_new_ns` = one group addition. Check against the README's
figures is in `tests/test_bench.py::test_model_reproduces_the_readme_figures`.

## What this does not claim

No ECDLP speed claim. Timings are single-core wall clock on one machine; the C kernel stops at `n = 63`, so the
`n = 131` per-target cost comes from the Rust record, not from this harness. The cost model's linear-algebra constant
is an assumption, and its answer for `m = 4, 5` depends on it (see `RESULTS.md`, sensitivity).

## Amendment 1: outcome of C1 and the successor run C1b (written before C1b was run)

C1 was run as specified (records `c1-calib.jsonl`, `c1-calib-k11-12.jsonl`, meta files alongside; 72 rows, all `ok`).
Against the predictions above:

| | outcome | verdict |
|---|---|---|
| P1 | enumeration `k >= 9`: ratios 1.17, 1.14, 1.19, 0.91. LRF: 1.97, 2.32, 1.56, 1.09 for `k = 9..12` | **enumeration confirmed, LRF falsified** (`k = 9, 10, 11` above 1.5) |
| P2 | slope over `k = 9..12`: enumeration 1.90 [1.86, 2.06]; LRF 0.66 [0.55, 0.76] | enumeration confirmed, **LRF falsified** (interval excludes 1.0) |
| P3 | A/A median 0.993; 5th to 95th percentile 0.777 to 1.287 for times above 1 ms (19 pairs) | median confirmed, **band falsified** (wider than 0.8 to 1.25) |
| P4 | 72 of 72 rows `ok` | confirmed for C1 (C2 pending) |
| P5 | peak RSS 74.8 to 75.2 MiB for both solvers across `k = 7..12` | confirmed |

Diagnosis of the LRF falsification (diagnostic run, not a record: the script was `/tmp/iso.py`, numbers quoted
here): the in-process, forked, pinned and GC-off variants of the same call all agree, so fork, pinning and GC are
not the cause. `Solver.lrf` allocates a `2^20`-row output buffer per call (48 MB at `m = 3`) with `np.zeros`. That
allocation alone takes 73 microseconds at its fastest and 2.3 ms at its median (glibc serves the first one by
`mmap`, later ones from the heap with a 48 MB memset), against a raw `rf_lrf` call of 233 microseconds at `k = 10`.
The recorded minimum of three runs caught the fast path (185 microseconds); a median does not. With `cap = 2^12`
the call takes 145 microseconds (min) and 151 (median). So at `k <= 11` the recorded and the C1 LRF timings both
contain a size-independent allocator artifact, and the speed-ups of `scaling-m3.jsonl` at `k = 9..11` are
therefore not clean: they are understated or noisy by up to 10x at the median. This does not touch the algorithm
or the claims at `k >= 12` (ratio 1.09 there).

Remedy in the harness only (the library is untouched): the solver adapters call with `cap = 2^12` and retry at
`2^16`, `2^20` on overflow (`csd/bench/solvers.py::_with_cap`); the cap used is in each row's counters.

### Run C1b (same grids, same seeds, same flags as C1; records `c1b-calib.jsonl`, `c1b-calib-k11-12.jsonl`)

Predictions, fixed before C1b (informed by the diagnosis above, so P6 is a check of the diagnosis, not a blind
test):

- **P1b.** `k = 9..12`: LRF median ratio harness/record in [0.4, 1.5]; enumeration in [0.67, 1.5].
- **P6.** Cell medians of `t_min(C1b) / t_min(C1)`: LRF at most 0.6 for `k = 7..10` and at least 0.8 at `k = 12`;
  enumeration at most 0.9 for `k = 7, 8` and in [0.9, 1.1] for `k >= 10`.
- **P2b.** LRF slope over `k = 9..12` in [0.85, 1.25]; enumeration unchanged (about 1.9).
- **P3b.** A/A band replicates C1: 5th percentile in [0.70, 0.85], 95th in [1.20, 1.40] for times above 1 ms.
- **P4b.** All rows `ok`.

## Amendment 2: outcome of C1b and C2, and a retraction (written after the runs; the predictions above are unchanged)

C1b (72 rows) and C2 (84 `ok` rows) are in `RESULTS.md` with every verdict. Two corrections to Amendment 1:

1. **Retraction.** Amendment 1 said the recorded speed-ups of `scaling-m3.jsonl` at `k = 9..11` are "understated or
   noisy by up to 10x at the median". The record is a minimum of three runs and C1b reproduces it
   (speed-ups 35, 66, 135 against 38, 66, 145). The allocator artifact biases measurements that are medians, or
   repeated calls in a long-lived process, or runs whose allocator state is not controlled; it does not
   invalidate the record. What C1 shows is only that the harness, run without the fix, did not reproduce the
   record at `k <= 11`.
2. **Two predictions of C1b were falsified** (P3b, and the enumeration part of P6); see `RESULTS.md`.
