# Results: benchmark harness calibration and cost model

Harness `csd/bench/`, driver `scripts/bench.py`, protocol `PROTOCOL.md` (predictions fixed before each run, see the
amendments). Machine: 2 vCPU Intel Xeon @ 2.1 GHz (PCLMULQDQ), gcc, Python 3.13, runs pinned to core 0 (`*.meta.json`
hold commit, hashes, load, governor). Single thread, minimum of 3 repetitions per row, 3 seeds per cell, arms A and B
for the noise floor. 228 `ok` rows (72 + 72 + 84) and 36 rows skipped as inapplicable; no failed row (`mismatch`, `invalid`, `crash`, `error`).

## What was established

1. The harness reproduces the recorded `m = 3` sweep once the solver wrappers use a small output buffer (C1b: ratio to
   the record 0.9 to 1.2 for both solvers at `k = 9..12`). Without that, LRF at `k <= 11` is 1.6 to 2.3 times the
   record because of an allocator artifact, not the algorithm (below).
2. Growth over `k = 9..12` (slope of `log2 t` per step of `k`): enumeration 2.09 [1.97, 2.14], LRF 1.01 [0.88, 1.11]
   (the record's own LRF slope on the same range: 0.95 [0.87, 1.09]). Speed-up of LRF over enumeration doubles with
   each step of `k`, as `ic-rootfind/README.md` says.
3. All solvers return the same decomposition sets on 12 of 12 planted instances (`m = 3, 4`; `py_ref`, `enum`, `lrf`,
   `mitm4`), and MutantXL closes at degree 7 on the three `n = 10, k = 4` instances (1.58 s median; the README has 1.73 s).
4. On this machine, differences below about 20% between two timings are not resolved (A/A, below).
5. Under the cost model, a relation finder beats Pollard rho on ECC2K-130 only if relations cost about 190 core-seconds
   or less each at `k` near 31 (table below). For `m = 3` that is 0.13 ns per target `R`, impossible; for `m = 4`
   it is 61 ms per target covering every 4-summand decomposition.

## Predictions and verdicts

| run | prediction | observed | verdict |
|---|---|---|---|
| C1 P1 | harness / record time in [0.67, 1.5] for `k >= 9` | enum 1.17, 1.14, 1.19, 0.91; LRF 1.97, 2.32, 1.56, 1.09 | enum confirmed; **LRF falsified** (`k = 9, 10, 11`) |
| C1 P2 | slopes enum 2.0, LRF 1.0, intervals contain them | enum 1.90 [1.86, 2.06]; LRF 0.66 [0.55, 0.76] | enum confirmed; **LRF falsified** |
| C1 P3 | A/A median in [0.95, 1.05]; 5 to 95% in [0.8, 1.25] above 1 ms | 0.988; 0.777 to 1.287 (19 pairs) | median confirmed; **band falsified** |
| C1 P4 | no failed rows | 72 of 72 `ok` | confirmed |
| C1 P5 | RSS varies by less than 10 MiB across `k` | 74.8 to 75.2 MiB | confirmed |
| C1b P1b | LRF in [0.4, 1.5], enum in [0.67, 1.5] (`k = 9..12`) | LRF 0.95, 1.03, 1.16, 1.02; enum 0.93, 0.91, 1.15, 1.13 | confirmed |
| C1b P6 | C1b/C1: LRF <= 0.6 at `k = 7..10`, >= 0.8 at 12 | 0.29, 0.35, 0.50, 0.43; 0.95 | confirmed |
| C1b P6 | C1b/C1: enum <= 0.9 at `k = 7, 8`; in [0.9, 1.1] at `k >= 10` | 0.72, 0.90; 0.80 (10), 0.90 (11), 1.24 (12) | `k = 7, 8` confirmed (8 at the limit); **`k = 10, 12` falsified** |
| C1b P2b | LRF slope in [0.85, 1.25]; enum about 1.9 | LRF 1.01 [0.88, 1.11]; enum 2.09 [1.97, 2.14] | confirmed (enum intervals overlap C1's) |
| C1b P3b | A/A 5th in [0.70, 0.85], 95th in [1.20, 1.40] above 1 ms | 0.888 to 1.081 (15 pairs) | **falsified** |
| C1b P4b, C2 | no failed rows; planted instances have a decomposition; sets agree | 72 of 72; 84 `ok` (36 skipped as inapplicable); 12 of 12; one `decs_sha` per instance | confirmed |

Reading the falsifications: the enumeration timings move by 20 to 25% between two runs of the same grid
(`k = 10` and `k = 12` above), which is the width of the A/A band, so the tight bands of P6 and P3b were not
supportable on this machine. The A/A band itself is estimated from 15 to 19 pairs and differs between two runs
(0.78 to 1.29, then 0.89 to 1.08). Pooled over C1 and C1b (34 pairs above 1 ms): median 0.991, 5 to 95% 0.83 to 1.19.
Use that as the resolution: timings closer than about 20% are not distinguished here.

## The allocator artifact (C1 versus C1b)

`csd.ic_rootfind.Solver.lrf/brute/mitm4` allocate `np.zeros(cap * 2 * m)` per call with `cap = 2^20`, 48 MB at
`m = 3`. At `k = 10` that allocation took 73 µs at its fastest and 2.3 ms at its median in a probe (the script is
not kept; numbers in `PROTOCOL.md`), against 233 µs for the raw kernel call; with `cap = 2^12` the whole call takes
145 µs (min) and 151 µs (median). Fork, CPU pinning and GC were ruled out by the same probe. The adapters now start at
`2^12` and retry at `2^16`, `2^20` on overflow; the library is unchanged. The effect is size-independent, so it
dominates the 0.1 to 1 ms searches at `k <= 11` and disappears by `k = 12`.

Calibration after the fix, speed-up of LRF over enumeration (paired by instance, geometric mean, 95% bootstrap CI;
the record column is the median of per-instance ratios in `scaling-m3.jsonl`):

| k | n | C1b | record |
|---|---|---|---|
| 7 | 19 | 6.2 [4.6, 8.2] | 6.4 |
| 8 | 22 | 16.1 [13.9, 18.3] | 12.9 |
| 9 | 25 | 35.1 [29.6, 39.5] | 37.5 |
| 10 | 28 | 65.7 [59.7, 71.6] | 65.7 |
| 11 | 31 | 134.8 [123.8, 144.3] | 145.2 |
| 12 | 34 | 313.5 [292.2, 357.3] | 292.8 |

LRF costs 182 to 418 ns per factor-base point at `k = 9..12` (median over seeds), the README's 0.2 to 0.6 µs.

## Cost model against Pollard rho (Run M; a model, not a measurement)

`python3 scripts/bench.py model --preset ecc2k130 --m M --from-records experiments/ic-rootfind/n131.jsonl`
(`model-ecc2k130-m{3,4,5}.md`). Per-target cost 1850 ns and 150 ns per addition are the `n = 131` record's; rho is
2^60.8 iterations (cofactor 4, negation and Frobenius), 9.6e3 core-years. The model reproduces the README's
1.0 core-year per scan of 2^44 targets, 3e13 core-years for the relation stage, and one third of a relation per
scan (`tests/test_bench.py`).

| m | best k | total | gap to rho | total with a free relation finder (floor) | break-even cost per relation at k = 31 | per target R, all decompositions |
|---|---|---|---|---|---|---|
| 3 | 44 | 2.9e13 core-years (2^69.6 s) | 2^31.5 | 2^43.0 s, **above** rho by 2^4.9 | 207 s | 0.13 ns (infeasible: forming `R` takes 150 ns) |
| 4 | 45 | 1.1e14 core-years (2^71.6 s) | 2^33.5 | 2^35.7 s, below rho by 2^2.5 | 187 s | 61 ms |
| 5 | 45 | 5.7e14 core-years (2^73.9 s) | 2^35.8 | 2^37.5 s, below rho by 2^0.7 | 167 s | 2.3e7 s (mu = 2^17) |

The last two columns are solver-agnostic: a Gröbner basis, a root finder or a SAT call is compared with them at its
own `(m, k)`. At `m = 4` the repository's LRF takes `|FB|^2` targets of 1.85 µs, 2^47 times over the 61 ms budget.
The floor is the linear-algebra stage (`2^(2k)` sparse operations, memory 2^49 bytes at `k = 44`) plus forming
targets; it is why `m = 3` loses even with a free relation finder.

Sensitivity (`model-sensitivity.txt`): the `m = 3` floor stays above rho for every tested `t_la` in [1, 100] ns and
`t_new` in [15, 1500] ns (it needs `t_la <= 0.008` ns); the `m = 4` floor is below rho for all of them; the `m = 5`
floor is below rho for `t_la <= 10` ns and above it at 100 ns. `t_la`, the cost of one sparse-matrix entry update, is an
assumption, not measured here.

## Limits

- One machine, 2 vCPUs, shared; the noise floor above is for this machine. Single thread only.
- The C kernel stops at `n = 63`; nothing here measures `n = 131` (the model takes it from the Rust record).
- No Gröbner (F4/F5, msolve, Magma), SAT or Rust-kernel adapter is registered yet; MutantXL is the only algebraic
  solver, and it applies to at most 12 Boolean variables. `PROTOCOL.md` and `csd/bench/solvers.py` show where to add one.
- The model assumes `|FB| = 2^k`, `mu = |FB|^m / (m! 2^n)`, the `V = <1, w, ..>` candidate count, and does not apply
  kernel limits (the Rust kernel needs `k >= 35` at `n = 131`).

## Reproduce

    python3 scripts/bench.py run --grid calib --seeds 3 --solvers enum,lrf --aa --pin 0 --out experiments/ic-bench/c1b-calib.jsonl
    python3 scripts/bench.py run --grid smoke --seeds 3 --solvers py_ref,enum,lrf,mitm4,mutantxl --aa --pin 0 --out experiments/ic-bench/c2-planted.jsonl
    python3 scripts/bench.py summarize experiments/ic-bench/c1b-calib.jsonl experiments/ic-bench/c1b-calib-k11-12.jsonl --baseline lrf
    python3 scripts/bench.py model --preset ecc2k130 --m 4 --from-records experiments/ic-rootfind/n131.jsonl
    python3 -m pytest -q tests/test_bench.py
