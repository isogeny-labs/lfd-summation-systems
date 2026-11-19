# A faster instrument for the HKY fall-degree closure (`hky_fp_tower_fast`, `--engine hky2`)

Date: 2026-11-17. Branch `wt/hky-speed` (local commits on origin/main `94479f7e3e`, not pushed).

Every number below is `[pcc:unverified]`: no pcc claim was filed, nothing here is a registered
study result, and no leaderboard or `benchmarks/README.md` row is claimed. Labels used:
**measured** (a command of this directory produced it on the machine of section 1),
**given** (taken from the task or from the other repository, not re-measured here),
**estimate** (a model, with its fit points), **proof sketch** (an argument, not machine-checked).

## 0. Summary

- The reference engine `src/cryptanalysis/hky_fp_tower.rs` is unchanged. The new engine
  `src/cryptanalysis/hky_fp_tower_fast.rs` computes the same closures `V_{F,c}` and the same
  `dims`, `one`, `falls`, first and last fall degree. Selected by `pkm_tower_pilot --engine hky2`.
- **Measured speed-up over the reference, same cells, same seeds, same machine (4 vCPU):**
  t = 8: 10.27 s → 0.304 s (34×); t = 9: 38.96 s → 1.507 s (26×); **t = 10 (60,460 columns,
  the largest depth run with both engines): 999.9 s → 42.9 s (23×)**. Section 4.
- **Runs beyond the reference's reach on this machine (measured, new engine only):**
  t = 11 in 241 s; t = 12 (190,051 columns) in 1,732 s, 1 enters at c = 6 as the reference found
  on rented pods in 6,857 s and 7,269 s on 16 vCPU (**given**); **t = 13 (313,912 columns) in 11,788 s (3 h 16 min) with 9.6 GB: 1 enters at c = 6**, the value the other repository's engine gives; the reference timed out at 14,400 s on 16 vCPU on one draw and aborted on the other (**given**). Section 8.
- Where the time went in the reference (section 3) is **not** what was assumed. The stored rows
  are 6 % to 8 % dense at t = 8..10, not dense. The three costs are: a schedule (batches of 4,096
  products, last in first out) that reduces 4 to 6 products per pivot found (t = 8, 9) and makes
  denser rows; a sequential second pass inside every batch (41 % of the wall time at t = 9 on 4
  threads, so most of it on 16); and a division per multiply-add.
- What paid, in the order measured (t = 9, section 4): delayed reduction 1.9×; no heap 1.6×;
  products formed in the parallel section 1.3×; dense rows with a SIMD loop 1.24×; the
  summand-swap split 1.15× under the reference's schedule and 2.15× under the new one; **the
  schedule (sparsest rows first, small batches) 2.5× alone**. What did not pay: dense rows and SIMD beyond that 1.24× (the
  threshold does not matter, AVX-512 against the portable loop is 6 %), lead-only reduction
  (3 times more work), low-degree-first order, finer priority buckets. Section 6.
- The degree lift of the other repository carries over to the Kummer tower ring of the study
  cell (same ring up to a diagonal change of basis; section 9). It is **not implemented** here.
- Not available in this container and therefore not used: Conductor, cairn, research-kb, ecdb,
  AWS / EKS / isolab, `perf` and any hardware counter, `tools/isolated_bench.py` isolation (another
  job ran on the same 4 vCPU throughout). No cloud job, pod or paid resource was started.

## 1. Machine, build, commands

| | |
|:--|:--|
| CPU | Intel Xeon (model name `Intel(R) Xeon(R) Processor @ 2.10GHz`), 4 vCPU, 1 thread per core, `avx2`, `avx512f/dq/ifma`, `bmi2`; L2 8 MiB total, L3 260 MiB |
| memory | 16 GiB (`MemTotal 16473764 kB`) |
| kernel, toolchain | Linux 6.18.44-fc-v553; rustc 1.97.0 (2d8144b78 2026-07-07), `cargo build --release` (opt-level 3, no `target-cpu` flag) |
| contention | another job (worktree `wt/optlevel`) ran on the same vCPUs during every measurement; the 1-minute load before and after each run is in `raw/*/conditions.log` (range 1.4 to 3.4 during the A/B of section 4; one vCPU-equivalent or more was often taken) |
| counters | none: no `perf`, no PMU. **Cycle and instruction counts are missing.** The deterministic unit is the engine's `muladds` (one per entry of a pivot row subtracted, plus one per term of a product) and `rows_reduced` |
| threads | rayon default (4) unless stated |

**Build.** The root crate does not compile at `94479f7e3e`: `crates/hashed-ecdl-tensor/src/lib.rs`
does not declare its `metrics` module, and `src/cryptanalysis/sat.rs`, `koblitz_bench.rs`,
`ecbench/workload.rs`, `src/ml_dataset/custody.rs` and `factor_bases.rs` do not parse
(`rustfmt --edition 2021 --check <file>`: unclosed delimiters, a duplicated match arm). This was
not repaired (out of scope). The pilot is built by `harness/`, a crate that compiles only the
modules the pilot names, **from the repository's own files** (`#[path]` and two symlinks), with
stubs for what those files name in modules that reach the broken ones:

- `jv_bitops` (opt-in counters, always off here as by default), `groebner_f4` (`Ordering`,
  `cmp_monomial`, verbatim lines 65-104), `glv_invariant_base` (one function, unreachable from the
  closure), `pq_groebner_f2::Mask`;
- `ic_boundary`: `PrimePoint`, `PrimeCurve`, `GroupOps`, `CountedGroup` with the function bodies
  copied verbatim (line ranges in the file); only the field-operation tally is dropped.

That the harness builds the study's systems was checked against committed values: at t = 12 it
draws `x_r = 150938225` and `1002906059`, the targets of `RESULTS_T12_BRACKET.md` section 2
(measured: `pilot_baseline … --t-min 12 --max-t 12 --hky-cmax 3`).

The sparse checkout lacks five source paths Cargo needs for the root manifest; for the first
(failed) build attempt they were written from `git archive HEAD` into the worktree (no
sparse-checkout setting was changed; they are tracked files with unchanged content).

```bash
export CARGO_TARGET_DIR=/home/user/crypto-wt/hky-speed/target          # 0.8 GB after all builds
cd research/hky_engine_speed_20261117/harness
cargo build --release --example pkm_tower_pilot --example hky_fp_tower_check --example hky_engine_equality
(cd unit && cargo test --release)                                      # 10 tests, both engines
# frozen binaries (bin/ is not committed; sha256 below)
#   bin/pilot_baseline      built from 94479f7e3e before any change            379e7ac164a527ec…
#   bin/pilot_candidate     built from this branch                              47d6c281ecf9cc85…
#   bin/hky_engine_equality                                                     dbc2926e89ced7f4…
cd .. && ./run_gate.sh raw/gate-default                                # section 5
./run_ab.sh raw/ab-t6-9 6 9 3                                          # section 4
./run_ab.sh raw/ab-t10 10 10 1 hky2 hkyA && ./run_ab.sh raw/ab-t10b 10 10 2 hky2
```

The study's cell is `--p 2013265921 --kinds kummer --m 2 --controls tower --planted 0 --random 2`
with `--hky-cmax 7`; `run_ab.sh` uses `--random 1` (draw 0), the gate both draws.

## 2. Baseline (reference engine, frozen binary) — measured

`bin/pilot_baseline --engine hky`, draw 0, 4 threads. Wall time of the whole cell (all closures
c = 2 … last): median of 6 runs for t = 6..9 (the two A/A arms `hkyA`, `hkyB` of `raw/ab-t6-9`,
3 rounds each, interleaved with the other arms), one run for t = 10 (`raw/ab-t10`).

| t | closures run | last c | columns | 1 enters | wall s median (min – max) | muladds of the last closure | rows reduced | stored entries | entries / (n²/2) |
|--:|:--|--:|--:|:--|--:|--:|--:|--:|--:|
| 6 | 2..5 | 5 | 1,586 | yes | 0.264 (0.244 – 0.278) | 95,455,207 | 9,517 | 132,179 | 0.105 |
| 7 | 2..5 | 5 | 3,473 | yes | 1.505 (1.453 – 1.518) | 882,098,928 | 20,581 | 504,805 | 0.084 |
| 8 | 2..5 | 5 | 6,885 | yes | 10.27 (9.93 – 10.49) | 7,233,276,712 | 40,193 | 1,786,892 | 0.075 |
| 9 | 2..5 | 5 | 12,616 | yes | 38.67 (37.88 – 39.15) | 25,900,021,077 | 55,789 | 6,237,999 | 0.078 |
| 10 | 2..6 | 6 | 60,460 | yes | 999.9 (one run) | 597,217,990,555 | 98,429 | 103,886,880 | 0.057 |

A/A spread (same binary, two arms): medians differ by 0.1 % to 2.2 % (t = 6..9). Ratios in section 4 use the arm `hkyA` (t = 9: 38.96 s). An earlier pair of
t = 10 runs made while builds and sweeps ran beside them took 1,391 s and 1,354 s
(`raw/baseline_contended_t9_t10.jsonl`); they are not used for any ratio. t = 11 was not run with
the reference (it would take hours here).

Dimensions the reference reports (`dims[d] = dim V_{F,c} ∩ R_{≤d}`), identical on both draws:

| t | c = 4 | c = 5 | c = 6 |
|--:|:--|:--|:--|
| 6 | 0, 0, 2, 111, 605 (of 794) | whole space, 1 enters | |
| 7 | 0, 0, 2, 115, 829 (of 1,471) | whole space, 1 enters | |
| 8 | 0, 0, 2, 119, 1061 (of 2,517) | whole space, 1 enters | |
| 9 | 0, 0, 2, 123, 1301 (of 4,048) | whole space, 1 enters | |
| 10 | 0, 0, 2, 127, 1549 (of 6,196) | 0, 0, 2, 138, 3163, 18569 (of 21,700) | whole space, 1 enters |

These are the values of the other repository's engine for its chain ring (`LIFT.md` section 5:
605, 829, 1061, 1586, 3473, 6885, 12616, 18569, 60460), which is what section 9 explains.

## 3. Profile of the reference — measured (wall-clock phase timers; no cycle counts)

Level 0 of the new engine is a port of the reference with phase timers. That it is the same
algorithm is shown by its counters: `muladds`, `rows_reduced` and `entries` equal the reference's
exactly on every cell (unit test `level_zero_counts_what_the_reference_counts`; the `L0` and
`hkyA` rows of `raw/ab-t6-9/times.md`), and its wall time is within the A/A spread (t = 9:
39.58 s against 38.96 s).

t = 9, closure c = 5, 4 threads (`raw/ab-t6-9/L0.r2.log`):

| phase | wall s | share |
|:--|--:|--:|
| parallel reduction of each batch against the echelon as it stood | 22.65 | 58 % |
| **sequential** reduction of each batch row against the pivots the batch added, products, insertion (one thread) | 16.10 | 41 % |
| of which forming the products (CPU s, one thread) | 3.26 | 8 % |

Findings:

1. **The sequential pass is 41 % of the wall time on 4 threads.** It does not shrink with more
   threads, so on the 16-vCPU pods it was the larger part (if the parallel part scales, 16.1 s of
   21.8 s). The reference was mostly single-threaded where it mattered.
2. **The rows are not dense.** Stored entries are 5.7 % to 10.5 % of the triangle n²/2 (table of
   section 2). With hybrid storage at t = 9 only 1,922 of 12,616 rows are more than one third
   dense. The dense-kernel idea is worth 1.24× in total (section 4), not an order of magnitude.
3. **Most reductions are redundant.** At t = 9 the reference reduces 55,789 products to find
   12,616 pivots. With batches of 4,096 taken from the end of the queue, thousands of products are
   reduced against an echelon that does not yet contain the rows of their own batch.
4. `1` is the **last** pivot found: when the reference stops, the echelon has all n pivots
   (`rank = [12616]`). The early exit on `1` exists and is sound but saves nothing with that
   order.
5. The division: one `% p` per multiply-add. Removing it alone is 1.9× (level 1).

## 4. Each change, with its A/B — measured

Arms (all on the same binary `bin/pilot_candidate` except the reference), each adding one thing:

| arm | what it adds |
|:--|:--|
| `hkyA`, `hkyB` | reference, frozen binary (A/A) |
| `L0` | port of the reference in the new module (`--hky2-level 0 --hky2-order lifo --hky2-sym off`) |
| `L1` | delayed reduction: the accumulator stays below p², one conditional subtraction of p² per update, one Barrett reduction per column visited |
| `L2` | no heap, no touched list: the accumulator is scanned upward; a row that meets no pivot added since it was reduced is not reduced twice |
| `L3` | products `x·g` formed when dequeued, in the parallel section; the queue holds (row, variable) pairs; multiplicative hash for the column index |
| `L4` | rows more than 1/3 dense stored as `u32` slices, subtracted by a vectorised loop (AVX-512 here by run-time detection) |
| `L4sym` | `L4` + the summand-swap split (proved per system, section 7) |
| `L4sched` | `L4` + the schedule: sparsest source rows first, batch 32 |
| `hky2` | the default: `L4` + split + schedule |

`./run_ab.sh raw/ab-t6-9 6 9 3`: 3 rounds, arm order rotated each round, draw 0, 4 threads.
Wall milliseconds of the whole cell, median (min – max); `muladds` and `rows reduced` of the last
closure (deterministic, equal in the three rounds). Full table: `raw/ab-t6-9/times.md`.

| arm | t = 7 ms | t = 8 ms | t = 9 ms | t = 9 muladds | t = 9 rows reduced | t = 9 step ratio | t = 9 against reference |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `hkyA` | 1477 (1453 – 1503) | 10273 (10130 – 10489) | 38963 (38381 – 39002) | 25,900,021,077 | 55,789 | | 1.00 |
| `hkyB` | 1508 (1507 – 1518) | 10259 (9934 – 10328) | 38097 (37881 – 39149) | 25,900,021,077 | 55,789 | | 1.02 |
| `L0` | 1508 (1469 – 1559) | 10193 (10180 – 10540) | 39580 (38793 – 40178) | 25,900,021,077 | 55,789 | 0.98 | 0.98 |
| `L1` | 912 (904 – 925) | 5563 (5561 – 5599) | 20499 (19414 – 20546) | 25,900,021,077 | 55,789 | 1.93 | 1.90 |
| `L2` | 606 (593 – 629) | 3500 (3417 – 3890) | 12738 (12613 – 12798) | 25,899,471,221 | 55,789 | 1.61 | 3.06 |
| `L3` | 390 (382 – 395) | 2702 (2692 – 2715) | 9939 (9367 – 9970) | 25,882,773,684 | 55,789 | 1.28 | 3.92 |
| `L4` | 320 (303 – 336) | 2160 (2142 – 2265) | 8024 (7868 – 8389) | 25,882,773,684 | 55,789 | 1.24 | 4.86 |
| `L4sym` | 239 (227 – 240) | 1454 (1422 – 1457) | 6962 (6837 – 6983) | 19,404,396,897 | 34,177 | 1.15 (on `L4`) | 5.60 |
| `L4sched` | 120 (117 – 124) | 539 (505 – 555) | 3238 (3182 – 3655) | 9,532,473,315 | 39,389 | 2.48 (on `L4`) | 12.0 |
| `hky2` | 82 (72 – 95) | 304 (279 – 331) | 1507 (1457 – 1522) | 3,378,346,764 | 19,522 | 2.15 (on `L4sched`) | **25.9** |

t = 10, `raw/ab-t10`, `raw/ab-t10b` (order: `hky2`, reference, `hky2`, `hky2`; load 2.7 to 5.7):

| arm | runs | wall s median (min – max) | muladds, closure c = 6 | rows reduced | ratio |
|:--|--:|--:|--:|--:|--:|
| reference | 1 | 999.9 | 597,217,990,555 | 98,429 | 1.00 |
| `hky2` | 3 | 42.9 (41.9 – 45.1) | 92,397,155,280 | 148,887 | **23.3** |

The t = 10 ratio has one reference run; with the contended pair (1,354 s, 1,391 s) it would read
32, which is why those are not used. Every arm of both tables passed the equality comparison
against `hkyA` on the same files (`hky_engine_equality compare`: 27 of 27 `EQUAL` for t = 6..9,
1 of 1 for t = 10) and the full gate of section 5 on the same binary.

**Threads** (`raw/scaling`, t = 9, `hky2`, 2 runs each): 1 thread 4.91 s, 2 threads 2.60 s
(1.89×), 4 threads 1.50 s (3.27×), with the other job on the machine. The sequential share is
6 % of the wall time at t = 9 (`par=1.35s seq=0.09s`) and 1.6 % at t = 11 (3.7 s of 231 s).
Nothing above 4 threads was measured.

## 5. Equality gate — measured

`run_gate.sh OUT [hky2 flags]` runs the frozen reference (`--engine hky`) and the candidate
(`--engine hky2`) on the same cells and compares the rows with `hky_engine_equality compare`
(Rust, `examples/hky_engine_equality.rs`): the same cells on both sides; the same system (`x_r`,
curve, tower, generator degrees, `planted_ok`); for every degree c the same `cols`, `dims`, `one`
and `censored`; the same `falls`, `first_fall_degree`, `last_fall_degree`, `censored_at`,
`timed_out`. Exit 1 on any difference.

| set | cells | what it exercises | reference outcome |
|:--|--:|:--|:--|
| `m2-tower-t1-8` | 16 | study cell, t = 1..8, both draws, c up to the last fall degree | d_F = 3, 4, 5 |
| `m2-tower-planted` | 12 | planted targets (consistent systems), t = 2..7, c = 2..5 | 1 in no closure |
| `m2-null-t2-7` | 18 | null control (random coefficients; **not symmetric**: split refused, 70 closures) | mixed |
| `m2-naive-t2-4` | 6 | naive control: no tower, free variables only | mixed |
| `m3-tower-t2-3` | 12 | m = 3: three towers and a free chain unknown; split accepted | mixed |
| `m4-tower-t2` | 2 | m = 4: four towers, two free unknowns | 1 in no closure up to c = 4 |
| `p0-kummer-isogeny` | 36 | p = 786433, Kummer and isogeny towers (general rule `y² = a·y·y′ + b·y′ + c·y + d`), t = 1..6 | mixed |
| `censor-max-cols` | 2 | `--hky-max-cols 1500`: the same closure is censored for `max_cols` in both | censored at c = 5 (t = 7) |
| t = 6..9 A/B | 4 × 27 files | section 4 | |
| t = 9, 10 (`raw/gate-t9-10`) | 4 + 1 | both draws at t = 9, 10 against the reference rows of `raw/baseline_contended_t9_t10.jsonl` (contention changes time, not rows); the clean t = 10 pair of section 4 | d_F = 5, 6 |

Result: **every set `EQUAL` (104 cells in the script, plus the A/B files) for the default and for
each of** `L0`, `L1`, `L2`, `L3`, `L4`, `L4sym`, low-degree-first with batch 64, lead-only, and the
portable loop with dense threshold 1/8 (`raw/gate-*`). m ≥ 3 cells are small (t ≤ 3); the generic
path with free variables is also covered by the unit test `free_variables_and_three_blocks`.

**What the comparison can and cannot see.** On a closure that contains 1 the reported `dims` are
those of the whole space whatever the arithmetic did before 1 was found. So a fault in the last,
largest closure of an inconsistent cell is invisible to `dims`; it can only show as 1 entering at
the wrong degree. The power of the gate comes from closures **without** 1: the planted set (1 in
no closure, every closure a proper subspace up to c = 5), the c = 4 closures at t ≥ 6 and c = 5 at
t ≥ 10. This was found the hard way: the first negative control (below) passed the comparison on
an inconsistent system.

**Negative controls** (test hooks `--hky2-fault`, never set in a measurement), compared against
the planted set with `compare --expect-mismatch`:

| fault | result |
|:--|:--|
| 1: products by variable 0 not formed | 60 differences over 12 cells: detected |
| 2: the delayed reduction subtracts p² − 1 instead of p² | 40 differences over 12 cells: detected |
| 2 at level 0 | not detected, correctly: level 0 has no delayed reduction, the hook is inert there (`raw/gate-L0`) |
| 2 on an inconsistent Kummer system (t = 4, unit test, first version) | **not detected**: the reason for the paragraph above |

The suggested control "perturb one generator coefficient" was not used as the gate's control: a
generic change of a coefficient gives another system with the same dimension profile.

**Censoring.** `max_cols` is checked on the monomial count before anything is built, identically.
`max_entries` and `budget` limit an instrument's own storage and time; the new engine stores fewer
entries and takes less time, so it is censored later than the reference or not at all. No gate
cell is censored by them in either engine.

## 6. Schedule, and what was tried and rejected — measured

The space does not depend on the order in which the queue is served. The work does, by a factor of
several, and not monotonically. `muladds` of one closure (deterministic), split, level 4:

| order (split unless noted) | batch | t = 8, c = 5 (1 enters) | t = 9, c = 5 (1 enters) | t = 10, c = 5 (fixpoint, every product reduced) |
|:--|--:|--:|--:|--:|
| reference schedule: last in first out, **unsplit** | 4096 | 7.23e9 | 2.59e10 | 6.29e10 |
| last in first out | 4096 | 3.34e9 | 1.94e10 | 2.67e10 |
| last in first out | 256 | 9.58e8 | 8.40e9 | 6.17e10 |
| last in first out | 64 | 1.51e9 | 6.65e9 | 5.48e10 |
| last in first out | 16 | 7.75e8 | 2.35e10 | 1.42e11 |
| lowest degree first | 64 | 1.00e9 | 1.08e10 | 7.98e10 |
| **sparsest source row first** (one bucket per doubling of the row's terms) — the default | 32 | **4.47e8** | **3.38e9** | **2.38e10** |
| sparsest first | 128 | 5.35e8 | 3.54e9 | 2.32e10 |
| sparsest first | 4096 | 1.81e9 | 9.08e9 | 2.36e10 |
| sparsest first, **unsplit** | 32 | 1.12e9 | 9.53e9 | 8.19e10 |
| sparsest first, one bucket per 4× / 8× / 16× | 32 | 8.86e8 / 7.54e8 / 1.59e9 | 5.81e9 / 8.16e9 / 7.22e9 | 3.50e10 / 4.35e10 / 4.72e10 |
| sparsest first, 8 buckets per doubling (interactive run, code since removed, rows not kept) | 32 | 6.21e8 | 4.47e9 | |
| last in first out, lead-only reduction | 64 | 1.93e9 | 2.27e10 | 1.86e11 |

`./run_schedule_sweep.sh raw/schedule` regenerates every row but one (`raw/schedule.out`, level 4).

Reading: a stored row is reduced against the pivots that existed when it entered, so the order of
insertion decides how dense the rows are (stored words at t = 10, c = 5: 1.01e7 to 2.62e7 across
the table), and stale large batches reduce products that turn out redundant. "Sparsest first" is
the only order that was at or near the best in all three columns and insensitive to the batch
size. It is not uniformly better: unsplit, on the fixpoint closure at t = 10, it does 1.3 times
more work than the reference's schedule (8.19e10 against 6.29e10). The default batch is `max(32, 2 × threads)`. This is a heuristic chosen on t = 8..10; nothing
proves it is good at t = 13.

Rejected, with numbers:

| idea | measurement | verdict |
|:--|:--|:--|
| lead-only reduction (`--hky2-lead-only on`: stop at the first non-pivot column, keep the tail unreduced, so a product with a new leading monomial enters as it is) | split, last in first out, batch 64: t = 8: 1.93e9 against 1.51e9 muladds; t = 9: 2.27e10 against 6.65e9 (stored words 2.97e7 against 7.1e6); t = 10, c = 5: 1.86e11 against 5.48e10 | 1.3 to 3.4 times more work: the unreduced tails are denser, not sparser. Kept as an option, off |
| dense rows as the main lever (idea 1/2 of the task) | threshold 1/2, 1/3, 1/5, 9/10 at t = 9: 1454, 1502, 1474, 1519 ms; AVX-512 1502 ms, AVX2 1541 ms, portable 1601 ms (`raw/scaling/times.md`) | inside the noise; SIMD is worth 6 %. 17 % of the rows are dense at t = 9, 6 % at t = 11 |
| bitmap of touched columns instead of the heap | single thread, t = 8, reference schedule, machine loaded (interactive runs, rows not kept): 16.4 s against 15.2 s for the heap with delayed reduction | no gain (probably the read-modify-write chain on the bitmap word; no counters to confirm); replaced by the plain scan of the accumulator, 10.2 s against 15.7 s in the next pair |
| low-degree-first order, to find 1 early | table above | worse than sparsest first |
| "rank equals the number of columns" exit | | cannot fire before 1 is a pivot (1 is a column); the exit on 1 covers it |

## 7. The summand-swap symmetry — proof at run time, measured effect

For m = 2 the system the pilot builds is `S₃(y_{1,0}, y_{2,0}, x_R)` in the quotient by two copies of
the same tower. The swap σ of the two towers is a symmetry of **this** system: checked, not
assumed. `check_symmetry` compares, with no sampling,

1. that σ is an involution of the tower variables that maps every tower rule to the rule of the
   image variable (`next` mapped, a, b, c, d equal), so σ is a degree-preserving automorphism of
   the quotient ring that permutes normal-form monomials;
2. that every generator in normal form is mapped to itself or to its negative, coefficient by
   coefficient.

If either fails the engine runs unsplit and records `symmetry: "rejected: <why>"`. In the gate:
accepted on all Kummer and isogeny tower cells with m = 2, 3, 4 (the first chain link `S₃(x₁, x₂, u₁)`
is symmetric in the first two summands); **refused on all 70 null-control closures** ("generator 0
is not mapped to ±itself": the null control has random coefficients), which then take the unsplit
path and agree with the reference. Unit test `a_false_symmetry_is_rejected_and_the_answer_is_unchanged`
covers a broken generator coefficient and unequal towers.

With σ proved and p odd, `V_{F,c} = V⁺ ⊕ V⁻` (the module documentation gives the argument: the
least closed space containing a σ-stable set is σ-stable; for `g ∈ V^s` and a swapped pair of
variables only one is applied and the product is projected to both parts). Columns per part at
t = 13, c = 6: 157,145 and 156,767 of 313,912.

Effect (section 4 and 6): the multiply-adds fall by 2.2× under the reference schedule at t = 8
(7.23e9 → 3.34e9), 1.3× at t = 9, and by 2.5× to 3.4× under the sparse-first schedule (t = 8:
1.12e9 → 4.47e8; t = 9: 9.53e9 → 3.38e9; t = 10, c = 5: 8.19e10 → 2.38e10). **Not the 4× of dense elimination**: the rows are
sparse, so the cost is in the number of terms, not in n². Stored entries fall by 1.6× (t = 9, reference → `hky2`:
6.24e6 → 3.83e6).

## 8. Beyond the reference: t = 11, 12, 13 — measured, new engine only

`bin/pilot_candidate --p 2013265921 --ladder-t none --planted 0 --random 1 --engine hky2
--hky-cmax 6 --budget … --hky-max-cols 2000000 --hky-max-entries 14000000000 --kinds kummer --m 2
--controls tower --t-min T --max-t T --hky-progress raw/tT/progress.jsonl --out raw/tT/hky2.jsonl`
(t = 11 with `--hky-cmax 7`), draw 0 (`x_r` = 150938225 at t = 12, the study's draw 0), 4 threads,
the neighbour job running. These are instrument runs, not study rows (section 11.6).

| t | c | columns | part columns | dims / rank | 1 ∈ V | rows reduced | muladds | stored words | wall s | parallel / sequential s |
|--:|--:|--:|:--|:--|:--|--:|--:|--:|--:|:--|
| 11 | 4 | 9,109 | 4,588 / 4,521 | 0, 0, 2, 131, 1805 | no | 1,442 | 3,718,432 | 112,273 | 0.014 | |
| 11 | 5 | 35,443 | 17,755 / 17,688 | dim 25,033 | no | 37,852 | 24,926,440,157 | 14,702,864 | 9.8 | 9.5 / 0.3 |
| 11 | 6 | 110,056 | 55,144 / 54,912 | whole space | **yes** | 262,484 | 539,254,650,848 | 111,433,956 | 231.1 | 227.4 / 3.7 |
| 12 | 4 | 12,951 | 6,515 / 6,436 | 0, 0, 2, 135, 2069 | no | 1,621 | 3,728,860 | 123,586 | 0.019 | |
| 12 | 5 | 55,455 | 27,767 / 27,688 | 0, 0, 2, 146, 3727, 32057 | no | 44,725 | 25,954,076,158 | 18,941,587 | 9.8 | 9.4 / 0.3 |
| 12 | 6 | 190,051 | 95,175 / 94,876 | whole space | **yes** | 436,789 | 3,690,966,694,706 | 352,689,051 | 1,722.4 | 1,712.6 / 9.7 |
| 13 | 4 | 17,902 | 8,997 / 8,905 | 0, 0, 2, 139, 2341 | no | 1,808 | 3,770,769 | 139,597 | 0.022 | |
| 13 | 5 | 83,682 | 41,887 / 41,795 | 0, 0, 2, 150, 4021, 39657 | no | 52,274 | 26,162,190,291 | 22,305,343 | 12.9 | 12.4 / 0.5 |
| 13 | 6 | 313,912 | 157,145 / 156,767 | whole space | **yes** | 683,496 | 24,217,016,013,719 | 1,070,913,677 | 11,774.7 | 11,743.8 / 30.8 |

- **t = 11** (`raw/hky2_first_t11.jsonl`): whole cell 241 s. Falls 4, 5, 6; last fall degree 6. The
  dimension 25,033 of `V_5` equals the other repository's `k11d4` value (**given**, `LIFT.md` section 5).
- **t = 12** (`raw/t12/`): whole cell **1,732 s** on 4 shared vCPU; peak resident memory 3.2 GB
  (`raw/t12/rss.log`). 1 enters at c = 6 with 190,051 columns, falls 4, 5, 6 — the reference's
  answer on the pods (6,857 s and 7,269 s on 16 vCPU, **given**). When the engine stopped it had
  95,161 + 94,863 of 190,051 pivots. This is the one cell above t = 10 where the new engine's
  answer can be compared with the reference's, and only through the given record: the two were
  not run on the same machine.
- **t = 13** (`raw/t13/`, draw 0, `x_r` = 1118191946): whole cell **11,788 s (3 h 16 min)** on 4
  shared vCPU (1-minute load 3.2 to 4.8 through the run), peak resident memory **9.6 GB**
  (`raw/t13/rss.log`). **1 enters at c = 6**; `V_5` has dimension 39,657 of 83,682 and does not
  contain 1; falls 4, 5, 6, first fall degree 4, last fall degree 6. This is the value expected
  from the other repository's engine (**given**: last fall degree 6 at t = 13). There is no
  reference run to compare with: on the pods the reference hit its 14,400 s limit on one draw and
  aborted on the other (**given**). `1` was again the last pivot: the echelon held all 313,912
  pivots when the engine stopped.
- **What this run is not.** One draw, one schedule, one machine, no second instrument at this
  size. "1 enters" rests on the arithmetic of this closure, which the dims comparison cannot check
  (section 5); the supporting evidence is the agreement at t ≤ 10 with the reference, at t = 11, 12
  with the given values, and the planted and fault controls. A second run with another schedule
  (`--hky2-order lifo --hky2-sym off`) or a capped-F4 bracket would be the independent check; it
  was not made (another 3 hours or more here).

## 9. The degree lift: does the theorem carry over?

`LIFT.md` of the other repository computes `V_{D+1} ∩ M_{≤D}` and whether `1 ∈ V_{D+1}` from a
closed `V_D`. Assessment for this engine's ring (**proof sketch**, nothing implemented or run):

**The rings.** The Kummer tower rule of the study cell is `y_j² = y_{j+1}` (a = c = d = 0, b = 1),
`y_{t−1}² = d_top`. A square-free tower monomial is `y_0^a`, `0 ≤ a < 2^t`, the degree is the
number of binary digits of `a`, and `y_j · y_0^a = y_0^{a + 2^j}` with `y_0^{2^t} = d_top`: the
quotient is `F_p[y]/(y^{2^t} − d_top)` per summand. Rescaling `y = g·x` with `g^{2^t} = d_top`
(the tower's own `g`) gives `F_p[(Z/2^t)²]`, the other engine's ring, and the rescaling is
diagonal on monomials, commutes with the degree and with multiplication by variables up to
scalars. So the closures correspond and have the same dimensions. Measured agreement: every
dimension in section 2 equals the other repository's value at the same depth and degree.

**Hypotheses of the theorem, and whether they hold here:**

| hypothesis (where it is used) | Kummer tower, m = 2 (study cell) | general tower rule (isogeny, Dickson: a, c ≠ 0) | m ≥ 3 (free chain unknowns) |
|:--|:--|:--|:--|
| H1 multiplication by a variable maps `span M_{≤w}` into `span M_{≤w+1}` (Lemma 1, 2) | holds (rewriting never raises the degree) | holds | holds |
| H2 `V_D` is the **closed** degree-D closure, given as a reduced echelon basis, pivot = leading monomial, rows `e_pivot +` free columns (Lemma 2, 4) | holds for the definition; **this engine does not produce that basis** (its rows are reduced only against earlier pivots) — a back-substitution is needed first | same | same |
| H3 σ is a symmetry, p odd (Lemma 3) | holds, proved per system (section 7) | holds when the two towers are equal | σ = swap of the first two towers holds |
| H4 a variable times a monomial of degree D has **at most one** monomial of degree D + 1, and it is that monomial times a variable (Lemma 4: "at most one regular top coordinate"; the definition of special columns) | holds: the product is one monomial | **holds by an argument not in LIFT.md**: if `y_j ∤ x` the product is the monomial `x·y_j`; if `y_j | x` the rule gives `a·x·y_{j+1} + (degree ≤ D)`, and if `y_{j+1} | x` the same step repeats up the tower, so there is at most one top monomial `x·y_{j′}`, with coefficient a product of a's | free variables multiply without rewriting: one top monomial |
| H5 the generators of a regular column are exactly its D + 1 divisor variables (the coverage assertion; not needed by Lemma 5) | holds | **fails**: the carry terms `a·x·y_{j+1}` add generators to a column. Lemma 5 needs only "at least one", which still holds; the assertion and the top coefficients (±1, 2 there) must be generalised | fails as stated for non-square-free free-variable monomials (a divisor variable can repeat); "at least one" holds |

Conclusion: for the Kummer m = 2 cell — the cell at t = 13 and t = 14 — the theorem applies as
stated (the rings are the same up to a diagonal change of basis and the summation generator is
the same polynomial). For general towers and m ≥ 3 the proof goes through with H5 weakened, as a
sketch only; it would need its own validation against direct closures before any use.

**Not implemented.** What an implementation here needs: (1) the reduced echelon basis of both
parts of `V_D` on its free columns (dense `u32`: 224,688 × 25,136 and 224,234 × 25,120 entries at
t = 14, D = 6, 45 GB; the back-substitution from this engine's echelon costs up to r²·f/2 ≈
6e14 multiply-adds per part in a dense modular product, which the engine does not have);
(2) the lift itself, about 900 lines in the compact Python reference
(`csd/ic_bounds_lift_independent.py`), with followers (up to 5.25 GB each, 56 of them, cached);
(3) validation against direct closures on every case that fits (t ≤ 10 here). The other
repository's measured cost at depth 14: 751 s (1 enters) and 4,679 s (1 does not) on one H200
GPU (**given**). Section 10 gives the alternative that needs no new mathematics: the direct
closure at c = 7 with this engine.

## 10. Memory, and the estimate for t = 13 and t = 14

### 10.1 Memory by representation

Columns `n = Σ_{k ≤ c} C(2t, k)`; split parts `n₊ = (n + s)/2`, `n₋ = (n − s)/2`,
`s = Σ_{k ≤ c/2} C(t, k)`.

| problem | n | n₊ / n₋ |
|:--|--:|:--|
| t = 13, c = 6 | 313,912 | 157,145 / 156,767 |
| t = 14, c = 6 (rank 448,922, **given**) | 499,178 | 249,824 / 249,354 |
| t = 14, c = 7 | 1,683,218 | 841,853 / 841,365 |

| representation | t = 13, c = 6 | t = 14, c = 6 | t = 14, c = 7 | status |
|:--|--:|--:|--:|:--|
| dense `u32`, n × n, unsplit | 394 GB | 997 GB | 11.3 TB | arithmetic |
| dense `u32`, rows from their leading column (triangle), split | 99 GB | 249 GB | 2.8 TB | arithmetic |
| reduced echelon on free columns, dense `u32`, split: peak `n₊²/4` entries per part during a closure; final r × f | 49 GB peak | 125 GB peak; **45 GB** final | 1.4 TB peak | arithmetic; final value **given** (LIFT.md: 21 + 21 GB) |
| reference: sparse rows, 8 bytes per entry, plus every pending product materialised in the queue | echelon about 3e9 entries = 24 GB at the density of section 2 (0.057 n²/2); queue unbounded — one draw aborted (exit 134, **given**) | | | estimate |
| **this engine** (hybrid rows, 8 bytes per sparse entry, 4 per dense slot; the queue holds 8-byte (row, variable) items; 8 n₊ bytes of scratch per thread) | 1.07e9 words, **9.6 GB resident** (peak) | 1.7e9 to 5.4e9 words: **14 to 43 GB** | 5.7e9 to 2.5e10 words: **45 to 200 GB** | t = 13: **measured**; t = 14: estimate, see 10.3 |

What fits (this engine): t = 13 in 16 GB (measured: 9.6 GB); t = 14, c = 6 in 128 GB if the estimate
holds, 256 GB with margin; t = 14, c = 7 directly needs 256 GB at the low end of the estimate and
does not fit 256 GB at the high end: 512 GB.

### 10.2 The model for the closure in which 1 enters (c = 6, t = 10..13)

`M(t) = 2 κ n₊³` multiply-adds, wall = `M × ρ`, stored words = `w n²`.

| fit point | n₊ | M (measured) | κ | wall s of the closure | ρ ns per muladd (4 threads, this machine) | w |
|:--|--:|--:|--:|--:|--:|--:|
| t = 10 | 30,318 | 92,397,155,280 | 0.001658 | 39.4 | 0.426 | 0.00942 |
| t = 11 | 55,144 | 539,254,650,848 | 0.001608 | 231.1 | 0.429 | 0.00920 |

The dominant term is the reduction of products against sparse pivot rows; n₊³ with a constant κ
fits the two points within 3 %. The prediction for t = 12 and t = 13 was written and committed
(`PREDICTION_t12_t13.md`, commit `64fb066690`) **before** either was run:

| t = 12, closure c = 6 | predicted | band | measured | verdict |
|:--|--:|:--|--:|:--|
| answer | 1 enters at c = 6 | | 1 enters at c = 6 | as predicted |
| multiply-adds | 2.81e12 | 2.0e12 – 3.7e12 | 3.69e12 | inside, at the upper edge |
| wall | 1,210 s | 850 – 1,700 s | 1,722 s | **outside by 1.3 %: refuted on time** |
| stored words | 3.4e8 | | 3.53e8 | 4 % above |

So κ is not constant: 0.00166, 0.00161, **0.00214** at t = 10, 11, 12, and ρ rose to 0.467 ns.
The model was revised before the t = 13 closure had run (`PREDICTION_t12_t13.md`, second part,
commit `08a4fb3b0c`): `M ∝ n₊^3.22` (least squares on three points; the last step alone has
exponent 3.52), 1.75e13 to 2.16e13 multiply-adds and 8,200 s to 10,800 s at t = 13 on this
machine.

| t = 13, closure c = 6 | revised prediction | measured | verdict |
|:--|--:|--:|:--|
| answer | 1 enters at c = 6 | 1 enters at c = 6 | as predicted |
| multiply-adds | 1.75e13 – 2.16e13 | 2.42e13 | 12 % above the upper value |
| wall | 8,200 – 10,800 s (refuted above 13,000 s) | 11,775 s | 9 % above the range, below the bound |
| stored words | 9.3e8 | 1.07e9 | 15 % above |
| resident memory | about 9 GB | 9.6 GB | |

**Both held-out predictions erred low.** κ is 0.00166, 0.00161, 0.00214, 0.00312 at t = 10..13;
the local exponent of M in n₊ is 2.95, 3.52, 3.75. ρ rose from 0.43 to 0.49 ns. A fixed-exponent
model underestimates; the estimates of 10.3 inherit that bias and should be read with a factor of
2 added on the high side.

### 10.3 t = 14 — estimate, two analogues only

No run at t = 14 was made. The closures needed are c = 6 (1 does not enter: a fixpoint closure,
every product of every row of degree < 6 is reduced) and c = 7 (1 expected to enter).

- **c = 6, fixpoint.** The fixpoint closures one degree below the last fall degree measured here:
  (t = 6, c = 4): κ = 0.0229, w = 0.064, rank/n = 0.76; (t = 10, c = 5): κ = 0.00925, w = 0.0215,
  rank/n = 0.856. (t = 14, c = 6) is the next "first depth of its degree", rank/n = 0.899. κ and w
  fall by 2.5× and 3× per degree; taking κ between 0.004 and 0.00925 gives
  `M = 2 κ (249,824)³` = 1.2e14 to 2.9e14 multiply-adds.
  A second reading, by items: a fixpoint closure reduces one item per (row of degree < c, variable
  of the first tower). At (t = 10, c = 5) that is 3,163 rows × 10 = 31,631 items at 7.5e5
  multiply-adds each, 7.4 % of the stored words per item. At (t = 14, c = 6) the rows of degree
  ≤ 5 are 73,587 (**given**), so 1.03e6 items, and 7.4 % of 1.7e9 to 5.4e9 stored words (10.1)
  each: 1.3e14 to 4.1e14. Taken together: **1.2e14 to 4e14 multiply-adds**, 5 to 17 times the
  measured t = 13 closure. (Inside one degree the fixpoint closure barely grows: c = 5 costs
  2.38e10, 2.49e10, 2.60e10, 2.62e10 at t = 10..13.)
- **c = 7, direct, 1 entering.** κ of the closure in which 1 enters at the first depth of its
  degree: c = 5 (t = 6): 0.0103; c = 6 (t = 10): 0.00166. With κ between 0.0003 and 0.0016:
  `M = 2 κ (841,853)³` = **3.6e14 to 1.9e15**.
- The other repository's engine did the c = 6 closure in 4.69e15 multiply-adds (dense, GPU,
  **given**); the estimate here is 12 to 40 times fewer. A ratio of that order between the two
  kinds of elimination is plausible but was not measured on a common case.

### 10.4 What machine, how long — estimate

Rate: ρ = 0.43 to 0.49 ns per multiply-add on this machine with 4 threads and a neighbour
(0.486 at t = 13), i.e. at most 1.94 core-ns, the value used below (the neighbour took part of
the 4 vCPUs, so the true per-core cost is lower). Per-core speed of another machine and scaling beyond 4 threads are **not
measured**. Two readings are given: linear scaling to the core count, and a bandwidth-limited one
(a sparse entry is 8 bytes read per multiply-add; at 1.9 core-ns that is 4 GB/s per core, so about
30 to 40 cores saturate a 150 GB/s memory system; 30 is used).

| problem | multiply-adds | 16 cores | 64 cores, linear | 64 cores, bandwidth-limited (30 effective) | RAM |
|:--|--:|--:|--:|--:|:--|
| t = 13, c ≤ 6 (1 enters at 6) | 2.42e13 (**measured**; 11,788 s here on 4 shared vCPU) | 49 min | 12 min | 26 min | 9.6 GB measured; 16 GB |
| t = 14, c = 6 | 1.2e14 – 4e14 | 4 – 13.5 h | 1 – 3.4 h | 2.2 – 7.2 h | 128 – 256 GB |
| t = 14, c = 7 direct | 3.6e14 – 1.9e15 | 12 – 64 h | 3 – 16 h | 6.5 – 34 h | 256 – 512 GB |
| t = 14, c = 7 by lift | not implemented here | | | | 256 – 512 GB (followers), **given** 751 – 4,679 s on one H200 GPU |

**Honest statement.** t = 13: **done on this 4-vCPU, 16 GB machine in 3 h 16 min** (one draw); on a 16-core machine about 50 minutes per draw if the rate scales with cores, which was not measured. t = 14: expected to finish on one 64-core machine with
512 GB in one to three days (c = 6 and c = 7 direct), **if** the two-point trends of 10.3 hold; they
are extrapolations across degrees with one analogue each, both held-out predictions inside one
degree erred low (10.2), and they can be wrong by a factor of 3 or more. The c = 6 closure alone
(which decides `1 ∉ V_6`, i.e. `d_F ≥ 7`) is the better-founded half: a 64-core, 256 GB machine,
an afternoon to a day. Nothing at t = 14 was run.

## 11. What a study amendment would have to say (not written, preregistration not edited)

The engine is an instrument. An amendment to `research/prime_fall_degree_20261027/PREREGISTRATION.md`
(and to the `pkm_falldegree_unblock_20261107` numbering) would have to state:

1. **Instrument.** Closures may be computed by `hky_fp_tower_fast` (`--engine hky2`, level 4,
   sparse-first, split when proved), commit and binary hashes named. Definition, prime
   (P1 = 2013265921), seeds, draws, `--hky-cmax`, `--hky-max-cols` unchanged.
2. **Equivalence evidence.** The gate of section 5 on the frozen binaries, with the negative
   controls, and the statement of its blind spot (closures that contain 1).
3. **Limits that change meaning.** `--hky-max-entries` and `--budget` are limits of the
   instrument: the fast engine counts stored terms of its own rows and is censored later. An
   amendment must fix them for the new instrument (e.g. `max_entries` 1.4e10 as before, which the
   new engine does not approach at t = 13) and say that "censored" rows of the two instruments are
   not comparable.
4. **Row schema.** Rows carry `engine: hky_fp_tower_fast` and an `hky2` block (options, per-closure
   timers, `symmetry`). Tools that select `engine == "hky_fp_tower"` (`prime_fall_degree_tools`)
   would skip them until told otherwise; G3 (capped F4 agrees with the closure) should be re-run
   with the new engine on its N ≤ 20 cells.
5. **Soundness of a positive result.** A row "1 enters at c" from the fast engine rests on the
   arithmetic of its last closure, which the dims comparison cannot check (section 5). The
   amendment should require, for a registered t = 13 or t = 14 row, either a second run with a
   different schedule (`--hky2-order lifo`, `--hky2-sym off`: different rows, same answer) or the
   bracket by capped F4 already used at t = 12.
6. **Unregistered runs.** The t = 12 and t = 13 runs of section 8 were made before any amendment;
   they are instrument checks and cannot be scored as study rows.

## 12. Not tried, and why

- **Lift** (section 9): needs a dense reduced-echelon basis and a modular matrix product this
  engine does not have; a week of work with its validation, not an afternoon.
- **Blocked elimination** (several candidates per pass over a pivot row, to cut memory traffic):
  the 4-vCPU machine cannot show a bandwidth limit, so its gain could not have been measured here.
  It is the first thing to do if a 64-core run scales badly.
- **Back-substituted (reduced) echelon rows**: the lead-only result says fuller reduction makes
  sparser rows; full back-substitution may reduce the work further, at the price of updating
  stored rows. Not tried.
- **Unchecked indexing, prefetch, Montgomery form** in the sparse loop: single-thread cost is
  about 1.45 ns per multiply-add at t = 9 including products; without cycle counters a 10 %
  change cannot be told from the neighbour's load.
- **A commutation (syzygy) criterion.** At t = 13 the engine reduced 1,366,992 products to find
  313,912 pivots; 77 % reduce to zero. `x_i·(x_j·g)` and `x_j·(x_i·g)` are the same element, and
  a signature-style rule could skip most of the repeats, as `sig_fp_tower` does for F4. It needs
  an ordering argument for soundness (a fall invalidates the naive rule) and was not attempted.
- **Per-batch scratch.** Each rayon job of each batch allocates and zeroes an accumulator of
  8 n₊ bytes; a few per cent at most at t = 12, not removed.
- **A second schedule at t = 13** as an independent check of "1 enters" (section 8).
- **Reference at t = 11, 12**: hours on this machine; the pods' t = 12 result is taken as given.
- **Isolation** (`tools/isolated_bench.py`), instruction counts, a `bench.record/v1` record, S3
  upload, leaderboard, `benchmarks/` copy, PDF and graphs: the tools need services not present
  here (S3, Conductor) or would publish a row the validation gate forbids without a replay
  certificate. To be done with the PR, not claimed here.
- **m ≥ 3 at depth**: the gate covers m = 3 (t ≤ 3) and m = 4 (t = 2) for equality only; no timing.

## 13. Files

| path | content |
|:--|:--|
| `src/cryptanalysis/hky_fp_tower_fast.rs` | the engine, `check_symmetry`, 10 unit tests |
| `examples/pkm_tower_pilot.rs` | `--engine hky2`, `--hky2-level`, `--hky2-sym`, `--hky2-order`, `--hky2-batch`, `--hky2-dense`, `--hky2-lead-only`, `--hky2-sparse-width`, `--hky2-fault`; reference rows unchanged |
| `examples/hky_engine_equality.rs` | `compare`, `compare --expect-mismatch`, `times` |
| `harness/`, `harness/unit/` | standalone build and unit-test crates |
| `run_gate.sh`, `run_ab.sh` | the gate and the interleaved A/B |
| `PREDICTION_t12_t13.md` | the held-out prediction, committed before the runs |
| `raw/` | every row, log and `conditions.log` behind the tables |
