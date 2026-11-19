# Optimisation round on the hky3 closure engine for the no-1 production closures

Date: 2026-11-18. Worktree `wt/engine-opt` (sparse checkout of isogeny-labs/crypto), local commits only.
Engine: `src/cryptanalysis/hky_fp_tower_fast2.rs` (`--engine hky3`), driven by `examples/pkm_tower_pilot.rs`,
built through `harness/` (a copy of `research/rand_vs_struct_fb_20261118/harness/`, which compiles only the modules
the pilot needs). Every number is `[pcc:unverified]` (no pcc claim filed; no study, leaderboard or
`benchmarks/README.md` row is claimed). Labels: **measured** (this machine, a command of this directory),
**given** (pod records of `research/hky_pod_runs_20261117/`, or the earlier READMEs), **estimate**.

The machine: Intel Xeon @ 2.1 GHz, 4 vCPU, L2 2 MiB/core, L3 260 MiB, 16 GiB, shared with five other workers'
builds and runs: the 1-minute load was **10 to 24 on 4 vCPU** throughout, so every process got about one core.
Wall times here are therefore noisy and roughly single-core; the deterministic counters are the primary evidence,
and CPU seconds of the process are reported beside wall.

## 0. Summary

- **What dominates a no-1 full closure** (§1, measured on the proxies): the parallel reduction of products
  (89–91 % of wall); within it, **reductions to zero** (51–57 % of all multiply-adds; their number,
  products − rank, is fixed by the final space, so no schedule can lower it) and, by rows, the **subtraction of rows
  whose pivot was already a pivot at the last interreduction** (57–65 % of all multiply-adds). The last 20 % of the
  rank costs 87–89 % of the multiply-adds. Interreduction 6–9 %, second pass 2–3 %, the calling thread alone
  0.1–0.4 %. Interreduction frequency moves the total by ±10 % only (rejected as a lever).
- **Two changes, adopted, neither skips any reduction** (§2): (1) *epoch coordinates*: stored row tails are written
  in the columns that were not pivots at the last interreduction; a product's terms at those pivots are final, so
  their rows are subtracted at once into an accumulator of `n − rank` entries instead of `n` (no scan, no cascade,
  many rows dense → SIMD); (2) *the sweep*: per batch, that subtraction is done by tiles of columns, each tile a task
  that reads each needed row **once** and applies it to every product of the batch that needs it (measured sharing:
  a row is needed by up to 19 products of a batch in the last rank decile).
- **Identical mathematical output**: same `dims`, `one`, falls, censoring — and also the same `muladds`,
  `rows_reduced`, `entries`, zero count and reducer count (the reduced form of a vector against a fixed echelon is
  unique, and so is the set of rows subtracted with non-zero factor) — in every run of this round (unit test
  `epochs_and_sweep_do_the_same_reductions`, every A/B run). **Gate PASS** (§4): 110 reference cells EQUAL, faults 1–5
  detected (4 and 5 new, attacking the epoch rewrite and the sweep), and the new **closedness check** run on every
  gate closure without 1: **314 closures, 708,997 reductions, 0 non-zero**.
- **Measured gain** (§3, 4 shared vCPU, interleaved, same counters): production-shaped closures (high rank fraction)
  **1.70× wall / 2.17× CPU** on the full no-1 closure planted t = 10 c = 6 (60,460 columns, 228 → 134 s); **1.35× wall / 1.48× CPU** at t = 11 c = 6 (110,056 columns, 306 → 226 s); **1.47×** at t = 10 c = 5 (85 % rank);
  1.35× and 1.16× on the low-rank proxies t = 12, 14 c = 5 (58 %, 39 % rank). Memory: stored row bytes unchanged
  (t = 14 c = 5: 283 MB both); RSS +6 % at t = 11 c = 6 from the sweep's accumulators (capped at 2 GiB).
- **Not measured here, argued** (§3.3): at k = 14..16 the base engine's per-thread accumulator is 2–4.6 MB
  (out of L2) and each product streams ~10³ reducer rows from DRAM; the new reduction works in a ~0.2–0.6 MB
  accumulator and reads each row once per batch. The pod's per-thread rate at k = 14 was 3.1× worse than this
  machine's at t = 11, which is the effect these changes target; `pod/scale_bench.sh` measures it.
- **Predictions** (`PREDICTIONS.md`, written before any pod data): V_6 on 32 vCPU, old → new:
  k = 14 4.6 h → 3.0 h; k = 15 30 h → 18 h; k = 16 177 h → 93 h (wide bands there), RSS 23 / 53 / 117 GB → +0..10 %.
- Frozen binary (no target-cpu flag, glibc ≥ 2.34): `research/hky_engine_opt_20261118/bin/pilot_opt2`,
  sha256 `c627ef4781bdb30c18e5d00a165a13d2f3dd79cff6691e05962970ee0557cf9f` (bin/ is git-ignored; rebuild:
  §6). Base for every A/B: `bin/pilot_base` (commit `c2afb6b334`), sha256
  `a8ff8af55b6d5faf6c66fffb2f00b4c20fd714e776545521f3570fcd66b039e1`.

## 1. Profile of the no-1 closures (base engine) — measured

Proxies: m = 2, Kummer tower, p = 2013265921, the pilot's draw 0 (`--planted 0 --random 1`), only the closure
c = 5 (`--hky-cmin 5 --hky-cmax 5`; LFD = 6 at t = 10..14, so V_5 has no 1); plus t = 11 c = 6 (1 enters only at
the last pivot: 99 % of the rank is built, the shape of a full closure). Runner `cell1.sh`, A/B runner `ab.sh`,
summaries `ab_summary.sh`. Traces `raw/profile/*.trace.jsonl.gz` (one line per batch and interreduction).

| cell | columns | rank (%) | muladds | products | to zero (all) | parallel | interreduce | second pass | calling thread alone |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| t = 10 c = 5 | 21,700 | 18,569 (85 %) | 1.603e10 | 63,261 | 44,692 (71 %) | 89 % | 9 % | 2 % | 0.1 % |
| t = 11 c = 5 | 35,443 | 25,033 (71 %) | 1.634e10 | 75,703 | 50,670 (67 %) | 90 % | 7 % | 3 % | 0.1 % |
| t = 12 c = 5 | 55,455 | 32,057 (58 %) | 1.744e10 | 89,449 | 57,392 (64 %) | 91 % | 6 % | 3 % | 0.1 % |
| t = 14 c = 5 | 122,438 | 47,849 (39 %) | 1.830e10 | 121,045 | 73,196 (60 %) | — | — | — | — |
| t = 11 c = 6 | 110,056 | 109,135 (99 %) | 3.215e11 | 515,081 | 405,946 (79 %) | 95 % | 3 % | 2 % | 0.05 % |

(t = 14 c = 5 dims `[0,0,2,154,4323,47849]` equal the pod's hky2 record.) Share of the reduction multiply-adds
(new counters, identical in both engines), as shares of all multiply-adds: against rows of epoch pivots (phase 0) **65 %** (t = 10 c = 5), 60 % (t = 12),
61 % (t = 14), **57 %** (t = 11 c = 6); spent on products that reduce to zero **53 %**, 57 %, 51 %, **52 %**.
The cost is concentrated at the end: the last 20 % of the rank does 87–89 % of the multiply-adds (traces), where a
product is reduced against ~1,000–1,300 rows of ~800 terms (`trace_summary.sh` of the speed2 round). Stored words
grow with the rank as before (`words/rank` 20 → 830 over the closure at t = 10).

Where the reduction's time goes, structurally: the accumulator has `n` (part) entries; an interreduced row has
terms only at non-pivot columns, and right after an interreduction the stored rows are **53 % dense on the non-pivot
columns after their lead** (trace field `shape`, t = 10 c = 5 at 96 % of the rank: 6.65e6 terms in 1.26e7 possible).
So most of the work is sparse-format scatter of near-dense data, into an accumulator that at k = 14..16
(250k–575k columns per part, 2–4.6 MB) does not fit in L2, re-reading the same rows for every product of a batch.

Scaling with columns: on the c = 5 proxies the work hardly grows (1.6e10 → 1.8e10 for 5.6× the columns) because the
rank fraction falls (85 % → 39 %); the production closures have rank fraction 90 % (k = 14), so the high-rank
proxies (t = 10 c = 5, t = 11 c = 6) are the relevant ones, and the t = 13 → 14 pod data give muladds ∝ n^4.15.

Interreduction frequency (`raw/sweep0/`, base engine): threshold 0.1 / 0.05 / 0.25 (default) / 0.5 gives
1.551e10 / 1.601e10 / 1.603e10 / 1.770e10 muladds at t = 10 and 1.722 / 1.725 / 1.744 / 1.960e10 at t = 12: ±10 %,
not a lever; default kept.

## 2. The changes

### 2.1 Epoch coordinates (`--hky3-np on`, default)

`Echelon` gains `np_cols` (the non-pivot columns at the last interreduction, ascending) and `of_col` (its inverse,
derived). Stored rows keep their leading column in part coordinates and write their tail in `np_cols` coordinates
(a dense row stores `v[1+i]` for epoch column `np_start(lead)+i`). After every interreduction every row's tail has no
pivot column, so `new_epoch` drops the new pivots from `np_cols` and rewrites every row (in parallel; it panics if a
tail met a pivot — an invariant check). Rows inserted later are fully reduced, hence also free of epoch pivots.

`eliminate` then has a phase 0: a term of the product at an epoch-pivot column is final (no stored tail contains
such a column), so that row is subtracted at once with the product's own coefficient; the other terms go to an
accumulator of `np_cols.len()` entries; then v2's cascade (bitmap while sparse, then the scan) runs in epoch
coordinates against the rows inserted since the epoch. Products are still formed from `row_full`, in part columns.

### 2.2 The sweep (`--hky3-sweep on`, default; `--hky3-sweep-ratio 4`, `--hky3-sweep-tile-bytes 1048576`, `--hky3-sweep-max-words 2^28`)

`reduce_products` plans phase 0 for every product of a batch (which epoch-pivot rows, with which factor), and when
the phase-0 work is at least 4 multiply-adds per accumulator entry and the accumulators fit in 2^28 words, it
gives each product its own epoch-coordinate accumulator, groups the (row, product, factor) triples by row, and runs
phase 0 **by tiles of epoch columns**: a tile is one task that owns those columns of every product's accumulator,
reads each needed row's segment once and applies it to all its users (row segment in L1, the tile of all
accumulators ≤ 1 MiB). The cascade then runs per product. Otherwise every product is reduced alone, as in 2.1.
Measured sharing (trace fields `core_sum`/`core_union`): in the last rank decile of t = 10 c = 5, a batch needs
5.5e5 distinct rows for 1.06e7 uses (19×); in blocks of 8 consecutive items only 1.9×, which is why the sweep is
over the whole batch.

### 2.3 Smaller items

- Storage forms of a batch's inserted rows computed in parallel (a dense form is as long as the epoch).
- The sweep's sort of triples is `par_sort_unstable` and its per-row metadata is built in parallel, so no
  per-batch section grows with the batch on one thread (§3.3).
- New counters in `Fast2Stats` and the pilot's JSON: `core_muladds`, `zero_muladds`, `swept_batches`,
  `final_row_bytes`, `sweep_bytes`, `verify_*`; trace fields `ops0`, `zero_ops`, `core_sum`, `core_union`,
  `core_union8`, interreduction `shape`.
- `--hky3-verify-closed on`: the closedness check (§4.2).
- `--hky3-fault 4` (a new epoch drops the last term of every row with ≥ 2 tail terms), `--hky3-fault 5` (the sweep
  skips its first tile): negative controls of the new code.

## 3. A/B — measured

Interleaved (base / np-only / opt in turn, 3 reps), 4 threads, other workers' load 7–14 on 4 vCPU (so each process
got ~1–1.5 cores; wall ≈ CPU). Counters identical across every run of a cell (`ab_summary.sh`: "counters identical").

| cell | arm | n | wall median (min–max) s | CPU s | peak RSS MB |
|:--|:--|--:|--:|--:|--:|
| t = 10 c = 5 | base | 3 | 27.7 (27.7–34.4) | 27.9 | 146 |
| | np only | 3 | 20.7 (18.4–23.2) | 19.1 | 138 |
| | **opt** (np + sweep) | 3 | **18.8** (18.1–21.8) | 18.2 | 185 |
| t = 12 c = 5 | base | 3 | 35.8 (33.7–40.3) | 33.3 | 245 |
| | np only | 3 | 31.0 (29.5–34.4) | 30.1 | 230 |
| | **opt** | 3 | **26.6** (21.6–28.0) | 26.4 | 302 |
| t = 14 c = 5 | base | 3 | 37.9 (24.5–40.6) | 36.0 | 345 |
| | np only | 3 | 35.9 (21.8–43.9) | 35.8 | 314 |
| | **opt** | 3 | **32.8** (20.0–36.3) | 31.0 | 416 |
| t = 11 c = 6 | base | 1 | 306.1 | 578 | 1,267 |
| | **opt** | 1 | **226.4** | 391 | 1,342 |

(`raw/ab_main_thr4/`, `raw/ab_t11c6/`; t = 11 c = 6 load 10 → 7.6 during the pair.) Single-thread CPU (`RAYON_NUM_THREADS=1`,
`raw/ab_np1_thr1/`, load ~20): t = 10 c = 5 base 30.0 s, np 21.1 s, np-off (same binary, A/A of the code paths) 27.9 s;
t = 12 c = 5 34.1 / 32.7 / 34.1 s.

**Full no-1 closure** (`raw/ab_planted/`, the consistent "planted" system, 1 in no closure, t = 10 c = 6: 60,460
columns, rank 60,458, 3.640e11 muladds, 433,961 products of which 350,398 reduce to zero in the first pass;
phase 0 177.6e9 = 49 %, zero products 320.0e9 = 88 % of muladds): base **228 s wall / 499 s CPU** (load 1.3), opt
**134 s / 230 s** (load 8.6): **1.70× wall, 2.17× CPU**, counters and dims identical; RSS 465 → 590 MB (row bytes
410 MB, sweep 41 MB).

### 3.2 Memory

Stored row bytes (new counter `final_row_bytes`): t = 14 c = 5 283 MB (opt) vs ≈ 282 MB (base: 35.29e6 words almost
all sparse at 8 B); t = 11 c = 6 1.19 GB vs ≈ 1.2 GB. The sweep's accumulators add 21 MB (t = 10), 100 MB (t = 14),
81 MB (t = 11 c = 6) and are capped by `--hky3-sweep-max-words` (2 GiB). Dense threshold 1/2 instead of 1/3 in
epoch coordinates (`raw/ab_dense/`): row bytes 81 vs 86 MB (t = 10), 282 vs 283 MB (t = 14), wall 12.6 vs 12.3 s and
22.5 vs 21.1 s — no memory to gain, slightly slower: **rejected**, default 1/3 kept.

### 3.3 32 threads (reasoning; not measurable here)

- The base engine's reduction scatters into an accumulator of `n` part columns (k = 14: 250k → 2 MB, k = 16: 575k →
  4.6 MB per thread), larger than a core's L2, and streams every reducer row from DRAM for every product: a late
  product reads ~10³ rows × ~10³–10⁴ terms × 8 B. On the pod the hky2 run at k = 14 did 5.7e9 muladd/s on 32 vCPU,
  0.18 G/s per vCPU, against 0.56 G/s per CPU-second on this machine at t = 11 c = 6 (accumulator 440 KB in L2): **3.1× slower per
  thread** — cache and bandwidth, not arithmetic.
- New: phase 0's accumulator is `n − rank` epoch columns late in the closure (k = 14: ~25k per part, 200 KB), the
  tile of all accumulators of a batch is ≤ 1 MiB, and each row is read once per batch (19× fewer row reads in the
  last decile at t = 10). This is the change most likely to pay more on the pod than here.
- Parallel structure per batch: form products (one task per item), plan (per product), sort (parallel), tiles
  (≥ 4 × threads tiles per part by construction, independent column ranges, no locks), cascade (one task per
  product), second pass (parallel rounds, as v2), storage forms (parallel), queue pushes (serial, O(inserted ×
  variables)). Interreduction: blocks of 256 rows (a barrier each) and the epoch rewrite (parallel). The calling
  thread alone measured 0.81 s of 226 s at t = 11 c = 6 (0.4 %).
- Determinism: tiles and tasks write disjoint ranges; results are reassembled in batch order; the schedule depends
  on counts only. Thread count changes neither the counters nor the result (unit test of v2 still passes; the
  new unit test compares with/without epoch and sweep).

## 4. Equality gate and soundness — measured

### 4.1 Gate

`run_gate.sh raw/gate-opt --hky3-verify-closed on` (copy of the speed2 gate: definitional reference `--engine hky` of
`bin/pilot_base` against `--engine hky3` of the candidate, compared by `bin/hky_engine_equality`):

| set | cells | result |
|:--|--:|:--|
| m2-tower-t1-8 | 16 | EQUAL |
| m2-tower-planted (1 in no closure) | 12 | EQUAL |
| m2-null-t2-7 | 18 | EQUAL |
| m2-naive-t2-4 | 6 | EQUAL |
| m3-tower-t2-3 | 12 | EQUAL |
| m4-tower-t2-3 | 2 | EQUAL |
| m5, m6, m7 (t = 2) | 2 + 2 + 2 | EQUAL |
| p0-kummer-isogeny (p = 786433) | 36 | EQUAL |
| censor-max-cols | 2 | EQUAL |
| neg-fault1 / 2 / 3 (v2's) | | detected: 60 / 40 / 40 differences |
| **neg-fault4** (epoch rewrite drops a term) | | **detected: 40 differences** |
| **neg-fault5** (sweep skips a tile) | | **detected: 38 differences** |

**GATE PASS** (`raw/gate-opt.log`, records gzipped in `raw/gate-opt/`). The gate binary was `bin/pilot_gate`
(commit `0ed3cc1a59`); the frozen binary differs from it only by two byte counters (`final_row_bytes`,
`sweep_bytes`); the gate was therefore **re-run on the frozen binary** `bin/pilot_opt2` (`raw/gate-opt2.log`, reference records reused): **GATE PASS** again, every set EQUAL, faults 1–5 detected (60 / 40 / 40 / 40 / 38), closedness check 314 closures, 708,997 reductions, 0 non-zero.

No-1 proxies: every A/B run above has identical `dims` (and counters) across arms; the t = 14 c = 5 profile equals the
pod's hky2 record `[0,0,2,154,4323,47849]`. Unit tests (`harness/unit`, `cargo test --release`): 18 pass, including
`epochs_and_sweep_do_the_same_reductions` (planted and random, t = 3, 4, split and unsplit, c = 3..5, every
combination of np/sweep against v2: equal dims, one, muladds, rows_reduced, entries, zero, reducers; asserts the
sweep and phase 0 actually ran) and `the_closedness_check_passes_and_sees_a_fault`.

### 4.2 Soundness

**No optimisation of this round skips work.** Every product the base engine reduces is reduced, against the same
echelon, in the same order of batches; only *how* a single reduction is carried out changes (coordinates, order of
subtraction within it, memory layout, which thread does which columns). Why the output is identical:

- *The reduced vector is unique.* For a fixed echelon E (rows with distinct leading columns, each row zero before its
  lead), the vector `red(r)` ∈ `r + span(E)` that is zero at every pivot column is unique (two such differ by an
  element of span(E) zero at every pivot column, and a non-zero element of span(E) is non-zero at the lead of its
  first row with non-zero coefficient). v2's `eliminate` returns it; so does the new one: phase 0 subtracts, for
  each epoch-pivot column j in the product, the row of j times the product's coefficient at j — final, because no
  stored tail has an epoch-pivot column (invariant: after an interreduction every tail is free of pivots; rows
  inserted later are fully reduced; `new_epoch` and `to_epoch` assert it) — and the cascade then clears the other
  pivot columns in increasing order as before. Lazy sums below p² taken in a different order give the same residue.
- *Same counters.* The coefficient of each echelon row in `r − red(r)` is unique (the rows are independent), so the
  set of rows subtracted with a non-zero factor, hence `muladds` (Σ of their term counts), `reducers`, the zero
  count and every downstream decision (yield, batch size, buckets by term count, interreductions) is the same.
  Measured equal in every run.
- *Storage.* A row in epoch coordinates is the same vector (`row_full` maps it back exactly); the interreduction's
  result rows are unique for the same reason, so stored vectors agree after every epoch.
- Closedness of the final space therefore follows exactly as for v2 (its README, "Why interreduction keeps the
  answer"): every inserted row of degree < c has all its products formed and reduced into the span.

**Closedness check** (kept available as asked, `--hky3-verify-closed on`): after a closure without 1, the echelon is
copied to part coordinates (one epoch of every column, so neither the epoch code nor the sweep is used) and every
generator, and `x_v · g` for **every** variable v (not only one of each pair the symmetry swaps; both projections)
and every final row g of degree < c, is reduced and must be zero. Measured: all 314 gate closures without 1,
**708,997 reductions, 0 non-zero**; t = 10 c = 5 proxy: 126,521 reductions, 0 non-zero (64 s, 2.5× the closure).
Negative controls: fault 1 (products by x₀ not formed) leaves 688 non-zero at t = 10 c = 5 and 2,104 over the gate's
fault-1 cells; it also flags faults 2–5 (44 / 96 / 220 / 200). Cost: one extra echelon copy (≈ 2× memory) and ≈ 2.5×
the closure's time — fine for gate cells and k = 14; it does not fit at k = 16 on 256 GB.

## 5. Rejected or not pursued — with numbers

| idea | evidence | verdict |
|:--|:--|:--|
| skip products that must reduce to zero (F5/syzygy-type criteria, x_a(x_b r) = x_b(x_a r)) | reductions to zero are 51–57 % of all muladds, 60–79 % of reductions; their **number** is products − rank, fixed by the final space, so only skipping can cut it | **not adopted**: it skips work, so the paper would need the closedness check, which costs ≥ the closure (measured 2.5×) — no net gain for a certified run; no proof written |
| probabilistic membership test (random functionals on W^⊥) | — | not sound with certainty; rejected outright |
| interreduction threshold 0.05 / 0.1 / 0.5 | ±10 % muladds (§1) | default 0.25 kept |
| dense threshold 1/2 in epoch coordinates | row bytes −0..6 %, wall +2..7 % (§3.2) | rejected |
| reducer sharing in blocks of 8 items (instead of the whole batch) | 1.9× sharing vs 19× | superseded by the batch-wide sweep |
| multi-RHS interleaved accumulators `acc[col][k]` over all columns | accumulators K × n × 8 B do not fit cache at k ≥ 14 | rejected by reasoning; the sweep tiles epoch columns instead |

## 6. Build, run, files

```bash
export CARGO_TARGET_DIR=/home/user/crypto-wt/engine-opt/target
cd research/hky_engine_opt_20261118/harness
cargo build --release --example pkm_tower_pilot --example hky_engine_equality
(cd unit && cargo test --release)          # 18 tests
cp $CARGO_TARGET_DIR/release/examples/pkm_tower_pilot ../bin/pilot_opt2
cd .. && ./run_gate.sh raw/gate-X --hky3-verify-closed on
./ab.sh raw/ab_X 3 "10:5 12:5 14:5" "base:bin/pilot_base:" "opt:bin/pilot_opt2:" && ./ab_summary.sh raw/ab_X/runs.jsonl
pod/scale_bench.sh bin/pilot_base bin/pilot_opt2 OUT     # on a 32-vCPU pod; DRY=1 prints the plan
```

`pod/scale_bench.sh` (bash + coreutils only): ladder cell planted t = 10 c = 6 (no 1, the closure runs to the end)
at RAYON_NUM_THREADS = 1, 2, 4, 8, 16, 32 (64 if present) and main cell the study draw t = 12 c = 6 (190,051 columns, 1 only as the last pivot;
estimate 5–15 min at 32 threads) at 16, 32 (64); both arms ABBA-interleaved; per run one JSON line: wall, user and
system CPU (bash `times`), utilisation, peak RSS (VmHWM polled from /proc), load, cap/timeout, and the pilot's own
record (dims, muladds, rows_reduced, entries, every hky3 stat) embedded. Cells and thread lists are environment
variables; ~70 min with the defaults at the predicted rates. Smoke-tested here on t = 8 c = 5.

**For the merge with `wt/paper-engine`** (checkpoint/resume and `hky.basis/v1` export): `Echelon` has a new field
`np_cols: Vec<u32>` that must be checkpointed with `pivot` and `rows` (rebuild with `Echelon::with_epoch(pivot,
rows, np_cols)`; `of_col` is derived). Stored rows' tails are in epoch coordinates, so the export must take each row
with `Echelon::row_full(ri)` (part columns), not the raw `RowS`. The sweep's buffer is scratch (all zero between
batches) and needs no checkpoint. `Fast2Options` gains `np_epochs`, `sweep`, `sweep_ratio`, `sweep_tile_bytes`,
`sweep_max_words`, `verify_closed` (none changes the result); `Fast2Stats` gains `core_muladds`, `zero_muladds`,
`swept_batches`, `verify_products`, `verify_nonzero`, `verify_s`, `final_row_bytes`, `sweep_bytes`.

Files: `README.md`, `PREDICTIONS.md`, `pod/scale_bench.sh`, `cell1.sh`, `ab.sh`, `ab_summary.sh`, `run_gate.sh`,
`harness/` (+ `unit/`), `raw/` (profile traces, A/B runs, gate log and records).
