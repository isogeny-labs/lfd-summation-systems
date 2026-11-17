# A GPU path for the split-presentation closures with m >= 4 summands: what was measured, what was built

Date: 2026-11-17. Branch `wt/gpu-split-m`. Labels: **measured** (a command of this directory
produced it on the machine of section 1), **given** (from the task or another report, not
re-measured), **estimate** (a model, stated), **argument** (a proof sketch, not machine-checked).
**No GPU was available: every statement about a real device is an estimate, written before any GPU
run.** No result here is a new last-fall-degree value; no cloud resource was created.

## 0. Verdict

1. **A dense GPU engine does not beat the sparse CPU engine on these cells; at the one target cell
   it can hold it is predicted to tie.** The sparse engine `hky2` of the sibling repository does
   the unit closures with 66 to 450 times fewer multiply-adds than the best dense formulation
   (measured, section 2.3), a 141 GB card does exact float64 products roughly 200 to 900 times faster
   than 32 cores do sparse row operations (estimate, section 8), and the dense basis of a closure that
   contains 1 needs about `n_q^2/4` cells per part. Result: (5,3) at D = 8 is predicted at 2 to 6
   minutes on an H200 against 257 s on 32 vCPU (given); (7,2) at D = 7 the same if 1 enters there;
   (7,2) at D = 8, (6,3) at D >= 8, (8,2) and (9,2) do not fit on any single device (260 GiB to
   2,500 GiB, section 4.2) while the sparse engine reaches them on a CPU host.
2. **What was built** is still the fastest engine in *this* repository and its first one that runs
   m >= 4 on a device: `csd/ic_bounds_split_gpu.py`, compact reduced echelon form + summand-swap
   split behind the numpy/CuPy/fake-device layer of `csd/ic_bounds_big.py`. On CPU it is 2.8 to 26
   times faster than `InPlaceClosure` on the shared cells (interleaved A/B, section 7), reproduces
   every certified value it was run on (section 6) and runs unchanged on the fake device.
3. **The degree lift carries over** to the split presentation with free auxiliaries (argument,
   section 5) **but saves little**: the space it eliminates in is 65 % to 86 % of the columns of
   the direct closure (measured on four cells), because V_D has pivots on only a third of the
   top-degree monomials. For m = 2 at depth 14 that fraction was about 7 %.
4. **Recommendation.** Decide the listed cells with the sparse CPU engine on a large-memory host.
   Use this engine as the independent second instrument for cells up to about 230,000 columns, and
   the GPU job of section 9 if a device time for (5,3), D = 8 is wanted as a measurement. The only
   GPU design with a chance of a real gain is the hybrid of section 4.1(a); its ceiling is the share
   of time `hky2` spends in dense rows, measured in section 2.4.

## 1. Machine, files, commands

Machine: 4 shared vCPU (Intel Xeon 2.1 GHz, AVX-512), 15 GiB, no GPU, no CuPy, numpy 2.5.3 with
OpenBLAS 0.3.34 (4 threads), Python 3.13. Load average is recorded with every timed run.

| file | role |
|---|---|
| `csd/ic_bounds_split_gpu.py` | the engine (`Coordinates`, `PartBasis`, `SplitClosure`, `closure`) |
| `scripts/ic_split_gpu_gate.py` | equality gate against `csd/ic_bounds_split_m.py`, negative controls |
| `scripts/ic_split_gpu_pod.py` | the pod job (section 9) |
| `tests/test_ic_bounds_split_gpu.py` | 12 tests (gate on small cells, controls, fake device, batch echelon, checkpoint) |
| `experiments/gpu-split-m/profile_reference.py` | profile of `InPlaceClosure` (section 2.1) |
| `experiments/gpu-split-m/symmetry_group.py` | which summand/auxiliary permutations are symmetries (section 3) |
| `experiments/gpu-split-m/lift_regular_count.py` | what the lift would save (section 5) |
| `experiments/gpu-split-m/ab_cpu.py` | interleaved CPU A/B (section 7) |
| `experiments/gpu-split-m/raw/` | every record quoted below |

```bash
python3 -m unittest tests.test_ic_bounds_split_gpu                      # or pytest tests/test_ic_bounds_split_gpu.py
python3 scripts/ic_split_gpu_gate.py --controls --out experiments/gpu-split-m/raw/gate-live.json
python3 scripts/ic_split_gpu_gate.py --live '' --recorded '4,3:7;3,5:6;6,2:6;4,4:6;5,3:3-6' --arms cpu:split \
        --out experiments/gpu-split-m/raw/gate-recorded.json
python3 experiments/gpu-split-m/ab_cpu.py --out experiments/gpu-split-m/raw/ab_cpu.jsonl
```

Existing engines are untouched: `csd/ic_bounds_split_m.py`, `csd/ic_bounds_m3_batched.py` and
`csd/ic_bounds_big.py` have no diff on this branch.

## 2. Profile of the problem (measured)

### 2.1 Two regimes

Instances: the certified ones of `experiments/ic-bounds-split-m/` (first line of each cell file).
`profile_reference.py` runs `InPlaceClosure` unchanged under a timing subclass
(`raw/profile_reference.jsonl`).

| cell, D | columns | 1 in V_D | rank when it stops | products per pivot | products reducing to 0 | density of a raw product | density of the stored basis, free columns | time in `_reduce` / `_absorb` |
|---|--:|:--|--:|--:|--:|--:|--:|--:|
| (4,2), 5 | 1,023 | no | 357 | 2.2 | 17 % | 1.8 % | 3.9 % | 28 % / 59 % |
| (4,2), 6 | 1,792 | yes | 1,792 (100 %) | 3.8 | 32 % | 11.8 % | (no free column left) | 18 % / 78 % |
| (5,2), 5 | 3,598 | no | 637 | 2.1 | 14 % | 0.5 % | 0.6 % | 29 % / 66 % |
| (5,2), 6 | 7,937 | yes | 7,563 (95 %) | 2.5 | 34 % | 8.2 % | 60 % | 10 % / 88 % |
| (4,3), 6 | 9,740 | no | 2,517 | 2.9 | 27 % | 0.4 % | 1.0 % | 24 % / 72 % |

* **A closure without 1 is sparse and small**: rank 18 % to 35 % of the columns, basis 0.6 % to 4 %
  dense. A dense engine wastes two orders of magnitude there, on any hardware.
* **A closure with 1 is the whole cost and is dense at the end**: 1 is found when 95 % to 100 % of
  the columns are pivots, and it is the closure that decides a cell. In the committed records the
  last closure takes 50 to 400 times longer than all the others of its cell together.
* `InPlaceClosure` spends 66 % to 88 % of its time in `_absorb` (full-width update of a
  `rank x n` basis, a Python back-substitution of the batch), not in the reduction products.

### 2.2 How the fill evolves in the closure with 1 (5,2), D = 6

Reduced echelon basis of `InPlaceClosure`, sampled every four batches (`trace` of the same record):

| rank | of which pivot degree < D | free columns | density on the free columns | nonzeros | cells `rank x free` |
|--:|--:|--:|--:|--:|--:|
| 952 | 288 | 6,985 | 0.2 % | 1.6e4 | 6.7e6 |
| 2,823 | 743 | 5,114 | 0.8 % | 1.2e5 | 1.44e7 |
| 3,695 | 880 | 4,242 | 5.7 % | 9.0e5 | 1.57e7 |
| 4,325 | 1,070 | 3,612 | 24 % | 3.8e6 | 1.56e7 |
| 5,962 | 1,854 | 1,975 | 55 % | 6.5e6 | 1.18e7 |
| 7,563 | 3,224 | 374 | 60 % | 1.7e6 | 2.8e6 |

The basis is below 1 % dense until about a third of the columns are pivots, then fills within a few
batches to 25 % to 60 %. The dense block is real in the second half; `rank x free` peaks at
`n^2/4` (1.57e7 = 0.25 n^2) when the fill is still 6 %.

### 2.3 The sparse CPU engine on the same cells

`hky2` (sibling repository, `src/cryptanalysis/hky_fp_tower_fast.rs`; frozen binary
`bin/pilot_candidate` of its report `research/hky_engine_speed_20261117/`, `--kinds kummer --controls tower --engine hky2`, p = 2013265921,
its own random draw, 4 threads on this machine). Its dimension profiles equal the certified ones of
this repository on every closure without 1 ((5,2): 51, 637; (4,3): 30, 424, 2,517; (6,2): 64, 995,
7,195; (5,3): 39, 711, 5,243, 27,838), and its part sizes and part ranks equal this engine's
((5,2), D = 5: [499, 138]; (4,3), D = 6: [1,607, 910]): two instruments, different primes, same split.

| cell, D | columns | 1 | `hky2` wall s (4 threads) | `hky2` multiply-adds | rows stored dense (> 1/3) / all | this engine, CPU, s | this engine, multiply-adds issued | ratio of multiply-adds |
|---|--:|:--|--:|--:|--:|--:|--:|--:|
| (5,2), 6 | 7,937 | yes | 0.19 | 3.07e8 | 1,732 / 7,107 | 3.9 | 2.45e10 | 80 |
| (4,3), 7 | 18,324 | yes | 1.7 | 6.81e9 | 7,119 / 16,752 | 40.4 | 4.52e11 | 66 |
| (6,2), 7 | 61,183 | yes | 9.1 | 2.32e10 | 6,219 / 39,937 | 803 | 1.05e13 | 450 |
| (5,3), 7 | 106,970 | no | 2.1 | 2.72e8 | 0 / 27,838 | not run (9.7 GiB, section 6) | | |
| (6,3), 7 | 440,772 | no | 14.4 | 3.68e8 | 0 / 60,140 | does not fit here | | |
| (5,3), 8 | 226,325 | yes | 257 on 32 vCPU (**given**) | not reported | | | | |

`hky2` does 4e8 to 9e8 multiply-adds per second and thread on the closures with 1 (measured,
multiply-adds / wall / 4). It stops (6,2), D = 7 at rank 39,937 of 61,183 (65 %); this engine and
`InPlaceClosure` find 1 at 95 % to 100 %.

Why the counts differ by two orders of magnitude: `hky2` keeps an unreduced echelon and subtracts
only the pivot rows a product actually meets (a product is 8 % to 12 % dense). The dense form
multiplies every product by every pivot row (`P[:, piv] @ B`) and keeps the basis reduced, which
costs a rank-s update of every stored row per batch. Skipping pivot rows that no row of the batch
meets (implemented, `hit_frac`) removes 65 % to 80 % of the dense count ((4,3), D = 7: 4.5e11 issued
of 1.28e12; (6,2), D = 7: 1.05e13 of 5.2e13), not 99 %.

### 2.4 Share of dense rows in the sparse engine

`hky2` stores a row densely once it is more than one third full and subtracts it with a SIMD loop.
Interleaved A/B on this machine, 3 rounds, same binary, `--hky2-dense 1/3` (default) against
`--hky2-dense 1000/1` (never dense), wall milliseconds of the closure with 1
(`raw/hky2_dense_rows_ab.txt`; the multiply-add counts are identical in both arms):

| cell, D | dense rows / all | default ms (3 runs) | never dense ms (3 runs) | gain from dense rows |
|---|--:|--:|--:|--:|
| (4,3), 7 | 7,119 / 16,752 | 1703, 1648, 1708 | 2595, 2562, 2517 | 1.50 |
| (6,2), 7 | 6,219 / 39,937 | 9186, 8783, 9060 | 10797, 9888, 10025 | 1.11 |

Dense rows save 0.86 s of 2.56 s and 0.97 s of 10.0 s. If the dense loop is f times faster than the
sparse one on those rows, the time still spent in dense rows is `saved / (f - 1)`: for f = 2 that is
50 % and 11 % of the run, for f = 4 it is 17 % and 4 %. **A hybrid that moved every dense row to a
GPU at zero cost would gain at most 2x on (4,3) and 1.1x on (6,2)** (f >= 2 assumed; f was not
measured). On the larger cell the time is in the sparse rows.

## 3. Symmetries (measured, and checked at run time)

Candidates: a permutation of the m summands together with a permutation of the m - 2 auxiliaries
(the degree-preserving monomial permutations that permute the shift variables). For each of the
48 / 720 / 48 candidates of (4,2), (5,2), (4,3) at D = 5, `symmetry_group.py` tests exactly
(a) whether it preserves the span of the inputs and (b) whether it maps the closed V_5 into itself
(`raw/symmetry_group.json`):

| cell | candidates | preserve the inputs | preserve V_5 |
|---|--:|--:|--:|
| (4,2) | 48 | 2 | 2 |
| (5,2) | 720 | 2 | 2 |
| (4,3) | 48 | 2 | 2 |

**The group is Z/2: the exchange of x_1 and x_2**, in all three cells, and V_D has no symmetry the
inputs do not show. Reason: g_1 = S_3(x_1, x_2, z_1) is symmetric in x_1, x_2; the chain
z_1 -> z_2 -> ... orders the other summands, and the last input contains the target, so the chain
cannot be reversed. The full S_m symmetry of the ideal in the x variables is broken by the
auxiliaries. Not examined: symmetries that are not monomial permutations.

Consequence: one split into two parts, of 56 % to 67 % and 44 % to 33 % of the columns (table of
section 4.2). The engine checks at run time, exactly, that sigma is an involution, preserves
degrees, intertwines the shift of u_{1,j} with that of u_{2,j} on every monomial, commutes with
every other shift, and fixes every input; otherwise it reports `symmetry: refused: <reason>` and
runs unsplit (test `test_asymmetric_input_refuses_the_split`). No isotypic decomposition beyond
this one exists for these candidates, so the dimension can be cut by 2, not more.

## 4. Designs

### 4.1 Comparison

Model for a closure that contains 1, per part of n_q columns, r = rank, rho = products reduced per
pivot (2.3 to 3.8 measured):

* dense reduced form, all columns stored (`InPlaceClosure`): reduce rho n^3/2 + update n^3/2;
  memory n^2 cells.
* dense compact form (this engine): reduce rho * int r (n - r) dr = rho n^3/6, update n^3/6;
  memory max r (n - r) = n^2/4 cells; with the split, sum over the two parts. Measured multiply-adds
  issued: c * sum n_q^3 with c = 0.33, 0.16, 0.28, 0.45, 0.14 for (4,2), (5,2), (4,3), (3,5), (6,2) (model
  with rho = 3 and no skipping: 0.67).
* sparse unreduced echelon (`hky2`): measured only, section 2.3.

| design | what it is | number that decides | verdict |
|---|---|---|---|
| (b) dense on device, compact + split | this engine | 66 to 450 times the multiply-adds of `hky2`; a GPU is 200 to 900 times faster per multiply-add than 32 cores (estimate, section 8); memory n_q^2/4 | **built**: ties at (5,3) D = 8 (estimate), cannot hold anything larger. Isotypic split: factor 2 only (section 3). Out-of-core tiling: not built, see below |
| (a) hybrid, sparse on CPU + dense block on GPU | Faugere-Lachartre layout: sparse pivot rows on the host, the rows that fill (section 2.2) as a dense reduced block on the device | ceiling = share of `hky2` time in dense rows: at most 2x on (4,3), 1.1x on (6,2) (section 2.4) | **not built**. It needs the sparse engine as its base; that engine is Rust in a read-only repository. The only design that could beat it, and by at most that factor |
| (c) sparse on device | CSR echelon, dense accumulator per batch row, gather/scatter | row operations are bandwidth-bound: about 20 bytes per multiply-add, so 1.5e10 to 3e10 per second at 10 % to 20 % of 3 TB/s (estimate) against 1.3e10 to 2.9e10 for 32 cores at the measured 4e8 to 9e8 per thread | **rejected**: parity at best, and the unreduced echelon is a data-dependent chain of subtractions per row, which a device cannot batch without reducing the echelon (which is design b) |
| (d) degree lift | V_{D+1} n M_<=D and "1 in V_{D+1}" from a closed V_D | eliminates in 65 % to 86 % of the columns (section 5) | hypotheses hold (argument); **saves at most a third**, not built |
| float32 / 32-bit integers | half the memory, twice the rate | exact float32 products need n p^2 < 2^24; unit instances need p of order 1e5 (protocol amendment 1) | **rejected** for float32. Montgomery on 32-bit integers has no BLAS on the device: a custom kernel, bandwidth-bound like (c) |

Out-of-core tiling for (b): the basis is already stored in tiles of 4,096 rows and every kernel
walks the tiles, so spilling tiles to host memory is a small change. It is not built because the
multiply-add count, not the memory, is what rules the next cells out: (7,2) at D = 8 would need
c * sum n_q^3 = 5e15 to 1.1e16 multiply-adds (7 to 30 minutes of pure GEMM at 6e12 to 12e12 per
second, estimate) plus 260 GiB moved through the device for every batch; (6,3) at D = 8 needs
5e16 to 1e17 and 1,421 GiB.

### 4.2 Memory per cell (upper bound when 1 enters) and what fits

`n_+`, `n_-` exact (counted). "compact + split" = 8 bytes * (n_+^2 + n_-^2)/4 * 1.14 (dead columns
are dropped at 12.5 %); add a working set of about 6 * batch * n_+ * 8 bytes (12 GiB for (5,3),
D = 8 at batch 2,048). Checked against a run: (4,3), D = 7 peak 0.37 GiB measured, bound 0.40.

| cell, D | columns | n_+ | n_- | `InPlaceClosure` n^2 | compact, no split n^2/4 | **compact + split** | sum n_q^3 | smallest device |
|---|--:|--:|--:|--:|--:|--:|--:|:--|
| (6,2), 7 | 61,183 | 40,127 | 21,056 | 28 GiB | 7.0 | **4.4 GiB** | 7.4e13 | 16 GB |
| (4,4), 7 | 91,276 | 49,708 | 41,568 | 62 | 15.5 | **8.9** | 1.9e14 | 16 GB (batch 512), 48 GB |
| (5,3), 7 (no 1: rank 27,838) | 106,970 | 62,832 | 44,138 | 85 | 21.3 | 12.5; **actual 9.7** | | 48 GB |
| **(5,3), 8** | 226,325 | 131,271 | 95,054 | 382 | 95 | **56** | 3.1e15 | 80 GB, 141 GB |
| **(7,2), 7** | 195,604 | 131,346 | 64,258 | 285 | 71 | **45** | 2.5e15 | 80 GB, 141 GB |
| (6,3), 7 (no 1: rank 60,140) | 440,772 | 266,898 | 173,874 | 1,447 | 362 | 216; **actual 104** | | 141 GB (pointless: 14 s on 4 vCPU with `hky2`) |
| (7,2), 8 | 471,041 | 311,297 | 159,744 | 1,653 | 413 | **260** | 3.4e16 | none |
| (8,2), 7 | 532,220 | 365,860 | 166,360 | 2,110 | 528 | **343** | 5.4e16 | none |
| (6,3), 8 | 1,136,799 | 675,969 | 460,830 | 9,628 | 2,407 | **1,421** | 4.1e17 | none |
| (9,2), 7 | 1,280,760 | 900,332 | 380,428 | 12,222 | 3,055 | **2,029** | 7.9e17 | none |
| (8,2), 8 | 1,466,345 | 989,163 | 477,182 | 16,020 | 4,005 | **2,561** | 1.1e18 | none |
| (7,3), 8 | 4,384,467 | 2,678,523 | 1,705,944 | 143,227 | 35,807 | **21,414** | 2.4e19 | none |
| (9,2), 8 | 3,981,645 | 2,741,023 | 1,240,622 | 118,118 | 29,529 | **19,222** | 2.2e19 | none |

16 GB: through (6,2) D = 7 and (4,4) D = 7. 48 GB: adds the closures without 1 up to (5,3) D = 7.
80 GB and 141 GB: add (5,3) D = 8 and (7,2) D = 7. Nothing else in the goal list fits on one device
in dense form, with or without the split.

## 5. The degree lift (argument; counts measured)

`experiments/ic-bounds-k14/LIFT.md` proves the lift for m = 2 in the chain ring. Its proof uses the
ring through four properties only; each holds for `SplitSpace` with free auxiliaries:

1. *Columns sorted by degree, a shift raises the degree by at most 1* (Lemmas 1, 2). True: a chain
   shift raises the weight by 1 or lowers it (carry), an auxiliary shift raises it by exactly 1.
2. *A shift maps a monomial to one monomial with coefficient 1, injectively* (Lemma 4: the image of
   a row with a regular pivot has one regular top coordinate). True: `SplitSpace.cols` is a table
   of single monomials; injectivity is asserted when the engine builds its tables.
3. *Every monomial of degree D + 1 is a shift of a monomial of degree D, and only of such* (Lemma
   5, s_r >= 1). True: remove a set bit of the mask or lower a positive auxiliary exponent. The
   number of such divisors is popcount(mask) + #{i : e_i > 0}, between 1 and D + 1, not always
   D + 1; the proof does not use the number, but the engine's assertion "every regular column is
   met with all its D + 1 generators" would have to count per column.
4. *sigma u_{1,j} sigma = u_{2,j}* (Lemma 3). True, and the remaining variables commute with sigma,
   which adds the trivial case pi_q v W^p = v W^p for p = q and 0 otherwise.

So statements 1 to 4 of the theorem hold for the split presentation (**argument**, not
machine-checked, not implemented, not tested numerically here).

What it would save (`lift_regular_count.py`, definitions of LIFT.md section 1, closed V_D of the
certified instance, `raw/lift_regular_count.jsonl`):

| cell, D -> D+1 | dim V_D | weight-D monomials | of them regular | weight-(D+1) monomials | of them regular | dim T (both parts) | columns of the direct closure at D+1 | T / columns |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| (4,2), 5 -> 6 | 357 | 522 | 246 | 769 | 268 (35 %) | 1,167 | 1,792 | 65 % |
| (5,2), 5 -> 6 | 637 | 2,223 | 444 | 4,339 | 447 (10 %) | 6,853 | 7,937 | 86 % |
| (4,3), 6 -> 7 | 2,517 | 5,282 | 1,914 | 8,584 | 1,766 (21 %) | 14,041 | 18,324 | 77 % |
| (6,2), 6 -> 7 | 7,195 | 16,508 | 5,300 | 34,832 | 6,201 (18 %) | 47,787 | 61,183 | 78 % |

The lift replaces arithmetic by a gather only on the regular top columns: 10 % to 35 % of them.
The elimination space T keeps 65 % to 86 % of the columns, so the cubic cost falls by a factor
1.6 to 3.6 at best, against an implementation as large as the engine. In the m = 2 run the same
ratio was about 7 % (LIFT.md section 3, correction (d)); LIFT.md section 10 already notes that the
lift loses when most top columns are special. **Not worth building for these cells.**

## 6. Correctness gate (measured)

`scripts/ic_split_gpu_gate.py`; candidate arms `cpu:split`, `cpu:plain`, `fake-gpu:split`,
`fake-gpu:plain` (`plain` = one part, no symmetry). A comparison is EQUAL when the number of
columns, dim V_D, the low-degree profile and "contains 1" agree, the candidate's own verification
pass succeeded on every closure without 1 (every operator image of every final row of degree < D
reduces to zero, every input in the span), and the exact-arithmetic guard never fired. On closures
without 1 of at most 5,000 columns the **spaces** are compared as well: every candidate basis row,
lifted from orbit coordinates to monomials, must lie in the reference V_D.

**Live reference** (`raw/gate-live.json`): `InPlaceClosure` run here, on (4,2) D = 3..6,
(5,2) D = 3..6, (4,3) D = 3..6, (3,5) D = 3..5, (6,2) D = 3..5, (4,4) D = 3..5, four arms each:
**84 of 84 comparisons EQUAL, 60 of them with the spaces compared, 42 on the fake device**. The reference also equals its committed record on each.

**Committed reference records** of the same instances (`raw/gate-recorded*.json`), where the
reference needs a pod:

| cell, D | columns | 1 in V_D | dim V_D | arm | candidate s (this machine, 4 vCPU) | reference s (its record) | result |
|---|--:|:--|--:|---|--:|--:|:--|
| (4,3), 7 | 18,324 | yes | 18,324 | cpu:split | 40.4 | 959.3 | EQUAL |
| (4,3), 7 | | | | fake-gpu:split | 39.9 | | EQUAL |
| (4,3), 7 | | | | cpu:plain | 108.7 | | EQUAL |
| (3,5), 6 | 17,548 | no | 7,007 | cpu:split | 15.1 | 67.4 | EQUAL |
| (3,5), 7 | 33,932 | yes | 33,932 | cpu:split | 248.8 | 2,112.2 | EQUAL |
| (6,2), 6 | 26,351 | no | 7,195 | cpu:split | 42.3 | 59.9 | EQUAL |
| (6,2), 7 | 61,183 | yes | 61,183 | cpu:split | 803.2 | 4,203.5 | EQUAL |
| (4,4), 6 | 39,796 | no | 3,672 | cpu:split | 19.0 | 24.8 | EQUAL |
| (5,3), 3..6 | 1,045 .. 44,366 | no | 1, 39, 711, 5,243 | cpu:split | 0.001, 0.02, 0.68, 35.2 | (pod log, dims only) | EQUAL |
@@G537@@

The reference times are from the pods of `experiments/ic-bounds-split-m/` and are not an A/B with
this machine; the A/B is section 7. The low-degree profile is compared wherever the record has one
(all rows except (5,3), whose record holds dimensions only).

Certified values reproduced by the candidate's own ladder (closed and verified without 1 below, 1
at the value): (4,2) = 6, (5,2) = 6, (4,3) = 7, (3,5) = 7, (6,2) = 7. Not reproduced here: (4,4) = 7 (the
closure with 1 has 91,276 columns: 8.9 GiB and an estimated 3 hours on 4 vCPU), and (5,3) "no 1 at
D = 7" (9.7 GiB of basis on a 15 GiB machine: not attempted); for (5,3) the closures at D = 3..6 equal
the record (1, 39, 711, 5,243). Both are in the pod job's default list.

**What the gate can and cannot see.** On a closure with 1 the reported dimensions are those of the
whole space whatever happened before; the power of the gate is in the closures without 1 and in
the degree at which 1 enters.

**Negative controls** ((4,2), D = 5, reference dim 357; `--controls`, `raw/gate-live.json`):

| control | arm | detected | by dimensions alone | what differs |
|---|---|:--|:--|---|
| one operator never applied in the closure loop (`fault='drop-operator'`) | split | yes | yes | dim 243, and the engine's own verification pass fails (it uses the full operator list) |
| one gather table rotated by one column (`fault='wrong-gather'`) | split, plain | yes | yes | 1 enters at D = 5 (dim 1,023) |
| one coefficient of g_1 changed by 1 (candidate only) | split, plain | yes | **no** | same dimensions 357, another space: only the comparison of spaces sees it |
| no fault | split | passes | | |

The changed coefficient is the control suggested in the task; it confirms what the sibling report
found: a generic change of one coefficient gives a system with the same dimension profile, so a
dimension-only gate cannot fail on it. That is why the gate compares spaces on small cells.

**Fake device.** Every `fake-gpu` arm runs the same code with arrays that raise on any mix of host
and device arrays in an operation, an index or an assignment (class `_Dev` of `ic_bounds_big`):
42 live comparisons and the (4,3), D = 7 closure pass in that mode. It cannot check CuPy semantics
that differ from numpy (section 10).

## 7. CPU A/B against `InPlaceClosure` (measured)

`ab_cpu.py`: both engines through their `closure()` (verification pass included), numpy backend,
arm order rotated every round, 3 rounds (2 for the cell whose reference run takes 105 s),
`raw/ab_cpu.jsonl`. Median seconds (min - max); load average before the runs in the last column
(4 shared vCPU; the load is the BLAS threads of the previous run and the other sessions on the
host).

| cell, D | 1 | rounds | `InPlaceClosure` s | this engine, split, s | this engine, no split, s | ratio (split) | ratio (no split) | load before |
|---|:--|--:|--:|--:|--:|--:|--:|--:|
| (4,2), 5 | no | 3 | 0.162 (0.154 - 0.179) | 0.057 (0.056 - 0.061) | 0.064 (0.054 - 0.071) | **2.8** | 2.5 | 2.9 - 2.9 |
| (4,2), 6 | yes | 3 | 3.14 (3.13 - 3.36) | 0.261 (0.24 - 0.268) | 0.389 (0.386 - 0.45) | **12.0** | 8.1 | 2.9 - 3.0 |
| (4,3), 5 | no | 3 | 0.689 (0.666 - 0.714) | 0.172 (0.141 - 0.206) | 0.191 (0.183 - 0.243) | **4.0** | 3.6 | 3.1 - 3.2 |
| (4,3), 6 | no | 3 | 20.7 (20.6 - 21.0) | 2.5 (2.35 - 2.59) | 5.6 (5.4 - 5.66) | **8.3** | 3.7 | 2.6 - 3.2 |
| (5,2), 5 | no | 3 | 1.01 (0.927 - 1.06) | 0.19 (0.177 - 0.239) | 0.275 (0.255 - 0.284) | **5.3** | 3.7 | 3.0 - 3.1 |
| (5,2), 6 | yes | 2 | 103 (99.9 - 106) | 3.94 (3.91 - 3.98) | 10.2 (9.95 - 10.4) | **26.1** | 10.1 | 1.6 - 2.9 |

Not interleaved, one run each, reference on another machine (so not a ratio of engines): the
closures with 1 of section 6, 40 s against 959 s, 249 s against 2,112 s, 803 s against 4,204 s.

The gain is 2.8x to 26x and grows with the size; it comes from the compact form (a third of the
multiply-adds, a quarter of the memory traffic), the split (2.2x to 2.6x on the larger cells,
compare the two candidate columns), skipping pivot rows no product meets, and a batch echelon
without a Python loop per row. It is a gain over this repository's dense engine only: `hky2` is
20 to 90 times faster again on the same machine ((5,2): 21x, (4,3): 24x, (6,2): 88x; section 2.3).

## 8. Predictions for the GPU (estimates, written before any GPU run)

Model: `seconds = c * sum n_q^3 / R + host + calls`, with

* `c` = multiply-adds issued per `sum n_q^3`, measured on CPU for the five closures with 1 that
  were run: 0.14 to 0.45 (section 4.1);
* `R` = exact float64 GEMM rate of the device. **Assumed** 6e12 to 12e12 multiply-adds per second
  for an H200/B200-class card (the assumption of `experiments/ic-bounds-k14/README.md` section 9;
  no rate of this engine was measured on a device). An A100-class 80 GB card: assumed 2.5e12 to
  4e12;
* `host + calls` = per echelon round about 3 ms of host elimination (64 rows by up to 320 columns
  in numpy, estimated from the CPU runs) and about 30 device calls at an **assumed** 100
  microseconds; 0.04 to 0.1 rounds per pivot (measured: 0.05, 0.10, 0.04 on (4,3), (3,5), (6,2)).

| cell, D | assumption | multiply-adds | GEMM s (H200) | host + calls s | **predicted wall, H200** | predicted, A100 80 GB | device memory, upper bound | CPU sparse engine |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| (6,2), 7 | 1 enters (certified) | 1.05e13 (counted on CPU) | 0.9 to 1.8 | 14 | **15 to 20 s** | 20 to 25 s | 4.4 + 2 GiB | 9.1 s on 4 vCPU (measured) |
| **(5,3), 8** | 1 enters (stage B reading) | 4.4e14 to 1.4e15 | 37 to 234 | 66 to 132 | **2 to 6 min** | 3.5 to 11 min | 56 + 12 GiB | 257 s on 32 vCPU (given) |
| **(7,2), 7** | if 1 enters at 7 | 3.5e14 to 1.1e15 | 30 to 190 | 60 to 115 | **1.5 to 5 min** | 3 to 9 min | 45 + 12 GiB | not measured; about 3 min on 32 vCPU by size |
| (7,2), 8 | if not | | | | **does not fit** (260 GiB) | | | 27 to 40 min on 32 vCPU (power law n^2.5..3 from the 257 s point; estimate) |
| (6,3), 7 | no 1 (measured with `hky2`) | | | | fits a 141 GB card (104 GiB); minutes | no | 104 + 20 GiB | 14.4 s on 4 vCPU (measured) |
| **(6,3), 8** | | | | | **does not fit** (1,421 GiB) | | | 4 to 9 h on 32 vCPU (estimate, same power law) |
| **(8,2), 7 / 8** | | | | | **does not fit** (343 / 2,561 GiB) | | | 0.6 to 0.9 h / 8 to 19 h (estimate) |
| **(9,2), 7 / 8** | | | | | **does not fit** (2,029 / 19,222 GiB) | | | 5.5 to 13 h at D = 7 (estimate) |

The CPU estimates in the last column are a power law through a single given point and are worth an
order of magnitude at most; the memory of the sparse engine at those sizes was not estimated.
Prediction to be scored by the first GPU run: (5,3), D = 8 on an H200 finishes with 1 in V_8 in
120 to 360 s of closure time with a device high-water mark below 70 GiB; (6,2), D = 7 in 15 to 20 s; the gate passes on the
device; `guard_hits` = 0.

## 9. The pod job

```bash
# image: any CUDA 12.x runtime image with Python >= 3.10, e.g. runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04
pip install 'numpy>=1.26,<2.1' 'cupy-cuda12x>=13.3,<14'   # CuPy 13 wheels for CUDA 12 (not tested here); cc/gcc present (csd/elim.py compiles a C kernel for the reference)
git clone <this repository> && cd lfd-summation-systems && git checkout wt/gpu-split-m
python3 scripts/ic_split_gpu_pod.py --device gpu --out /workspace/gpu-split-m \
    --cells '6,2:3-7;5,3:3-8;7,2:3-8' --budget-s 3600 2>&1 | tee /workspace/gpu-split-m.log
```

Card: 80 GB or more (H100/H200/A100-80); with 48 GB replace the cell list by `'6,2:3-7;4,4:3-7;5,3:3-7'`.
Host: 8 vCPU and 64 GB RAM are enough (the reference closures of the gate and the pivot search run
on the host; a checkpoint is at most `--checkpoint-gb` = 40 GiB).

What it does, in order, all under `--out`: `env.json` (versions, device, memory, `nvidia-smi`, git
commit, source hashes); `selftest.json` (exact arithmetic on the device against Python integers; the
job stops if it fails); `gate.json` (equality gate on the device, arms `gpu:split` and `gpu:plain`
against the CPU reference on (4,2), (5,2), (4,3) D = 3..6 and the recorded (4,3) D = 7, with the
negative controls; **no target cell is run if it fails**); `closures.jsonl` (one record per closure:
status, dim, profile, 1, verification, seconds per phase, device high-water mark from the CuPy
memory pool, counted multiply-adds; after each closure an `{"schema":"optlevel.ops/v1",...}` line on
stdout); `job.log` (every line, UTC time). A closure that passes its time budget or its memory
budget stops with `status: "stopped: ..."` and a checkpoint when it is small enough
(`--resume` continues from it and skips recorded closures). An exception is recorded with its
traceback and the job goes on with the next cell. Exit status: 0 all closures decided, 1 a ladder
stopped, 2 to 4 environment / self-test / gate failure.

Expected per default cell (estimates of section 8): (6,2) through D = 7 under a minute, 6 GiB;
(5,3) through D = 7 a few minutes (the D = 7 closure without 1 is dense-wasteful: 10 GiB), D = 8
2 to 6 minutes, below 70 GiB; (7,2) through D = 6 seconds, D = 7 either 1 enters (1.5 to 5 minutes,
below 60 GiB) or it closes without 1 and D = 8 stops at the memory budget (recorded as a lower bound
D >= 8 for that instance, nothing more). Whole job: under 30 minutes on an H200 if the estimates
hold. The same command with `--device fake-gpu --cells '4,2:3-6'` was run here end to end.

## 10. Unverified, and limits

* **Everything about a real GPU**: no CuPy import was executed. The engine uses only array calls
  that `ic_bounds_big` already used on an H200 (fancy-index gather and assignment, `@`, `any`,
  `argmax` on a boolean array, `concatenate`, in-place arithmetic with `out=`), but this file was
  never run on a device. First suspects if the gate fails there: `argmax` over a boolean array,
  assignment through two index arrays in `PartBasis.expand`, a zero-size tile.
* The device rates, the call latency and therefore every time in section 8 are assumptions.
* The engine is dense: closures without 1 are 100 to 100,000 times more work than in the sparse
  engine (section 2.3). It is not the tool for those.
* `peak_device_gib` is the CuPy pool size on a device and the tile bytes on CPU; the fake device
  reports live bytes.
* Checkpoints are whole-basis `.npz` files; resume was tested at (4,2), D = 5 only.
* The lift argument of section 5 is a reading of the m = 2 proof, not a proof checked by a machine
  or by a run.
* The symmetry statement of section 3 is measured on three cells at D = 5 and one instance each.
* `hky2` numbers of section 2.3 are from its frozen binary with its own instances (prime
  2013265921); they are not registered results of either repository, and the 257 s point is given.
* No pcc claim, benchmark record or leaderboard row is made: this is an engineering comparison in
  a repository that has none of those tools.
