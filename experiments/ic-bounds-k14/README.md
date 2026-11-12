# Closure engine for V_D at k = 14, D = 6 (CPU numpy / GPU cupy)

Status of every statement below is labelled **proof**, **exact** (certified computation, command named),
**empirical**, **estimate** or **unverified**.  A closure that stops growing is not a bound unless the full
sweep passed (`closed: true`).

| file | role |
|---|---|
| `csd/ic_bounds_big.py` | the engine; one self-contained file (Python >= 3.10, numpy; cupy only for `--device gpu`; threadpoolctl only if importable) |
| `scripts/ic_bounds_big_run.py` | driver: `closure`, `plan`, `selftest` run remotely without repository imports; `instance` runs locally |
| `scripts/ic_bounds_big_validate.py` | local validation harness against the reference engine (imports `csd.ic_bounds_d5_verify`) |
| `experiments/ic-bounds-k14/instances/*.json` | instance files consumed by `closure --instance` |
| `experiments/ic-bounds-k14/validation/*.json` | validation records (each carries the sha256 of engine, driver, harness) |
| `experiments/ic-bounds-k14/estimate.py` | the extrapolation behind the estimate table |

## 1. What is computed

Same mathematics and canonical form as the reference engine `csd/ic_bounds_d5_verify.py`: modulo the chains the
ring is F_p[(Z/2^k)^2]; a squarefree chain monomial is an exponent pair (a1, a2) of weight popcount(a1) +
popcount(a2); u_{i,j} adds 2^j to a_i mod 2^k.  V_D is the least subspace of span(M_<=D) containing
NF(g) (and NF(m g) for weight(m) <= D - deg g, `--gens all`, the default; `--gens g` takes NF(g) only) and closed
under f -> u_{i,j} f for every f supported on weight < D.  Columns: weight descending, then key a1 | a2 << k.
With g symmetric, V_D = V+ (+) V- (swap split, orbit coordinates, only the k block-1 shifts, each shifted row
projected to both parts).  Each part is a reduced row echelon basis (RREF) for that column order: pivot columns plus
coefficients on the free columns, stored as float32 tiles, products exact in float64.  `Part.sha()` is the
reference definition (sha256 over pivots sorted, free columns, rows sorted by pivot as uint16 for p < 65536 and
uint32 otherwise), computed from host copies, so hashes are directly comparable with the reference records.

## 2. Engine

* One array module `xp` (numpy or cupy).  Index bookkeeping (pivot lists, column maps, layouts, processed flags)
  is host numpy; only B, batch blocks and residual blocks are on the device.  Device-specific code is the class
  `Device` (transfer, pool, memory info) and the two-line `use_device`.  `--device fake-gpu` runs the identical
  device code path with a numpy stand-in that raises on any mixing of host and device arrays (in operators, ufuncs,
  numpy functions, indices and assignments) and reports the peak of live "device" bytes; it found one slip (a 0-d
  reduction result) during development.
* B is stored in row tiles; a rank-t update touches only the rows that meet the new pivots and only columns in
  [min lead, last nonzero column of L), in place, tile by tile and column block by column block.  New pivot columns
  are deleted lazily (kept exactly zero, checked by `check_rref`) and compacted tile by tile when more than
  `--compact-frac` (0.05) of the physical columns are dead; the device pool is freed after compaction.
  No operation copies a part's whole matrix.
* Sparse-aware reduction.  All rows of one batch (same source part, same shift, same target part) share their support T.
  R = Y[:, free] - Y[:, piv] @ B uses only the pivot columns of T on which the batch is nonzero and only the
  columns of B at or after the smallest such pivot (B[r] vanishes before piv[r]).  `madds` counts what is done,
  `madds_dense` what the dense method (n x rank x free per image, plus update and RREF) would do.
* Batches: every low pivot of one source part not yet processed, up to `min(--batch, budget)` rows, one shift
  and one target part at a time.  `Engine.batch_rows` sizes them from `--work-gb` (R, its filtered copy and the
  in-place block RREF; S, Y, filtered Y, the selected inner block and gathers; 128 MiB for the `modp_` chunks;
  25 % kept for the column-block temporaries of `Engine.col_block`).  The size uses a bound on the physical columns
  derived from the rank, never the lazy layout, so it is the same after a resume.
* Threads on the CPU: BLAS threads, `--threads` (default: the minimum of `len(os.sched_getaffinity(0))` and the cgroup CPU quota,
  `/sys/fs/cgroup/cpu.max` (v2) or the v1 cfs quota, falling back silently; never `os.cpu_count()`);
  the driver sets the BLAS environment before numpy is imported (`--threads` overrides an existing
  `OPENBLAS_NUM_THREADS`; without the flag an existing variable wins).  Elementwise passes are single threaded.
* **Spilling to the host (`--dev-gb G`, default: no limit).**  Tiles live either on the device or, spilled, in host numpy
  arrays.  The tile budget is `G - work - one compaction tile - 1.5 GiB context slack` (`--tile-budget-mb` sets it
  directly, for tests).  Before a part is used as the reduction/update target (`Engine.activate`, called when the
  target loop starts and again in `residual`, `add_residuals` and `combine`) all of ITS tiles are moved to the device,
  after tiles of the OTHER part (highest tile index first, just enough) were moved to the host with `cupy.asnumpy`;
  `cupy.asarray` moves them back, and the pool is freed after each spill.  The same `ensure_room` runs before a new tile
  is allocated and before a compaction replaces a tile.  Reading the rows of an inactive source part (batch source
  block, sampled checks, the hash) works from either side (`Part._rows`, `Part._rows_host`).  Bytes and seconds moved in
  each direction, the number of spills, fetches and target switches are logged (`spill` in every log line and in the
  record).  Pinned memory is not used.  A checkpoint stores tile contents; with a budget the tiles of a resumed run
  start on the host and are fetched when their part becomes the target.  `fake-gpu` raises when its live device bytes
  exceed the budget (`dev_limit`), and `check_term` asserts between batches that the tiles on the device fit the tile
  budget and that the fake device counted them.
* **Target-major order.**  Within a batch all images into part 0 are reduced and added first, then all images into part 1
  (previously shift-major).  The full sweep checks all source rows of both parts against target part 0 first, then against
  part 1 (`sweep.q`, `sweep.pi`, `sweep.pos` in the checkpoint; `rows_checked` counts each source row once per target
  part).  This is always on: one code path.  Why the invariant of section 3.2 is unchanged: the source block S of a batch
  is captured before any image is added, the images are fixed linear functions of S, and adding the residual of each image
  to the basis gives span(basis + images) whatever the order of the images; each low pivot is still processed once, at its
  state when its batch started.  In the sweep the basis does not change while a pass is clean, and any nonzero residual ends the
  pass, is added, and restarts the sweep.
* Random pre-filter (`--compress` of the inherited file): **removed** (see 7).  Plain mode is the only mode.

### Kept, changed, removed from the inherited `csd/ic_bounds_big.py`
The inherited file (about 1,440 lines, never run) was not repaired; the engine was rewritten.  I read its first
500 lines (header, exact arithmetic, instance helpers, layout, thread pool, the start of `Part.residual`) and not
the rest.  Kept in spirit and re-derived/re-checked: the per-weight-pair layout, the composite sparse
lift-shift-project map (`Layout.shift_map`, checked against the reference by hash), the idea of a floor-form
modular reduction with a proof, lazy deletion of pivot columns.  Changed: reduction kernel (now guarded by a
correction step, tested against Python integers), storage (device tiles), the block RREF (in place in R), the
certificate and the checkpoint/resume layer.  Removed: the thread pool (replaced by BLAS threads, the
GPU path has none), the random pre-filter, everything not exercised by a test.

## 3. Soundness

1. **Hash comparability (proof).**  The RREF of a subspace for a fixed column order is unique, so `basis_sha256`
   depends on the subspace only.  Any complete closure (any batch order, tiling, device) that reaches V_D gives the
   reference hashes; equality of both hashes is therefore a test of the whole computation.
2. **Closure invariant (proof, as in the reference).**  A *low* row has pivot weight < D; because columns are
   ordered by weight descending its support lies on weight < D columns, so all its shifts stay within weight <= D.
   Adding rows only modifies an old row by multiples of rows whose pivot is later in the column order (an old row has
   zeros before its pivot and a new row has zeros before its leading column), and such a multiple can only appear
   on a column of weight < D, hence only low rows are added to low rows.  Each low pivot is processed once, at its
   state when its batch starts, and every processed state's shifts are reduced into the basis, so the final space
   contains all shifts of all processed states; the final low rows lie in the span of the processed states (a final
   low row = the state it had when processed plus later-added low rows, which are processed in their turn).  Hence
   the final space is closed.  The certificate below does not depend on this argument.
3. **Exact arithmetic (proof + exact test).**  p < 2^24, so float32 storage of entries in [0, p) is exact.  Every
   float64 product multiplies entries in [0, p); the inner dimension accumulated between two reductions is at most
   `max_inner(p) = floor((2^50 - 4p) / (p-1)^2)` and the accumulator starts in [0, p), so every partial sum is an
   integer of magnitude <= 2^50 < 2^53 and is exact for any summation order, blocking or FMA (for p = 65537:
   262,143).  Larger inner dimensions are split with a reduction in between (`matmul_sub`, `Part.accumulate`).  The
   remaining elementwise products (scaling a row by a modular inverse, the outer product in the small RREF) have both
   factors below p, so |value| < p^2 < 2^48, and the shift coefficients (+-1, +-2) give |value| < 4p; both are reduced
   with `modp_` (|X| <= 2^50) before they are used again.
   `modp_` (floor of (X + 0.5)/p, then X - qp, then a guard correction) is exact for |X| <= 2^50; the proof is in its
   docstring.  `selftest_arith` compares it with Python integers on 3.3 million values per run (random up to
   +-2^50, multiples of p +- 1, +-p, +-2^50, six primes up to 2^24 - 3) with zero failures and zero guard hits, compares
   products whose exact result is 9.0e15 (within 2^-5 of 2^53) with int64 host arithmetic, and checks the split matmul.
   (Hardware caveat, **unverified here**: on the GPU this assumes IEEE double DGEMM (no emulated or tensor-core
   reduced-precision path); `selftest --device gpu` tests exactly that with the near-2^53 product.)
4. **Certificate (proof).**  Let V have RREF basis rows b_i.  An element f of V supported on weight < D equals
   sum f[piv_i] b_i (RREF), and f[piv_i] = 0 for every high pivot, so f is a combination of the low rows.  Hence V is
   closed under every u_{i,j} on weight < D iff every shift of every low row lies in V, i.e. reduces to zero against the
   basis.  `closed: true` is set only when a full sweep (every shift of every final low row, both projections) over a
   basis that did not change during the sweep gave zero residuals; a nonzero residual is added (`repairs`,
   `repair_rows`), the closure continues, and a new sweep starts.  Under the swap split V is sigma-stable
   (g symmetric and the generator rows m g are sigma-stable; p odd), so V = V+ (+) V-, and closure under the block-1
   shifts implies closure under the block-2 shifts u_{1,j} = sigma u_{0,j} sigma.  The randomized sweep (t random
   combinations of the low rows of each part; a non-closed space passes with probability at most (number of maps)
   p^-t, number of maps = (number of parts)^2 x (number of shifts), by a union bound over the linear maps) is only a first check.
5. **What `closed` does not cover.**  It certifies that the computed space is closed and contains the generators; its
   minimality (that it is V_D and not larger) rests on every added row being a consequence of an earlier one: each row
   is a reduction residual of a shift of a basis row, so it lies in V_D.  Rank, profile and `contains_one` are exact
   statements about V_D when `closed` is true; with `--stop-on-one`, `contains_one: true` is exact (1 is a pivot, so 1 lies
   in the space built so far, which is inside V_D) while `closure_complete` is false.

## 4. Instances (`instance` subcommand, local; exact)

`scripts/ic_bounds_big_run.py instance` draws the target with `random.Random(seed)` and
`ic_curve_sweep.random_rational_target` on `ic_dreg.PrimeCurve`, enumerates the full solution list Z
(`solutions`, copied from the reference), cross-checks the S_3 coefficients with `ic_dreg.s3_prime_coeffs`, and writes
k, p, a, b, seed, target, generator dict, Z, attempts, column counts and source hashes.

Reproduction of the reference generator (exact): `instance` reproduces `generator`, `target`, `points` and
`attempts` of the saved reference files `instance-k10-p12289-s10`, `instance-k10-p13313-s11`,
`instance-k10-p15361-s12`, `instance-k13-p40961-s12` (all four identical), and `eval_rank_below_D` (0, 98, 60) of the
first three (commands in section 6; the k = 13 evaluation rank was not recomputed).

* `unit-p65537-a1-s10-k{10,11,12,13,14}.json`: a = 1, b = 0, seed 10, `--unit-nonresidue`; first draw R = (4271, 11720),
  x(R) a quadratic non-residue, attempts = 1; same a, b, target at every k.  Z is empty at every k (enumerated, exact).
* `proper-p65537-s12-k{10,13,14}.json`: a = 4021, b = 9917, seed 12, `--want`; target (1424, 826); |Z| = 14, 1028, 4068
  (k = 10, 13, 14), full lists stored.

**Proof (unit instance: Z is empty for every k <= 15).**  E: y^2 = x^3 + x over F_p, p = 65537 = 2^16 + 1.  (0,0) is a
rational 2-torsion point and E has the form y^2 = x(x^2 + B), B = 1.  By the standard 2-descent (Silverman, AEC X.4.9)
alpha: E(F_p) -> F_p^*/F_p^*2, O -> 1, (0,0) -> B, (x,y) -> x is a group homomorphism.  The subgroup H of order 2^k of the
cyclic group F_p^* (order 2^16) lies in the squares (the subgroup of order 2^15) iff k <= 15, i.e. iff 2^(k+1) | p - 1.
Suppose S_3(x1, x2, X) = 0 with x1, x2 in H and X = x(R), R in E(F_p), x(R) a non-residue.  S_3 = 0 means X = x(P1 + P2)
or x(P1 - P2) for points P1, P2 over F_p-bar with x(P_i) = x_i.  (a) P1, P2 in E(F_p): then
alpha(R) = alpha(P1) alpha(P2)^{+-1} = x1 x2 mod squares = 1 (x1, x2 in H are squares), contradicting alpha(R) = X.
(b) neither point rational: both lie on the quadratic twist E^d(F_p) (d a non-residue), and so does P1 +- P2, so X would
be the abscissa of a point of E^d(F_p); but X is the abscissa of R in E(F_p) with X^3 + X = Y^2 a nonzero square, and an
abscissa with nonzero right-hand side is rational on exactly one of E, E^d; contradiction.  (c) one rational, one
not: Frobenius sends P1 + P2 to P1 - P2, and x(P1 + P2) in F_p forces P1 - P2 = +-(P1 + P2), i.e. 2 P1 = O or 2 P2 = O,
so that point has y = 0 and is rational on both E and E^d, which reduces to (a) or (b).  Hence Z = {}.  The enumeration
at k = 10..14 (instance files, `points: 0`) confirms it.  (The proper instance has no such obstruction: |Z| > 0.)

## 5. Validation (local, `ulimit -v 1500000`, one thread; oracle = reference engine hashes)

| case | mode | result |
|---|---|---|
| 18 closures of `small.json` (k = 6, 7; D = 3..5; 3 instances; `--gens g`) | plain cpu, fake-gpu, tiny tiles/batches (tile 5 rows, batch 3, column block 33, compact 0.02), tiny + fake-gpu | 18/18 cases x 4 modes: rank, profile, part ranks, part profiles, contains_one, both hashes equal the reference: YES; the reference also agrees with the saved dense records (rank, profile, contains_one): yes |
| `k9-p12289-s10-D5.json` | plain | rank 12616, hashes equal the saved reference; closed True; 25.1 s (reference 59.5 s); sweeps clean |
| `k9-p12289-s10-D5.json` | fake-gpu | rank 12616, hashes equal the saved reference; closed True; 28.1 s (reference 59.5 s); sweeps clean |
| `k10-p12289-s10-D5.json` (k10: 21,700 columns, rank 18,569, profile [0,0,2,138,3163,18569]) | plain | rank 18569, hashes equal the saved reference; closed True; 102.4 s (reference 277.0 s); sweeps clean |
| `k10-p12289-s10-D5.json` (k10: 21,700 columns, rank 18,569, profile [0,0,2,138,3163,18569]) | fake-gpu | rank 18569, hashes equal the saved reference; closed True; 101.8 s (reference 277.0 s); sweeps clean |
| p = 65537 vs the reference engine (generator identical; hash dtype uint32): unit-65537-k6-D5, unit-65537-k7-D5, unit-65537-k8-D5, unit-65537-k6-D6, unit-65537-k7-D6, proper-65537-k6-D5, proper-65537-k7-D5, proper-65537-k8-D5, proper-65537-k6-D6, proper-65537-k7-D6 | plain, fake-gpu, tiny | all hashes equal: YES; contains_one per case: unit-k6-D5=True, unit-k7-D5=True, unit-k8-D5=True, unit-k6-D6=True, unit-k7-D6=True, proper-k6-D5=False, proper-k7-D5=False, proper-k8-D5=False, proper-k6-D6=False, proper-k7-D6=False |
| `--nosplit` vs the reference with `split=False` (8 cases) | plain, tiny + fake-gpu | all equal: YES |
| p = 15361, k = 7, D = 5 | `--threads 1` and `--threads 2` (BLAS env set by the driver) | both equal the reference: YES |
| k8 D5 p = 15361 (device cpu), tiles 40 rows, batch 16 | SIGTERM mid-run twice then `--resume`; SIGTERM in the generator phase and in the full sweep then `--resume`; kill -9 during the second checkpoint then `--resume` | hashes equal the uninterrupted run and the reference: YES; interrupted exit codes [75, 75] (75 = SIGTERM); same number of batches as the uninterrupted run: True; kill -9 left ['.tmp-00000002', 'ck-00000001'] and resume used ck-00000001 |
| k8 D5 spill p = 15361 (device fake-gpu, --work-gb 0.2 --dev-gb 0.2142), tiles 40 rows, batch 16 | SIGTERM mid-run twice then `--resume`; SIGTERM in the generator phase and in the full sweep then `--resume`; kill -9 during the second checkpoint then `--resume` | hashes equal the uninterrupted run and the reference: YES; interrupted exit codes [75, 75] (75 = SIGTERM); same number of batches as the uninterrupted run: True; kill -9 left ['.tmp-00000002', 'ck-00000001'] and resume used ck-00000001 |
| k7 D5 p = 15361 (device fake-gpu), tiles 40 rows, batch 16 | SIGTERM mid-run twice then `--resume`; SIGTERM in the generator phase and in the full sweep then `--resume`; kill -9 during the second checkpoint then `--resume` | hashes equal the uninterrupted run and the reference: YES; interrupted exit codes [75, 75] (75 = SIGTERM); same number of batches as the uninterrupted run: True; kill -9 left ['.tmp-00000002', 'ck-00000001'] and resume used ck-00000001 |
| spill: `unit-12289-k9-D5` (`--dev-gb 0.25673 --work-gb 0.2 --tile-rows 128`, tile budget 55 MB, fake-gpu enforcing the budget) | fake-gpu with spilling | hashes equal the saved reference record: YES; closed True; spills 282, fetches 282, target switches 84, moved 416 MB to the host and 416 MB back; peak live device 72 MB under a 262 MB limit; 26.8 s |
| spill: `unit-12289-k10-D5` (`--dev-gb 0.401671 --work-gb 0.25 --tile-rows 128`, tile budget 150 MB, fake-gpu enforcing the budget) | fake-gpu with spilling | hashes equal the saved reference record: YES; closed True; spills 787, fetches 787, target switches 78, moved 1927 MB to the host and 1927 MB back; peak live device 179 MB under a 411 MB limit; 128.1 s |
| spill: `p15361-k7-D6` (`--dev-gb 0.213482 --work-gb 0.2 --tile-rows 64`, tile budget 13 MB, fake-gpu enforcing the budget) | fake-gpu with spilling | hashes equal the reference engine: YES; closed True; spills 224, fetches 224, target switches 60, moved 84 MB to the host and 84 MB back; peak live device 33 MB under a 218 MB limit; 10.1 s |
| spill: `proper-65537-k8-D5` (`--dev-gb 0.214494 --work-gb 0.2 --tile-rows 64`, tile budget 14 MB, fake-gpu enforcing the budget) | fake-gpu with spilling | hashes equal the reference engine: YES; closed True; spills 366, fetches 366, target switches 52, moved 143 MB to the host and 143 MB back; peak live device 32 MB under a 219 MB limit; 10.8 s |
| spill: `unit-65537-k8-D5` (`--dev-gb 0.214494 --work-gb 0.2 --tile-rows 64`, tile budget 14 MB, fake-gpu enforcing the budget) | fake-gpu with spilling | hashes equal the reference engine: YES; closed True; spills 366, fetches 366, target switches 52, moved 143 MB to the host and 143 MB back; peak live device 32 MB under a 219 MB limit; 10.9 s |
| `--stop-on-one`, unit p = 12289, (k, D) = (6,5), (7,5), (6,6) (1 is in V_D) | plain cpu | exit 2, `contains_one: true`, `closure_complete: false`, `closed: false`; the full run equals the reference: YES; rank at stop 1584/1586, 3407/3473, 2510/2510 |
| `selftest --device fake-gpu` (8 built-in cases, cpu vs device, plus embedded reference rank and hash prefixes; arithmetic tests) | fake-gpu | passed: True |

## 6. Commands

All local commands: `cd` to the repository root, `ulimit -v 1500000`, `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

```
# instances (reproduction check against the saved reference files, then the new ones)
python3 scripts/ic_bounds_big_run.py instance --k 10 --prime 12289 --a 1 --b 0 --seed 10 --D 5 --eval-rank --out W/i1.json   # = instance-k10-p12289-s10.json
python3 scripts/ic_bounds_big_run.py instance --k 10 --prime 13313 --a 1 --b 0 --seed 11 --want --D 5 --eval-rank --out W/i2.json
python3 scripts/ic_bounds_big_run.py instance --k 10 --prime 15361 --a 4021 --b 9917 --seed 12 --want --D 5 --eval-rank --out W/i3.json
python3 scripts/ic_bounds_big_run.py instance --k 13 --prime 40961 --a 4021 --b 9917 --seed 12 --want --D 5 --out W/i4.json
for k in 10 11 12 13 14; do python3 scripts/ic_bounds_big_run.py instance --k $k --prime 65537 --a 1 --b 0 --seed 10 --unit-nonresidue --D 6 --eval-rank --out experiments/ic-bounds-k14/instances/unit-p65537-a1-s10-k$k.json; done
python3 scripts/ic_bounds_big_run.py instance --k 10 --prime 65537 --a 4021 --b 9917 --seed 12 --want --D 6 --eval-rank --out experiments/ic-bounds-k14/instances/proper-p65537-s12-k10.json
for k in 13 14; do python3 scripts/ic_bounds_big_run.py instance --k $k --prime 65537 --a 4021 --b 9917 --seed 12 --want --out experiments/ic-bounds-k14/instances/proper-p65537-s12-k$k.json; done

# validation (oracle: the reference engine; W = scratch directory)
python3 scripts/ic_bounds_big_validate.py small18  --work W
python3 scripts/ic_bounds_big_validate.py refs     --work W --modes plain,fake-gpu
python3 scripts/ic_bounds_big_validate.py p65537   --work W
python3 scripts/ic_bounds_big_validate.py nosplit  --work W
python3 scripts/ic_bounds_big_validate.py threads  --work W
python3 scripts/ic_bounds_big_validate.py resume   --work W --k 8 --D 5 --delays 2,3
python3 scripts/ic_bounds_big_validate.py resume   --work W --k 7 --D 5 --delays 1.5,1.5 --device fake-gpu
python3 scripts/ic_bounds_big_validate.py negative --work W --k 7 --D 5 --gens g
python3 scripts/ic_bounds_big_validate.py spill    --work W
python3 scripts/ic_bounds_big_validate.py resume   --work W --k 8 --D 5 --delays 2,3 --device fake-gpu --extra "--work-gb 0.2 --dev-gb 0.2142" --tag-suffix -spill
python3 scripts/ic_bounds_big_validate.py negative --work W --k 7 --D 5 --gens g --device fake-gpu --extra "--work-gb 0.2 --dev-gb 0.20465" --tag-suffix -spill
python3 scripts/ic_bounds_big_validate.py scan     --work W --k 6 --D 4 --gens g --n 40
python3 scripts/ic_bounds_big_validate.py stopone  --work W
python3 scripts/ic_bounds_big_validate.py rates    --work W --cases 5:6,5:7,5:8,5:9,6:6,6:7,6:8,6:9,5:10,5:11 --tag d5d6
python3 scripts/ic_bounds_big_run.py selftest --device fake-gpu --out experiments/ic-bounds-k14/validation/selftest-fake-gpu.json
python3 experiments/ic-bounds-k14/estimate.py
```

Single runs: `python3 scripts/ic_bounds_big_run.py closure --instance FILE --D 5 --device cpu|fake-gpu|gpu [--tile-rows 5 --batch 3 --col-block 33 --compact-frac 0.02] --out REC.json`.
Test-only hooks: `--withhold part:index` (skip the shifts of one low row), environment `ICBIG_TEST_CKPT_SLEEP` / `ICBIG_TEST_CKPT_SLEEP_FROM`
(slow a checkpoint down so that it can be killed in the middle).  Exit codes of `closure`: 0 closed, 2 stopped on 1 in V_D, 3 sweep failed with `--no-repair`, 75 SIGTERM after a checkpoint.

## 7. Pre-filter, negative controls

**Pre-filter.**  The inherited `--compress` was never exercised and was not validated by the required negative control (margin 0, tiny T), so it was removed.  Nothing here can change which space is computed except the closure and the sweeps.

**Negative control (a): withheld shifts** (`--withhold part:index`, test hook; p = 15361, k = 7, D = 5, `--gens g`, `--no-repair`; `experiments/ic-bounds-k14/validation/negative-withheld-k7-D5-gensg.json`).  The shifts of one pending low row are skipped (the row is marked processed).  15 runs; 2 changed the final basis (withholding the first low row leaves rank 1 instead of 3469): all 2 failed the certificate (exit 3, `closed: false`; once with the randomized sweep, once with `--t-random 0` so that the FULL sweep alone caught it), and with repair enabled the run ended at the correct hashes with `closed: true` and `repairs` = 1.  The other withheld rows did not change the result (their shifts lie in the span of the others) and gave clean sweeps and correct hashes, as they must.

**Scan of every withheld row** (k = 6, D = 4, `--gens g`, parts 0 and 1, positions 0..39; randomized sweep off, no repair): 80 runs, 1 changed the final basis, 1 of those were caught by the full sweep, 0 undetected, 0 false alarms, 79 harmless omissions with a clean sweep.

The same 15 runs with spilling (fake-gpu, `--work-gb 0.2 --dev-gb 0.20465`, so the sweep runs target-major with tiles moving between host and device): 2 changed the final basis, all caught (exit 3) and repaired to the correct hashes: True.

**Negative control (b)** does not apply (pre-filter removed).  Not done: a control in which the sweep is run on a corrupted *stored* matrix entry (the sweep would fail by construction; `check_rref` checks the structure of sampled rows).

## 8. Rates

One thread (`OPENBLAS_NUM_THREADS=1`), engine defaults with `--work-gb 0.4 --batch 256`, unit instance p = 12289, `experiments/ic-bounds-k14/validation/rates-*.json`.  "closure madds" counts the closure phase (reduce + update + block RREF, multiply-adds actually done); "sparse / dense" divides it by the dense method (every image reduced against all pivots: n x rank x free, plus the same update and RREF); rates are madds divided by the seconds spent in that kernel (the reduce rate includes the certificate sweeps).  Rows with "(full)" are the regime 1 in V_D (all columns pivots), which is not the regime of the target run.

| D | k | columns | rank | closure madds | sparse / dense | reduce G madds/s | update G madds/s | RREF G madds/s | seconds | max RSS MB |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 6 | 1586 | 1586 (full) | 3.079e+08 | 0.521 | 2.69 | 1.96 | 0.62 | 0.3 | 49 |
| 5 | 7 | 3473 | 3473 (full) | 2.220e+09 | 0.410 | 5.57 | 2.27 | 0.80 | 1.2 | 77 |
| 5 | 8 | 6885 | 6885 (full) | 1.290e+10 | 0.342 | 6.94 | 2.53 | 0.96 | 4.9 | 177 |
| 5 | 9 | 12616 | 12616 (full) | 6.347e+10 | 0.307 | 7.82 | 2.13 | 1.04 | 23.3 | 491 |
| 6 | 6 | 2510 | 2510 (full) | 2.147e+09 | 0.630 | 6.44 | 3.00 | 0.99 | 1.0 | 62 |
| 6 | 7 | 6476 | 6476 (full) | 2.351e+10 | 0.452 | 10.53 | 3.46 | 2.03 | 5.3 | 174 |
| 6 | 8 | 14893 | 14893 (full) | 1.757e+11 | 0.336 | 12.42 | 4.01 | 2.44 | 33.0 | 657 |
| 6 | 9 | 31180 | 31180 (full) | 1.022e+12 | 0.262 | 13.57 | 4.53 | 1.79 | 177.3 | 1018 |
| 5 | 10 | 21700 | 18569 | 3.196e+11 | 0.247 | 10.67 | 1.85 | 0.99 | 113.3 | 804 |
| 5 | 11 | 35443 | 25033 | 1.018e+12 | 0.198 | 10.52 | 1.89 | 0.91 | 333.1 | 1084 |

Largest local case by time: D = 5, k = 11 (35443 columns, 333.1 s).  Single-thread rates: reduction about 10 G madds/s (a large float64 GEMM per call), update about 1.7 G madds/s (rank-t updates with small t are memory bound: gather, cast, tiny GEMM, reduce, scatter), block RREF about 1 G madds/s (latency bound).  The reference engine took 277 s in total for k = 10, D = 5 (229 s closure; counted madds 1.575e12, all of them dense); this engine takes 113.3 s in the rates run (94.8 s in the validation run), closure and full sweep included, on a shared machine.

## 9. Estimate (ESTIMATE, not a measurement)

**Estimate**, from `estimate.py`: closure-phase multiply-adds (one engine, counted as in section 8), power law in the number of columns with the local slope between consecutive depths; the band uses the smallest and largest of the last three slopes.  Assumptions: the local runs (k <= 9 for D = 6, k <= 11 for D = 5) are in the regime 1 in V_D (full rank) for k <= 9 and not for k >= 10, the target run is in the non-full regime and has a different p, so the extrapolation from k <= 9 to k = 14 is more than a factor 30 in the number of columns and **can be off by a large factor**; backtest on the measured D = 5 points from k <= 9 data only: k = 10: measured/mid = 1.25; k = 11: measured/mid = 1.12.  The certificate sweep multiplies the count by 1.71 (measured on D = 5, k = 10 and 11, used for the "with sweep" column).  The hours column is a **lower bound**: it divides the madds by a sustained 6 to 12 T exact-FP64 madds/s on a B200 (an **assumption**: 30 to 60 % of a nominal 20 T madds/s) and ignores the memory-bound update, the latency-bound block RREF and the host-side bookkeeping, which take most of the single-thread time (in the D = 5, k = 10 run the reduction GEMMs are 36 % of the kernel time).  The single-thread column divides by the measured blended rate of the non-full runs (total madds over total seconds, all kernels, closure and sweep: 5.3 G madds/s).  The count itself is very likely an underestimate (the backtest underestimates by 12 to 25 %, and a rough hand count for the target, images x selected pivots x free columns, gives a few times 10^15 for the closure phase; this is not a measurement).

| D | k | columns | closure madds (low) | (mid) | (high) | with sweep (mid) | one thread at the measured blended rate, days | GEMM-bound GPU hours at 12 .. 6 T madds/s |
|---|---|---|---|---|---|---|---|---|
| 5 | 12 | 55455 | 2.93e+12 | 3.35e+12 | 3.86e+12 | 5.73e+12 | 0.01 | 0.000 .. 0.000 |
| 5 | 13 | 83682 | 7.74e+12 | 9.99e+12 | 1.32e+13 | 1.71e+13 | 0.04 | 0.000 .. 0.001 |
| 5 | 14 | 122438 | 1.90e+13 | 2.75e+13 | 4.10e+13 | 4.71e+13 | 0.10 | 0.001 .. 0.003 |
| 6 | 10 | 60460 | 4.95e+12 | 5.15e+12 | 5.44e+12 | 8.82e+12 | 0.02 | 0.000 .. 0.000 |
| 6 | 11 | 110056 | 2.06e+13 | 2.22e+13 | 2.47e+13 | 3.81e+13 | 0.08 | 0.001 .. 0.002 |
| 6 | 12 | 190051 | 7.59e+13 | 8.43e+13 | 9.81e+13 | 1.44e+14 | 0.32 | 0.003 .. 0.008 |
| 6 | 13 | 313912 | 2.51e+14 | 2.87e+14 | 3.48e+14 | 4.92e+14 | 1.08 | 0.010 .. 0.028 |
| 6 | 14 | 499178 | 7.58e+14 | 8.90e+14 | 1.12e+15 | 1.53e+15 | 3.35 | 0.030 .. 0.089 |


## 10. Memory plan (`plan`, gpu)

D = 6, `plan --device gpu` (matrices: max over the rank of r x (free + lazily deleted columns) x 4 bytes per part, plus one tile per part; plus the work budget, one tile for compaction, 1.5 GiB allowance for context and cuBLAS):

| k | work-gb | tile-mb | batch | matrices worst case GB | predicted peak GB |
|---|---|---|---|---|---|
| 10 | 8 | 256 | 1024 | 2.3 | 12.0 |
| 11 | 12 | 512 | 1024 | 6.9 | 20.9 |
| 12 | 20 | 1024 | 2048 | 19.7 | 42.2 |
| 13 | 25 | 2048 | 2048 | 52.2 | 80.7 |
| 14 | 30 | 2048 | 2048 | 125.8 | 159.3 |

At k = 14: columns 499,178, parts 249,824 and 249,354, low columns 61,272 and 61,166, matrix worst case 65.5 GB and 65.3 GB per part (at rank about 125,000 of each part), 125.8 GB with the tile slack; predicted peak **159.3 GB** (D = 5: 44.8 GB).  The fake-gpu run measures the peak of live arrays: k = 9, D = 5, `--work-gb 0.5`: 125 MB against 100 MB of matrices plus a 500 MB budget (the budget was not binding on a case this small).  Fragmentation of the cupy pool is not modelled (cupy retries after freeing cached blocks on an allocation failure).

**Spilling (`--dev-gb`)**, `plan --dev-gb G` reports the tile budget, the worst case of one part (the target), the predicted device peak and the predicted host memory for the spilled tiles (both parts at their worst rank at the same time):

| k, D | `--dev-gb` | `--work-gb` | `--tile-mb` | tile budget GB | one part (target) worst GB | predicted device peak GB | predicted host spill (worst) GB |
|---|---|---|---|---|---|---|---|
| 14, 6 (H200, 141 GB card, about 139 GB usable) | 137 | 28 | 2048 | 105.5 | 63.0 | 137 (within budget by construction) | 20.3 |
| 14, 6 (H200) | 139 | 28 | 2048 | 107.5 | 63.0 | 139 | 18.3 |
| 11, 6 (forced-spill check on a real card) | 18.5 | 12 | 512 | 4.5 | 3.5 | 18.5 | 2.4 |
| 14, 5 | 139 | 20 | 1024 | 116.5 | 4.7 | 31.8 (no spill) | 0 |

At k = 14, D = 6 on the H200 the host holds at most about 20 GB of tiles, and a switch of the target part moves at most that amount in each direction (about 1 to 10 s over PCIe from pageable memory; pinned memory is not used).  With target-major order there are two switches per batch and two per sweep.  Host RAM (377 GB) is far above this.

## 11. Remote command lines

Remote layout: code in `/workspace/code` (`ic_bounds_big.py` and `ic_bounds_big_run.py` side by side), instances in `/workspace/in`, outputs in `/workspace/out`,
checkpoints in `/workspace/ckpt`.  Run `selftest` first; do not start a closure unless it exits 0.

```
cd /workspace/code
python3 ic_bounds_big_run.py selftest --device gpu --out /workspace/out/selftest-gpu.json
python3 ic_bounds_big_run.py plan --k 14 --D 6 --device gpu --work-gb 30 --tile-mb 2048 --batch 2048

# D = 6, unit instance (restart any of them with the same command: --resume continues the newest complete checkpoint)
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k10.json --D 6 --device gpu --work-gb 8  --tile-mb 256  --batch 1024 --ckpt /workspace/ckpt/k10-D6 --ckpt-minutes 10 --resume --status /workspace/out/k10-D6.status.json --out /workspace/out/k10-D6.json
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k11.json --D 6 --device gpu --work-gb 12 --tile-mb 512  --batch 1024 --ckpt /workspace/ckpt/k11-D6 --ckpt-minutes 15 --resume --status /workspace/out/k11-D6.status.json --out /workspace/out/k11-D6.json
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k12.json --D 6 --device gpu --work-gb 20 --tile-mb 1024 --batch 2048 --ckpt /workspace/ckpt/k12-D6 --ckpt-minutes 20 --resume --status /workspace/out/k12-D6.status.json --out /workspace/out/k12-D6.json
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k13.json --D 6 --device gpu --work-gb 25 --tile-mb 2048 --batch 2048 --ckpt /workspace/ckpt/k13-D6 --ckpt-minutes 30 --resume --status /workspace/out/k13-D6.status.json --out /workspace/out/k13-D6.json
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k14.json --D 6 --device gpu --work-gb 30 --tile-mb 2048 --batch 2048 --ckpt /workspace/ckpt/k14-D6 --ckpt-minutes 30 --resume --stop-on-one --status /workspace/out/k14-D6.status.json --out /workspace/out/k14-D6.json

# D = 5, k = 14
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k14.json --D 5 --device gpu --work-gb 20 --tile-mb 1024 --batch 2048 --ckpt /workspace/ckpt/k14-D5 --ckpt-minutes 20 --resume --stop-on-one --status /workspace/out/k14-D5.status.json --out /workspace/out/k14-D5.json
```

**H200 (141 GB, about 139 GB usable, host 377 GB): k = 14, D = 6 with spilling** (predicted device peak 137 GB, host about 20 GB):

```
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k14.json --D 6 --device gpu --dev-gb 137 --work-gb 28 --tile-mb 2048 --batch 2048 --ckpt /workspace/ckpt/k14-D6 --ckpt-minutes 30 --resume --stop-on-one --status /workspace/out/k14-D6.status.json --out /workspace/out/k14-D6.json
```

Quick forced-spill check at k = 11, D = 6 on a real card (tile budget 4.5 GB against 6.9 GB of matrices, so about 2.4 GB spills; the basis hashes in the record, `basis_sha256`, must equal those of the no-spill k = 11 run; the log lines show a nonzero `spill`):

```
python3 ic_bounds_big_run.py closure --instance /workspace/in/unit-p65537-a1-s10-k11.json --D 6 --device gpu --dev-gb 18.5 --work-gb 12 --tile-mb 512 --batch 1024 --status /workspace/out/k11-D6-spill.status.json --out /workspace/out/k11-D6-spill.json
python3 ic_bounds_big_run.py selftest --device gpu --out /workspace/out/selftest-gpu.json     # includes two forced-spill cases (cpu vs device)
```

Settings: `--work-gb` is the working-buffer budget on the device (R, block RREF, gathers; see section 2); the matrices come on top (section 10), so
(without `--dev-gb`) `work-gb` + matrices + about 3 GB must stay under the device memory (k = 14: 30 + 125.8 + 2 + 1.5 = 159.3 GB of 180 GB; if the pool fragments, lower `--work-gb` to 24).
`--tile-mb` is the tile of B (rows per tile = tile bytes / (4 x part columns)); `--batch` caps the rows per batch (the budget lowers it).  Checkpoints are written
with one slot: the new slot is complete before the old one is deleted, so `/workspace/ckpt` needs room for two copies of B (about 2 x 125 GB at the worst rank for
k = 14, far less for the smaller k, and less than the matrix size once the rank passes half the columns); `--ckpt-slots 2` keeps two.  `--stop-on-one` is recommended for
k = 14 (stop when the weight-0 column becomes a pivot of the symmetric part; the record then has `contains_one: true`, `closure_complete: false`); leave it off to
complete the closure anyway.  SIGTERM (spot preemption) finishes the current batch, checkpoints and exits 75; start the same command again.  A resumed run follows the
same sequence of batches and ends at the same hashes as an uninterrupted one (validated).  `closed: true` and `sweep_full.complete: true` in the record are the
certificate; a run with `closed: false` is not a bound.

## 12. Cupy calls that could not be checked here

Everything with `cupy` was written against the documented behaviour and exercised only through the `fake-gpu` stand-in (numpy semantics, host/device mixing forbidden).
Not checked on a real device:

* `cupy.asarray` / `cupy.asnumpy` for int64, float64, float32 and empty arrays; `cupy.zeros`, `cupy.array(view)` (copy), `cupy.outer`, `cupy.floor(x, out=x)`, `cupy.any(x, axis=0/1)`, `cupy.argmax` on an int8 array (`_rref_rows`), `x.astype`, `x.min()`, `x.max()` on views, `float(0-d array)`.
* Indexing: `tile[idx_dev, a:b]` (array index mixed with a slice, get and set), `tile[:nr][:, idx_dev]`, `A[:, idx_dev]`, `R[:, idx_dev] = block`, `Yt[:, pos_dev] += block` (augmented assignment with an index array, unique indices), `S[arange_dev, idx_dev] = 1.0`, `Sphys[idx_dev] = block`, `Mf[:, idx_dev] = M`, `R[t:t+tn] = Ln`, in-place `-=`, `+=`, `*=` on non-contiguous views.
* `@` on float64 2-D arrays including non-contiguous operands (cupy copies them to contiguous; the copies are small compared with the budgeted blocks but are not counted in it) via cuBLAS DGEMM, and that its result is exact for integers below 2^53 (`selftest` part (b) tests this; a cuBLAS float64-emulation mode would show up there).
* Runtime and bookkeeping: `cupy.cuda.runtime.memGetInfo`, `runtimeGetVersion`, `driverGetVersion`, `getDeviceProperties` (the `name` field may be bytes or str, both handled), `cupy.cuda.cublas.getVersion(cupy.cuda.device.get_cublas_handle())` (wrapped in try/except in `selftest`), `cupy.get_default_memory_pool().free_all_blocks()` and `.total_bytes()`, `cupy.get_default_pinned_memory_pool().free_all_blocks()`, `cupy.cuda.Stream.null.synchronize()`.
* Verified by the coordinator on an H200 (CuPy 14.2.0, cuBLAS 12.9) for the engine version before spilling: `selftest --device gpu` passed (exact arithmetic, 0 guard hits, 8 cases equal to cpu and reference) and k = 10, D = 6 closed with a clean sweep, so the calls above worked as written.  Not yet checked on a real device: the spilling path (`cupy.asnumpy` of 0.5 to 2 GB tiles, `cupy.asarray` back, pool freeing after each spill) and the target-major loop.
* Behaviour under memory pressure (pool fragmentation, `OutOfMemoryError` after the automatic free-and-retry) is not modelled.

## 13. Limits and unverified items

* The version of this engine before spilling ran on an H200 (selftest passed; k = 10, D = 6 closed, see section 12); the version in this directory (spilling, target-major order) has been validated only through `fake-gpu` (same code path, numpy semantics, budget enforced).  Real-device speed of the spill path, host-to-device transfer rates and cuBLAS behaviour at k >= 11 are unverified.
* **Largest local cases**: D = 5, k = 11 (35,443 columns) and D = 6, k = 9 (31,180 columns, full rank, 1 in V_D); D = 5, k = 10 (rank 18,569, 1 not in V_D) matches the reference hashes.
  Nothing at the size of the target (499,178 columns) has been run; the rank, profile, memory and time at k = 12..14 are **unverified**.  The estimate in section 9 is an extrapolation and is labelled as such.
* The local instances with small k are in the regime 1 in V_D (full rank): they exercise the closure and sweep code, but not the late phase of the target run (pivots in weight-6 columns with free columns left).  The k = 10 case (rank 18,569 of 21,700) and the proper instances (1 not in V_D) do.
* Update and block-RREF kernels are memory- and latency-bound in the single-thread measurements (about 1.7 and 1 G madds/s against 10 G madds/s for the reduction); on a GPU the block RREF is a Python loop of about ten device calls per row
  (batch of n rows: about 10 n calls per image), which may dominate; moving the 32-row echelon to the host or batching it is the first thing to try if the measured rate is far below the DGEMM rate.
* Hash comparability covers the two hashes (both parts); `rref_ok` samples 64 rows per part; `Part.sha()` compacts the physical columns first.
* The certificate (`closed`) is only as good as the exactness of the arithmetic (section 3.3) on the device that ran it.
* `--threads` sets BLAS threads only; elementwise passes are single threaded on the CPU path.  threadpoolctl was not installed here, so the in-process limit was not exercised (the driver sets the environment instead).
* `source_sha256` in the instance files is the hash of the engine and driver at the time they were written; the instance content does not depend on them (generator code `s3_coeffs`, `chain_generator`, `solutions` is copied unchanged from the reference).
* The validation records carry two harness hashes (`7f5dbbbb...` and `23eb1271...`): the harness gained the `scan` command (and its `--n` flag) while the first runs were in progress; the engine and driver hashes are identical in every record.
* Hypotheses of the mathematics (that the reference's chain model is the right object for the bound, the conjectured final rank 450,000) are not examined here.
