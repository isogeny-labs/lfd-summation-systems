# CERTCHECK: an independent checker for the "closed containing space" certificate

Evidence labels (AGENT_RULES): **exact** = certified computation with the command and file named; **proof** =
argument written out here; **empirical** = pattern in data. Nothing below is **kernel**.

Files (all new):

| file | role |
|---|---|
| `scripts/ic_bounds_cert_check.py` | the checker, one self-contained file (numpy; cupy only for `--device gpu`) |
| `scripts/ic_bounds_cert_check_test.py` | test harness; MAY import the engines to produce inputs |
| `experiments/ic-bounds-k14/validation/certcheck-*.json` | validation records (arith, spec, positive, reference, negative, devices) |

## 1. What is claimed and what is checked

Setting: F_p[(Z/2^k)^2], monomials (a1, a2), weight popcount(a1) + popcount(a2), shifts u_{i,j} = multiplication by
x_i^(2^j) (add 2^j to a_i mod 2^k), g = NF(S_3(x1, x2, X)). V_D is the least subspace U of span(M_<=D) with g in U and
u_{i,j} f in U for every f in U cap span(M_<D) and every i, j.

**Lemma (proof).** Let W be a subspace of span(M_<=D) with (a) g in W, (b) u_{i,j} f in W for all f in W cap span(M_<D)
and all i in {1, 2}, j < k. Then V_D <= W. In particular, if 1 is not in W, then 1 is not in V_D.
*Proof.* V_D is the intersection of all subspaces satisfying the two defining conditions (the family is closed under
intersection and contains span(M_<=D)); W is one of them. QED. (A shift of a monomial of weight < D has weight <= D:
adding 2^j to a with bit j clear gives weight + 1, with bit j set the carry only lowers the weight, and the wrap
modulo 2^k drops a bit; so every shifted vector lies in the column space.)

The checker establishes (a), (b) and (c) 1 not in W from the stored basis of W (the engine's checkpoint at the end of
the closure phase) and reports `certificate_valid` = format ok and (a) and (b); `one_excluded` = (c);
`all_true` = both. Equal to the three facts of the request.

### 1.1 Description of W
Part "+" has columns {a1 <= a2}, part "-" has {a1 < a2} (orbit coordinates, inherited column order: weight descending,
then key a1 | a2 << k). lift(Y)[h] = Y[pos h] for a representative h and +-Y[pos sigma h] otherwise (0 on the
antisymmetric diagonal); project(H)[m] = H[rep m] +- H[sigma rep m] (2 H on the symmetric diagonal).
Each part is stored as a basis in the form row i = e_{piv[i]} + B[i] on the free columns (free = complement of the
pivots, increasing). W = lift(W+) (+) lift(W-). Because lift/project are linear and project o lift is the identity up
to the unit 2 on each part, **v in W iff project_+(v) in W+ and project_-(v) in W-**, and for a part
**y in W_part iff y[free] - y[piv] @ B = 0 (mod p)**. (The row space of [e_piv | B] is exactly the graph
{y : y[free] = y[piv] B}; no echelon property is used for this.)

### 1.2 Check 4 (1 not in W)
1 is the monomial (0, 0), the weight-0 column of the "+" part (the last column). project_+(1) = 2 e_{(0,0)} and
project_-(1) = 0. If (0, 0) is a free column of "+", then every y in W+ with y[piv] = 0 has y[free] = 0, but 2 e_{(0,0)}
has zero pivot coordinates and a nonzero free coordinate, so it is not in W+: 1 is not in W. This uses only that the
column is free; it does not depend on the echelon order. The checker does both: reports `weight0_column_free` and
runs the membership test of project(1) (`one.in_W_per_part`).

### 1.3 W cap span(M_<D) is spanned by the rows whose pivot has weight < D (proof, tested)
Let the "low" rows be those with pivot weight < D. (>=) A low row is e_piv + B[i]; it lies in M_<D iff B[i] vanishes on
the free columns of weight >= D (its pivot column has weight < D, and the other pivot columns carry 0 by the form).
The checker verifies exactly this for every low row (`low_rows_have_no_support_on_weight_ge_D`). With the column order
weight-descending this is implied by the leading-column property (all entries of row i left of its pivot are 0, hence
all entries lie on columns of weight <= weight(pivot)), which the checker also tests
(`rref_leading_column_ok`, reported separately: it is needed for the profile = dim(W cap M_<=j) statement and for
comparing hashes with the engine, not for the proof). (<=) Any v in W is sum_i v[piv_i] row_i. If v is supported on
weight < D then v[piv_i] = 0 for every row whose pivot has weight >= D, so v is in the span of the low rows.
Dense test (`certcheck-spec.json`): k = 6, D = 4, reference basis, dim(W cap M_<D) computed by dense rank
(rank W - rank of W on the weight >= D columns) equals the number of low rows.

For the split: M_<D is sigma-stable, so v in W cap M_<D implies (v +- sigma v)/2 in W cap M_<D, hence
W cap M_<D = lift(low rows of W+) (+) lift(low rows of W-); by linearity it suffices to test the lifted low rows of both
parts. Checked count: **rows** = number of low rows of both parts (73,587 at the target).

### 1.4 Shifts: all 2k directly (default) or block 1 only
sigma W = W (W = lift W+ (+) lift W-, and lift of "+" is sigma-fixed, of "-" is negated). For f in W cap M_<D:
sigma f is in W cap M_<D and u_{2,j} f = sigma u_{1,j} sigma f. So closure under the block-1 shifts implies closure
under the block-2 shifts. **Choice: the default `--shifts all` tests all 2k shifts u_{1,j}, u_{2,j} directly in full
coordinates** (lift the low row to the full column space, shift, project to both parts, test membership). Reason:
it makes the conclusion independent of the sigma argument and of the engine's choice to apply only block 1; it costs
twice the block-1 time (section 5). `--shifts block1` is available and, on every case tested here, gives the same
verdict (validation records contain both). Triples checked = rows x shifts x 2 target parts
(4,120,872 for `all`, 2,060,436 for `block1` at the target).

## 2. How the membership of every shifted low row is tested (exact arithmetic)

For a source part pi, target part q and shift (i, j), the map "lift, shift, project to q, reduce modulo W_q" is linear
in the source coordinates s (the low columns of pi). The checker builds for each (pi, q, shift) the matrix
T[s, f] over the free columns f of q:
* a low source column s lifts to the full columns rep(s) (coefficient 1) and sigma rep(s) (coefficient sign of pi);
  each is shifted (a translation of the monomial, an injective map), and projected onto the column m of q with
  coefficient c in {1, -1, 2, 0};
* if m is a free column of q: T[s, f(m)] += c; if m is the pivot of row r: T[s, :] -= c B[r, :].
Then the residual of a low row y (vector on the low columns) is y @ T, and the triple (row, shift, q) passes iff the
residual is 0 modulo p. This is the same statement as `v[free] - v[piv] @ B = 0` for the shifted-and-projected v,
evaluated for all low rows of pi at once as one dense float64 matrix product (rows x Lp times Lp x F).

Exactness (proof). All entries are integers in [0, p) stored in float64 (the float32 tiles are exact for p < 2^24 and
are checked to hold integers in [0, p)). A float64 product of two such integers is exact; every partial sum of a
product block is a non-negative integer, so it is exact as long as it is <= 2^53 whatever the summation order, BLAS
blocking or FMA contraction (all integers up to 2^53 are representable). `max_inner(p)` is the largest m with
m (p-1)^2 + (p-1) <= 2^53 - 2p; a block with inner dimension <= m plus a carried residue < p is <= 2^53 - 2p; longer
inner dimensions are cut into pieces of m with a reduction between pieces. For p = 65537, m = 2,096,639 against an
inner dimension of at most 61,272 here and 224,688 in the membership test of g. `modp(X)` for |X| <= 2^53 - 2p:
q = floor(fl(X/p)): the quotient is below 2^53/p and its rounding error is at most 2^-53 of it, i.e. below 1/p < 1, so q is within 1 of floor(X/p);
q p is an integer of magnitude <= |X| + p <= 2^53, so q p and
X - q p are exact, the result lies in [-p, 2p) and two selections bring it into [0, p). Tested against Python
integers at the bound (`certcheck-arith.json`: p = 12289, 65537, 114689, 16777213; values within 40 of +-(2^53 - 2p),
every residue boundary near q_max p, 20,000 random values; products with all entries p - 1 at inner dimension
max_inner(p) + 7 for p = 16777213, i.e. across the chunk boundary, and with a forced inner dimension of 5 and 7).

## 3. Independence statement

The checker file imports no repository module (it imports argparse, hashlib, io, json, os, resource, sys, tempfile,
time, numpy; cupy only for `--device gpu`).

READ for the specification, as permitted: the module docstring and the `Layout` class of `csd/ic_bounds_d5_verify.py`,
its `s3_coeffs`, `chain_generator`, `Part.sha`, `Part` (to know the form of the basis), and in
`csd/ic_bounds_big.py` the checkpoint code (`save_checkpoint`, `_state_json`, `_config_json`, `find_checkpoint`,
`load_checkpoint`).

**READ BEYOND the permitted scope (disclosure).** While locating the checkpoint code in `csd/ic_bounds_big.py` I also
read, in that file: the module docstring (lines 1 to 75), the `Layout` class (lines 425 to 540, including
`shift_map` and `lift_full`), the head of `class Part` (lines 600 to 660, `tile_rows_filled`) and `Part.compact` and
`Part.source` (lines 828 to 880), the first lines of `Engine.run_closure` and `Engine.run` (1255 to 1290, 1395 to
1445) and `Engine.record` (1445 to 1475). I read `scripts/ic_bounds_big_run.py` (the driver) in full and the
`experiments/ic-bounds-k14/runs/*.json` records. I did not read the sweep code, the residual or add code, or the
shift application in the big engine beyond `Layout.shift_map`. Consequence: the layer decomposition used here for
"lift, shift, project" (a source column contributes through its representative and its sigma partner) has the same
shape as the engine's `shift_map`; it is forced by the mathematics (the lift has two preimage columns), but I cannot
claim I would have arrived at it without having seen it. The conclusion does not depend on it, for two reasons that
are tested: (i) the layout, shifts and projection are compared with the reference engine's `Layout` (equal, `spec`
group), which the request allowed me to read, and (ii) `certcheck-spec.json` contains a dense full-coordinate closure
test (no parts, no projection, no layers; rows reduced against the full-coordinate RREF of lift(W+) + lift(W-)) that
gives the same verdicts as the checker on the untouched and on a corrupted basis, and the built-in self-test
(`--selftest`) runs the whole pipeline with `--shifts all` on bases produced by a naive full-coordinate closure
(`naive_basis`) and on an unsplit single part. Nothing was copied from the engines: the checker's closure test is a
product-matrix formulation (one matrix per (source part, target part, shift)) that does not occur in them.

Shared with the engines by necessity (the specification): column order, orbit coordinates and lift/project
conventions, the on-disk checkpoint format, and the definition of `Part.sha()`. A bug in the *specification* (for
example a wrong normal form of S_3) would affect both; it is covered by: S_3 expansion checked against the reference
`s3_coeffs` and numerically against the addition law (S_3(x1, x2, x(P1 +- P2)) = 0 at 120 random point pairs per curve,
three curves); the instance's `generator` dict is recomputed from (a, b, target[0]) and the verdict records whether they
agree.

## 4. Validation (all runs under `ulimit -v 1500000`, one thread; records in `validation/certcheck-*.json`)

### 4.1 Positive: bases written by `scripts/ic_bounds_big_run.py closure --ckpt` at the end of the closure phase

Checker `--shifts all` (every u_{i,j}), compared with the engine record (`--engine-record`): rank, part ranks, profile, basis hashes, contains_one. `rows / triples` = low rows / (row, shift, target) triples checked; fail = triples leaving W. Same verdicts with `--shifts block1` (half the triples).

| case | engine rank | profile (pivots of weight <= j) | closed | rows / triples (all) | fail | g in W | m g in W | 1 excluded | hashes = engine | outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| unit-p12289-k7-D4 | 829 | [0, 0, 2, 115, 829] | True (block1: True) | 115 / 3220 | 0 (block1: 0) | True | True | True | True | pass |
| proper-p15361-k7-D4 | 829 | [0, 0, 2, 115, 829] | True (block1: True) | 115 / 3220 | 0 (block1: 0) | True | True | True | True | pass |
| unit-p12289-k8-D4 | 1061 | [0, 0, 2, 119, 1061] | True (block1: True) | 119 / 3808 | 0 (block1: 0) | True | True | True | True | pass |
| proper-p65537-k8-D4 | 1061 | [0, 0, 2, 119, 1061] | True (block1: True) | 119 / 3808 | 0 (block1: 0) | True | True | True | True | pass |
| unit-p65537-k8-D4 | 1061 | [0, 0, 2, 119, 1061] | True (block1: True) | 119 / 3808 | 0 (block1: 0) | True | True | True | True | pass |
| proper-p15361-k7-D5 | 3469 | [0, 11, 102, 466, 1467, 3469] | True (block1: True) | 1467 / 41076 | 0 (block1: 0) | True | True | True | True | pass |
| unit-p12289-k6-D5 | 1586 | [1, 13, 79, 299, 794, 1586] | True (block1: True) | 794 / 19056 | 0 (block1: 0) | True | True | False | True | pass |
| unit-p65537-k10-D5 | 18569 | [0, 0, 2, 138, 3163, 18569] | True (block1: True) | 3163 / 126520 | 0 (block1: 0) | True | True | True | True | pass |
| proper-p65537-k10-D5 | 18569 | [0, 0, 2, 138, 3163, 18569] | True (block1: True) | 3163 / 126520 | 0 (block1: 0) | True | True | True | True | pass |

The case `unit-p12289-k6-D5` is the case with 1 in V_D: the checker reports `one_excluded: false` (and `closed: true`; W is the whole space there, no free column, so `closed` is vacuous and `weight0_column_free` is false).

### 4.2 Bases of the reference engine `csd/ic_bounds_d5_verify.py`, converted into the checkpoint format by the harness

| case | reference rank | profile | closed | fail | g in W | 1 excluded | hashes = `Part.sha()` of the reference | outcome |
|---|---|---|---|---|---|---|---|---|
| unit-p12289-k7-D4 | 829 | [0, 0, 2, 115, 829] | True | 0 | True | True | True | pass |
| proper-p15361-k7-D4 | 829 | [0, 0, 2, 115, 829] | True | 0 | True | True | True | pass |
| proper-p15361-k7-D5 | 3469 | [0, 11, 102, 466, 1467, 3469] | True | 0 | True | True | True | pass |
| unit-p12289-k6-D5 | 1586 | [1, 13, 79, 299, 794, 1586] | True | 0 | True | False | True | pass |

### 4.3 Negative controls (unit, p = 12289, k = 8, D = 4; engine basis edited, both `all` and `block1` shown)

| control | what was done | checker result | reason reported | outcome |
|---|---|---|---|---|
| untouched-control | unmodified basis from the engine is accepted | valid=True, format_ok=True, closed=True, g in W=True, fail=0 (block1: valid=True, fail=0) |  | pass |
| low-row-coefficient-changed | part 0, row 170 (pivot weight 3), free column 978: +1 mod p | valid=False, format_ok=True, closed=False, g in W=True, fail=640 (block1: valid=False, fail=320) | closure: 640 (row, shift, target) triples leave W | pass |
| high-row-coefficient-changed-closure-only | part 0, row 1 (pivot column 6, weight 4), free column 224: +1 mod p; g and m*g only involve weights <= 3 < D so they stay in W; only the shift of some low row lands on this pivot. rows tried before this one: 0 | valid=False, format_ok=True, closed=False, g in W=True, fail=16 (block1: valid=False, fail=8) | closure: 16 (row, shift, target) triples leave W | pass |
| high-row-coefficient-harmless | no harmless high row among the first 200: every high row of part 0 is hit | - | - | pass |
| low-row-removed-all-tested | each low row of the first 80 of each part removed in turn; an accepted removal means the removed row was redundant (the smaller space is still closed and contains g), which is a valid smaller certificate | tested 119, rejected 119 (closure only 117, g only 0, both 2), format invalid 0, accepted 0 | closure / g | pass |
| low-row-removed-closure-only | part 0, row 0 removed; format valid, g still in W, rejected by the closure check | valid=False, format_ok=True, closed=False, g in W=True, fail=214 (block1: valid=False, fail=107) | some m*g (weight m <= 1) not in W; closure: 214 (row, shift, target) triples leave W | pass |
| interrupted-closure | 61 checkpoints kept (35 taken inside the closure phase); first, middle, last shown; all 35 checked | rejected 33 of 35; accepted 2 (each accepted checkpoint has the final basis hash: True) | closure failures | pass |
| engine-withheld-row-closure-done | engine test hook --withhold 0:3 --gens g --no-repair; the checkpoint is the closure_done one | valid=True, closed=True | none: the withheld row's shifts were redundant (the basis is the true closure; hash equal to a full run is not claimed) | pass |
| g-coefficient-changed | instance generator dict replaced; the basis was computed for the original g, so g is not in W and the dict disagrees with S_3(a, b, target[0]) | valid=False, format_ok=True, closed=True, g in W=False, fail=0 | generator:instance dict differs from the formula S_3(a, b, target[0]); g not in W; some m*g (weight m <= 1) not in W | pass |
| g-other-target | instance generator dict replaced; the basis was computed for the original g, so g is not in W and the dict disagrees with S_3(a, b, target[0]) | valid=False, format_ok=True, closed=True, g in W=False, fail=0 | generator:instance dict differs from the formula S_3(a, b, target[0]); g not in W; some m*g (weight m <= 1) not in W | pass |
| tile-byte-flipped | sha256 of a tile differs from the manifest | valid=False, format_ok=False, closed=None, g in W=None, fail=None | format:p0_t00000.npy: sha256 mismatch | pass |
| duplicate-pivot | two rows with the same pivot column | valid=False, format_ok=False, closed=None, g in W=None, fail=None | format:part 0: pivot columns not distinct | pass |
| coefficient-equal-p | an entry equal to p (outside [0, p)) | valid=False, format_ok=False, closed=None, g in W=None, fail=None | format:entries_in_0_p_and_integral | pass |

A corrupted high row (pivot weight D) that is hit by the shift of some low row is rejected by the closure check alone (g and m g only involve columns of weight <= 3 < D, so they stay in W): see `high-row-coefficient-changed-closure-only`. A high row that no shifted low row touches could not be found (every high row of the "+" part is hit), recorded in `high-row-coefficient-harmless`. In principle such a corruption is harmless: the edited space would still contain g and be closed, i.e. another valid certificate.

### 4.4 Specification, exactness and devices

* `certcheck-spec.json` (passed=True): own S_3 expansion = reference `s3_coeffs` and vanishes at x(P1 +- P2) (three curves); own layout, shifts and projection = reference `Layout` (k, D = 6,4 / 7,5 / 8,4); dim(W cap M_<D) = number of low rows (dense); dense full-coordinate closure test agrees with the checker on an untouched and a corrupted basis.
* `certcheck-arith.json` (passed=True): `modp` and block products against Python integers near the bound (cpu and fake-gpu).
* `certcheck-devices.json` (passed=True): identical verdicts, counts, profile and hashes for cpu, fake-gpu, tiny column/row blocks and a forced inner dimension of 7 (chunked reduction path); CLI `--selftest --device fake-gpu` passes.


## 5. Predicted cost at k = 14, D = 6, p = 65537

Data (from `runs/proper-p65537-k14-D6.json`): columns 499,178; parts 249,824 ("+") and 249,354 ("-");
ranks 224,688 and 224,234; free columns F = 25,136 and 25,120; low rows 36,845 and 36,742 (73,587 in total);
low part columns Lp = 61,272 and 61,166.

Work of the closure check (dense, float64): 2 n_pi Lp_pi F_q per (pi, q, shift) = 4.53e14 flop per shift, summed
over the four (pi, q); **1.27e16 flop for `--shifts all` (28 shifts), 6.34e15 for `block1`**.

| | GPU, 1 x 141 GB (H200 class) | CPU, 24 vCPU |
|---|---|---|
| sustained double-precision rate assumed | 30 to 45 Tflop/s (peak about 60, tensor cores) | 0.3 to 0.8 Tflop/s (local 4-core measurement: 48 Gflop/s per thread) |
| closure `all` | 5 to 11 min of products | 4.4 to 12 h |
| closure `block1` | 2.5 to 6 min | 2.2 to 6 h |
| T-matrix builds (112 x gathers of 5 GB rows, modp passes over 12 GB) | about 1 to 3 min | 20 to 40 min |
| file reading and hashing, twice (about 90 GB read, sha256 twice) | 3 to 8 min | same |
| whole run `all` (estimate) | **15 to 40 min** | **5 to 14 h** |
| memory | device: B of one part 22.6 GB + S of both parts 18.0 GB + T block 12.3 GB + row temporaries about 15 GB, peak about 70 GB (`--work-gb 45`); host: a few GB (one tile at a time, 0.1 to 0.3 GB) | RAM: same arrays, about 75 GB resident |

These are predictions from flop counts and published or measured rates, not measurements at the target: I could
not run at the target size or on a GPU. The engine's own sweep on the H200 took about 25 minutes, so a GPU time
within the same order for `all` (and about half for `block1`) is plausible, but the dense formulation does not skip
the leading zeros of the low rows (about a factor 2 in flop that is left unused: a safe overestimate).

## 6. Remote command

On the machine (files: `/workspace/code/ic_bounds_cert_check.py`, checkpoint `/workspace/ckpt/<name>/ck-*/`, instance in
`/workspace/in`, engine record in `/workspace/out`):

```
cd /workspace/code
python3 ic_bounds_cert_check.py --selftest --device gpu --out /workspace/out/certcheck-selftest-gpu.json   # first: exactness of modp/products and tiny pipeline on the GPU
python3 ic_bounds_cert_check.py --device gpu --shifts all \
  --instance /workspace/in/proper-p65537-s12-k14.json --ckpt /workspace/ckpt/<name> --D 6 \
  --engine-record /workspace/out/<engine-result>.json \
  --work-gb 45 --row-chunk 4096 --progress-seconds 30 \
  --status /workspace/out/certcheck-status.json --out /workspace/out/certcheck-<name>.json \
  2>&1 | tee /workspace/out/certcheck-<name>.log
```
CPU: replace `--device gpu` by `--device cpu --threads 24 --work-gb 40`. `--shifts block1` halves the time.
Progress lines (JSON) go to stdout; `--status` holds the latest progress (phase, done units of total, failures,
ETA, Tflop done). Exit codes: 0 certificate valid, 2 format ok but a fact fails, 3 format failure. The verdict JSON
holds all counts, `format_ok`, `g_in_W`, `closed`, `one_excluded`, `certificate_valid`, `all_true`, the dimension and
profile of W, the part hashes, per-file sha256 (every tile and npz), the sha256 of the checker source, the instance,
the manifest, timings, and the comparison with the engine record (rank, part ranks, profile, hashes, contains_one).

## 7. cupy calls that could not be exercised here (no GPU, no cupy)

`--device fake-gpu` runs the identical code path on numpy arrays that refuse to mix with host arrays and refuse
Python sequences or host arrays as indices (it caught real slips during development), but it cannot validate these
cupy behaviours: `cupy.asarray`, `cupy.asnumpy`, `cupy.zeros`, `cupy.empty`, `cupy.floor`, `cupy.where`;
slice reads and slice assignments on device arrays; advanced indexing with device int64 index arrays: row gather
`Bd[rows]`, `Bd[rows, f0:f1]` (index array combined with a slice), column gather `Yd[:, idx]`, scatter-add
`T[i, j] += v` and row update `T[rows] -= X` (indices are unique inside each call, so the read-modify-write is
well defined); in-place `*=`; `astype`; `.any(axis=1)`; boolean `|`; the matrix product `@` (cuBLAS dgemm) on
non-contiguous slices `A[:, a:b]`; `cupy.get_default_memory_pool().free_all_blocks()`;
`cupy.cuda.runtime.memGetInfo()`. The remote self-test (`--selftest --device gpu`) exercises all of these on tiny data and
the exactness of `modp` and the products at the inner-dimension bound on the device; run it first.

## 8. Limits

* The proof is about the stored basis: it shows the lemma's hypotheses for W; it does not show that the engine's W is
  the *least* such space (not needed). If `closed` is false for a basis taken before the closure ended, that says
  nothing against the final one.
* The dense formulation costs about twice what a sparse one could (leading zeros of low rows are not skipped).
* No resume: a failure after an hour restarts from the beginning (the status file shows how far it got).
* Not verified by me: GPU behaviour (section 7), timing at the target (section 5), and the dimension/profile of W
  at k = 14 (they come out of the run).
* A modified checkpoint that changes only the files and the manifest consistently passes the file-hash check by
  construction; the hashes protect against transport errors, the mathematical checks against wrong content.
