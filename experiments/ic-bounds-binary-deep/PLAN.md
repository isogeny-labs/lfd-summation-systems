# Binary deep closure: engine, validation and preregistered run plan

Status 2026-11-12. Nothing in this directory has been run at k >= 11. Labels: **exact** (certified computation named),
**proof**, **empirical**, **conjecture**, **unverified** (not run anywhere).

## 1. What was built

* `csd/ic_bounds_bin_big.py` (engine), `scripts/ic_bounds_bin_big_run.py` (driver: `closure`, `ladder`, `plan`,
  `selftest`, `instance`, `zcheck`).  The engine needs `csd/ic_bounds_big.py` only for its device layer (cpu / cupy /
  `fake-gpu`); on the pod both files are uploaded next to each other.
* **Definition matched** (exact; the same text is at the top of the engine): `csd/ic_boolean_closure.py::boolean_degree_closure`
  and `csd/ic_local_certificate.py::local_profile`.  Columns = squarefree monomials of degree <= D, order
  `sorted(key=(popcount, mask), reverse=True)` (`monomial_masks`); V_D = least subspace containing the inputs of degree
  <= D and closed under f -> x_i f for every variable and every element of degree < D (S -> S u {i}, colliding sources
  add mod 2: `_products`); fall at D iff dim(V_D cap M_<=D-1) > dim V_(D-1) (`low_degree_dimension`); unit iff the
  constant column is a pivot; certificate: evaluation rank of M_<D on Z = |Z| = codim V_D, LFD = last strict fall.
  No swap split (there is none in characteristic 2).
* **Storage: bit-packed uint64 words**, 1 bit per entry, row tiles, lazy deletion of pivot columns and compaction.
  Reason: the reduced echelon basis has at worst rank x (N - rank) entries, N^2/4 bits = 41 GB at N = 1,149,017
  (2k = 32, D = 6) against 1.3 TB as float32 and 165 GB as uint8.  float32 appears only in transient blocks: a block of
  the basis is unpacked to 0/1 float32, multiplied (exact: every partial sum is an integer below 2^24, asserted
  `N < 2^24`; the selftest checks a 200,000-term sum on the device), the parity is packed back and XORed in.
* Algorithm: RREF of a residual block by panels of 64 columns (host Gauss-Jordan on the panel word, one product to apply
  it); one product `C_old @ New` per tile to add the new rows; sparse-aware residuals `Y_free + Y_piv @ B` using only
  the pivot columns on which the batch is nonzero; every low row is shifted once at its state when its batch starts;
  deterministic full sweep (every variable x every row of pivot degree < D, plus every input, reduces to zero) as the
  certificate of closedness; a nonzero residual is absorbed (repair) and a new sweep starts.  Checkpoint: one, at the
  end of a closure (`--ckpt`; `--resume` loads it and runs the certificate sweep).  Canonical output: `rref_sha256`
  (depends on the subspace and the column order only), `profile` p_j = dim(V_D cap M_<=j), `quotient_ranks`, falls,
  `local_checks`.

## 2. Validation (exact on the data named; `validation/`)

`validation/validate.py`: 65 systems, all equal (`validation/records.jsonl`, `VALIDATION_TABLE.md`, `summary.json`).
Drawn from `experiments/ic-curve-sweep` (k = 4..7), `experiments/ic-bounds-binary` (k = 8..10) and the subfield
construction (n = 8, 12, 16; V = F_(2^k) in F_(2^2k), order-two target, a = 0 and 1): classes affine 21, order-two 19,
unit (|Z| = 0, not const) 13, const 12; k = 4..10; Koblitz 36, ordinary 29; Frobenius-invariant bases 35; scaled
subfields 6.  For every D up to the certificate degree (268 degrees in total): dim V_D equal to a live run of the
existing engine and to the stored certified record (59 systems have one), fall profile, LFD, certificate degree, local
checks and |Z| equal, and the canonical RREF hash of the new basis equal to the hash of the existing engine's basis
re-reduced independently (65 of 65); cpu, and `fake-gpu` with small tiles and batches (batch 3..200, tiles 8..512 rows,
column blocks 64..1024, lazy-column compaction every few pivots, basis tiles on the host for k = 6..8): all hashes
equal between the two modes.  Instances are rebuilt from the stored specifications by the new independent builders and
their `system_hash` equals the stored one; |Z| from the new field enumerator equals the stored count (59 systems with a
stored record) and the certificate itself (codim V_D = |Z| = evaluation rank) certifies completeness of Z in all 65; the
enumerator equals `binary_solutions` set-wise on 6 further systems and on the k = 10 calibration instance (`zcheck`).
**Negative controls** (`validation/negative.json`): (a) withholding the shifts of one row: with `--no-repair` the sweep
fails (`closed` false, rank deficit up to 164) whenever the withheld row matters (5 of 16 trials; the others are
redundant rows and the closure is still complete, correctly); with repair the sweep repairs (up to 24 repairs) and the
hash equals the clean hash in 16 of 16; (b) a basis row zeroed in a finished closure: the sweep reports it in 4 of 4.
Also exact: `closure --ckpt` then `--resume` reproduces the hash; chain (below) end to end.
Not covered by this validation: k >= 11 and D >= 6 against the existing engine (the existing engine cannot reach them);
the cupy device itself.

## 3. Cost model (empirical calibration, then extrapolation = conjecture)

Measured on this machine (one thread, `closure`, `flops_update + flops_residual` counted by the engine, total / N^3):
random k = 10, D = 5 (N = 21,700): 0.283 (112 s); subfield k = 10, D = 5: 0.309; random k = 8, n = 17, D = 6
(N = 14,893): 0.175.  The plan uses **0.28 N^3 float32 flops** per closure at the certificate degree (the lower degrees
cost under 1% of it), at an assumed **25 TFLOP/s** effective on the H200 (conservative against 67 TFLOP/s peak
float32; **unverified**).  A pessimistic column multiplies by 4 (kernel launch overheads, the device pool off,
residual cost above the measured range).  The `plan` subcommand prints memory exactly from the buffers named in
`csd/ic_bounds_bin_big.plan` (basis worst case at rank N/2, transient blocks capped at `--work-gb`, tile slack, CUDA
context 1.5 GB): `plan-numbers.json`.

| 2k | D | columns | basis worst case (GB, bit-packed) | predicted device peak (GB) | N^3/3 update flops | total estimate at 25 TFLOP/s (s) |
|---|---|---|---|---|---|---|
| 22 | 5 | 35443 | 0.05 | 8.1 | 1.48e+13 | 0 |
| 22 | 6 | 110056 | 0.44 | 23.8 | 4.44e+14 | 12 |
| 22 | 7 | 280600 | 2.87 | 26.8 | 7.36e+15 | 206 |
| 24 | 5 | 55455 | 0.11 | 15.1 | 5.68e+13 | 2 |
| 24 | 6 | 190051 | 1.32 | 24.9 | 2.29e+15 | 64 |
| 24 | 7 | 536155 | 10.46 | 35.1 | 5.14e+16 | 1438 |
| 26 | 5 | 83682 | 0.26 | 23.5 | 1.95e+14 | 6 |
| 26 | 6 | 313912 | 3.59 | 27.5 | 1.03e+16 | 289 |
| 28 | 5 | 122438 | 0.55 | 23.8 | 6.12e+14 | 17 |
| 28 | 6 | 499178 | 9.07 | 33.4 | 4.15e+16 | 1161 |
| 30 | 5 | 174437 | 1.11 | 24.5 | 1.77e+15 | 50 |
| 30 | 6 | 768212 | 21.47 | 46.5 | 1.51e+17 | 4231 |
| 32 | 5 | 242825 | 2.15 | 25.7 | 4.77e+15 | 134 |
| 32 | 6 | 1149017 | 48.04 | 73.9 | 5.06e+17 | 14158 |

D = 6 at 2k = 32 (4 h) and 2k = 30 (1.2 h) are what the budget cannot carry for more than one system: the cells below
cap the degree at 5 for k = 15, 16 (the LFD is expected to be 4; if a certified LFD 5 shows up at k >= 15 and its
certificate needs D = 6, that single closure is a decision for the main session).  TF32 tensor cores (`CUPY_TF32=1`) would
be exact for 0/1 operands and fp32 accumulation below 2^24 and could be several times faster; **not used, not verified**:
the self-test (which includes the exactness cases) would be the gate.

## 4. Preregistered cells (frozen file `predictions-2026-11-12.json`)

`make_cells.py` built every instance (equations by `binary_descent`, Z by the field enumerator, structure e = r_eff -
m_eff by linear algebra) and wrote `predictions-2026-11-12.json` **before any closure at k >= 11** (none was run; the only
closures run are the k <= 10 validation and calibration ones).  The file is write-once.

**sha256 of `predictions-2026-11-12.json` (file bytes): `7feee5dbafbde72b4b8ddf94b310ee4b8b56a45446b886bf61b8f98c038e33e1`.**

Laws on the table (affine class, per system; hit iff the certified LFD equals the prediction): B = 2 + [k+e >= 3] +
[k+e >= 9] + [k+e >= 19]; post-hoc = max(2 + [k+e >= 3] + [k+e >= 9], A(e)) (conjecture of the law README, fitted on the
data it was scored on); A(e) = 2/3/4/5 for e <= -1 / <= 3 / <= 8 / >= 9; A_Z = 5 if |Z| > 470 else 4; for the scaled
subfield with order-two target the written inverse-family value max(3, k/2 + 1) (6 at k = 10, 7 at k = 12).  Const
systems: LFD 2 (S2).  Order-two systems and units (|Z| = 0, not const) are recorded without a scored prediction (B and
post-hoc are stored as `*_unscored`; the earlier P2/P3 scoring did include units in the affine class, and the draws at
n = 2k - 1 are mostly units: they are scored as affine here too, and reported separately).

What each cell decides (the numbers are in the table):

* C1 (14,11): B 5 vs post-hoc/A/A_Z 4; the discriminating cell of the README.  C6a (16,12), C6b (18,13), R1, R2: the
  same discrimination at k + e = 19..21 (e = 8 or 7), where B says 5 and the others 4.
* C3 (n = 2k - 1, k = 11..16): all small-|Z| draws; B, post-hoc and A_Z say 4, A says 3.  A certified 5 anywhere refutes
  all four (the k-driven floor rising to 5) and the largest k with LFD 4 brackets where it rises.
* C2 (trace-zero hyperplane n = k + 1, k = 11..13, e = k): every law says 5.  An LFD 6 (certificate at D = 6) is the
  first 6 on a proper non-subfield space and answers "where does LFD 6 first occur".
* C4 (Koblitz, random and Frobenius-invariant bases at (15,11) and (21,12)): the invariant bases have r_eff = 2k more
  often (V inside ker Tr); (21,12) with e = 3..4 sits on the A thresholds (3 vs 4).
* C5 (subfield k = 10 and 12): 6 and 7.

Largest degree that may be needed = largest predicted + 1 (certificate degree is LFD or, in 8 of 1486 old systems,
LFD + 1), capped as in the table (`max_degree` in `cells.json`).

Cells in run order (cheap and decisive first; `cells.json` is what the chain reads).  Time = sum of the per-system
nominal estimates (0.28 N^3 / 25 TFLOP/s at the cap degree); device peak from `plan` with `--work-gb 20`:

| order | cell | (n,k) | systems | bases / curves | classes | e | k+e | \|Z\| range | laws (B / post-hoc / A / inverse) | D cap | columns at cap | device peak (GB) | nominal time | pessimistic (x4) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | C1 | (14,11) | 8 | invariant+random / koblitz0+koblitz1+ordinary | {'affine': 7, 'order-two': 1} | 8 | 19 | 210-335 | 5 / 4 / 4 | 6 | 110056 | 24 | 616 s | 2464 s |
| 2 | C2a | (12,11) | 3 | trace-zero / koblitz0+koblitz1 | {'affine': 2, 'order-two': 1} | 11 | 22 | 1035-2070 | 5 / 5 / 5 | 6 | 110056 | 24 | 231 s | 924 s |
| 3 | C3a | (21,11) | 3 | random / koblitz0+ordinary | {'affine': 3} | 1 | 12 | 0-0 | 4 / 4 / 3 | 5 | 35443 | 8 | 183 s | 732 s |
| 4 | C3b | (23,12) | 2 | random / ordinary | {'affine': 2} | 1 | 13 | 0-4 | 4 / 4 / 3 | 5 | 55455 | 15 | 124 s | 496 s |
| 5 | C2b | (13,12) | 2 | trace-zero / koblitz0 | {'affine': 2} | 12 | 24 | 4004-4004 | 5 / 5 / 5 | 6 | 190051 | 25 | 296 s | 1184 s |
| 6 | C5a | (20,10) | 2 | subfield / koblitz0+koblitz1 | {'order-two': 2} | 10 | 20 | 1023-1023 | inverse family 6 | 7 | 137980 | 24 | 188 s | 752 s |
| 7 | C6a | (16,12) | 4 | random / koblitz0+ordinary | {'affine': 4} | 8 | 20 | 238-284 | 5 / 4 / 4 | 6 | 190051 | 25 | 592 s | 2368 s |
| 8 | C3c | (25,13) | 2 | random / ordinary | {'affine': 2} | 1 | 14 | 0-2 | 4 / 4 / 3 | 5 | 83682 | 24 | 136 s | 544 s |
| 9 | C4a | (15,11) | 4 | invariant+random / koblitz0 | {'affine': 4} | 7 | 18 | 88-124 | 4 / 4 / 4 | 5 | 35443 | 8 | 244 s | 976 s |
| 10 | C2c | (14,13) | 2 | trace-zero / koblitz1 | {'affine': 2} | 13 | 26 | 8234-8234 | 5 / 5 / 5 | 6 | 313912 | 28 | 916 s | 3664 s |
| 11 | C4b | (21,12) | 4 | invariant+random / koblitz1 | {'affine': 4} | 3 | 15 | 2-20 | 4 / 4 / 3 | 5 | 55455 | 15 | 248 s | 992 s |
| 12 | C6b | (18,13) | 3 | random / ordinary | {'affine': 3} | 8 | 21 | 226-286 | 5 / 4 / 4 | 6 | 313912 | 28 | 1374 s | 5496 s |
| 13 | C3d | (27,14) | 2 | random / ordinary | {'affine': 2} | 1 | 15 | 2-4 | 4 / 4 / 3 | 5 | 122438 | 24 | 168 s | 672 s |
| 14 | C5b | (24,12) | 1 | subfield / koblitz0 | {'order-two': 1} | 12 | 24 | 4095-4095 | inverse family 7 | 7 | 536155 | 35 | 2045 s | 8180 s |
| 15 | C3e | (29,15) | 1 | random / ordinary | {'affine': 1} | 1 | 16 | 0-0 | 4 / 4 / 3 | 5 | 174437 | 24 | 128 s | 512 s |
| 16 | C3f | (31,16) | 1 | random / ordinary | {'affine': 1} | 1 | 17 | 0-0 | 4 / 4 / 3 | 5 | 242825 | 26 | 244 s | 976 s |
| 17 (reserve) | R1 | (14,11) | 9 | random / koblitz0+ordinary | {'affine': 9} | 8 | 19 | 234-292 | 5 / 4 / 4 | 6 | 110056 | 24 | 693 s | 2772 s |
| 18 (reserve) | R2 | (17,12) | 4 | random / ordinary | {'affine': 4} | 7 | 19 | 104-140 | 5 / 4 / 4 | 6 | 190051 | 25 | 592 s | 2368 s |

Non-reserve total: 7,733 s nominal (2.1 h, $11), 8.6 h pessimistic ($45).  Reserve cells R1, R2: 1,285 s nominal.
Setup (pod creation, pip install of cupy, self-test): about 15 to 30 minutes.  **Expected cost $12 to $20, pessimistic
$55** at $5.29/h; the pod is terminated by the copy loop when the chain is finished.

**What to cut first**, in this order: R1/R2 (reserve, already last); C5b (one system, 34 min nominal, the most
expensive); C3e, C3f (k = 15, 16 units: the largest columns, little new information beyond C3a-d); C2c and C6b
(k = 13 at D = 6, 5 to 23 min each); then the Koblitz cells C4a/C4b.  Keep C1, C2a, C3a, C3b, C6a, C5a: they carry the
three README questions at k = 11, 12.  The chain skips (and records) a system whose estimate no longer fits before
`--deadline` and kills one that exceeds its own limit (`timeout`).

## 5. Pod side

`pod/chain_bin.py` (self-test through `nopool_run.py` = `ic-bounds-k14/relaunch/nopool_dev_only.py`, device pool off,
pinned pool on; then the cells one system after another; status JSON after every step; every record in
`/workspace/out` as soon as it exists: per-degree `<system>.json.D<d>.json`, final `<system>.json`, `results.jsonl`;
logs in `/workspace/logs`) and `pod/launch.sh` (copy loop `ic-bounds-k14/relaunch/mirror.sh` started first and detached,
logs and every record in its copy list; `--datacenters US-NC-1,EU-FR-1` polled once a minute for up to three hours;
no checkpoints, so the disk is 60 GB; lifetime cap = argument; the chain gets `--deadline` = cap - 1800 s).
Tested locally end to end with the simulated device on two tiny cells and one budget-skipped cell (`WORKSPACE`
override, `--device fake-gpu --tiny`); `bash -n` on launch.sh; **launch.sh itself was not run** (it calls Runpod).

Launch (main session only):

```
export RUNPOD_API_KEY=...   # from the secret store
export POD_STATE_DIR=<a directory that survives the session>
experiments/ic-bounds-binary-deep/pod/launch.sh 36000 60
```

## 6. Unverified / uncertain

* Everything on the card: the 25 TFLOP/s rate, memory use, the device pool being off (many transient allocations per
  product: if the self-test timing on the pod is several times slower than the local per-flop figure, re-run the chain
  with the default pool: `chain_bin.py --no-nopool`; the lifetime lessons say the pool was a problem for the k = 14 lift,
  not necessarily here).
* cupy calls never executed (the `fake-gpu` stand-in only checks that no host array is mixed into device operations and
  runs numpy semantics): `cp.asarray`, `cp.asnumpy`, `cp.zeros`, `cp.ascontiguousarray`, `ndarray.view(uint8/uint64)`
  on a contiguous 2-d array, `reshape`, uint8 `>>`, `<<`, `&`, `^=` on slice views, `sum(axis=2, dtype=uint8)`,
  `astype`, float32 `@` (cuBLAS SGEMM) at sizes up to 8k x 8k x 64k, `flatnonzero`, `any(axis=...)`, `concatenate`, fancy
  indexing and assignment with int64 device index arrays including the mixed forms `a[idx, w0:w1]`, `A[:, idx] = B`,
  `A[i, j] = A[i, j] ^ b`, `cp.cuda.runtime.memGetInfo`, `get_default_memory_pool().total_bytes()` (zero with the pool off).
* Cost model: three calibration points only (0.175 to 0.309 N^3); k >= 11 structure (rank growth, residual density) is
  extrapolated; the residual part could be larger on structured (subfield, trace-zero) systems at large N.
* The predictions are about LFD; a closure that is not certified (no certificate degree up to the cap) is recorded as
  `capped` with a lower bound and counts as neither hit nor miss.
