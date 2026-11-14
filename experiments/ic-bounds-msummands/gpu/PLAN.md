# GPU plan for the m-summand cells that were not run locally

Status 2026-11-12. **Nothing in this directory has been run on a GPU** (and no rented machine was used). `chain_msum.py` was only exercised
with `--dry` and, on small instances of the same formats, locally: `ic_bounds_msum_run.py instance` -> `ic_bounds_bin_big_run.py ladder --device cpu`
and `ic_bounds_msum_run.py lfd-file` reproduce the recorded LFD of `B-m4n6k2-ordA-ran-pla` (8) and `P-m4k2-p97-pl` (9). Labels: exact (measured here),
model (extrapolation, unverified).

## Which engine runs what

* **Binary** (any m): `csd/ic_bounds_bin_big.py` through `scripts/ic_bounds_bin_big_run.py ladder` (bit-packed GF(2), cupy). It is generic in the equation list
  and in r = number of variables (columns = squarefree monomials of degree <= D in r variables); nothing in it assumes m = 2 (the instance builder
  `ic_bounds_bin_big_run.py instance` does: it builds `[basis, basis]` systems; this plan builds instances with `ic_bounds_msum_run.py instance`, format
  `ic-bin-instance/1` plus `m`, `msum_system_hash`). Verified here on the cpu device against the python engine on 10 systems (`validation/engines-binary.json`).
* **Prime, m >= 3**: `csd/ic_bounds_big.py` is specific to m = 2: its `Layout(k, D, split)` holds exactly two exponent blocks (`a1`, `a2`, key `a1 | a2 << k`),
  the swap split (`split=True`) pairs the two blocks, the generator is a dict keyed by `(e1, e2)`, `chain_generator`/`solutions`/`evaluate_full`/`eval_rank` are the S_3 chain of
  two summands, and `count_columns(k, D)` is `sum C(2k, w)`. With `split=False` the swap projection disappears, but the two-block layout, the shift list
  (`i in (0, 1)`) and the generator/evaluation code remain. **Missing**: an m-block `Layout` (columns = squarefree masks of m k variables, shifts for all m k variables, one part),
  generator rows keyed by m-tuples, `evaluate_full`/`eval_rank` over m exponent vectors, and the grid-based generator of `csd/ic_bounds_msum.py`. That is a new file
  importing `ic_bounds_big` (Part, Device, arithmetic kernels are reusable); it was NOT written. Until it exists the prime cells below run on the CPU engine
  `CompactClosure` (generic in m; used here for m = 3, 4 at k <= 4), which has no checkpointing and keeps the basis in host memory (uint16): the chain script runs it as is.
  The prime cost rows are therefore for the CPU engine's operation count, converted to time at an assumed fp64 rate (model).

## Cost model (stated assumptions)

* Binary: engine's own `plan()` (memory exact bookkeeping) and `0.28 N^3` float32 flops at 30 TFLOP/s (`experiments/ic-bounds-binary-deep/PLAN.md`); pessimistic column x4.
  Calibration on this machine (cpu, 3 threads): `B-m5n8k3-koblitz0-ran-ord` N = 32,768 columns at D = 15, 16: 958 s for the whole ladder = `0.28 N^3 / 958 s` = 10 GFLOP/s effective,
  i.e. the GPU model (30 TFLOP/s) is a 3000x extrapolation of this one measured point and the earlier single-thread calibration (26 GFLOP/s) of the binary-deep plan.
* Prime (CompactClosure operation count): `madds = (m k) n^3 / 6` for a closure that fills to rank ~ n (every absorbed row is shifted by each of the m k variables and reduced
  against rank x free columns); time = `2 madds / 20 TFLOP/s` (fp64 on an H200 would need the missing engine); memory = stored residues `n^2 / 4` entries at rank n/2
  (float32 on the GPU engine, uint16 on CompactClosure). Calibration: m = 3, k = 5, D = 7, n = 16,384: model 1.1e13 madds, measured 210-450 s on one CPU thread (recorded).
  The m = 3, k = 6, D = 7 run (63,004 columns, 1381 s) closed with rank 19,958 and is not comparable to a full-rank closure.

## Cells (all instances are fixed by the specs in `specs/`; seeds as in `scripts/ic_bounds_msum_run.py`)

### Binary (the selected curves/bases/targets: koblitz0 random base order-two target; koblitz1 invariant-if-exists base random target; ordA random base planted target)

| cell | m | n | k | r | d_in (source) | \|Z\| | cap D | columns at cap | basis worst case GB | device peak GB (plan) | fits 141 GB | model s | pessimistic s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G-m4n9k4-koblitz0-ran-ord | 4 | 9 | 4 | 16 | 8 (built) | 66 | 10 | 58651 | 0.13 | 16.7 | True | 1.9 | 8 |
| G-m4n9k4-koblitz1-inv-ran | 4 | 9 | 4 | 16 | 12 (built) | 115 | 14 | 65519 | 0.16 | 20.0 | True | 2.6 | 10 |
| G-m4n9k4-ordA-ran-pla | 4 | 9 | 4 | 16 | 12 (built) | 312 | 14 | 65519 | 0.16 | 20.0 | True | 2.6 | 10 |
| G-m3n10k6-koblitz0-ran-ord | 3 | 10 | 6 | 18 | 3 (built) | 412 | 5 | 12616 | 0.01 | 3.7 | True | 0.0 | 0 |
| G-m3n10k6-koblitz1-inv-ran | 3 | 10 | 6 | 18 | 6 (built) | 1227 | 8 | 106762 | 0.42 | 23.8 | True | 11.4 | 46 |
| G-m3n10k6-ordA-ran-pla | 3 | 10 | 6 | 18 | 6 (built) | 375 | 8 | 106762 | 0.42 | 23.8 | True | 11.4 | 46 |
| G-m4n10k5-koblitz0-ran-ord | 4 | 10 | 5 | 20 | 8 (evaluated (equations not built)) | 2086 | 10 | 616666 | 13.84 | 39.4 | True | 2188.7 | 8755 |
| G-m4n10k5-koblitz1-inv-ran | 4 | 10 | 5 | 20 | 12 (evaluated (equations not built)) | 25016 | 14 | 1026876 | 38.37 | 66.3 | True | 10106.3 | 40425 |
| G-m4n10k5-ordA-ran-pla | 4 | 10 | 5 | 20 | 12 (evaluated (equations not built)) | 564 | 14 | 1026876 | 38.37 | 66.3 | True | 10106.3 | 40425 |
| G-m5n10k4-koblitz0-ran-ord | 5 | 10 | 4 | 20 | 15 (evaluated (equations not built)) | 4806 | 17 | 1048365 | 39.99 | 68.0 | True | 10754.1 | 43016 |
| G-m5n10k4-koblitz1-inv-ran | 5 | 10 | 4 | 20 | 20 (evaluated (equations not built)) | 0 | 20 | 1048576 | 40.01 | 68.1 | True | 10760.6 | 43042 |
| G-m5n10k4-ordA-ran-pla | 5 | 10 | 4 | 20 | 20 (evaluated (equations not built)) | 125 | 20 | 1048576 | 40.01 | 68.1 | True | 10760.6 | 43042 |
| G-m3n11k7-koblitz0-ran-ord | 3 | 11 | 7 | 21 | 6 (bound, not computed) | None | 8 | 401930 | 5.88 | 30.3 | True | 606.0 | 2424 |
| G-m3n11k7-koblitz1-inv-ran | 3 | 11 | 7 | 21 | 6 (bound, not computed) | None | 8 | 401930 | 5.88 | 30.3 | True | 606.0 | 2424 |
| G-m3n11k7-ordA-ran-pla | 3 | 11 | 7 | 21 | 6 (bound, not computed) | None | 8 | 401930 | 5.88 | 30.3 | True | 606.0 | 2424 |
| G-m3n12k8-koblitz0-ran-ord | 3 | 12 | 8 | 24 | 6 (bound, not computed) | None | 8 | 1271626 | 58.84 | 86.1 | True | 19191.8 | 76767 |
| G-m3n12k8-koblitz1-inv-ran | 3 | 12 | 8 | 24 | 6 (bound, not computed) | None | 8 | 1271626 | 58.84 | 86.1 | True | 19191.8 | 76767 |
| G-m3n12k8-ordA-ran-pla | 3 | 12 | 8 | 24 | 6 (bound, not computed) | None | 8 | 1271626 | 58.84 | 86.1 | True | 19191.8 | 76767 |
| G-m4n12k6-koblitz0-ran-ord | 4 | 12 | 6 | 24 | 12 (bound, not computed) | None | 14 | 14198086 | 7333.74 | 7410.2 | False | 26713216.5 | 106852866 |
| G-m4n12k6-koblitz1-inv-ran | 4 | 12 | 6 | 24 | 12 (bound, not computed) | None | 14 | 14198086 | 7333.74 | 7410.2 | False | 26713216.5 | 106852866 |
| G-m4n12k6-ordA-ran-pla | 4 | 12 | 6 | 24 | 12 (bound, not computed) | None | 14 | 14198086 | 7333.74 | 7410.2 | False | 26713216.5 | 106852866 |

Notes. r = 16 (`G-m4n9k4-*`) are the three r = 16 cases that the local column budget (40,000) excluded: the CPU measurement above predicts about 50 min to 2 h each on this machine
(0.28 N^3 = 7.8e13 flops at 10-26 GFLOP/s). d_in is exact for the cells marked built/evaluated and a theoretical upper bound (6 for m = 3, 12 for m = 4: popcount(e) <= 3 for e <= 8 per block) for r > 20; for
the order-two targets the true d_in is smaller (3 for m = 3 here), so those cap degrees are over-estimates. `G-m4n12k6-*` (r = 24, 14.2 M columns, 7.4 TB basis) does not fit any single device and is listed only
to size the question: not plannable. m = 5 at r = 20 (`G-m5n10k4-*`): d_in = 20 = r for two of the three instances (the whole space, 1,048,576 columns, D = 20), so the closure is the full-space one.
The unit instance (`G-m5n10k4-koblitz1-inv-ran`, |Z| = 0) needs D up to the LFD; its cap is r.

### Prime (engine: CompactClosure on cpu; see above)

| cell | m | k | p | deg g | \|Z\| | Prop. 1 variables | D | columns | madds model | model s (20 TFLOP/s fp64) | pessimistic x4 | stored GB (float32 / uint16 at rank n/2) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GP-m4k4-p97-pl | 4 | 4 | 97 | 12 | 416 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p97-pl | 4 | 4 | 97 | 12 | 416 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m4k4-p193-r | 4 | 4 | 193 | 12 | 316 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p193-r | 4 | 4 | 193 | 12 | 316 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m4k4-p257-pl | 4 | 4 | 257 | 12 | 553 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p257-pl | 4 | 4 | 257 | 12 | 553 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m4k4-p449-r | 4 | 4 | 449 | 12 | 60 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p449-r | 4 | 4 | 449 | 12 | 60 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m4k4-p641-pl | 4 | 4 | 641 | 12 | 163 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p641-pl | 4 | 4 | 641 | 12 | 163 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m4k4-p7681-r | 4 | 4 | 7681 | 12 | 0 | 12 | 12 | 64839 | 7.27e+14 | 72.7 | 290.8 | 3.92 / 1.96 |
| GP-m4k4-p7681-r | 4 | 4 | 7681 | 12 | 0 | 12 | 13 | 65399 | 7.46e+14 | 74.6 | 298.4 | 3.98 / 1.99 |
| GP-m5k3-p97-pl | 5 | 3 | 97 | 15 | 96 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p97-pl | 5 | 3 | 97 | 15 | 96 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p193-r | 5 | 3 | 193 | 15 | 305 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p193-r | 5 | 3 | 193 | 15 | 305 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p257-pl | 5 | 3 | 257 | 15 | 1 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p257-pl | 5 | 3 | 257 | 15 | 1 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p449-r | 5 | 3 | 449 | 15 | 0 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p449-r | 5 | 3 | 449 | 15 | 0 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p641-pl | 5 | 3 | 641 | 15 | 736 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p641-pl | 5 | 3 | 641 | 15 | 736 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p7681-r | 5 | 3 | 7681 | 15 | 0 | 15 | 15 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m5k3-p7681-r | 5 | 3 | 7681 | 15 | 0 | 15 | 16 | 32768 | 8.80e+13 | 8.8 | 35.2 | 1.0 / 0.5 |
| GP-m3k6-p7681-r | 3 | 6 | 7681 | None | None | None | 8 | 106762 | 3.65e+15 | 365.1 | 1460.3 | 10.62 / 5.31 |
| GP-m3k6-p7681-pl | 3 | 6 | 7681 | None | None | None | 8 | 106762 | 3.65e+15 | 365.1 | 1460.3 | 10.62 / 5.31 |
| GP-m3k7-p7681-r | 3 | 7 | 7681 | None | None | None | 7 | 198440 | 2.73e+16 | 2735.0 | 10940.0 | 36.67 / 18.34 |
| GP-m3k7-p7681-r | 3 | 7 | 7681 | None | None | None | 8 | 401930 | 2.27e+17 | 22725.8 | 90903.2 | 150.45 / 75.23 |
| GP-m3k7-p7681-pl | 3 | 7 | 7681 | None | None | None | 7 | 198440 | 2.73e+16 | 2735.0 | 10940.0 | 36.67 / 18.34 |
| GP-m3k7-p7681-pl | 3 | 7 | 7681 | None | None | None | 8 | 401930 | 2.27e+17 | 22725.8 | 90903.2 | 150.45 / 75.23 |

Notes. `GP-m4k4-*` (65,536 monomials): the discriminating cell; frozen predictions in `../predictions-2026-11-12.json` (G7: LFD in [13, 17], point 13 = d + 1; Proposition 1 applies with 12 of 16 variables).
Only D = 12 and 13 need a closure (D < 12 has no input). `GP-m5k3-*` (32,768): d = 15 = m k, LFD = 16 forced (G8). Closures here stop at D = 16 where the certificate is trivial (all monomials). `GP-m3k6-*` D = 8 and `GP-m3k7-*`
are extensions of the recorded m = 3 table (k = 6 D = 7: one instance, recorded; k = 7 D = 7: 198,440 columns); their generators were not built here (m = 3 builders exist; the instance step builds them on the pod:
2^18 and 2^21 grid evaluations). m = 3, k = 7, D = 8 (401,930 columns) needs 150 GB in float32: it does not fit one 141 GB device without host spill (the engine's spill mode exists in `ic_bounds_big.py`, not for these layouts).
Run time on this machine is dominated by the CPU: the model above says minutes to hours on a fp64 GPU only if the missing m-block engine is written; with CompactClosure on a 16-core host expect
`madds / (3e11 madds/s)`: about 40 min (m = 4, k = 4) and 10 min (m = 5, k = 3) per instance (model).

## How to run (everything reproducible from the repository; nothing is uploaded from here)

```
# pod layout: code/{csd/, scripts/, gpu/{cells.json, cost.json, specs/}, chain_msum.py}, WORKSPACE=/workspace
python3 chain_msum.py --device gpu --family binary --deadline 14400          # selftest, then every binary cell in cells.json order
python3 chain_msum.py --device cpu --family prime --skip-selftest --only GP-m5k3-p97-pl GP-m4k4-p97-pl
# each step builds its instance (ic_bounds_msum_run.py instance), compares the system hash / generator hash with cost.json, runs the ladder,
# writes out/<cell>.json, out/results.jsonl and out/chain-status.json
```

Order recommended if the budget forces a cut: `G-m4n9k4-*` (3), `GP-m5k3-p97-pl`, `GP-m4k4-p97-pl`, `GP-m4k4-p257-pl`, `G-m3n10k6-*`, `G-m4n10k5-koblitz0`, `G-m5n10k4-ordA`, then the rest. Cut first: r = 24 cells, m = 3 k = 7.
Reading decided in advance: a certified LFD for `G-m4n9k4-koblitz1-inv-ran` / `ordA` (d_in = 12, r = 16) above `min(r, d_in + 2) = 14`, or equal to 17, refutes prediction G4; `GP-m4k4` LFD 17 on any instance refutes (b), 13 supports it.
