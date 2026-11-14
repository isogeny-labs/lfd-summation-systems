# Binary deep closures: results (2026-11-12)

Two-summand binary systems (Weil descent of S_3 over F_(2^n), factor base of dimension k, 2k Boolean variables) at
k = 11 to 16, and the scaled-subfield family at k = 10 and 12. Plan, cells and frozen predictions: `PLAN.md`,
`cells.json`, `predictions-2026-11-12.json`. Labels: **exact** (certified computation), **empirical**, **conjecture**.

## What was run

- 57 systems in 18 cells, all of the plan including the two reserve cells, on one rented NVIDIA H200
  (pod `icb-bin-deep2`, data centre US-NC-1), engine `csd/ic_bounds_bin_big.py` (sha256 `3ba6e423cbe04c8c…`,
  the same in every record), chain `pod/chain_bin.py`. On-card self-test: `ok`, 7 cases (`fetched/selftest-gpu.json`).
- Chain start 2026-11-12T01:50:13Z, end 03:07:51Z (4,658 s); 57 of 57 systems `done`, none timed out, none skipped,
  0 repairs in every closure. Sum of the closure times: 4,601 s.
- Every record is certified by the local certificate of the existing binary engine: at the certificate degree D the
  closure is closed under the deterministic full sweep and codim V_D = |Z| = evaluation rank of M_<D on Z
  (status `exact-proper-local`), or 1 is in V_D for the systems without solutions (`exact-unit`). The last fall
  degree is the last degree with a strict fall.
- Records: `fetched/<name>.json` (final) and `fetched/<name>.json.D<d>.json` (per degree), `fetched/results.jsonl`,
  `fetched/chain-status.json`, logs `fetched/chain-<name>.log`.

## Scoring against the frozen predictions (**exact** values; the scoring rule is the frozen one)

`python3 score.py` → `outcomes-2026-11-12.json`. It checks the sha256 of the predictions file against the value
frozen in `PLAN.md` (equal), and for every case that the record's `system_hash`, `points_sha256`, n, k and |Z| equal
the predicted case's, that the engine hash equals the repository file, and that the record is certified:
57 of 57 records, 0 problems.

Cell entries "p: h/s" mean: prediction p, h hits among s systems.

| cell | (n, k) | e | k+e | base / curve | class | systems | \|Z\| | certified LFD | cert. degree | B | post-hoc | A | A_Z | inverse |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | (14, 11) | 8 | 19 | invariant+random / koblitz0+koblitz1+ordinary | affine | 7 | 210-335 | 4 (x7) | 4 | 5: 0/7 | 4: 7/7 | 4: 7/7 | 4: 7/7 | - |
| C1 | (14, 11) | 8 | 19 | random / ordinary | order-two | 1 | 268 | 4 (x1) | 4 | - | - | - | - | - |
| C2a | (12, 11) | 11 | 22 | trace-zero / koblitz0 | affine | 2 | 2070 | 5 (x2) | 5 | 5: 2/2 | 5: 2/2 | 5: 2/2 | 5: 2/2 | - |
| C2a | (12, 11) | 10 | 21 | trace-zero / koblitz1 | order-two | 1 | 1035 | 5 (x1) | 6 | - | - | - | - | - |
| C3a | (21, 11) | 1 | 12 | random / koblitz0+ordinary | affine | 3 | 0 | 4 (x3) | 4 | 4: 3/3 | 4: 3/3 | 3: 0/3 | 4: 3/3 | - |
| C3b | (23, 12) | 1 | 13 | random / ordinary | affine | 2 | 0-4 | 4 (x2) | 4 | 4: 2/2 | 4: 2/2 | 3: 0/2 | 4: 2/2 | - |
| C2b | (13, 12) | 12 | 24 | trace-zero / koblitz0 | affine | 2 | 4004 | 6 (x2) | 6 | 5: 0/2 | 5: 0/2 | 5: 0/2 | 5: 0/2 | - |
| C5a | (20, 10) | 10 | 20 | subfield / koblitz0+koblitz1 | order-two | 2 | 1023 | 6 (x2) | 6 | - | - | - | - | 6: 2/2 |
| C6a | (16, 12) | 8 | 20 | random / koblitz0+ordinary | affine | 4 | 238-284 | 4 (x4) | 4 | 5: 0/4 | 4: 4/4 | 4: 4/4 | 4: 4/4 | - |
| C3c | (25, 13) | 1 | 14 | random / ordinary | affine | 2 | 0-2 | 4 (x2) | 4 | 4: 2/2 | 4: 2/2 | 3: 0/2 | 4: 2/2 | - |
| C4a | (15, 11) | 7 | 18 | invariant+random / koblitz0 | affine | 4 | 88-124 | 4 (x4) | 4 | 4: 4/4 | 4: 4/4 | 4: 4/4 | 4: 4/4 | - |
| C2c | (14, 13) | 13 | 26 | trace-zero / koblitz1 | affine | 2 | 8234 | 6 (x2) | 6 | 5: 0/2 | 5: 0/2 | 5: 0/2 | 5: 0/2 | - |
| C4b | (21, 12) | 3-4 | 15-16 | invariant+random / koblitz1 | affine | 4 | 2-20 | 4 (x4) | 4 | 4: 4/4 | 4: 4/4 | 3-4: 1/4 | 4: 4/4 | - |
| C6b | (18, 13) | 8 | 21 | random / ordinary | affine | 3 | 226-286 | 5 (x3) | 5 | 5: 3/3 | 4: 0/3 | 4: 0/3 | 4: 0/3 | - |
| C3d | (27, 14) | 1 | 15 | random / ordinary | affine | 2 | 2-4 | 4 (x2) | 4 | 4: 2/2 | 4: 2/2 | 3: 0/2 | 4: 2/2 | - |
| C5b | (24, 12) | 12 | 24 | subfield / koblitz0 | order-two | 1 | 4095 | 7 (x1) | 7 | - | - | - | - | 7: 1/1 |
| C3e | (29, 15) | 1 | 16 | random / ordinary | affine | 1 | 0 | 4 (x1) | 4 | 4: 1/1 | 4: 1/1 | 3: 0/1 | 4: 1/1 | - |
| C3f | (31, 16) | 1 | 17 | random / ordinary | affine | 1 | 0 | 4 (x1) | 4 | 4: 1/1 | 4: 1/1 | 3: 0/1 | 4: 1/1 | - |
| R1 | (14, 11) | 8 | 19 | random / koblitz0+ordinary | affine | 9 | 234-292 | 4 (x9) | 4 | 5: 0/9 | 4: 9/9 | 4: 9/9 | 4: 9/9 | - |
| R2 | (17, 12) | 7 | 19 | random / ordinary | affine | 4 | 104-140 | 4 (x4) | 4 | 5: 0/4 | 4: 4/4 | 4: 4/4 | 4: 4/4 | - |

| class | systems | B | post-hoc | A | A_Z |
|---|---|---|---|---|---|
| affine, all | 52 | 24 | 45 | 31 | 45 |
| affine, \|Z\| > 0 | 45 | 17 | 38 | 31 | 38 |
| affine, \|Z\| = 0 (units) | 7 | 7 | 7 | 0 | 7 |

Scaled subfield, inverse-family value max(3, k/2 + 1): 3 of 3 (6 at k = 10 on two curves, 7 at k = 12).

## What the predictions decided

- **Law B is refuted at its last threshold.** B predicts 5 for k + e >= 19. All 24 affine systems at k + e = 19 and
  20 (C1, R1, R2, C6a; k = 11 and 12) have LFD 4.
- **The post-hoc law, A and A_Z are refuted at (18, 13).** They predict 4; all 3 systems of C6b (k + e = 21) have
  LFD 5.
- **Every law is refuted on the trace-zero hyperplane from k = 12.** All predict 5; the 4 systems of C2b (13, 12)
  and C2c (14, 13) have LFD 6. At k = 11 (C2a) the value is 5, as predicted.
- **Law A is refuted on small e.** For n = 2k - 1 (e = 1) A predicts 3; all 11 systems of C3a to C3f (k = 11 to 16)
  have LFD 4, as B, post-hoc and A_Z predict.
- **The written inverse-family value holds** on the three scaled-subfield systems: 6 at k = 10, 7 at k = 12. At
  k = 12 the only strict falls are at degrees 3 and 7 (`fetched/C5b_subf_koblitz0_ord_0.json`).
- No frozen law is right on all cells: the best two (post-hoc and A_Z) score 45 of 52 on the affine class.

## Readings (empirical; not laws, fitted after the fact)

- Largest certified values by k on the systems computed here: trace-zero hyperplane 5, 6, 6 at k = 11, 12, 13;
  scaled subfield 6, 7 at k = 10, 12; all other bases 4 for k = 11, 12, 14, 15, 16 and 4 or 5 at k = 13.
- On the non-hyperplane affine systems the value is 4 for k + e = 12 to 20 and 5 at k + e = 21. A threshold at 21
  instead of 19 fits these cells, but (18, 13) is the only cell with k + e >= 21 that is not a hyperplane, so this
  cannot be told apart from an effect of k = 13 at e = 8. It was not predicted.
- At k = 15 and 16 only one system each was run, both with n = 2k - 1, e = 1 and no solution; they say nothing about
  larger e at those depths.

## Limits

- One engine, one run per system, one card. The engine was validated against the existing engine on 65 systems at
  k <= 10 (`validation/`), which is below every size run here; at k >= 11 there is no second implementation and no
  independent checker of the stored certificate. The bases were not kept (hashes `rref_sha256` only).
- Degree cap 5 at k = 15, 16 by plan; both systems certified at degree 4, so the cap did not bind. No system ended
  `capped`.
- The instance draws are the plan's (seed 20261112); cells have 1 to 9 systems. C2b and C2c have two systems each
  on one curve family each.
- Predictions were frozen before any closure at k >= 11; the file hash is checked by `score.py`.

## Machine and cost

- First launch (pod `icb-bin-deep`, `z5wxz84ligol2q`, created 01:45:34Z): one upload got HTTP 502 from the proxy
  and the launch script terminated the pod; nothing ran. `pod/launch.sh` now retries uploads and the start call.
- Second launch (pod `icb-bin-deep2`, `dr8hmd8fjr5msp`, created 01:48:25Z): terminated by the copy loop at
  03:19:06Z after its final copy (`fetched/mirror-icb-bin-deep2.log`); `ls` reports the pod not found.
- Cost computed from these times at $5.29/h: about $8.0 for the second pod and about $0.3 for the first. The account
  balance went from $100.22 (01:44Z) to $89.60 (03:19Z), which also includes pods of other work running in the same
  account, so the billed amount for this run was not read separately.
- The plan estimated 2.1 h nominal for the 16 main cells; the whole chain with the reserve cells took 1.3 h.

## Replication (2026-11-12, plan frozen in `REPLICATION.md` before the run)

Pod `icb-bin-repl` (`33m82rz6vc9i5f`, H200, US-NC-1), created 03:36:59Z, terminated by the copy loop at 05:21:27Z
after its final copy. Wrapper `pod/repl_bin.py`; both on-card self-tests `ok`. Archives `fetched-repl/repl-r2.tgz`,
`repl-r3.tgz`, unpacked to `fetched-repl/r2/`, `r3/`. Comparison: `python3 compare_repl.py --tar … --run r2-gpu=fetched-repl/r2/out --run r3-cpu=fetched-repl/r3/out --out replication-2026-11-12.json`.

| run | device and settings | systems | equal | degree records equal | basis hashes equal on proper subspaces | problems |
|---|---|---|---|---|---|---|
| r2 | GPU, batch 768, tile 512 MB, work 12 GB (first run: 1024, 1024 MB, 20 GB); tile rows and tile counts differ on every large system | 57 | 57 | 307 of 307 | 186 of 186 | 0 |
| r3 | CPU (numpy, OpenBLAS, 16 threads) | 46 | 46 | 232 of 232 | 133 of 133 | 0 |

**Pass** by the frozen criterion: every system of both runs is equal to the first run on every compared field at every
degree, including the canonical basis hash. Every value of the first run, including LFD 6 on the trace-zero
hyperplanes at k = 12, 13, LFD 5 at (18, 13) and LFD 7 on the subfield base at k = 12, is reproduced by a second GPU
execution with different batching and tiling; the 46 smaller systems are also reproduced on the CPU device. The 11
larger systems were not run on the CPU (by plan).

What this does not change: both executions run the same engine source, so an error in the algorithm that does not
depend on device, batch or tiling would replicate; the guard against that is still the validation against the
existing engine at k <= 10. There is still no second implementation at k >= 11.

Cost from the pod's lifetime at $5.29/h: 1 h 44 min, about $9.2.
