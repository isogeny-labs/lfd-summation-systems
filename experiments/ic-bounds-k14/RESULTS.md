# Results on the rented GPU (unit instance over p = 65537)

Instance: y^2 = x^3 + x over F_65537, target R = (4271, 11720) (x(R) a non-residue, so the solution set is
empty at every depth k <= 15; proof and enumeration in README.md). Same (a, b, R) at every depth.
NVIDIA H200 (141 GB), CuPy 14.2.0. Engine `csd/ic_bounds_big.py`: sha256 ba99c529… for every row except
k = 14, D = 6, which used e4f15a08… (host spill, `--dev-gb 105`, memory pool disabled); both versions passed
`selftest --device gpu` on the card (`runs/selftest-gpu-h200*.json`) and reproduce the reference hashes.
Label of every row: **exact** (closure run to the end, deterministic full sweep clean, 0 repairs).

| D | k | columns | dim V_D | profile | 1 in V_D | closed | seconds | multiply-adds |
|---|---|---|---|---|---|---|---|---|
| 5 | 10 | 21700 | 18569 | (0, 0, 2, 138, 3163, 18569) | no | yes | 121 | 5.75e+11 |
| 5 | 11 | 35443 | 25033 | (0, 0, 2, 142, 3441, 25033) | no | yes | 151 | 2.10e+12 |
| 5 | 12 | 55455 | 32057 | (0, 0, 2, 146, 3727, 32057) | no | yes | 208 | 5.53e+12 |
| 5 | 13 | 83682 | 39657 | (0, 0, 2, 150, 4021, 39657) | no | yes | 242 | 1.22e+13 |
| 5 | 14 | 122438 | 47849 | (0, 0, 2, 154, 4323, 47849) | no | yes | 142 | 2.79e+13 |
| 6 | 10 | 60460 | 60460 | (1, 21, 211, 1351, 6196, 21700, 60460) | yes | yes | 73 | 9.34e+12 |
| 6 | 11 | 110056 | 110056 | (1, 23, 254, 1794, 9109, 35443, 110056) | yes | yes | 197 | 4.68e+13 |
| 6 | 12 | 190051 | 190051 | (1, 25, 301, 2325, 12951, 55455, 190051) | yes | yes | 414 | 2.04e+14 |
| 6 | 13 | 313912 | 313912 | (1, 27, 352, 2952, 17902, 83682, 313912) | yes | yes | 1291 | 7.83e+14 |
| 6 | 14 | 499178 | 448922 | (0, 0, 2, 154, 4451, 73587, 448922) | no | yes | 8077 | 4.69e+15 |

## What these rows establish (for this instance)

- **exact** 1 not in V_5(k) for k = 10..14, with closed V_5 of dimension 18,569 / 25,033 / 32,057 / 39,657 /
  47,849: the dimensions and profiles of the earlier instances over p = 12289 and p = 40961
  (`experiments/ic-bounds-d5/`), now from a third engine. The k = 14 value equals the prediction frozen in
  `experiments/ic-bounds-freeext/predictions-2026-11-10.json` (`d5-14`: 47,849, profile (0,0,2,154,4323,47849)).
- **exact** V_6(k) = M_<=6 and 1 in V_6(k) for k = 10, 11, 12, 13 (prediction `d6-full-through-13`): the last
  fall degree of this unit-ideal system is exactly 6 at k = 10, 11, 12, 13.
- **exact** V_6(14) is closed, has dimension 448,922 < 499,178 (codimension 50,256) and does not contain 1
  (prediction `d6-falls-at-14`): the last fall degree at k = 14 is at least 7. With the kernel
  depth-restriction theorem, 1 is not in V_6(k) for every k >= 14 on this instance (given 2^k | p - 1, i.e.
  k <= 16 over this prime; the formal statement is in `DepthRestrictionM2.not_one_of_le`).
- Codimension of V_D in M_<=D at the first deficient depth 4D - 10: 10, 189, 3131, 50,256 for D = 3, 4, 5, 6.
- The low profile of V_6(14) is (0, 0, 2, 154, 4451, 73587): it agrees with V_5(14) through weight 3 and is
  larger at weights 4 and 5 (4451 against 4323; 73,587 against 47,849).

## Limits

- One instance (j = 1728 curve, unit ideal), one engine at k = 14, one run. The engine was validated by
  canonical-basis hash against the reference engine up to 21,700 columns locally; on the card it reproduced
  the known dimensions at k = 10..13 (D = 5) and the known k = 10, D = 6 outcome (full rank, parts
  30,318 + 30,142, as in `experiments/ic-bounds-d5/verify/k10-p12289-s10-D6.json`; the hashes are not
  comparable across primes).
- The k = 14 basis (about 45 GB) was not kept; the record holds its hashes. A replication on a proper
  instance (random curve) is recorded below.
- The exact value at k = 14 was later obtained by the degree lift (section "Degree lift at k = 14" below): 7, in one run
  without an independent certificate. A direct closure in M_<=7 (1,683,218 columns) is out of reach here.

## Replication on a proper instance (random curve)

Instance: y^2 = x^3 + 4021 x + 9917 over F_65537, seed 12, target drawn until the solution set is nonempty;
|Z| = 1028 at k = 13 and 4068 at k = 14; evaluation at Z is onto, so dim(I cap M_<=D) = columns - |Z|.
Engine e4f15a08…, same card. Every row **exact** (closed, clean full sweep, 0 repairs; random elements of the
final space vanish on Z).

| D | k | columns | dim(I cap M_<=D) | dim V_D | V_D = ideal part | profile | seconds |
|---|---|---|---|---|---|---|---|
| 5 | 13 | 83682 | 82654 | 39657 | no | (0, 0, 2, 150, 4021, 39657) | 147 |
| 5 | 14 | 122438 | 118370 | 47849 | no | (0, 0, 2, 154, 4323, 47849) | 226 |
| 6 | 13 | 313912 | 312884 | 312884 | yes | (0, 0, 2, 1924, 16874, 82654, 312884) | 2157 |
| 6 | 14 | 499178 | 495110 | 448922 | no | (0, 0, 2, 154, 4451, 73587, 448922) | 7430 |

- At k = 13: V_5 is short of the ideal part and V_6 equals it, so the last fall degree is exactly 6
  (local-to-global criterion; the evaluation of M_<6 on Z has rank |Z|).
- At k = 14: V_6 is closed with dimension 448,922, the same dimension, profile and part ranks as on the unit
  instance, and 46,188 short of the ideal part: last fall degree at least 7.
- The dimensions of V_5 and V_6 agree between the j = 1728 unit instance and this curve at both depths.

## Later runs on the second pod (2026-11-10, 14:26 UTC onwards)

Engine e4f15a08…, driver 6c8bb10e…, same card type (H200). Every closure row is **exact** (closed, clean full sweep,
0 repairs). Records in `runs/`.

### Cross-checks on the unit instance over p = 65537

| check | record | outcome |
|---|---|---|
| Reference engine `csd/ic_bounds_d5_verify.py`, k = 14, D = 5 (CPU, 16,333 s) | `reference-engine-unit-p65537-k14-D5.json` | rank 47,849, same profile, both basis hashes equal to the GPU engine's (cfbdb672…, fceae22f…) |
| Reference engine, k = 10, D = 6 (CPU, 3,129 s) | `reference-engine-unit-p65537-k10-D6.json` | full rank 60,460, 1 in V_6, hashes equal (the basis of a full space is the identity, so this confirms full rank) |
| Second run of k = 14, D = 6, basis kept (6,573 s) | `rerun-unit-p65537-k14-D6.json` | rank 448,922, same profile and part ranks, both basis hashes equal to the first run (eeac52ef…, a58c6c7d…) |
| Independent checker on the stored basis of that run (731 s) | `certcheck-unit-p65537-k14-D6.json` | format and echelon form ok; g in W; all 28 shifts of all 73,587 low rows in both parts reduce to zero (4,120,872 checks, 0 failures); 1 not in W; rank, profile, hashes equal the engine record |

### More depths for D = 5

| instance | k | columns | dim V_5 | profile | 1 in V_5 | seconds |
|---|---|---|---|---|---|---|
| unit, p = 65537 | 15 | 174437 | 56649 | (0, 0, 2, 158, 4633, 56649) | no | 193 |
| random curve, p = 65537 (seed 12) | 15 | 174437 | 56649 | (0, 0, 2, 158, 4633, 56649) | no | 332 |
| random curve, p = 65537 (seed 12) | 16 | 242825 | 66073 | (0, 0, 2, 162, 4951, 66073) | no | 994 |

These are the values of the closed form frozen as prediction `d5-closed-form` (56,649 and 66,073). The profile
recursion p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k) therefore holds for k = 10..15.

### Two more primes

Unit instance y^2 = x^3 + 3x over p = 163841 = 5 * 2^15 + 1, target (73730, 28291), seed 21. X = 73730 is a non-residue,
X^3 + 3X = 28291^2, and the subgroup of order 2^k lies in the squares for k <= 14, so the argument of the paper's unit
lemma applies; the enumeration gives 0 points at k = 13 and 14. Random curve y^2 = x^3 + 31415x + 92653 over
p = 114689 = 7 * 2^14 + 1, target (80359, 998), seed 22: 594 solutions at k = 13 and 2,370 at k = 14.

| instance | D | k | dim V_D | profile | reading | seconds |
|---|---|---|---|---|---|---|
| unit, p = 163841 | 5 | 13 | 39657 | (0, 0, 2, 150, 4021, 39657) | 1 not in V_5 | 90 |
| unit, p = 163841 | 6 | 13 | 313912 | (1, 27, 352, 2952, 17902, 83682, 313912) | full, 1 in V_6: **LFD = 6** | 1533 |
| unit, p = 163841 | 5 | 14 | 47849 | (0, 0, 2, 154, 4323, 47849) | 1 not in V_5 | 131 |
| unit, p = 163841 | 6 | 14 | 448922 | (0, 0, 2, 154, 4451, 73587, 448922) | closed, 1 not in V_6: **LFD >= 7** | 6393 |
| random curve, p = 114689 | 5 | 13 | 39657 | (0, 0, 2, 150, 4021, 39657) | short of the ideal part (83,088) | 95 |
| random curve, p = 114689 | 6 | 13 | 313318 | (0, 0, 2, 2358, 17308, 83088, 313318) | 313,912 - 594, equals the ideal part: **LFD = 6** | 1639 |
| random curve, p = 114689 | 5 | 14 | 47849 | (0, 0, 2, 154, 4323, 47849) | short of the ideal part | 153 |

Part ranks at k = 14, D = 6 over p = 163841: 224,688 + 224,234, as over p = 65537.

### Degree lift at k = 14 (D = 6 -> 7, unit instance over p = 65537)

Method, proof and validation: `LIFT.md`. Record: `lift-unit-p65537-k14-D6to7.json` (751 s; lift 93e8d19b…).
Input: the checked basis of the second k = 14, D = 6 run (dim 448,922; 25,136 + 25,120 free columns).

| quantity | + part | - part |
|---|---|---|
| regular weight-7 columns | 586,501 | 586,501 |
| non-regular weight-7 columns | 5,519 | 5,519 |
| dim T | 30,655 | 30,639 |
| dim K | 30,655 | 30,639 |
| falls (rounds 1, 2, 3) | 494, 12,260, 12,382 | 494, 12,260, 12,366 |

- Fixed point after 3 rounds; every regular weight-7 column covered exactly once; `guard_hits` 0.
- dim V_7(14) = 448,922 + 2 * 586,501 + 30,655 + 30,639 = 1,683,218 = |M_<=7|, and 1 is in V_7(14).
- With 1 not in V_6(14): **the last fall degree is exactly 7 at k = 14 on this instance.**
- On the card the lift also reproduced V_6(10) = M_<=6 from V_5(10) (`lift-unit-p65537-k10-D5to6-gpu.json`, 8.4 s).
- Control added on 2026-11-11 (CPU, same lift and engine files): D = 5 -> 6 at k = 10 on the random curve over
  p = 65537 with 14 solutions. Lift: dim V_6 = 60,446 = 60,460 - 14, 1 not in V_6, dim K = 1,958 / 1,947 against
  dim T = 1,965 / 1,954, 2 rounds, 37 s; direct closure: the same dimension, profile (0, 7, 197, 1337, 6182, 21686, 60446)
  and both low-basis hashes, 761 s (`validation/lift-k10d5-proper.json`).
- Label: exact computation by a method with a written proof; **no independent certificate** (no sweep in M_<=7, one
  run, one implementation, one instance).

### Lost with the pod, and rerun

The second pod was deleted by the provider at about 04:30 UTC on 2026-11-11, when the account balance ran out. Two
jobs were on it and nothing of them was fetched: the soundness control (k = 14, D = 6 on the random curve over
p = 114689, then the lift) and the reference engine at k = 11, D = 6. Both were rerun the same day; see the next
section.

## Relaunch (2026-11-11, 06:15 UTC onwards; `relaunch/`, plan in `RELAUNCH.md`)

Same engine, lift, driver and checker files as the recorded runs (hashes verified at upload). H200 pod `icb-k14b`
with the CuPy pools disabled through `relaunch/nopool_run.py`; both self-tests passed on the card
(`runs/selftest-gpu-h200-relaunch.json`, `runs/lift-selftest-gpu-h200-relaunch.json`).

### Soundness control and exact value on the random curve over p = 114689 (k = 14, 2,370 solutions)

| step | record | outcome |
|---|---|---|
| closure D = 6, basis kept (11,678 s) | `proper-p114689-k14-D6.json` | **exact**: closed, clean full sweep, 0 repairs; dim 448,922, profile (0, 0, 2, 154, 4451, 73587, 448922), part ranks 224,688 + 224,234; 1 not in V_6; 47,886 short of the ideal part (496,808): **LFD >= 7** |
| lift D = 6 -> 7 (4,679 s) | `lift-proper-p114689-k14-D6to7.json` | fixed point after 3 rounds; falls 494 + 494, 12,260 + 12,260, 11,197 + 11,181; dim K = 29,470 / 29,454 against dim T = 30,655 / 30,639 (codimension 1,185 per part, 2,370 in all); every regular weight-7 column covered once; `guard_hits` 0 |

- dim V_7(14) = 1,680,848 = 1,683,218 - 2,370: **exactly the ideal bound**. A lift that produced anything outside the
  ideal would have exceeded it. 1 is not in V_7.
- dim (V_7 cap M_<=6) = 496,808 = 499,178 - 2,370, so the evaluation of M_<7 on Z has rank 2,370, V_7 is the degree-7
  part of the ideal, and by the local-to-global criterion **LFD = 7 exactly on this curve**.
- Label: as for the unit instance (written proof, no independent certificate), now on two instances.
- The independent checker on the kept basis of this instance (`certcheck-proper-p114689-k14-D6.json`, 650 s on the
  card): g in W, W closed under all 28 shifts of its 73,587 low rows (4,120,872 triples, 0 failures), 1 not in W,
  dimension 448,922 and both basis hashes equal to the engine's record. So the exclusion of 1 from V_6(14) is now
  checked independently on two instances over two primes (65537 unit, 114689 random curve).

### Reference engine at k = 11, D = 6 (unit instance over p = 65537, CPU pod, 6,909 s, 23.4 GB)

`reference-engine-unit-p65537-k11-D6.json`: rank 110,056 (full), 1 in V_6, full sweep of 35,443 low rows clean, profile,
part ranks (55,144 + 54,912) and both basis hashes (b9b3695c…, 77fef68d…) equal to the GPU engine's
`unit-p65537-k11-D6.json`. Driver: `relaunch/ref_k11.py`, checked against `scripts/ic_bounds_d5_verify.py` at
k = 7, D = 6 (all fields equal).

## Job B (2026-11-11/12): second unit instance, y^2 = x^3 + 3x over p = 163841, k = 14

Pod `icb-k14e` (H200, US-NC-1), chain `relaunch/chain_u163841.py`, device memory pool disabled through
`relaunch/nopool_dev_only.py`; both self-tests passed on the card (`runs/selftest-gpu-h200-jobB.json`,
`runs/lift-selftest-gpu-h200-jobB.json`). Raw copies in `relaunch/fetched-e/`.

| step | record | outcome |
|---|---|---|
| closure D = 6, basis kept (6634 s) | `rerun-unit-p163841-k14-D6.json` | **exact**: closed, clean full sweep (147174 row checks), 0 repairs; dim 448922, profile (0, 0, 2, 154, 4451, 73587, 448922); 1 not in V_6; both basis hashes (e0c4afd7…, 44a89923…) equal to the first run of this instance (`unit-p163841-k14-D6.json`, batch 2048; this run batch 1024) |
| lift D = 6 -> 7 (677 s) | `lift-unit-p163841-k14-D6to7.json` | fixed point after 3 rounds; falls per round [[494, 494], [12260, 12260], [12382, 12366]]; dim K = [30655, 30639] = dim T; every regular weight-7 column covered once; `guard_hits` 0; `contains_one: true`; dim V_7 = 1683218 = |M_<=7| |
| independent checker (528 s) | `certcheck-unit-p163841-k14-D6.json` | g in W, W closed under all 28 shifts (4120872 triples, 0 failures), 1 not in W; dimension and both hashes equal to the engine's record |

- **LFD = 7 exactly on this instance**: 1 not in V_6(14) (two runs with equal hashes, independent check), 1 in V_7(14)
  (lift). The lift's counts (dim T, dim K, falls per round) are the same as on the unit instance over p = 65537.
- The exact value 7 at k = 14 now holds on three instances over three primes: unit over 65537, unit over 163841,
  random curve over 114689. Label unchanged: the lift has a written proof and passed its ideal-bound controls,
  but there is no independent certificate of V_7.

Two things went wrong on the way and are recorded in `RELAUNCH.md`: a first pod was abandoned after nine minutes
(self-test slow because the wrapper `nopool_run.py` also disabled the pinned-memory pool, which the recorded runs
of 2026-11-10 did not), and a first attempt of the closure stopped at about 4,100 s (inferred: disk full, two
30-minute checkpoints near peak size on a 150 GB disk; its log was not copied).

## Machine and cost

One NVIDIA H200 pod (Runpod, $5.29/h), created 06:41 UTC and terminated 12:42 UTC on 2026-11-10 (about 6.0 h,
about $32), plus two smoke-test pods (a few cents). All three pods were terminated and confirmed gone.

A second H200 pod (`icb-more1`, $5.29/h) ran from 14:26 UTC on 2026-11-10 until the provider deleted it at about
04:30 UTC on 2026-11-11 (about 14 h, about $74). It produced the records of the section "Later runs on the second pod".

Relaunch on 2026-11-11: H200 pod `icb-k14b` ($5.29/h) from 06:15 UTC, and an 8-vCPU pod `icb-ref11` ($0.44/h) from
06:14 to 08:16 UTC (about $0.90). The GPU pod was billed $27.28 (06:15 to about
11:24 UTC) and was terminated after the checker. No pod of this work remains (checked 2026-11-11, 19:30 UTC).
Billed in total for the three GPU pods: about $134 ($32 + $74.96 + $27.28), plus under $2 for the CPU and
smoke-test pods.

Job B (2026-11-11, 20:13 UTC to 2026-11-12, 00:27 UTC): three H200 pods, $0.80 (abandoned after the slow
self-test), $7.02 (closure attempt that stopped) and $7.33 billed so far for the successful run of 2 h 21 min
(about $12.5 expected when the last hour is billed). All terminated; no pod of this work exists (checked 00:27 UTC).
