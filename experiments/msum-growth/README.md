# Growth of the last fall degree with the depth, three and four summands

Every number in this file is `[pcc:unverified]`: one engine, one host, no pcc claim, no
S3 record, no independent replay. Preregistration (committed before any closure of this
study): [`PREREGISTRATION.md`](PREREGISTRATION.md).

## Question

How does the last fall degree (LFD) of the Petit–Kosters–Messeng two-power subgroup
summation system grow with the depth `k` when there are `m >= 3` summands? On a unit
ideal the LFD is the Caminata–Gorla solving degree for every degree-compatible order
(paper, Section `sec:solving`). The paper proves `LFD >= 7` for `m = 3` (`thm:carry`) and
`LFD <= mk + 1` (`cor:allm`), and had exact values `7` for `m = 3` at `k <= 5`, one lower
bound `LFD >= 8` at `k = 6` (one proper instance), and `13` for `m = 4` at `k = 4`.

## Method

* **System.** The direct normalized chain presentation
  `{NF(S_{m+1}(x_1..x_m, t))} ∪ G_{m,k}` over `F_p`, `H = mu_{2^k}` (the system of the
  paper's `m >= 3` theorems; why not the split or System (4) presentations: see the
  preregistration). Closures `V_D` are computed modulo the chains (chain bridge,
  `thm:bridge`).
* **Generator.** `S_{m+1}` evaluated on `H^m` through Semaev's resultant recursion with
  numbers (`S_4 = Res_X(S_3(x_1,x_2,X), S_3(x_3,t,X))`, `S_5 = Res_X(S_3(x_1,x_2,X),
  S_4(x_3,x_4,t,X))`, formal-degree Sylvester determinants), then the inverse DFT over
  `(Z/2^k)^m`. Checks run on every instance: `S_{m+1}` vanishes at `(x(P_1),..,x(P_m),
  x(P_1+..+P_m))` for 8 random point tuples and not at a perturbed last coordinate; the
  normal form evaluates to the computed values on 64 random grid points and on up to 64
  zeros; for `m = 3`, `k >= 3` the coefficient of `x_1^3x_2^3x_3^3` equals
  `-8(5t^3+3at+6b)` exactly (ratio 1 on all nine `m = 3` instances built), as in the proof of
  `thm:carry`; `deg g = 6` (`m = 3`), `8, 12, 12, 12` (`m = 4`, `k = 2..5`).
* **Engine.** The paper's second engine, `second-engine/engine` (`hky_fp_tower_fast2`,
  "hky3", files unchanged), called by a new Rust driver [`rust/`](rust/) on a ring of
  `mk` Kummer chain variables (`y_{i,j}^2 = y_{i,j+1}`, `y_{i,k-1}^2 = 1`), with the swap
  of the first two blocks as symmetry (verified by the engine; `split` in every record)
  and the engine's closedness check `verify_closed` on every closure that ends without
  `1` (every generator and every variable times every final basis row of degree `< D`
  reduces to zero; `verify_nonzero = 0` in every record). Driver binary SHA-256
  `f8ba1a3ede8964b51d5a407b902b34db0a78e11eabf646a47b63e8d341e3920d` (rustc 1.97.0,
  release), the same for every record. One closure per process ([`run.sh`](run.sh),
  address space capped at 14 GB), so `maxrss_mb` is the peak of that closure.
* **Reading.** Unit ideal (`|Z| = 0`, counted on all of `H^m`): `LFD` = least `D` with
  `1 ∈ V_D`, with every smaller degree certified without `1` (`V_D = 0` for `D < deg g`,
  `V_{deg g} = span(g)`, larger `D` by a verified closed closure). Proper ideal: the
  local-to-global criterion (`rank ev(M_{<D}) = |Z|`, `dim V_D = |M_{<=D}| - |Z|`), the
  value is then the last strict fall.
* **Host.** 4 vCPU (AVX-512), 15 GB, shared at times with another session's job (load
  average up to 6.7 during the first `k = 6, D = 8` closure); wall times are indicative
  only. A container restart killed two queued closures (`records/*.log` says which); they were rerun unchanged.

### Calibration on recorded cells (before any new cell)

The recorded proper instance `p7681-r` of `experiments/ic-bounds-m3-depth`
(`p = 7681`, `a = 1876`, `b = 2678`, `t = 6223`):

| k | D | recorded `|Z|`, profile | this engine | match |
|---|---|---|---|:-:|
| 2, 3, 4 | 7 | `|Z| = 0`, LFD 7 | `|Z| = 0`, `1 ∈ V_7` | yes |
| 5 | 7 | `12`; `0 6 109 564 1929 4932 9937 16372` | identical, closed, verified | yes |
| 6 | 7 | `48`; `0 0 0 0 0 3 2375 19958`, closed, no `1` | identical, closed, verified, 176 s (Python engine: 1381 s) | yes |

For `m = 4` the recorded `GP-m4k4-*` instances could not be rebuilt (their curve
parameters are not in the committed spec files); the recorded unit instance
`P-m4k3-p193-r` (`a = b = 167`, `t = 116`) gives `deg g = 12`, 4084 generator terms and
`|Z| = 0`, as recorded.

## Instances

Seeded (`--seed 0`, rule of the preregistration; seed string `msum-<p>-0`):

| name | p | a | b | t | `|Z|` (m = 3) at k = 2..7 |
|---|---|---|---|---|---|
| U1 | 2013265921 | 137549545 | 257513864 | 171427433 | 0 0 0 0 0 0 |
| U2 | 998244353 | 538530846 | 574242499 | 522220447 | 0 0 0 0 0 0 |
| U3 | 469762049 | 64397718 | 173339033 | 415975955 | 0 0 0 0 0 0 |
| P1 | 7681 | 1876 | 2678 | 6223 (`p7681-r`) | 0 0 0 12 48 – |
| (65537, seed 0) | 65537 | 54337 | 33232 | 33146 | 0 0 0 0 6 36 |

U1 and U2 are also unit for `m = 4` at `k = 2..5`. All targets have
`5t^3 + 3at + 6b != 0`. Instance data: [`records/instances-m3.jsonl`](records/instances-m3.jsonl),
[`records/instances-m4.jsonl`](records/instances-m4.jsonl).

## Results

### Three summands: LFD against the depth (unit instances U1, U2, U3)

| k | `3k` vars | `deg g` | `V_6` | `V_7` (dim, profile `p_0..p_7`) | `V_8` | LFD | instances | status |
|---|---|---|---|---|---|---|---|---|
| 2 | 6 | 6 | `span(g)` | contains 1 | – | **7** | U1 U2 U3 | exact (1 in `V_7`, `V_6 = span(g)` verified closed) |
| 3 | 9 | 6 | `span(g)` | contains 1 | – | **7** | U1 U2 U3 | exact |
| 4 | 12 | 6 | `span(g)` | contains 1 | – | **7** | U1 U2 U3 | exact |
| 5 | 15 | 6 | `span(g)` | contains 1 | – | **7** | U1 U2 U3 | exact |
| 6 | 18 | 6 | `span(g)` | 19,958; `0 0 0 0 0 3 2375 19958`; closed, verified, no 1 | contains 1 (U1, U2, U3) | **8** | U1 U2 U3 | exact on all three |
| 7 | 21 | 6 | `span(g)` | 27,092; `0 0 0 0 0 3 2384 27092`; closed, verified, no 1 | not computed (memory) | **`>= 8`** | U1 | lower bound certified |

Proper instances (reported separately):

| instance | k | `|Z|` | `V_7` | `V_8` | LFD | status |
|---|---|---|---|---|---|---|
| P1 `p7681-r` | 5 | 12 | `= I ∩ M_{<=7}` (recorded profile reproduced) | – | 7 | exact (recorded; reproduced) |
| P1 `p7681-r` | 6 | 48 | 19,958, closed, verified, deficient (`dim(I∩M_{<=7}) = 62,956`) | 106,714 `= 106,762 - 48`, closed, verified; `rank ev(M_{<8}) = 48 = |Z|`; profile `0 0 124 940 4000 12568 31132 62956 106714` | **8** | exact (local-to-global at 8; the strict fall at 8 is the last) |
| 65537 seed 0 | 6 | 6 | 19,958, closed, verified, deficient (`dim(I∩M_{<=7}) = 62,998`) | 106,756 `= 106,762 - 6`, closed, verified; `rank ev(M_{<8}) = 6 = |Z|`; profile `0 14 166 982 4042 12610 31174 62998 106756` | **8** | exact (local-to-global at 8; last strict fall at 8) |

Observations (exact on these instances):

* The closures do not depend on the instance at all in these cells: at `k = 6` the
  degree-7 closure has the same profile `0 0 0 0 0 3 2375 19958` on all five instances
  (three unit, two proper with `|Z| = 48` and `6`), and the degree-8 closures of U1, U2 and U3
  reach `1` after the same number of rows (`182,972` rows reduced; part ranks
  `47,560 + 47,387 = 94,947` of `106,762` columns when `1` enters) with the same peak
  storage. The same holds for `m = 4`, `k = 4` (part ranks `34,726 + 30,672` on both
  unit instances). At `k <= 5` the evaluation of `M_{<=5}` already has full rank on
  `Z`, so proper and unit instances differ only in `|Z|`.
* The step from 7 to 8 happens at `k = 6`, the same depth as the two-summand step from 4
  to 5 (`floor((k+14)/4)`). At `k = 2..6` the three-summand value is the two-summand
  value plus three.
* Depth restriction: the substitution `u_{i,k} -> 1` of `thm:restrict` (kernel-checked
  for the two-summand system) maps the depth-`(k+1)` system to the depth-`k` system for
  `m = 3` as well (the normal form of `S_4` has exponents `<= 4 < 2^k` for `k >= 3`, so it
  is fixed), so the written argument gives `1 ∉ V_7` at every `k >= 6` on U1, U2, U3 (and
  `LFD` nondecreasing in `k` while the ideal stays unit). This is a written extension,
  not a kernel statement; the `k = 7` closure checks it directly at one depth.

### Four summands (unit instances U1, U2)

| k | `4k` vars | `deg g` | LFD | columns `|M_{<=LFD}|` | status |
|---|---|---|---|---|---|
| 2 | 8 | 8 | 9 | 256 | paper (`tab:msum`, forced: `deg g = mk`) |
| 3 | 12 | 12 | 13 | 4,096 | paper; recorded `P-m4k3-p193-r` reproduced (1 in `V_13`) |
| 4 | 16 | 12 | **13** | 65,399 | exact here on U1 and U2 (`V_12 = span(g)`, 1 in `V_13`); paper: 13 on six proper instances |
| 5 | 20 | 12 | **`>= 13`** | 988,116 at `D = 13` | lower bound only (`V_12 = span(g)` misses 1); `V_13` not computed (about `10^6` columns) |

### Macaulay sizes against `m = 2`

`|M_{<=D}|` is the number of squarefree monomials of degree `<= D` in the `mk` chain
variables (the columns of the degree-`D` closure); `C(mk, D)` is its top layer.

| k | m = 2: LFD, `|M_{<=LFD}|`, `C(2k,LFD)` | m = 3: LFD, `|M_{<=LFD}|`, `C(3k,LFD)` | m = 4: LFD, `|M_{<=LFD}|`, `C(4k,LFD)` |
|---|---|---|---|
| 2 | 4, 16, 1 | 7, 64, 0 | 9, 256, 0 |
| 3 | 4, 57, 15 | 7, 502, 36 | 13, 4,096, 0 |
| 4 | 4, 163, 70 | 7, 3,302, 792 | 13, 65,399, 560 |
| 5 | 4, 386, 210 | 7, 16,384, 6,435 | `>= 13`, `>=` 988,116, 77,520 |
| 6 | 5, 1,586, 792 | 8, 106,762, 43,758 | – |
| 7 | 5, 3,473, 2,002 | `>= 8`, `>=` 401,930, 203,490 | – |
| 8 | 5, 6,885, 4,368 | (P4: 8) 1,271,626, 735,471 | – |

`m = 2` values are the paper's (`floor((k+14)/4)`, exact on its instances through
`k = 13`). At equal `k` the three-summand Macaulay space is 4 to 67 times larger
(`k = 2..6`) than the two-summand one, and the gap widens with `k` because the degree
is higher in more variables.

### Timings and memory (this host, 4 vCPU; seconds include the closedness check)

| cell | columns | rank at end | closure s | of which verify s | peak RSS | `F_p` multiply-adds |
|---|---|---|---|---|---|---|
| m=3 k=5 D=7, unit | 16,384 | 16,384 (1 in) | 12–15 | – | 0.2 GB | 6.2e10 |
| m=3 k=6 D=7, unit (x3) | 63,004 | 19,958 | 312–329 | 217–228 | 2.0 GB | 1.6e11 |
| m=3 k=6 D=8, unit (U1, U2, U3) | 106,762 | 94,947 when 1 enters | 1,425 / 1,257 / 1,277 | – | 5.2 GB | 5.44e12 |
| m=3 k=6 D=8, proper P1 | 106,762 | 106,714 | 2,176 | 6 | 5.1 GB | 9.10e12 |
| m=3 k=6 D=8, proper 65537 s0 | 106,762 | 106,756 | 1,667 | 7 | 5.0 GB | 7.41e12 |
| m=3 k=7 D=7, unit U1 | 198,440 | 27,092 | 205 | 128 | 2.8 GB | 1.67e11 |
| m=3 k=7 D=8, unit U1 | 401,930 | – | killed at 633 s | – | 13.5 GB and rising | – |
| m=4 k=4 D=13, unit (U1, U2) | 65,399 | 65,399 (1 in) | 905 / 913 | – | 2.3–2.4 GB | 5.27e12 |

The first `k = 6, D = 8` closure shared the host with another session's job (load
average up to 6.7); the other numbers were taken at load average about 4.

## Prediction outcomes

| | prediction | outcome |
|---|---|---|
| P1 | `LFD = 7` at `k = 2..5`, unit | **held** on U1, U2, U3 (12 systems) |
| P2 | `1 ∉ V_7` at `k = 6`, unit | **held** on U1, U2, U3 (closed, verified) |
| P3 | `LFD = 8` at `k = 6`, unit | **held** on U1, U2, U3 |
| P4 | `LFD(3,k) = floor((k+14)/4) + 3`; testable: `LFD(3,7) = 8` | **not decided**: `LFD(3,7) >= 8` certified (consistent with P4 and with the linear alternative `k + 2`); `V_8` at `k = 7` needs about 73 GB (below). The `k = 2..6` values agree with P4 and are also consistent with `k + 2` only from `k = 5` on (`k + 2 = 7, 8` at `k = 5, 6`, but 4, 5, 6 at `k = 2..4` against measured 7), so a plain linear law `k + 2` is already ruled out at `k <= 4`; the alternative still open is a slope-one growth starting at `k = 5`. |
| P5 | proper = unit value at equal `k` | **held** for P1 at `k = 5` (7) and `k = 6` (8); 65537 seed 0 at `k = 6` (8): **held** |
| P6 | `m = 4`: 13 at `k = 4` and `k = 5` | `k = 4`: **held** (U1, U2). `k = 5`: **not tested** (`>= 13` only) |

Nothing was falsified. The new facts: the three-summand LFD grows with the depth
(`7 -> 8` between `k = 5` and `k = 6`, exact on three unit and two proper instances, a
lower bound on all five), at the same depth as the two-summand step, and stays at least
`8` at `k = 7`. Two data points (`k = 6, 7`) cannot fix a rate; they are consistent with
one step per four depths (P4) and with the slower end of anything faster.

## What deciding the slope would cost (estimate, no spend)

Model from the measured closures: a closure that runs until `1` enters near full rank
costs about `c n^3` multiply-adds with `n` the column count (`c = 4.5e-3` from
`k = 6, D = 8`; `7.5e-3` for a proper closure that runs to the end), and its peak memory
is about `0.45 n^2` bytes (`5.2 GB` at `n = 106,762`; the killed `k = 7` run had reached
13.5 GB after 10 minutes). Local throughput was about `4e9` multiply-adds per second on
4 vCPU; a 32-vCPU host is assumed to give 6 to 8 times that (`2.5–3.3e10`/s), shared
hosts less. Low-rank closures that end without 1 (all `D = 7` cells) are far cheaper
than this model (`k = 7, D = 7`: 1.7e11 multiply-adds, 2.8 GB).

| cell | columns | multiply-adds | memory | 32 vCPU / 256 GB | decides |
|---|---|---|---|---|---|
| m=3 k=7 D=8 | 401,930 | 2.9e14 (to 4.9e14 if it runs to the end) | about 73 GB | 2.5–5.5 h | P4 (`1 ∈ V_8`: LFD(3,7) = 8) against slope-one growth from `k = 5` (`1 ∉ V_8`) |
| m=3 k=7 D=9 (only if `1 ∉ V_8`) | 695,860 | 1.5–2.5e15 | about 220 GB | 13–28 h, at the memory limit | LFD(3,7) = 9 |
| m=3 k=8 D=8 | 1,271,626 | 9e15 | about 730 GB | does not fit | whether 8 persists to `k = 8` |
| m=3 k=10 D=9 (P4's second step) | 22,964,087 | – | – | out of reach of this engine | the second step |
| m=4 k=5 D=13 | 988,116 | 4e15–2e16 | about 440 GB | does not fit | P6 at `k = 5` |

So one pod job, `m = 3`, `k = 7`, `D = 8` on U1 (about 3 to 6 hours of a 32-vCPU,
256 GB host, roughly US$6–12 at the US$2.08/h of the 2026-11-17 pods), decides between
P4 and faster growth; with the instance independence observed above, one instance is
likely to speak for all, but the reading rule asks for two. Beyond `k = 7` the dense
closure does not fit 256 GB; a second step (P4: `k = 10`) is not reachable with this
engine, and the slope beyond one step needs either a structural argument (the closures
are instance independent, so the obstruction is combinatorial) or a different
instrument.

### The pod runs (2026-11-20): undecided

Both attempts at `m = 3`, `k = 7`, `D = 8` on U1 ran on one 32-vCPU, 256 GB RunPod CPU
host (pod `yn1hjuubx6pstv`, US$1.76/h, same binary, SHA-256 `f8ba1a3e…3920d`), and
neither produced a closure.

* Run 1 (16:20 UTC) stopped at the driver's default `max_entries` of 1.2e9 and was
  recorded as censored (`records/pod-m3k7/`). The job script then ran `D = 9`
  unconditionally, which was also censored; the host idled for about two hours.
* Run 2 (18:12 UTC, `--max-entries 1.8e10`, `--budget 64800`) was still in the closure
  at 20:19 UTC (resident set 74 GB, only the instance header written: `deg g = 6`, 125
  terms, `|Z| = 0`, the leading coefficient check at ratio 1). The pod was terminated
  by the provider at about 22:35 UTC; every other pod on the account stopped in the
  same hour, consistent with the account balance running out (not confirmed). The
  output lived on the pod's container disk, so the run left no record. It had run
  about 4 h 20 min, inside the 2.5–5.5 h estimate above.

`LFD(3,7) >= 8` therefore stands as certified and `LFD(3,7) = 8` (P4) against `>= 9`
remains open. A rerun should write its output to a network volume or copy it off the
host periodically, and the job script must stop on a censored record.

## Files

* [`PREREGISTRATION.md`](PREREGISTRATION.md): instances, quantity, predictions, written
  before the runs (commit `cc7bffb`).
* [`rust/`](rust/): the driver (`msum_direct`), on the unchanged engine
  `second-engine/engine`. Build: `CARGO_TARGET_DIR=$PWD/target cargo build --release
  --manifest-path experiments/msum-growth/rust/Cargo.toml` from the repository root.
* [`run.sh`](run.sh): one closure per process under a 14 GB address-space limit.
* `records/*.jsonl`: one JSON record per run (instance, generator checks, `|Z|`,
  evaluation ranks, closure profile, `contains_one`, closedness check, timings, peak
  RSS, multiply-adds, binary SHA-256); `records/*.log`: the run logs, including the
  killed `k = 7, D = 8` run and the runs lost to a container restart and rerun.
* [`lfd-vs-depth.svg`](lfd-vs-depth.svg): LFD against `k` for `m = 2` (paper), `m = 3`
  and `m = 4` (this study), with P4 and the linear alternative.

Reproduce one cell (k = 6, D = 8, U1; about 25 minutes, 5.2 GB):

```sh
experiments/msum-growth/run.sh out.jsonl out.log --m 3 --k 6 --p 2013265921 --seed 0 --D 8
```
