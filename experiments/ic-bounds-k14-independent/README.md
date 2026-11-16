# Independent degree lift V_D -> V_{D+1}: implementation and local validation

Purpose: a second implementation of the lift of `../ic-bounds-k14/LIFT.md`, for the paper's open item
"the value at k = 14 by a second method (... or an independent lift)". **Nothing here has run at
k = 14 or on a GPU.** The k = 14 run is planned and costed in `PLAN.md`.

Labels: **measured** = a command run here (CPU, numpy 2.5.3, 2 BLAS threads), with the file it left;
**read** = taken from a named file; **proof** = argument written in the module docstring.

| file | role |
|---|---|
| `csd/ic_bounds_lift_independent.py` | the lift (numpy only; `xp` = numpy or cupy) |
| `scripts/run_ic_lift_independent.py` | driver: `lift`, `plan`, `arith` |
| `validate.py` | local validation harness (this directory) |
| `validation.json` | its record: 37 cases, 138 lift runs, 69 negative controls |
| `chain_independent.py` | pod-side chain for `PLAN.md` (dry run and CPU run of its small stage only) |
| `local/` | records of that CPU run (k = 10, D = 5 -> 6, the repository's instance files) |
| `PLAN.md` | plan, memory, time and cost for k = 14 |

## 1. What the implementation does

Input: the engine's checkpoint of V_D (both sigma-parts) and the instance (p, a, b, X, k). Output:
whether 1 is in V_{D+1}, dim(V_{D+1} n M_<=D) in total and per sigma-part, dim V_{D+1}.

Algorithm and proof: module docstring. In short, it is the fixed point of Lemma 1 of `LIFT.md`
(W_{i+1} = (W_i + sum_u u W_i) n L) carried out in **ordinary monomial coordinates**:

* the two stored parts are read as vectors l_x = e_x + phi(x) for every monomial x whose orbit is a
  pivot of both parts (set Pd), plus "single" rows supported on the other monomials N;
* red : span(M_<=D+1) -> F^Reg (+) F^Sp (+) F^N replaces e_x by -phi(x) (kernel: span of the l_x);
  a top is special if it is u n for some n in N; a regular top t has exactly D + 1 generators
  red(u l_{t/u}), each with coefficient 1 at t, so their differences span the part without regular
  coordinate (in orbit coordinates the first lift has coefficients 1, -1, 2 and (D+1)/2 generators on
  the diagonal; none of these cases exists here);
* one echelon form K on T = Sp (+) N; its rows with pivot in N are the new elements of V_{D+1} n L,
  and each is multiplied by every R_u (e_n -> red(e_{un})) until nothing is pending; at the end all
  of them are shifted once more and must reduce to zero (explicit fixed-point check);
* the regular tops are taken in chunks of a fixed pseudo-random order (not by the first lift's
  J-groups), each chunk against a frozen snapshot of K, with one reference vector per top.

Two properties used by the plan (**proof**, docstring "Soundness and completeness"):

1. every vector added to K is red of an element of V_{D+1}, so `contains_one: true` is proved as soon
   as the column of 1 becomes a pivot, whatever part of the schedule has run;
2. if K = T (`saturated: true`) then V_{D+1} = span(M_<=D+1), again by soundness alone.

Both need the stored space to be contained in V_D. That is an assumption on the input (section 6).

Checks made before any lift step; a failure raises `Refuse` and the driver exits with code 4:
manifest file sizes and sha256; entries are integers in [0, p); the engine's canonical hash of both
parts equals the closure record (`--closure-record`), which must say `closed: true`; the generator
rebuilt from (p, a, b, X, k) equals the instance's and lies in the stored space; rows with pivot of
weight < D have no support in weight D; u times 2 random combinations of those rows lies in the stored
space for all 2k shifts (this is exactly the statement that the generators skipped by Lemma 2 of
`LIFT.md` are redundant; a violation survives with probability at most 1/p^2 per shift);
N is swap-stable; every regular top has exactly D + 1 generators.

## 2. What is shared with the first implementation

Imports: none from the repository (`grep -n import csd/ic_bounds_lift_independent.py`: copy, hashlib,
json, os, time, numpy). The driver imports only the new module.

I did **not open** `csd/ic_bounds_lift.py`, `scripts/ic_bounds_lift_run.py` or
`scripts/ic_bounds_lift_validate.py` (only their sizes, with `wc`). The task allowed reading them for
the file format; the format was taken from the engine instead.

Read, and therefore possibly shared in shape:

* `../ic-bounds-k14/LIFT.md`, all of it: statement, proof and the algorithm of its section 3. The
  split into regular and special tops, the elimination of regular coordinates by differences, the
  reduction of the maps R_u modulo K ("followers") and "shift each fall row once" are the same
  mathematics, re-derived here in ordinary coordinates. A flaw in the theorem itself would be common
  to both; the comparison with direct closures (section 3) is the test of that.
* `csd/ic_bounds_big.py`: `Layout.__init__`, the head of `Part`, `Part.compact`, `Part.sha`,
  `_state_json`, `_config_json`, `save_checkpoint`, `find_checkpoint`, `load_checkpoint` (file format,
  column order, hash definition). `scripts/ic_bounds_big_run.py`: argument list and `cmd_closure`.
* `csd/ic_bounds_d5_verify.py` (the reference engine), all of it. It is the oracle of the validation,
  not a lift. My `Echelon`/`rref` are written from scratch but have the textbook block structure that
  the reference also has.
* `../ic-bounds-k14/CERTCHECK.md` sections 1 to 3 and 8 (lift/project convention), the checker's
  argument list, `README.md` section 3 of that directory, `RESULTS.md`, `POD.md`, `RELAUNCH.md`,
  `relaunch/chain_k14.py`, `chain_u163841.py`, `launch.sh`, `nopool_dev_only.py`.

Shared by necessity (specification): the checkpoint format; the column order (weight descending,
then a1 | a2 << k) and the orbit convention of the two parts, used only in the reader and in the
final split of the result into sigma-parts; the definition of the canonical hash; the formula of
S_3 (built here by polynomial arithmetic and compared with the reference `chain_generator`).

## 3. Validation (measured)

```
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 experiments/ic-bounds-k14-independent/validate.py --sets enter,stay,proper,k9,k10d5 \
  --work <scratch>/work --out experiments/ic-bounds-k14-independent/validation.json     # 1,717 s
```

Per case: the engine (subprocess, `scripts/ic_bounds_big_run.py closure --ckpt`) writes V_D; the
independent lift runs on that checkpoint in up to four settings; the results are compared with

* **oracle A**: the direct closure in degree D + 1 by the reference engine `csd/ic_bounds_d5_verify.py`
  (`Closure(..., split=True)`, full sweep clean in every case): dim V_{D+1}, 1 in or not,
  dim(V_{D+1} n M_<=D) per part. For the two k = 10, D = 5 -> 6 cases the direct closure is too large
  to recompute here and is **read** from records: unit, `runs/reference-engine-unit-p65537-k10-D6.json`
  and `runs/unit-p65537-k10-D6.json` (two engines, equal); proper,
  `validation/lift-k10d5-proper.json`, field `direct` (big engine on CPU; total only, no per-part
  profile, so for that case the per-part numbers [10,871, 10,815] are the lift's own, equal in both
  settings and summing to the recorded 21,686);
* **oracle B** (cases with at most 1,700 columns in M_<=D+1): a naive closure written in the harness
  in ordinary coordinates, from g alone, with no file and no sigma-split;
* **the first lift's recorded table** (`LIFT.md` section 5, parsed from the file, 37 rows found):
  dim V_D, dim V_{D+1}, 1 in, profile[D], and dim T, dim K summed over the parts (equal to
  `n_T - n_single` and `rank_K - n_single` here).

Settings: `plain` (defaults), `tiny` (batches of 7 rows, first chunk of 5 tops, reference table
budget 20 KB so that there are many chunks, products split at inner dimension 5), `low`
(`include_low`: no use of Lemma 2, basis held as an `xp` array), `seed9` (another order of the tops,
first chunk 64). `tiny` only up to 1,500 columns of M_<=D; above 3,000 columns only `plain` and `low`.
Also per case: the canonical hashes recomputed by the reader equal the engine's record, and phi(x)
computed by two routes (row gathers; products with whole parts) agrees on 40 random x.

| set | case | p | solutions | dim V_D | oracle A: dim V_{D+1}, 1 in, dim(V_{D+1} n M_<=D) per part | independent lift (all settings) | settings | K = T | oracle B | first lift (LIFT.md) | lift s / oracle s (`plain`) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| enter | unit-65537-k3-D3to4 | 65537 | 0 | 7 | 57, yes, [23, 19] | 57, yes, [23, 19]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.0 / 0.0 |
| enter | unit-65537-k4-D3to4 | 65537 | 0 | 9 | 163, yes, [49, 44] | 163, yes, [49, 44]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.0 / 0.0 |
| enter | unit-65537-k5-D3to4 | 65537 | 0 | 11 | 386, yes, [91, 85] | 386, yes, [91, 85]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.1 / 0.0 |
| enter | unit-65537-k6-D4to5 | 65537 | 0 | 605 | 1586, yes, [408, 386] | 1586, yes, [408, 386]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.1 / 0.6 |
| enter | unit-65537-k7-D4to5 | 65537 | 0 | 829 | 3473, yes, [750, 721] | 3473, yes, [750, 721]: equal | plain, tiny, low, seed9 | yes | - | equal | 2.1 / 3.0 |
| enter | unit-65537-k8-D4to5 | 65537 | 0 | 1061 | 6885, yes, [1277, 1240] | 6885, yes, [1277, 1240]: equal | plain, low, seed9 | yes | - | equal | 19.5 / 14.2 |
| enter | unit-12289-k4-D3to4 | 12289 | 0 | 9 | 163, yes, [49, 44] | 163, yes, [49, 44]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.0 / 0.0 |
| enter | unit-12289-k5-D3to4 | 12289 | 0 | 11 | 386, yes, [91, 85] | 386, yes, [91, 85]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.1 / 0.0 |
| enter | unit-12289-k6-D4to5 | 12289 | 0 | 605 | 1586, yes, [408, 386] | 1586, yes, [408, 386]: equal | plain, tiny, low, seed9 | yes | equal | equal | 0.1 / 0.6 |
| enter | unit-12289-k7-D4to5 | 12289 | 0 | 829 | 3473, yes, [750, 721] | 3473, yes, [750, 721]: equal | plain, tiny, low, seed9 | yes | - | - | 2.0 / 3.0 |
| enter | unit-163841-k5-D3to4 | 163841 | 0 | 11 | 386, yes, [91, 85] | 386, yes, [91, 85]: equal | plain, tiny, low, seed9 | yes | equal | - | 0.1 / 0.0 |
| enter | unit-163841-k6-D4to5 | 163841 | 0 | 605 | 1586, yes, [408, 386] | 1586, yes, [408, 386]: equal | plain, tiny, low, seed9 | yes | equal | - | 0.1 / 0.6 |
| enter | unit-163841-k7-D4to5 | 163841 | 0 | 829 | 3473, yes, [750, 721] | 3473, yes, [750, 721]: equal | plain, tiny, low, seed9 | yes | - | - | 2.1 / 3.0 |
| stay | unit-65537-k2-D2to3 | 65537 | 0 | 1 | 5, no, [1, 0] | 5, no, [1, 0]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.0 / 0.0 |
| stay | unit-65537-k3-D2to3 | 65537 | 0 | 1 | 7, no, [1, 0] | 7, no, [1, 0]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.0 / 0.0 |
| stay | unit-65537-k4-D2to3 | 65537 | 0 | 1 | 9, no, [1, 0] | 9, no, [1, 0]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.0 / 0.0 |
| stay | unit-65537-k6-D3to4 | 65537 | 0 | 13 | 605, no, [58, 53] | 605, no, [58, 53]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.2 / 0.1 |
| stay | unit-65537-k7-D3to4 | 65537 | 0 | 15 | 829, no, [60, 55] | 829, no, [60, 55]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.9 / 0.4 |
| stay | unit-65537-k8-D3to4 | 65537 | 0 | 17 | 1061, no, [62, 57] | 1061, no, [62, 57]: equal | plain, tiny, low, seed9 | no | - | equal | 3.0 / 0.9 |
| stay | unit-12289-k3-D2to3 | 12289 | 0 | 1 | 7, no, [1, 0] | 7, no, [1, 0]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.0 / 0.0 |
| stay | unit-12289-k6-D3to4 | 12289 | 0 | 13 | 605, no, [58, 53] | 605, no, [58, 53]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.2 / 0.1 |
| stay | unit-12289-k7-D3to4 | 12289 | 0 | 15 | 829, no, [60, 55] | 829, no, [60, 55]: equal | plain, tiny, low, seed9 | no | equal | - | 0.9 / 0.4 |
| stay | unit-163841-k6-D3to4 | 163841 | 0 | 13 | 605, no, [58, 53] | 605, no, [58, 53]: equal | plain, tiny, low, seed9 | no | equal | - | 0.2 / 0.1 |
| stay | unit-163841-k7-D3to4 | 163841 | 0 | 15 | 829, no, [60, 55] | 829, no, [60, 55]: equal | plain, tiny, low, seed9 | no | equal | - | 0.9 / 0.4 |
| proper | proper-65537-X2-k6-D3to4 | 65537 | 2 | 13 | 605, no, [58, 53] | 605, no, [58, 53]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.3 / 0.1 |
| proper | proper-65537-X2-k6-D4to5 | 65537 | 2 | 605 | 1584, no, [407, 385] | 1584, no, [407, 385]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.1 / 0.6 |
| proper | proper-65537-X2-k7-D3to4 | 65537 | 2 | 15 | 829, no, [60, 55] | 829, no, [60, 55]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.9 / 0.4 |
| proper | proper-65537-X2-k7-D4to5 | 65537 | 2 | 829 | 3471, no, [749, 720] | 3471, no, [749, 720]: equal | plain, tiny, low, seed9 | no | - | equal | 1.8 / 3.0 |
| proper | proper-65537-X2-k8-D3to4 | 65537 | 4 | 17 | 1061, no, [62, 57] | 1061, no, [62, 57]: equal | plain, tiny, low, seed9 | no | - | equal | 3.1 / 0.9 |
| proper | proper-65537-X2-k8-D4to5 | 65537 | 4 | 1061 | 6881, no, [1275, 1238] | 6881, no, [1275, 1238]: equal | plain, low, seed9 | no | - | equal | 19.4 / 14.7 |
| proper | proper-15361-k6-D4to5 | 15361 | 2 | 605 | 1584, no, [407, 385] | 1584, no, [407, 385]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.1 / 0.6 |
| proper | proper-15361-k7-D4to5 | 15361 | 4 | 829 | 3469, no, [748, 719] | 3469, no, [748, 719]: equal | plain, tiny, low, seed9 | no | - | equal | 1.9 / 3.0 |
| proper | proper-15361-k7-D3to4 | 15361 | 4 | 15 | 829, no, [60, 55] | 829, no, [60, 55]: equal | plain, tiny, low, seed9 | no | equal | equal | 0.8 / 0.4 |
| k9 | unit-65537-k9-D4to5 | 65537 | 0 | 1301 | 12616, yes, [2047, 2001] | 12616, yes, [2047, 2001]: equal | plain, low | yes | - | equal | 150.0 / 56.4 |
| k9 | proper-65537-X1424-k9-D4to5 | 65537 | 2 | 1301 | 12614, no, [2046, 2000] | 12614, no, [2046, 2000]: equal | plain, low | no | - | equal | 146.7 / 58.3 |
| k10d5 | unit-65537-k10-D5to6 | 65537 | 0 | 18569 | 60460, yes, [10878, 10822] (read from records) | 60460, yes, [10878, 10822]: equal | plain, low | yes | - | equal | 51.7 / - |
| k10d5 | proper-65537-X1424-k10-D5to6 | 65537 | 14 | 18569 | 60446, no, total 21686 (read from a record) | 60446, no, [10871, 10815]: equal | plain, low | no | - | - | 151.1 / - |

Totals (**measured**, `validation.json`, `totals`): 37 cases, 138 lift runs, **0 mismatches** with
oracle A; 15 cases where 1 enters at D + 1 and 22 where it stays out (11 unit, 11 proper); primes
12289, 15361, 65537, 163841; oracle B agrees in 24 of 24 cases; the first lift's recorded numbers are
equal in 29 of 29 rows that the two tables share. Every run ended `fixed_point: true` and
`sigma_stable: true`; the non-saturated runs passed the final fixed-point check (`verify_new_rows: 0`).

"lift s / oracle s" compares the lift with the direct closure of the *reference* engine on the same
machine; the lift is not faster at these sizes and is not meant to be.

### Negative controls (**measured**, `validation.json`, `negative`; 4 cases, 69 controls)

Cases: unit-65537 k = 7 D = 4 -> 5 (1 enters), unit-65537 k = 7 D = 3 -> 4 (stays), proper-65537-X2
k = 7 D = 4 -> 5, unit-163841 k = 6 D = 4 -> 5.

| control | runs | outcome |
|---|---|---|
| generator: coefficient of (0,0) or (1,1) + 1 | 8 | 8 refused: "the generator g is not in the stored space" |
| one entry + 1 in a stored row of pivot weight < D (in memory, file hashes bypassed) | 21 | 21 refused: 18 "not closed" (12 to 28 of the 24 or 28 shifted random elements leave the space), 3 "generator not in the stored space" |
| one entry + 1 in a stored row of pivot weight D (same) | 24 | 24 refused: "not closed" (4 to 16 shifted elements leave); the canonical hash also differs from the engine record in all 45 corrupted bases |
| one byte of a tile file flipped on disk | 4 | 4 refused by the reader: sha256 differs from the manifest |
| test hook in the harness: fall rows never shifted | 4 | 4 results changed: `contains_one` flips to false in the 2 cases where 1 enters (dimension 2,844 instead of 3,473; 1,523 instead of 1,586), dimension falls in the other 2 |
| test hook: direct generators of one shift dropped | 4 | 0 changed |
| test hook: all direct generators dropped | 4 | 1 changed (unit k = 7, D = 3 -> 4: 810 instead of 829) |

62 of 69 detected. The 7 undetected ones are the two "dropped generators" hooks in cases where the
dropped generators are redundant, so the correct result is unchanged; they show that the comparison
with the oracle does **not** test the completeness of the direct generators in most cases (the first
lift's controls found the same: `LIFT.md`, "the dropped rows were redundant"). A wrongly *added*
vector would be seen on the proper instances, where the result must stay inside the ideal: all 11
proper cases give the oracle's dimension exactly.

The driver's refusals were also exercised from the command line (exit code 4 for a changed generator
in the instance file and for a closure record with another hash; exit code 3 with `stopped: deadline`
for `--max-hours`; `--stop-on-one` stops with `contains_one: true`, `fixed_point: false`): run by
hand, not recorded in a file.

### Chain test (**measured**, `local/`)

`chain_independent.py --stage small --device cpu` in a flat copy of the files, with the repository's
instance files `unit-p65537-a1-s10-k10.json` and `proper-p65537-s12-k10.json`: arithmetic test,
engine closure D = 5 (106 s and 101 s), driver `lift` with `--closure-record`: both equal the expected
values; `local/chain-independent-status.json`. Records: unit 54 s, 863 MB peak, K = T after one chunk
(2,048 of 37,972 regular tops) and 2 fall rounds, 2.33e11 multiply-adds in products; proper 146 s,
894 MB, all 308,120 generator rows, 5 chunks, rank K = 3,921 = 3,935 - 14, final check 0 new rows.
The checker's command line of the chain was run by hand on the unit checkpoint (`all_true`, 23 s).
Scratch paths in the records were replaced by `SCRATCH` after the runs.

## 4. Not verified

* **GPU**: cupy is not installed here; no line has run under cupy. The code path is the same (`xp`),
  and `PLAN.md` puts an arithmetic test and the k = 10 cases on the card before the target.
* **k = 14**: memory and time are estimates (`PLAN.md` sections 4 and 5).
* The regime of the target (few special tops: V_D almost everything) is covered locally by the two
  k = 10, D = 5 -> 6 cases only (788 special against 37,972 regular tops), and there the oracle is
  read from records. All other cases have many special tops.
* A k = 11, D = 5 -> 6 run was started and stopped: it is not in the target's regime (22,020 special
  tops, T has 32,451 columns) and needs more memory than the limit of this session.
* p = 163841 only up to k = 7; p = 16777213 (product splitting forced by the prime) only in the
  arithmetic test, not in a lift. Splitting inside a lift is exercised by the `tiny` setting (inner
  dimension 5).
* `--basis-on device` with a real device; `--no-pool`.

## 5. Cost of the implementation (measured locally)

k = 10, D = 5 -> 6, unit: 2.33e11 multiply-adds in products, 2.7e10 in follower reductions, 1.7e10 in
K (`local/indep-lift-unit-...json`, `counts`); the first lift's record for the same case has 2.05e11
in all (`LIFT.md` section 6). On CPU most of the time is elementwise (`modp`) and row gathers, not
products (profile run by hand: 31 of 67 s in `modp`).

## 6. Remarks on the existing material

1. **The upper bound at k = 14 depends on a containment that no check covers.** Both lifts compute
   the degree-(D+1) closure of the *stored* space W. The independent checker of `CERTCHECK.md` proves
   W contains V_6 (its section 8 says it does not show W is the least such space, "not needed" -
   true for 1 not in V_6). For 1 in V_7 the other inclusion, W inside V_6, is what is needed, and it
   rests on the engine alone. `LIFT.md` section 10 states the hypothesis ("V_D must be the closed
   degree-D closure"); the paper's limitation paragraph on the value 7 lists the lift's proof,
   assertions and comparisons but not this dependence. A corrupted single entry is caught (controls
   above: a changed top row breaks closedness, because shifted low rows are reduced with the top
   rows), but a consistently larger closed space would not be.
2. `LIFT.md` section 7: "the split 703 + 702 of the 1,405 free weight-6 columns is a guess". The
   closure record gives 709 + 696 (`runs/rerun-unit-p65537-k14-D6.json`, `part_profiles`: 187,843 and
   187,492 weight-6 pivots of 188,552 and 188,188 columns). Only the sum enters the estimate.
3. `LIFT.md` section 7 lists 187,849 regular weight-6 orbits (an estimate written before the run);
   the lift records have 187,842. Not an error, but the table is not marked as superseded.
4. The numbers of the first lift's validation table that this work could compare (29 rows) are all
   reproduced. I found no error in the statement or proof of `LIFT.md` sections 1 and 2.
