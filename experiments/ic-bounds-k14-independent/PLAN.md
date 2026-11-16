# Plan: the value at k = 14 by an independent lift (NOT RUN)

Status: **plan only.** Nothing below has been run on a GPU or at k = 14. No pod was created and no
money was spent in preparing it. Labels: **read** = taken from a named file of this repository;
**measured** = a command I ran locally (CPU, numpy); **estimate** = derived, not measured;
**unverified** = an assumption I could not check. Implementation and local validation: `README.md`.

## 1. Question and predicted outcome

At k = 14 the paper has 1 not in V_6 (independent checker) and 1 in V_7 only by the first lift
(`csd/ic_bounds_lift.py`). This run repeats the step V_6 -> V_7 with the second implementation
(`csd/ic_bounds_lift_independent.py`, ordinary monomial coordinates, no shared code).

Instances (read: `../ic-bounds-k14/instances/`, `RESULTS.md`):

| tag | instance file | curve, target | role |
|---|---|---|---|
| `unit-p65537` | `unit-p65537-a1-s10-k14.json` | y^2 = x^3 + x over F_65537, X = 4271, 0 points | the paper's main unit instance |
| `unit-p163841` | `unit-p163841-a3-s21-k14.json` | y^2 = x^3 + 3x over F_163841, X = 73730, 0 points | second unit instance |
| `proper-p114689` (optional) | `proper-p114689-s22-k14.json` | y^2 = x^3 + 31415x + 92653 over F_114689, X = 80359, 2,370 points | soundness control of the new implementation at the target size |

**Prediction (fixed here, before any run).** For both unit instances the independent lift ends with
`fixed_point: true`, `saturated: true`, `contains_one: true`, and

| field of the record | predicted | source of the number |
|---|---|---|
| `dim_VD` | 448,922 | read: `runs/rerun-unit-*-k14-D6.json` (`rank`) |
| `dim_low` = dim(V_7 n M_<=6) | 499,178 | read: `runs/lift-unit-*-k14-D6to7.json` (`profile[6]`) |
| `dim_low_parts` | [249,824, 249,354] | read: same records, `low_part_ranks` |
| `dim_VD1` = dim V_7 | 1,683,218 | read: same records, `dim_VD1` |
| `n_special_tops` | 11,038 = 5,519 + 5,519 | read: `special_top_orbits` (weight 7 is odd, so no top is self-paired and each orbit column is two monomials in the + and - part together) |
| `n_regular_tops` | 1,173,002 = 586,501 + 586,501 | read: `regular_top_orbits` |
| `n_N_weight_D` | 1,406 | read: `free_type_monomials_weight_D` |
| `n_T - n_single`, `rank_K - n_single` | 61,294 = 30,655 + 30,639 | read: `dim_T`, `dim_K` |
| `falls_new` | 50,256 = 25,136 + 25,120 | read: `falls` |

`n_single` (stored rows whose orbit is a pivot of exactly one part) is **not predicted**: the records do
not contain it. From the per-part profiles of `rerun-unit-p65537-k14-D6.json` (weight-6 pivots 187,843
and 187,492 of 188,552 and 188,188 columns, so 709 + 696 free weight-6 columns) and the 1,406 free-type
weight-6 monomials, exactly one single row has weight 6 (derivation: 14 self-paired and 696 paired
free-type orbits; all 696 are free in the - part, 695 in the + part); the number at lower weight is
unknown. It was 4 at k = 6, D = 4 and 16 at k = 10, D = 5 (measured, `validation.json`).

For `proper-p114689` (optional): `contains_one: false`, `saturated: false`, `schedule_complete: true`,
`verify_new_rows: 0`, `dim_VD1` = 1,680,848 = 1,683,218 - 2,370, `dim_low` = 496,808,
`dim_low_parts` = [248,639, 248,169] (read: `runs/lift-proper-p114689-k14-D6to7.json`).

## 2. Reading rule

Per instance, decided from `out/indep-lift-<tag>-k14-D6to7.json` and the chain status:

* **CONFIRM**: preconditions hold (closure record `closed: true`, `rank` 448,922, `contains_one: false`;
  independent checker `all_true`; the lift's `basis.canonical_sha256` equal to the closure record's
  `basis_sha256`, which the driver enforces) and the record has `status: done`, `fixed_point: true`,
  `sigma_stable: true`, `closed_check.nonzero: 0` and every predicted field of section 1 equal.
  For a unit instance this says: 1 in V_7(14) and V_7(14) = M_<=7, by a second implementation.
* **REFUTE**: `fixed_point: true` and `contains_one` or any predicted dimension differs. Then one of
  the two lifts (or their common input) is wrong, and the paper's upper bound LFD <= 7 at k = 14 has
  no support until the difference is explained. Do not rerun with other settings before recording it.
* **UNDECIDED**: `status: refused` (reason in the record), exit by `--max-hours` without
  `contains_one: true`, a crash, a failed small check on the device, or a closure whose hashes differ
  from the recorded ones. One partial case is still evidence: a run stopped by the deadline with
  `contains_one: true` proves 1 in V_7 (every vector the lift adds lies in V_7; `README.md`,
  "Soundness and completeness"), without the dimensions.

What a CONFIRM does **not** show (section 8): that the stored degree-6 space is contained in V_6.
Both lifts start from the same kind of input, produced by the same engine.

## 3. Steps on the pod

One H200 pod per instance (141 GB card; the launcher line of `../ic-bounds-k14/relaunch/launch.sh`,
with the disk of the successful run of 2026-11-11: read `RELAUNCH.md`, "Second attempt"):

```
T=scripts/ic_bounds_big_pod.py
python3 $T prices                                   # check the H200 price before creating (unverified today)
python3 $T create --name icb-indep1 --gpu "NVIDIA H200" --min-ram-per-gpu 200 --min-vcpu-per-gpu 16 \
  --pip "numpy threadpoolctl cupy-cuda12x[ctk]" --volume-gb 0 --disk-gb 200 --max-seconds 36000 --yes-create
for f in csd/ic_bounds_big.py scripts/ic_bounds_big_run.py scripts/ic_bounds_cert_check.py \
         csd/ic_bounds_lift_independent.py scripts/run_ic_lift_independent.py \
         experiments/ic-bounds-k14-independent/chain_independent.py \
         csd/ic_bounds_lift.py scripts/ic_bounds_lift_run.py \
         experiments/ic-bounds-k14/instances/unit-p65537-a1-s10-k10.json \
         experiments/ic-bounds-k14/instances/proper-p65537-s12-k10.json \
         experiments/ic-bounds-k14/instances/unit-p65537-a1-s10-k14.json; do python3 $T push icb-indep1 "$f"; done
python3 $T push icb-indep1 experiments/ic-bounds-k14/relaunch/nopool_dev_only.py nopool_run.py   # device pool off, pinned pool on
python3 $T run icb-indep1 chain -- python3 chain_independent.py --stage all --tags unit-p65537 --first-lift
```

Start the copy loop (`../ic-bounds-k14/relaunch/mirror.sh`) **before** the chain, with
`out/chain-independent-status.json`, every `out/*.json` named below and the `logs/chain-*.log` files in
its list: two earlier results were lost with their pods (read: `RELAUNCH.md`). The exact syntax of
`push SRC DEST-NAME` and of `mirror.sh` should be checked against those two files before use
(**unverified**: I read `launch.sh`, I did not run the launcher).

`chain_independent.py` (in this directory; dry run and CPU run of its `small` stage: `README.md`)
runs, and stops at the first failure:

1. `run_ic_lift_independent.py arith --device gpu`: `modp` and the chunked product against Python
   integers at the bounds, for p = 12289, 65537, 163841, 16777213.
2. **Small check on the card** (k = 10, D = 5 -> 6, unit and proper instance over p = 65537): engine
   closure of V_5 (121 s on the card, read `RESULTS.md`), then the independent lift three times
   (device; device with `--basis-on device`, another seed and chunk size; CPU). Each must give
   60,460 / 21,700 / [10,878, 10,822] / 1 in, resp. 60,446 / 21,686 / [10,871, 10,815] / 1 not in
   (measured locally, `validation.json`, set `k10d5`; oracles there). This is the only test of the
   cupy code path before the target; nothing of the module has run under cupy (**unverified**).
3. **Closure V_6(14), basis kept**: the recorded command of the second k = 14 run,
   `ic_bounds_big_run.py closure --instance <inst> --D 6 --device gpu --dev-gb 105 --work-gb 20
   --tile-mb 1024 --batch 1024 --threads 8 --ckpt /workspace/ckpt/<tag>-k14-D6 --ckpt-minutes 100000
   --resume --status ... --out /workspace/out/<tag>-k14-D6.json`, through `nopool_run.py`
   (read: `relaunch/chain_k14.py`, `chain_u163841.py`). One checkpoint, written at the end of the
   closure phase: two timed checkpoints near the peak filled a 150 GB disk on 2026-11-11.
   Recorded: 6,573 s (p = 65537) and 6,634 s (p = 163841); host peak 82.5 GB (`maxrss_mb`), device
   budget 105 GB; checkpoint about 45 GB (224,688 x 25,136 + 224,234 x 25,120 float32 = 42.0 GiB).
   The chain requires `closed: true`, rank 448,922, 1 not in, and records whether both basis hashes
   equal the recorded ones (65537: eeac52ef..., a58c6c7d...; 163841: e0c4afd7..., 44a89923...). The
   hash depends on the subspace only, so equality identifies the space with the one the first lift
   and the checker used.
4. **Independent checker** on that checkpoint: `ic_bounds_cert_check.py --device gpu --shifts all
   --instance <inst> --ckpt <ckpt> --D 6 --engine-record <closure record> --work-gb 45 --row-chunk 4096`
   (read: `chain_k14.py`). Recorded 731 s and 528 s. Must end `all_true`.
5. **Independent lift**: `run_ic_lift_independent.py lift --ckpt <ckpt> --instance <inst>
   --closure-record <closure record> --device gpu --no-pool --max-hours 6 --table-gb 16 --batch 2048
   --verbose --out /workspace/out/indep-lift-<tag>-k14-D6to7.json`. The driver verifies every file
   hash of the manifest, recomputes the engine's canonical hash of both parts and refuses unless they
   equal the closure record, rebuilds g from (p, a, b, X, k) and refuses unless it equals the
   instance's generator and lies in the stored space, and tests closedness on random combinations.
6. Optional (`--first-lift`, about 12 min): the first lift on the same checkpoint, so that both lifts
   are on record for one set of files.
7. The chain never deletes the checkpoint. Fetch all records, then decide about the basis (section 6)
   and terminate the pod; confirm with `ls`.

The second unit instance: the same with `--tags unit-p163841` and its instance file, on a new pod,
only after the first reading is recorded.

## 4. Memory at k = 14 (estimate from recorded counts)

Counts in ordinary coordinates: |M_<=6| = 499,178 and 1,184,040 tops (C(28, w) monomials of weight w);
N has 50,256 + n_single monomials, T has 61,294 + n_single columns (section 1); 375,334 = 376,740 - 1,406
weight-6 sources; 10,509,352 = 375,334 x 28 generators, of which 8,211,014 = 1,173,002 x 7 have a regular
top and 2,298,338 are direct; 73,587 low rows are not shifted (their generators are tested by the
closedness check instead).

`python3 scripts/run_ic_lift_independent.py plan --n-free 25136,25120 --ranks 224688,224234
--n-special 11038 --regular 1173002 --n-single 100 --batch 2048 --table-gb 16` (measured output,
n_single = 100 assumed):

| item (float64 on the device unless said) | GB |
|---|---|
| one follower R_u modulo K, 50,356 x 61,394 at most | 23.0 |
| K, worst case over its rank (r = f = n_T / 2), four copies during an update | 28.1 |
| fall rows of the snapshot, worst case | 5.8 |
| reference table of a chunk (`--table-gb`) | 16.0 |
| batch buffers (2,048 rows) | 6.8 |
| echelon forms for the sigma-parts at the end (the follower is gone by then) | 14.2 |
| **device, sum of the upper bounds** | **93.9** of 141 |
| host: stored basis as float32 (42.0) + tables | about 45 |

With `--basis-on device` the basis moves to the card (sum of bounds 135.9 GB): too tight for a first
run. The first lift's host peak on the same input was 44.4 GB (read: `maxrss_mb`). The estimate is a
formula; the first log lines of the run (`plan` in the `tables` event) print it for the true counts.
If the card runs out of memory: rerun step 5 with `--table-gb 4 --batch 1024` (the result does not
depend on these; tested locally with other settings, `README.md`).

## 5. Time at k = 14 (estimate; the largest uncertainty of this plan)

Unit instance. The run stops as soon as K = T. Locally at k = 10, D = 5 -> 6 that happened after the
first chunk (2,048 of 37,972 regular tops, 5.4 %) and its fall rounds, with 2.3e11 multiply-adds in the
products (measured, `local/indep-lift-unit-p65537-a1-s10-k10-D5to6-cpu.json`, `counts`); the first lift
needed 2.05e11 in all there and 7.15e14 in 751 s at
k = 14, where it reached K = T after 5.8 % of its schedule (read: `LIFT.md` section 6, lift record).

* Products of the generator chunks: rows x 50,356 x f, f <= 61,394. If one chunk of 2,048 tops is
  enough: about 18,000 rows, 6e13. If 5.8 % to 11 % of the tops are needed (five or six doubling
  chunks): 0.57 to 1.16 million rows, 0.9e15 to 3.6e15.
* Fall images: at most 50,356 x 28 = 1.41 million rows x 50,356 x f; with f between 300 (the k = 10
  average, derived from the same counts: (2.33e11 - 16,618 x 3,147 x 3,919) / (28,016 x 3,147)) and
  30,000: 2e13 to 2.1e15.
* Reductions of followers and K: at most about 1e15 (bound 28 followers per round, 24,000 dense rows
  x fall rows x f each; the dense-row count is a guess).
* Total 1e14 to 7e15 multiply-adds. Recorded overall rates on the H200: 0.58e12 /s (engine closure),
  0.95e12 /s (first lift), 8.7e12 /s (checker, almost pure products). So **5 minutes to 3.5 hours**,
  plus reading and hashing 45 GB three times (file hashes, canonical hash, load: about 10 minutes) and
  the host-side row gathers for the followers (about 5 GB per follower, 28 per chunk and per fall
  round: estimate 10 to 60 minutes in all; `--basis-on device` removes them when the memory allows).
* **Planning figure: 0.5 to 4 hours per unit instance; deadline `--max-hours 6`.**

Proper instance (optional): all 10.5 million generator rows are computed (the first lift took 4,679 s
on it). Once K is near its final rank (codimension 2,370) a row costs 50,356 x 2,400 multiply-adds:
1.3e15 in all, plus the host gathers of 10.5 million rows (5,130 batches). Estimate 2 to 5 hours.

## 6. Keeping the basis

The k = 14 bases of the earlier runs were not archived (read: `paper/curve-lfd-bounds.tex`, "not
archived"; `RESULTS.md`, "not kept"). This plan keeps the checkpoint on the pod until the checker, the
independent lift and (optionally) the first lift have all run on the **same files**; the records carry
the manifest's file hashes (`basis.manifest_files_id`) and the canonical hashes.

Keeping it after the pod is gone is **not solved by the existing tooling** (unverified options, owner
decision): the launcher's `get` is for small files (default 20 MB, proxy limit about 100 MB: read
`POD.md`); the pod holds no cloud credentials by design; `--volume-gb` is a pod volume that
`terminate` removes. Possibilities: a Runpod network volume (not supported by the launcher today),
or 45 GB uploaded from the pod through presigned URLs. Since the closure is deterministic up to the
canonical hash (two runs per instance with equal hashes are on record), regenerating it costs 1.85 h
(about $10) each time; archiving matters only if further passes over V_6(14) are planned.

## 7. Cost

H200 at $5.29/h (read: `RESULTS.md`, 2026-11-10/11; **today's price unverified**).

| step, one unit instance | hours | cost |
|---|---|---|
| boot, push, arithmetic and small checks | 0.3 | $1.6 |
| closure V_6(14) | 1.85 (recorded 6,573 s, 6,634 s) | $9.8 |
| independent checker | 0.2 (recorded 731 s, 528 s) | $1.1 |
| independent lift | 0.5 to 4 (estimate), deadline 6 | $2.6 to $21.2, at most $31.7 |
| first lift again (optional) | 0.2 (recorded 751 s, 677 s) | $1.1 |
| **expected** | **3.1 to 6.6** | **$16 to $35** |
| **ceiling: pod lifetime `--max-seconds 36000`** | 10 | **$52.90** |

Both unit instances: expected $32 to $70, **ceiling $106** (two pods of 10 h), plus under $2 if a
smoke pod is needed. The optional proper control: closure 2.1 to 3.2 h (recorded 7,430 s on the
p = 65537 curve, 11,678 s on a slow pod), checker 0.2 h, lift 2 to 5 h: $23 to $45, ceiling $52.90
with the same cap. **Overall ceiling with the control: $160.** Stop rule: if the first instance ends
UNDECIDED by deadline or memory, do not start the second pod before the cause is recorded.

Cheapest version that is still independent evidence for the paper's sentence: `unit-p65537` only,
with `--stop-on-one` added to step 5 (stops when 1 enters K; proves 1 in V_7 but not the dimensions):
expected $14 to $25.

## 8. Assumptions and what this run cannot show

1. **Input containment (not checkable here).** Both lifts compute the degree-7 closure of the STORED
   degree-6 space W. The independent checker proves W contains V_6 (g in W, closed), which is the
   direction needed for "1 not in V_6". "1 in V_7" needs the other direction, W contained in V_6, i.e.
   that the engine added nothing but consequences of g. That rests on the engine (two runs with equal
   hashes, agreement with the reference engine at k = 14, D = 5 and at k <= 11, D = 6, and the
   dimension 448,922 shared with the proper instances), not on a certificate. A second lift does not
   change this. A lift from the doubly computed V_5(14) would, but it is out of reach: V_5(14) has
   74,589 free columns and almost every weight-6 monomial is special, so T has about 450,000 columns
   and K is a dense matrix of that size (the direct closure in another basis).
2. cupy: no line of the new module has run under cupy (not installed here). Used: advanced indexing
   with index arrays on both sides (`A[idx] = B`, `A[i, j] = 1`, `A[idx] += B`, `A[:, idx]`), boolean
   row masks, `@` on float64 (also 1-D @ 2-D), in-place arithmetic on slices (views), `argmax`,
   `concatenate`, `float()` / `int()` of 0-d arrays, `cupy.cuda.set_allocator(None)`. Step 2 exercises
   all of them at k = 10 with the same code path (`xp`), including the non-saturating schedule.
3. Exactness on the device: integer-valued float64 with partial sums below 2^53 (proof in the module
   docstring); it assumes cuBLAS `dgemm` and the elementwise kernels are IEEE-exact on such inputs,
   as the engine and the checker assume. Step 1 tests it at the bounds.
4. Time (section 5): the number of chunks to saturation at k = 14 is extrapolated from one local case
   and from the first lift's schedule.
5. Memory (section 4): formula, not measurement; n_single unknown.
6. Price, availability of an H200 with 200 GB RAM, the 200 GB disk, and the launcher's behaviour today.
7. The closure reproduces the recorded hashes (it did, twice per instance). If it does not, the chain
   records it and continues; the reading is then UNDECIDED until the difference is understood.
8. `chain_independent.py` was run here only as a dry run and, for its `small` stage, on CPU.
