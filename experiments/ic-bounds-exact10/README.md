# ic-bounds-exact10: preparing the exact LFD at depth 10, and the next step (item 4)

> **Note (2026-11-12).** Some runs listed below as not run were run later and are reported in `paper/curve-lfd-bounds.tex`, including the degree-6 closures at `k = 10..14` (`experiments/ic-bounds-k14/`). The lists below record the state when this file was written.

Two-summand prime chain systems: subgroup (or coset) chains of depth `k` over `F_p` with
`2^k | p-1`. The normalized generator is `g = NF(S_3(c x1, c x2, x_R))` modulo the chains. LFD
means the ordinary iterated last fall degree. Every claim below carries one of these labels:
**kernel** (Lean-checked), **proof** (argued in full here), **exact** (a certified computation,
with the certificate named), **empirical**, or **conjecture**.

## Results at a glance

| item | outcome | label |
|---|---|---|
| (a) block engine `csd/ic_chain_block_closure.py` vs the dense engine, k in {5,6,7}, D in {4,5,6}, 8 instances (unit and proper, j=0, j=1728, random, coset); 72 rows (33 unit, 39 proper) | 1 in V_D agrees in 72/72 rows. In all 45 rows where it stops without 1, sym_rank + anti_rank equals the dense dim V_D (completeness observed, not certified) | exact (dense closures, exclusion-certified) |
| (a) 1 NOT in V_D | 45 rows: 6 unit rows (k = 6 and 7, D = 4) and all 39 proper rows. The block engine reports `contains_one = False` in each | exact |
| (b) new module `csd/ic_bounds_exact10.py` (unit **and** proper, deterministic closure verification), same grid | 72/72 rows from scratch and 72/72 from V_{D-1}: exact dimension, dim(V_D cap M_<D) and 1 in V_D equal the dense values. Negative control rejected in 56/56 applicable rows. 24/24 profiles reproduce the exact LFD and strict falls of `prime_local_profile` (LFD 4 at k=5, LFD 5 at k=6,7) | exact |
| (b) mid-scale | k = 8 (3 stored systems) and k = 9 (1 stored system): V_4 verified closed with dim 4k^2+172k-571 (1061, 1301); V_5 reaches the stored shift-engine certificate (6633, 6885, 6803; 12592); exact LFD 5 reproduced | exact |
| runner `scripts/ic_bounds_exact10_degree6.py` at k = 8, D = 6 (proper, 252 points), interrupted after one round and resumed | `reached`: dim V_6 = 14641 = 14893 - 252 | exact (target reached, ev(M_<6) onto Z) |
| depth-10 runner memory | the original `run_ic_chain_depth10_degree6.py` initialisation needed more than 15 GB; the main session's 22:36 patch addresses it (see below) | exact (arithmetic on array shapes) |
| (c) LFD 7 | design, cost table and commands below. Dense split closures can certify each depth only up to k = 13 (113 GB). k = 14 needs at least 0.5 TB, or a new engine | model / conjecture |

No exact LFD at depth 10 is claimed here: that needs the heavy run listed under
"Heavy runs needed".

## Files

| path | what |
|---|---|
| `csd/ic_bounds_exact10.py` | `SplitClosure`: the symmetry-split degree-D closure for unit and proper ideals. Adds inputs, random rounds (dense or sparse mixing), a deterministic verification pass, certificate classification, `split_profile`, and `memory_model` |
| `scripts/ic_bounds_exact10_validate.py` | parts (a) and (b) at small k; writes `validation-*.jsonl` and `profiles-*.jsonl` |
| `scripts/ic_bounds_exact10_midscale.py` | stored depth-8 and depth-9 systems; writes `midscale-k8.jsonl` and `midscale-k9.jsonl` |
| `scripts/ic_bounds_exact10_degree6.py` | general memory-lean runner (any k, D; unit or proper; V_{D-1} input optional; checkpoint and resume; `--dry-run`) |
| `scripts/ic_bounds_exact10_plan.py` | part (c): per-depth primes, sizes, memory model, time extrapolation, commands; writes `plan.json` |
| `scripts/ic_bounds_exact10_summarize.py` | `summary.json`, `calibration.json`, `hashes.json`, `cost.svg`, `cost.png` |
| `degree6-k8-D6-p257-x190.json` | runner result at k = 8 |
| `logs/` | stdout of every run (`dryrun-k10.log` reconstructs the three depth-10 instances) |

## The engine and why its answers are certificates (proof)

Notation: `M_<=D` is the set of squarefree chain masks of weight `<= D`. Multiplication by a
chain variable is the exponent shift. For `D >= 2` the closure modulo the chains equals the
ordinary closure V_D (docstring of `csd/ic_fall_certificates.py`, as in all existing engines).

1. **Split.** The block swap sigma is a ring automorphism of the chain quotient. It permutes the
   variables, preserves weight, fixes the chain ideal and fixes `g`, because `S_3` is symmetric
   and both summands use the same scale. So sigma(V_D) is the least closed space containing
   sigma(inputs) = inputs, which gives sigma(V_D) = V_D. For any subspace W and char p != 2,
   `(1+sigma)W (+) (1-sigma)W = W + sigma W`, since `w = (1+sigma)w/2 + (1-sigma)w/2`. The
   `sym` and `anti` coordinate maps are injective on symmetric and antisymmetric vectors
   respectively. Hence `rank(sym part) + rank(anti part) = dim(W' + sigma W')` for the span W' of
   the inserted rows. Every inserted row lies in V_D, so this rank is at most dim V_D.
2. **Target.** The chain quotient is the group algebra `F_p[(Z/2^k)^2]`. It is semisimple and
   split because `2^k | p-1`, so every ideal is radical and is the vanishing ideal of its zero set
   Z, whose points are all rational. Therefore
   `dim(I cap M_<=D) = |M_<=D| - rank ev(M_<=D -> F_p^Z)` (the `target`), and
   `V_D subset I cap M_<=D`.
3. **reached.** If rank = target, then V_D = I cap M_<=D. If, in addition, `ev(M_<D)` is onto
   `F_p^Z`, the local-to-global step (proof in the docstring of
   `experiments/ic-curve-validation/code/lfd_engine.py`) gives LFD <= D. For the unit ideal this
   case is exactly `1 in V_D`.
4. **closed-deficient.** The verification pass reduces all 2k shifts of every low row (pivot
   weight `< D`) of both parts against both parts. If a full pass adds nothing, the bases were
   unchanged during the pass. Columns are sorted by weight in descending order and sigma
   preserves weight, so the low pivot rows of the two parts span
   `(W cap M_<D) = ((1+sigma)W cap M_<D) (+) ((1-sigma)W cap M_<D)`. A shifted vector lies in W
   exactly when both of its projections reduce to zero. So W contains every input (the rows
   `NF(mu g)`, `g` among them), is closed under every admissible product, and lies inside V_D.
   Hence W = V_D. If rank < target, some h in I of degree <= D is not in V_D, so
   **LFD >= D + 1**. For the unit ideal this is "1 not in V_D".
5. **Exact value.** Unit ideal: LFD is the least D with 1 in V_D. Proper ideal: LFD is the last
   strict fall. A fall at D occurs exactly when `dim(V_D cap M_<D) > dim V_{D-1}`, and both
   dimensions are exact once each degree is `reached` or `closed-deficient`.

Random combinations (dense or sparse mixing) only change how fast the rank grows. A stall is
never reported as a bound unless the clean verification pass has run.

## (a) Block-engine validation

`block_closure_contains_one`, the code imported by `run_ic_chain_depth10_degree6.py`, was
tested in two versions. The main session edited `csd/ic_chain_block_closure.py` at 22:36
(chunked V_{D-1} insertion) and again at 23:31:40, when only the default float64 tile of the
exact `_modmat` changed (same file size).
- Batches k56 and k7a-d imported sha256 `865032b874b2...`.
- Batch k7e and a full k=5,6 repeat (`validation-k56-v2.jsonl`, 48/48 rows; 16/16 profiles)
  imported the current `de5e39705f88...`.
`summary.json` maps each batch to the version it imported (`validation_batches`). The engine was
started
from the dense, exclusion-certified V_{D-1} at every (instance, k, D) with k in {5,6,7}, D in
{4,5,6}. The reference is the dense closure `chain_degree_closure`. At k = 7, D = 6 the dense
closure exceeds its product budget, so the reference there is the local-to-global theorem
applied to the dense certificate at D = 5 (proof; labelled `ref_source` in each row). Counts are
in `summary.json`:

- `block_one`: 1 in V_D agrees with the reference in every row with a block-engine run.
- In every row where 1 is not in V_D (all proper rows, and unit rows at k = 6, 7 with D = 4),
  the engine returns `contains_one = False`. Its sym_rank + anti_rank there equals the exact
  dense dim V_D. In these cases the random rounds found the whole closure (empirical
  completeness). By itself, though, the block engine's stall certifies nothing: only the new
  module's verification pass does.

## (b) Proper ideals

`SplitClosure` adds the proper-ideal target `|M_<=D| - rank ev(M_<=D)`, the `ev(M_<D)`
onto-check, the deterministic verification pass, and the exact-LFD bookkeeping. On the grid
above it was run from the generator rows alone (dense mixing) and from the dense V_{D-1} (sparse
mixing, 4 terms). Every exact dimension, `dim(V_D cap M_<D)` and `1 in V_D` equals the dense
value. A negative control presents V_{D-1} plus the generator rows as if they were V_D. The
verification pass must reject it exactly when that space is smaller than V_D, and it does so in
every applicable row (`checks.negative_control`; rows where V_{D-1} already contains 1 are
excluded). `split_profile` reproduces the exact LFD and the strict falls (from degree 3) of
`prime_local_profile` on every (instance, k).

Mid-scale results (exact): see the table above, `midscale-k8.jsonl` and `midscale-k9.jsonl`. At
k = 9 the split engine reached the D = 5 certificate in 432 s, against 1318 s for the stored
shift-engine run (wall times on a shared, loaded host).

## Depth-10 runner: memory (exact arithmetic on array shapes)

As first written, `block_closure_contains_one` built `np.zeros((18569, 60460), int64)` for V_5
(9.0 GB), plus an int64 copy of V_5 (3.2 GB) and int64 `sym`/`anti` projections (4.5 GB
each). That is more than the 15 GB host. Four processes of `run_ic_chain_depth10_degree6.py`,
each with about 16.6 GB of virtual memory and 12.7-13.3 GB anon-rss, were OOM-killed by the
cgroup (dmesg, about 22:34). The main session then patched the engine to insert V_5 in
2000-row chunks (22:36), and that run (pid 11870) was still going at the time of writing. Model
of its remaining peak, with combos 6000, low_batch 3000 and tile 4096: about 3.7 GB of uint16
parts at full rank, plus 2.9 GB of float64 Y, a transient 2.9 GB int64 copy, two 1.45 GB
projections and about 3.5 GB of reduction temporaries. That is roughly 10 GB at late rounds:
it fits on this host only with little else running. `scripts/ic_bounds_exact10_degree6.py`
needs about 6.6 GB under the same model (combos 2000, sparse mixing, tile 1024, memory-mapped
V_5). It also handles the two proper depth-10 instances.

## (c) Measuring the depth at which LFD reaches 7

**What must be measured (proof).** Fix a depth k where `ev(M_<5)` is onto Z (true in every
instance so far).
- LFD >= 7 at k iff V_6 != I cap M_<=6, and for the unit ideal iff 1 is not in V_6. If
  V_6 = I cap M_<=6 and ev(M_<6) is onto Z, local-to-global gives LFD <= 6. Otherwise some h in I
  of degree <= 6 is missing from V_6 = V_max(6, deg h).
- Reduction: V_6 = I cap M_<=6 iff V_6 cap M_<=5 = I cap M_<=5. For the forward direction,
  intersect. For the reverse direction, take h in I cap M_<=6. Every weight-6 mask is
  `m = x_v m'` with m' of weight 5 and no carry (clear one set bit). By onto-ness pick
  q in M_<5 with m' - q in I cap M_<=5 subset V_6. Then x_v(m' - q) is in V_6 (admissible), so
  h minus the corresponding combination lies in I cap M_<=5 subset V_6.
- So the quantity to measure is the fall at 6, `f_6(k) = dim(V_6 cap M_<=5) - dim V_5`, against
  the V_5 deficiency `delta_5(k) = |M_<=5| - |Z| - dim V_5`. LFD >= 7 iff `f_6(k) < delta_5(k)`.
  At depth 10 (unit, p = 12289): `delta_5 = 21700 - 18569 = 3131` (exact, verified-closed V_5).
  For comparison, the previous step: `f_5(9) = delta_4(9) = 2747`, while
  `f_5(10) = 3163 - 1549 = 1614 < delta_4(10) = 4647` (exact).

**Protocol per depth k.** k_7 is the least k whose closed V_6 is deficient.
1. `--degree 5`: a closed-deficient V_5 certifies LFD >= 6 at k and records dim V_5 and delta_5.
2. `--degree 6` starting from that V_5: `reached` gives exact LFD = 6 at k (strict fall at 6).
   `closed-deficient` gives LFD >= 7 at k, so k_7 <= k.
3. Only at the first deficient k: the exact value 7 needs `--degree 7` (unit: 1 in V_7).
Order of depths: k = 10 (unit, then the two proper instances), 11, 12, 13. Stop at the first
deficient depth. Instances: the probe's unit j=1728 instance (seed 10) plus one proper
`--want-solutions` instance per depth. At k = 10 three instances gave identical dim V_4 and
dim V_5 (empirical genericity), so one instance per depth locates the step; a second instance
checks it.

**Cost (model; `plan.json`, `cost.svg`).** The dense split engine is at full rank at the
`reached` end. Its peak memory model, single machine:

| k | p | |M_<=6| | D=6 peak | D=5 peak | D=7 peak |
|---|---|---|---|---|---|
| 10 | 12289 | 60,460 | 6.6 GB | 1.5 GB | 26 GB |
| 11 | 12289 | 110,056 | 17 GB | 2.9 GB | 92 GB |
| 12 | 12289 | 190,051 | 45 GB | 5.5 GB | 313 GB |
| 13 | 40961 | 313,912 | 113 GB | 11 GB | 989 GB |
| 14 | 65537 | 499,178 | 520 GB (uint32: no prime = 1 mod 2^14 below 2^16) | 35 GB | 5.7 TB |

Time is extrapolated from the measured reached runs (`calibration.json`: 21 runs from 3473 to
14893 columns; shared host with load about 5 on 4 cores, single BLAS thread). The measured runs
scale like cols^2 at these sizes (k=8 to k=9, D=5: 130 s to 432 s for 1.83x the columns). Dense
elimination tends to cols^3. The two laws bracket the estimate. Single-thread hours, an
order-of-magnitude estimate rather than a measurement:

| k | D=6, cols^2 law | D=6, cols^3 law (`plan.json` `hours_est`) |
|---|---|---|
| 10 | 1.8 h | 9.8 h |
| 11 | 6 h | 59 h |
| 12 | 18 h | 306 h |
| 13 | 49 h | 1,377 h |
| 14 | 123 h | 5,537 h |

Multi-threaded BLAS divides these by roughly the core count.

**Conjecture (labelled).** The onsets so far are k = 6 (LFD 5) and k = 10 (LFD >= 6). Petit,
Kosters and Messeng report Magma degree of regularity 6 at k = 10 and 11, which is consistent
with LFD 6 there. An equal gap would put k_7 at 14. We put k_7 in 12..17. This rests on two data
points and is not a bound. If k = 11, 12 and 13 all reach V_6, the dense method stops, because
k = 14 needs at least 0.5 TB. Options then:
- out-of-core blocked elimination;
- a non-split prime below 2^16 with a generic symmetric generator of the same support. This
  gives uint16 storage and about 260 GB, but the result is evidence about generic systems, not
  about the ECDLP instance;
- a top-degree syzygy formulation that computes f_6 directly in `M_=6` coordinates. This has the
  same column count, but could be done with sparse or Wiedemann methods.

## Reproduce

Source and data hashes are in `hashes.json`. Every Python run below used
`ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

```sh
python3 scripts/ic_bounds_exact10_validate.py --ks 5 6 --tag=-k56
python3 scripts/ic_bounds_exact10_validate.py --ks 7 --dense-max-cells 120000000 --instances j1728-p7681-x936 rand-p257-x58 --tag=-k7a
python3 scripts/ic_bounds_exact10_validate.py --ks 7 --dense-max-cells 120000000 --instances rand-p7681-x2244 rand-p12289-x11294 --tag=-k7b
python3 scripts/ic_bounds_exact10_validate.py --ks 7 --dense-max-cells 120000000 --instances j0-p7681-seed41 coset5-p12289-seed42 --tag=-k7c
python3 scripts/ic_bounds_exact10_validate.py --ks 7 --dense-max-cells 120000000 --instances proper-p12289-seed43 --tag=-k7d
python3 scripts/ic_bounds_exact10_validate.py --ks 7 --dense-max-cells 120000000 --instances proper-coset7-p7681-seed44 --tag=-k7e
python3 scripts/ic_bounds_exact10_validate.py --ks 5 6 --tag=-k56-v2   # repeat against the 23:31 block engine
python3 scripts/ic_bounds_exact10_midscale.py --save-v5 $SCRATCH/v5-k8
python3 scripts/ic_bounds_exact10_midscale.py --depth-dir experiments/ic-chain-depth9 --systems 0 --out experiments/ic-bounds-exact10/midscale-k9.jsonl
H=$(sha256sum $SCRATCH/v5-k8/basis.npy | cut -d' ' -f1)   # 636d184a...9f32
python3 scripts/ic_bounds_exact10_degree6.py --k 8 --degree 6 --prime 257 --a 1 --b 0 --target-x 190 --v5 $SCRATCH/v5-k8 --v5-sha256 $H --tag k8-D6-p257-x190 --ckpt $SCRATCH/ck-k8 --max-rounds 1 --out-dir $SCRATCH
python3 scripts/ic_bounds_exact10_degree6.py --k 8 --degree 6 --prime 257 --a 1 --b 0 --target-x 190 --v5 $SCRATCH/v5-k8 --v5-sha256 $H --tag k8-D6-p257-x190 --ckpt $SCRATCH/ck-k8   # resumes
python3 scripts/ic_bounds_exact10_degree6.py --dry-run --k 10 --prime 12289 --a 1 --b 0 --seed 10
python3 scripts/ic_bounds_exact10_degree6.py --dry-run --k 10 --prime 13313 --a 1 --b 0 --seed 11 --want-solutions
python3 scripts/ic_bounds_exact10_degree6.py --dry-run --k 10 --prime 15361 --a 4021 --b 9917 --seed 12 --want-solutions
python3 scripts/ic_bounds_exact10_summarize.py && python3 scripts/ic_bounds_exact10_plan.py && python3 scripts/ic_bounds_exact10_summarize.py
```

Run notes:
- The k = 7 grid was split into batches of two instances, to stay under 25 minutes per run and
  3 GB per process.
- The first `-k7d` batch completed `proper-p12289-seed43`, then hit the 3 GB virtual-memory cap
  inside the dense reference for the next instance. That instance was rerun alone as `-k7e`.
- Batches k56 and k7a loaded `csd/ic_bounds_exact10.py` before its last edit. That edit only
  relaxed the prime guard from p < 2^16 to p < 2^25; every prime used here is below 2^16.
- Seeds: `--seed 20261107`. Round seeds are `seed + 10 D` for scratch runs and `seed + 10 D + 1`
  for runs started from V_{D-1}; block-engine seeds are `seed + D`.
- Instances are listed in `INSTANCES` in the validate script. Targets drawn from seeds are
  recorded in every row (`target_x`), together with a generator hash.

## Heavy runs needed (not run here)

See the structured result. In brief:
1. Exact LFD at k = 10 for the unit instance: either the main session's patched runner, or
   `scripts/ic_bounds_exact10_degree6.py --k 10 --degree 6 --prime 12289 --a 1 --b 0 --target-x 9362 --v5 ~/.cache/ic-chain-depth10/rerun --v5-sha256 d4fe133411ff396d4b607165035615cde9e796f112f4a26bf23db1ef2612ee6b`.
   About 7 GB.
2. Exact LFD at k = 10 for the proper instances: 13313 (x = 7412, 98 points) and 15361
   (a = 4021, b = 9917, x = 7776, 60 points). A `reached` V_6 together with the verified
   deficient V_5 gives exact LFD 6.
3. Depths 11, 12, 13 at degrees 5 and 6, with the commands in `plan.json`: 17 GB, 45 GB and
   113 GB for degree 6.

The ecdb curve database was not consulted: these are toy-prime closure dimensions, not curve
facts, and no AWS secret was used.
