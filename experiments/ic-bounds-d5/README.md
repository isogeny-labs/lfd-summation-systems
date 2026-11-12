# ic-bounds-d5: the degree-5 closure dimension d_5(k) at depths 10-13

Item 3 of the bounds programme: prepare (and, where it fits, take) the measurement of the
degree-5 closure dimension `d_5(k) = dim V_5` of two-summand prime chain systems for
`k = 10, 11, 12, 13`, starting from the depth-10 probe (`scripts/run_ic_chain_depth10_probe.py`,
`dim V_5 = 18,569` at `k = 10`, 2 h and about 8 GB per instance).

Evidence labels: **exact** = certified exact computation (certificate named), **proof** = written
argument given here in full, **empirical** = pattern in data, **conjecture**, **model** = a
projection from measured scaling (not a measurement). Nothing below is kernel-checked.

## Results

| claim | evidence |
|---|---|
| A new engine (`csd/ic_bounds_d5.py`, split + free-column RREF + exact tiled float64 products) computes the same closures as the dense engine and the depth-10 probe's shift engine: 18 closures at `k = 6, 7` (`D = 3, 4, 5`, three instance families) agree with the dense engine in dimension, in every filtration dimension `dim(V_D cap M_<=w)` and in `1 in V_D`; 12 agree with the shift engine; the 20 saved depth-8 systems reproduce `dim V_4 = 1061` and the saved `D = 5` rank certificates and low dimensions. | exact (validation.json; dense engine with exclusion certificates, shift engine verified closed, split engine verified closed) |
| `k = 10`: the three depth-10 instances (p = 12289 / 13313 / 15361) give `dim V_5 = 18,569`, `dim(V_5 cap M_<=4) = 3163`, identical to the shift engine's results, in about 250 s and 0.84 GB each (instead of 2 h, 8 GB). A paused-and-resumed run gives the bit-identical basis (same SHA-256). | exact (deficient-closed certificate: verified closure, rank below `|M_<=5| - |Z|`) |
| **`k = 11`: `d_5(11) = dim V_5 = 25,033`** on two instances (p = 12289; `y^2 = x^3 + x`, seed 10, `|Z| = 0`; `y^2 = x^3 + 4021x + 9917`, seed 12, `|Z| = 348`), with identical filtration profiles `[0, 0, 2, 142, 3441, 25033]` and identical part ranks `12546 + 12487`. `|M_<=5| = 35,443`, so the deficiency is 10,410 (unit ideal) and 10,062 (proper). **LFD >= 6 for both systems.** | exact (deficient-closed, verified closure; one run each, about 15.5 min and 0.86 GB peak RSS on one shared core) |
| At `k = 11` every system whose `dim V_5` is 25,033 has LFD >= 6, whatever its solution set, because `|Z| <= 2 * 2^11 = 4096 < 10,410`. | proof (below) |
| The degree-4 identity `d_4(k) = 4k^2 + 172k - 571` holds at `k = 11, 12, 13` (`dim V_4` = 1805, 2069, 2341) on two instances each, and again at `k = 6..10`; `dim(V_4 cap M_<=3) = 4k + 87` at every `k = 6..13` (16 instances). | exact per instance (verified closures, d4-check.json); the identities as functions of k stay empirical |
| `d_5(k) <= 1 + 2k * lambda(k)` with `lambda(k) = dim(V_5 cap M_<=4)`, for every instance and depth. | proof (below) |
| `d_5(10) = 18,569`, `d_5(11) = 25,033` are instance-independent on the instances computed (3 and 2 instances, three primes, two curves, unit and proper ideals). | empirical |
| `k = 12`: about 0.9 h on one core (3.2 h if lambda were as large as `|M_<=4|`), peak RAM <= 2.3 GB, checkpoint <= 0.77 GB; `k = 13` (p = 40961): about 2.4 h (10.9 h), peak RAM <= 4.2 GB, checkpoint <= 1.75 GB; `k = 14` (p = 65537, uint32): about 6 h (34 h), <= 13.5 GB. | model (measured throughput at k = 10, 11; storage bound `c_part^2/4`; scenario A for lambda) |


### Filtration profiles `dim(V_5 cap M_<=w)`, `w = 0..5`

| k | instances | dim V_4 | profile of V_5 | `|M_<=5|` | top-weight leading terms missing (`C(2k,5) - (d_5 - lambda)`) |
|---|---|---|---|---|---|
| 8 | probe (closed-full) | 1061 | [1, 17, 137, 697, 2517, 6885] | 6885 | 0 |
| 9 | probe (closed-full) | 1301 | [1, 19, 172, 988, 4048, 12616] | 12616 | 0 |
| 10 | 3 | 1549 | [0, 0, 2, 138, 3163, 18569] | 21700 | 98 |
| 11 | 2 | 1805 | [0, 0, 2, 142, 3441, 25033] | 35443 | 4742 |

(`k <= 9`: V_5 contains 1 for the unit instance, i.e. it is all of `M_<=5`; the closure has not separated
from the ideal, so these are not values of a "generic" `d_5` below the ideal.)
Symmetric/antisymmetric parts at k = 10: `[0,0,2,70,1536,7703]` / `[0,0,0,66,1489,7703]`; at k = 11:
`[0,0,2,72,1674,10798]` / `[0,0,0,68,1625,10794]` (pivot counts per weight).

## The engine (csd/ic_bounds_d5.py)

Modulo the chains the two-summand system lives in the group algebra of `G = (Z/2^k)^2`; a squarefree
monomial is an exponent pair, its degree the Hamming weight, and `u_{i,j}` adds `2^j` to `a_i`. Three
changes against the depth-10 probe:

1. **Symmetry split.** The summand swap `sigma` fixes `g`; `V_D = V^+ (+) V^-` (p odd), each part kept on
   the `sigma`-orbit coordinates (`~|M_<=5|/2` columns). Only the `k` block-1 shifts are applied; each
   shifted row is projected to both parts.
2. **Free-column RREF.** Each part is stored as its pivot columns plus the coefficients on the
   non-pivot columns: `r x (c - r)` uint16 entries instead of `r x c` int64. Reducing a row is the single
   product `y[free] - y[piv] @ B` (the reduced row vanishes on the pivots), cost `r(c - r) <= c^2/4`.
   At the depth-10 point this stores 29 M entries (58 MB) instead of 400 M int64 (3.2 GB).
3. **Exact tiled products.** float64 GEMMs whose inner dimension is cut so that every partial sum is
   below `2^53`, reduced mod p after every tile (exact for residues; the `ic_chain_block_closure`
   technique).

**Certificate (proof).** Let W be the final space. (i) `W subset V_D`: every row added is a generator
row `NF(mu*g)` (`wt mu + 2 <= D`), a block-1 shift of a row of W of weight `< D`, or a `sigma`-projection
`(f +- sigma f)/2` of such a row; `V_D` contains the first two by definition and the third because it is
`sigma`-invariant (the inputs are, and the closure rule commutes with `sigma`). (ii) `V_D subset W`: W
contains the generator rows. `W = W^+ (+) W^-` is `sigma`-invariant, and `W cap M_<D` is spanned by the
part rows whose pivot has weight `< D` (columns are ordered by weight descending, so such a row is
supported on `M_<D`). The verification sweep re-reduces the `k` block-1 shifts of every final low row
of both parts and adds nothing, so `u_{1,j} f in W` for every `f in W cap M_<D`; for block 2,
`u_{2,j} f = +- sigma(u_{1,j} sigma f)` with `sigma f = +- f`, which lies in W because W is
`sigma`-invariant. So W is a subspace containing the inputs and closed under every admissible product,
hence contains the least one, `V_D`. Thus `W = V_D` exactly. If a sweep adds a row, the closure phase
resumes and a full new sweep follows; only a sweep with zero additions certifies.

**Consequences (proof).** If `dim V_5 < |M_<=5| - |Z|`, then since the chain ideal is radical with all
points rational, `I` modulo the chains is the vanishing ideal of `Z` and
`dim(I cap M_<=5) >= |M_<=5| - |Z| > dim V_5`; some `h in I` of degree `<= 5` is outside `V_5`, and since
`V_{max(d, deg h)} subset V_5` for `d <= 5`, `LFD >= 6`. `|Z| <= 2 * 2^k`: for fixed `x1`, `S_3(x1, x2, X)`
has degree `<= 2` in `x2` and is not the zero polynomial (for fixed points `P1, P3` only finitely many
`P2` satisfy `+-P1 +- P2 +- P3 = O`), so each of the `2^k` values of `E1` has at most two partners.

**`d_5 <= 1 + 2k lambda` (proof).** Let `L = V_5 cap M_<5` and `W' = span(g) + sum_v u_v L`. Then
`W' subset V_5`, and `W' cap M_<5 subset L`, so `u_v (W' cap M_<5) subset W'`: `W'` is closed and contains
`g`, hence `W' = V_5`, and `dim V_5 <= 1 + 2k dim L`. So `deg d_5 <= 1 + deg lambda` if both are
polynomials, and the growth question reduces to `lambda(k) = dim(V_5 cap M_<=4)` (3163, 3441 at
k = 10, 11; `|M_<=4|` = 6196, 9109).

## What the values of d_5 decide, and what they cannot

Decides, per instance and exactly:
- `LFD >= 6` at each depth where the verified closure is deficient (`k = 10, 11` done). At `k = 12, 13`
  it suffices that `d_5 < |M_<=5| - |Z|` (`|M_<=5|` = 55,455 / 83,682; `|Z|` = 1352 / 1484 on the
  prepared instances; at most `2^{k+1}` always).
- With `lambda` and the profile, whether the top-weight layer keeps losing ground: missing top leading
  terms grew from 98 (k = 10) to 4742 (k = 11).

Can be decided only with extra structure:
- **Degree of d_5.** Values at `k = 10..13` (four points) determine a cubic interpolant and test
  nothing; `k = 10..14` (five points, `k = 14` needs p = 65537 and uint32 storage) allow one check of a
  cubic and none of a quartic. A quintic with leading coefficient `4/15` (the order of `|M_<=5|`)
  cannot be excluded by any finite set of values. `k = 10` may still be transitional (only 98 top
  leading terms missing, against 4742 at k = 11): a fit that includes it may describe no regime.
  Polynomiality itself is not established for any `k`-range (`d_4` is quadratic on `k = 6..13`
  empirically, not by proof).
- The useful decision is through `lambda`: by the proved inequality, `lambda = O(k^2)` would give
  `d_5 = O(k^3)`, so `|M_<=5| - d_5 -> infinity` and, for systems that attain the generic values, the
  step to `LFD >= 6` is permanent (conjecture until a genericity theorem exists). Two values of
  `lambda` (3163, 3441) fix nothing yet; four (`k = 10..13`) give one check of a quadratic `lambda`.
- **Predicted next step.** `d_5` says nothing about the step `6 -> 7`, which is decided by
  `d_6(k)` against `|M_<=6| - |Z|`. Steps so far: `LFD = 4 -> 5` at `k = 6`, `5 -> >= 6` at `k = 10`
  (empirical; `4 -> 5` on the computed instances). A step every four depths would put `>= 7` near
  `k = 14` (conjecture, no evidence beyond two step locations).

Cannot decide:
- Any upper bound: `LFD <= 6` needs `V_6 = I cap M_<=6` (degree-6 closure; the same engine with
  `--degree 6`, see heavy runs).
- Genericity: equal values on 2-3 instances per depth are empirical; a statement for all instances
  outside an exceptional locus needs the parametric argument of manuscript Section 7, item 1.
- Behaviour as `k -> infinity`, the `m >= 3` systems, binary/Koblitz systems.

## Cost model and the plan for k = 12, 13 (model)

Measured (one contended core of a 4-core host shared with three CPU-bound jobs): reduction
throughput 1.7e10 multiply-adds/s (verification sweep), 7.5e9/s (closure phase, which also runs the
scalar C echelon kernel and the RREF updates). Work at `k = 11`: 4.8e12 (closure) + 4.9e12 (sweep)
= 0.80-0.83 of the bound `k * lambda * (c_sym^2 + c_anti^2)/4` per phase (k = 10: 0.79 and 0.49). Peak free-column storage
156.6 M entries against the bound 157.0 M (`c^2/4` per part, attained when `r = c/2`).

| k | p | `|M_<=5|` | part columns | storage/checkpoint bound | peak RAM (conservative) | 1-core hours, lambda linear | 1-core hours, lambda <= `|M_<=4|` |
|---|---|---|---|---|---|---|---|
| 12 | 12289 | 55,455 | 27,767 / 27,688 | 0.77 GB | 2.3 GB | 0.92 | 3.2 |
| 13 | 40961 | 83,682 | 41,887 / 41,795 | 1.75 GB | 4.2 GB | 2.4 | 10.9 |
| 14 | 65537 (uint32) | 122,438 | 61,272 / 61,166 | 7.5 GB | 13.5 GB | 6.0 | 34 |

(table.md and summary.json carry the exact model outputs.) The probe was 2 h / 8 GB at k = 10 with the
dense int64 shift engine; its storage scales like `d_5 * |M_<=5|` int64 (up to 56 GB at k = 13), which
is what the split + free-column design removes. With 4 BLAS threads on an idle host the sweep phase
should run several times faster; the closure phase is partly single-threaded (C echelon kernel).

What makes `k = 12, 13` feasible: (1) halving the columns and the shifts by the swap symmetry,
(2) storing only the free-column block `r(c - r) <= c^2/4` in uint16 (p < 2^16: p = 12289 for k <= 12,
p = 40961 = 5 * 2^13 + 1 for k = 13; k = 14 needs p = 65537 and uint32), (3) exact tiled float64
GEMMs instead of the scalar elimination kernel, (4) checkpoints after every round (or every
`--save-every` seconds) with a `--max-minutes` pause, so a run on a preemptible host can be split.
Implemented and validated here at k = 6-11.

## Prepared instances for the heavy runs (constructed here, not closed)

| k | p | curve | seed | target | `|Z|` | `|M_<=5|` | exact `dim V_4` |
|---|---|---|---|---|---|---|---|
| 12 | 12289 | a=1, b=0 | 10 | (9362, 1126) | 1352 | 55,455 | 2069 |
| 12 | 12289 | a=4021, b=9917 | 12 (--want-solutions) | (4408, 4441) | 1399 | 55,455 | 2069 |
| 13 | 40961 | a=1, b=0 | 10 | (2136, 2641) | 1484 | 83,682 | 2341 |
| 13 | 40961 | a=4021, b=9917 | 12 (--want-solutions) | (31102, 3136) | 1582 | 83,682 | 2341 |

At k = 12, 13 the subgroup has order `2^k` comparable to p, so these ideals are proper (no unit ideal
is expected); LFD >= 6 follows if `dim V_5 < |M_<=5| - |Z|`.

## Commands

```sh
# all runs: ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/ic_bounds_d5_validate.py                    # validation.json (about 8 min)
python3 scripts/ic_bounds_d5_probe.py --k 8 --no-early-stop --tag k8-p12289-s10-D5-split-full
python3 scripts/ic_bounds_d5_probe.py --k 9 --no-early-stop --tag k9-p12289-s10-D5-split-full
python3 scripts/ic_bounds_d5_probe.py --k 10                                     # probe instance
python3 scripts/ic_bounds_d5_probe.py --k 10 --prime 13313 --seed 11 --want-solutions
python3 scripts/ic_bounds_d5_probe.py --k 10 --prime 15361 --a 4021 --b 9917 --seed 12 --want-solutions
python3 scripts/ic_bounds_d5_probe.py --k 11                                     # d_5(11), unit ideal
python3 scripts/ic_bounds_d5_probe.py --k 11 --a 4021 --b 9917 --seed 12 --want-solutions
python3 scripts/ic_bounds_d5_probe.py --k 6 --degree 4 --engine shift            # (and split; D = 5)
python3 scripts/ic_bounds_d5_d4check.py                    # d4-check.json
python3 scripts/ic_bounds_d5_analyze.py                    # summary.json, table.md, d5.svg/png
```

### Heavy runs (not run here: beyond the 25-minute / 300 MB limits of this session)

```sh
# k = 12 (two instances); about 1 h each on one core, <= 2.3 GB RAM, <= 0.8 GB checkpoint
ulimit -v 6000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
python3 scripts/ic_bounds_d5_probe.py --k 12 --save-every 1200
python3 scripts/ic_bounds_d5_probe.py --k 12 --a 4021 --b 9917 --seed 12 --want-solutions --save-every 1200
# k = 13 (two instances); about 2.5 h each on one core, <= 4.2 GB RAM, <= 1.8 GB checkpoint
ulimit -v 9000000
python3 scripts/ic_bounds_d5_probe.py --k 13 --prime 40961 --save-every 1800
python3 scripts/ic_bounds_d5_probe.py --k 13 --prime 40961 --a 4021 --b 9917 --seed 12 --want-solutions --save-every 1800
# optional fifth value, k = 14 (p = 65537, uint32): about 6 h, <= 13.5 GB RAM, <= 7.5 GB checkpoint
ulimit -v 20000000
python3 scripts/ic_bounds_d5_probe.py --k 14 --prime 65537 --save-every 1800
# then
python3 scripts/ic_bounds_d5_analyze.py
# optional, for item 4 (exact LFD at k = 10): the degree-6 closure of the probe instance (unit ideal);
# 'reached' means 1 in V_6, so LFD = 6 exactly. Model: storage <= 0.91 GB, peak RAM <= 2.5 GB,
# closure <= 3.7 h on one core (+1.6 h verification sweep if deficient; then LFD >= 7)
ulimit -v 6000000
python3 scripts/ic_bounds_d5_probe.py --k 10 --degree 6 --prev split --save-every 1800
```

On a preemptible host add `--max-minutes N` (exit status 3 after a checkpoint; rerun the same
command). Each finished run writes `runs/<tag>.json`; compare `dim V_4` with the table above and
`dim V_5` across the two instances of each depth.

Checkpoints default to `~/.cache/ic-bounds-d5/<tag>` (`--ckpt` to move them); a paused run exits with
status 3 and the same command resumes it. The runs here used checkpoints in the session scratchpad and
deleted them afterwards (no checkpoint is committed; every record carries `basis_sha256`, a hash of the
pivots and free-column blocks of both parts).

## Files

- `runs/*.json`: one record per probe run (instance, target, `|Z|`, `dim V_{D-1}` with its engine,
  result and certificate status, profiles, part ranks, rank history, per-phase work and seconds,
  peak storage, peak RSS, source SHA-256s, basis SHA-256).
- `validation.json`: the validation of the split engine (dense, shift, saved depth-8 results).
- `d4-check.json`: exact `dim V_4` at `k = 6..13`.
- `summary.json`, `table.md`, `d5.svg`, `d5.png`: aggregation, cost model, projections.
- `logs/`: JSON-lines logs of every run.
- `verification.json`: SHA-256 of the sources and of every data file here.

Provenance note: `csd/ic_chain_block_closure.py` (imported by the probe only for the shift engine's
verification) was changed by another session while these runs were in progress; the records
`k11-...-s12-...-proper`, `k8-...-full`, `k9-...-full` carry its earlier hash. The split engine does not
use that module.
