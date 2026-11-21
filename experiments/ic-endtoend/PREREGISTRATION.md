# Preregistration: relation yield and linear algebra of PKM prime-field index calculus, end to end, against Pollard rho

Frozen 2026-11-20, before any data was taken. Every number that later appears
in `README.md` is `[pcc:unverified]` until a pcc claim is filed for it.

## Question

The paper (`paper/curve-lfd-bounds.tex`) bounds the cost of one
relation-search attempt by a Macaulay-matrix solver. A fully charged cost
comparison with Pollard rho also needs (a) the **relation yield**, the
probability that a random group element decomposes over the factor base, and
(b) the **linear-algebra stage**. This experiment measures both, runs the whole
attack end to end on small curves, verifies every recovered logarithm, and runs
Pollard rho on the same curves and targets.

## Setting (Petit, Kosters and Messeng, Section 3.3)

* Prime `p` with `2^k | p-1`. `H_k` is the subgroup of order `2^k` of `F_p^*`.
* `E : y^2 = x^3 + a x + b` over `F_p` of **prime order** `N` (so every point
  other than `O` generates the group, and the factor base lies in it).
* Factor base `F = {P in E(F_p) : x(P) in t*H_k}`, with `t = 1` (subgroup, the
  PKM choice) and, as a second arm, one coset `t*H_k` with `t` not in `H_k`.
  `F` is closed under negation. `|F|` counts **points** (both signs); the number
  of `+-` classes is `f = |F|/2`. `|F|` is measured exactly per curve and `k`.
* A relation is `R = e_1 P_1 + ... + e_m P_m` with `e_i in {+1,-1}` and
  `P_i` class representatives (repetition allowed), `m in {2,3}`.

## Primes and curves (fixed rule, no choice after data)

For each bit size `b in {16, 20, 24, 28, 32, 36, 40}`: `K = min(b-4, 24)` and
`p` is the smallest prime `>= 2^b` of the form `c*2^K + 1` (for `b = 16` this is
`65537`). Curves are drawn from a ChaCha-free splitmix64 stream seeded with
`20261120 ^ b` as `(a, b)` uniform in `F_p^2`; a curve is kept when the order
found by baby-step giant-step on a random point is prime (`n` prime and
`nP = O` with `n` in the Hasse interval proves `#E = n`, since `n > 4 sqrt p`).
The first two such curves per `p` are used (curve 0 and curve 1). Targets are
uniform random points (random `x`, square root, random sign), from a stream
seeded per cell.

## 1. Yield model (tested form)

Expected number of unordered `+-` decompositions of a uniform target:

    lambda_m = |F|^m / (m! * N)          (|F| = points, both signs)

which equals `2^m C(f, m) / N` to leading order (distinct classes; repetitions
are lower order). The tested prediction is

    Pr[R decomposes] = 1 - exp(-lambda_m)      (Poisson form; = lambda_m for small lambda)

Primary statistic per cell: `r = p_hat / (1 - exp(-lambda_m))`, with a 99%
Wilson interval for `p_hat` carried to `r`. Secondary: mean number of
decompositions per target against `lambda_m`.

**Exact methods only, never Groebner.**
* `m = 2`, exact census: when `2 f^2 <= 2^26`, compute the set
  `S = {x(P_i + P_j), x(P_i - P_j) : i <= j}`; then `Pr = 2|S| / (N - 1)`
  exactly (each `x` in `S` is two points of the sum set; `O` excluded).
* `m = 2`, sampling (all cells, and the only method when the census is too big):
  for each target `R` and every `Q in F`, test `x(R - Q) in t*H_k` by hash
  lookup (`R - Q in F` iff its `x` lies in `t*H_k`).
* `m = 3`, sampling with meet in the middle: build `T2 = {x(P_i +- P_j)}`
  (capped at `2^24` entries), then for every `Q in F` test `x(R - Q) in T2`.

**Cells.** `m = 2`: every `k <= K` with `lambda_2 in [2^-10, 2]`. `m = 3`:
every `k <= K` with `lambda_3 in [2^-10, 2]` and `|T2| <= 2^24`. Both curves
for the subgroup arm; curve 0 for the coset arm. **Sample size:** targets until
400 decomposing targets or a work cap of `2^30` element tests per cell or
`10^6` targets, at least 200 targets.

**Controls (m = 2, census and sampling).** Negative control: a factor base of
`f` uniformly random classes of the same size (no subgroup structure); it must
match the model. Positive control: an arithmetic-progression base
`{+-jG : 1 <= j <= f}`, whose sums collide (`|S| ~ 4f`, not `~2f^2`); the
detector must report `r << 1` there. A/A: the two curves at the same `(p, k)`
estimate the noise floor.

## 2. Linear-algebra model

Matrix `M`: one row per relation, one column per factor-base class (`f`
unknowns), entries `+-1` (or `+-2`, `+-3` on repetition), at most `m` nonzeros
per row, modulo the prime `N`. Relations come from `R_j = a_j G + b_j Q`
(`Q` the target), so one left-kernel vector `lambda` (`lambda^T M = 0`) with
`sum lambda_j b_j != 0` gives `log_G Q = -sum lambda_j a_j / sum lambda_j b_j`.
Only **one** dependency is needed for one logarithm.

* Generic budget (the model the task states): `f + slack` relations, scalar
  Wiedemann on the full matrix, cost `~ m |F|^2` multiply-adds mod `N`.
* What is run: relations are collected until a dependency is certain, the
  matrix is peeled to its 2-core (structured Gaussian elimination's first
  step), and scalar Wiedemann finds a kernel vector of the core.
  Predicted multiply-adds `C_LA in [1, 10] * m * n_core^2`.
* **Relations needed (prediction).** `m = 2` rows are edges of a signed graph
  on the `f` classes; a dependency is a balanced cycle or two unbalanced cycles
  in one component, so the first dependency arrives near the giant-component
  threshold: median `r*/f in [0.40, 0.55]`. `m = 3`: the 3-XORSAT / 2-core
  dependency threshold `0.918`: median `r*/f in [0.85, 0.97]`. The generic
  budget `f + slack` therefore overstates `m = 2` by about `2x`.

## 3. Relation search, end to end

* Walk `R_{j+1} = R_j + S_{h(R_j)}` over 32 precomputed steps
  `S_t = c_t G + d_t Q` (one group addition per attempt).
* `m = 2`: pair-sum table `T2` (`~2 f^2` group additions to build) and one
  lookup of `x(R_j)` per attempt; a hit is resolved by recomputing the pair
  and comparing points.
* `m = 3`: the same `T2` and, per attempt, `|F|` group additions and `|F|`
  lookups (meet in the middle).
* **Groebner route, recorded per attempt (not run):** columns of the degree-`D`
  Macaulay matrix in `2k` chain variables, `M_D = sum_{d<=D} C(2k, d)`, at
  `D` = certified LFD (`4` for `k <= 5`, `5` for `6..9`, `6` for `10..13`,
  `7` for `14`), the conjecture `floor((k+14)/4)` beyond, and the ceiling
  `k + 2`; and the paper's per-failed-attempt lower bound
  `C(2k, floor((k-10)/14) + 1)`. For `m = 3` only the ceiling `3k + 1`.
* **Cells.** `b in {20, 24, 28, 32, 36, 40}` (40 if runtime allows),
  `m in {2, 3}`, curve 0, `k in {k*-1, k*, k*+1}` where `k*` minimises the
  predicted total (`2f^2 + N/(4f)` for `m = 2`; `2f^2 + 1.4 N / f` for
  `m = 3`) subject to `|T2| <= 2^25` (`m = 2`) / `2^24` (`m = 3`).
  Three verified targets per cell. A run counts only if `xG == Q` holds.

## 4. Pollard rho baseline (measured, not assumed)

r-adding walk, `r = 32`, distinguished points (`theta = 2^-8`), restart after
`20/theta` steps without a distinguished point, coefficients recovered by
replaying the two colliding walks (replay steps are charged). Two variants:
**with the negation map** (canonical `+-` representative, one-step look-ahead
against fruitless 2-cycles; expected `sqrt(pi N / 4)`) — the **primary
baseline** — and without (`sqrt(pi N / 2)`). 400 trials per variant per
end-to-end curve; the e2e targets are the first three trials. Prediction:
mean / `sqrt(pi N/4)` in `[0.9, 1.25]` with the negation map and mean /
`sqrt(pi N/2)` in `[0.9, 1.15]` without.

## 5. Charging and the GAE band

Counted separately in both arms: group additions (doublings included),
hash lookups and inserts, multiply-adds mod `N`, and bookkeeping operations
(walk index and distinguished-point tests, peeling steps, union-find steps).
A group addition is 1 GAE. Every other counted operation is priced at
`0`, `0.1` and `1` GAE. Field multiplications are recorded inside the group
additions and not charged twice. Two arms are ranked only where the order
holds at all three prices.

## 6. Predictions and falsification criteria

* **Y1 (yield model).** A cell falsifies the model when the 99% interval for
  `r` lies entirely outside `[0.9, 1.1]`. The model fails for PKM bases if
  this happens in at least 10% of PKM cells (subgroup or coset) with at least
  100 expected hits; otherwise it holds.
* **Y2 (no trend).** Least-squares slope of `log r` on `k` (and on `log2 p`)
  across PKM cells: prediction `|slope| < 0.02` per unit.
* **Y3 (controls).** Random-base cells satisfy Y1; the arithmetic-progression
  base gives `r < 0.1` at every cell with `f >= 16`. If the positive control
  is not detected, the detector is broken and Y1 is void.
* **L1 (relations needed).** Medians of `r*/f` in the ranges of Section 2.
* **L2 (LA cost).** Measured `C_LA / (m n_core^2) in [1, 10]`.
* **R1 (rho).** Rho ranges of Section 4.
* **C1 (headline).** No end-to-end cell, `m = 2` or `3`, beats rho at any GAE
  price; IC/rho grows with `N` roughly like `N^{1/6}`
  (relation search `~N^{2/3}` against `~N^{1/2}`). A cell with IC/rho `< 1` at
  all three prices would falsify C1.

## What this does not decide

Small sizes only (`N < 2^41`); no Groebner relation search is run; `m >= 4`
is out of scope; the solving-degree growth for `m >= 3` is handled elsewhere.
