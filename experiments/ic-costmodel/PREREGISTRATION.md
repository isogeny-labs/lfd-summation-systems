# Preregistration: a fully charged cost model of PKM prime-field index calculus against Pollard rho

Frozen 2026-11-20, before any model output was computed. Every number that
later appears in `README.md` is `[pcc:unverified]` until a pcc claim is filed
for it. This is a model study: no curve is attacked, nothing is run on a pod,
nothing is spent.

## Question

The paper (`paper/curve-lfd-bounds.tex`, and the per-failed-attempt reading of
its Section `sec:meaning` in the PR #40 draft) bounds the degree a
Macaulay-type solver needs for one relation-search attempt of the
Petit–Kosters–Messeng (PKM) two-power-subgroup index calculus. Study A
(`experiments/ic-endtoend/`) measured the relation yield, the relations
needed, the linear algebra and the end-to-end cost at `N <= 2^40`. Study B
(`/home/user/lfd-wt/msum-growth/experiments/msum-growth/`, commit `9142020`)
measured the `m = 3` and `m = 4` last fall degree against the depth and the
cost of a closure. Combined, where does the minimum fully charged cost of
the attack stand against Pollard rho for `log2 N` from 40 to 256, under the
proved bounds, under the measured law, and under the proved ceilings?

## Parameters

* `n = log2 N` from 40 to 256 (grid: every integer; reported at
  `n in {40, 64, 80, 128, 160, 256}`). `p ~ N` (Hasse).
* Summands `m = 2..8`.
* Depth `k` (factor base `F = {P : x(P) in H_k}`, `|H_k| = 2^k`) free in
  `2 <= k <= n - 2`. The model assumes `2^k | p - 1` for every `k` used, i.e.
  the prime is PKM-friendly; for a generic prime (`v2(p - 1)` small, e.g.
  `1` for P-256 and secp256k1) the attack does not apply at all. `k <= n - 2`
  keeps `p > 2^{k+1}`, the hypothesis of the growth theorem.
* `|F| = 2^k` points (Study A: measured `|F|` within a few percent of
  `2^k`), `f = |F|/2` classes.

## The model

All costs are in group-addition equivalents (GAE). A group addition (or
doubling) is 1 GAE. Every other counted operation ("native": hash inserts
and lookups, bookkeeping, linear-algebra multiply-adds, Groebner/Macaulay
`F_p` multiply-adds) is priced at `pi in {0, 0.1, 1}` GAE (AGENTS.md §5a band).
Totals are `GA + pi * native`. Work is computed in `log2` throughout.

### Yield and attempts (measured, Study A)

* `lambda = |F|^m / (m! N)`, `Pr[decompose] = 1 - exp(-lambda)` (Study A:
  `r ~ 1`, no trend, `m = 2, 3`, `p` to `2^40`; for `m >= 4` this is an
  **extrapolation**).
* Relations needed `R_m = rho_m f`: `rho_2 = 0.53` (measured; one signed
  cycle suffices), `rho_3 = 0.92` (measured, the 3-XORSAT 2-core threshold);
  `m >= 4`: the m-XORSAT thresholds `0.977, 0.992, 0.997, 0.999, 0.9996`
  (`m = 4..8`), an **assumption** carried over from the `m = 3` match.
* Attempts `A = R_m / Pr`. Each attempt: 1 walk addition and 3 bookkeeping
  operations (as charged in Study A). Setup: `94 n` additions (Study A's 66
  scalar multiplications; measured `92 n` to `95 n`).

### Linear algebra (measured, Study A)

* `m = 2`: one balanced cycle by signed union-find; `15 R_2` native
  operations (measured bookkeeping plus multiply-adds per relation, `<= 10`).
* `m >= 3`: `4.73 m (0.65 f)^2` multiply-adds (measured `LA/(m n_core^2) =
  4.73`, `n_core / f ~ 0.65` at `m = 3`) plus `115 f` peeling bookkeeping.
  For `m >= 4` the same constants are an **assumption**.
* Plus `m R_m` additions to resolve relations.

### Per-attempt decomposition: the minimum of two routes

**(a) Exact enumeration / meet in the middle (an algorithm, no degree
assumption).** For `a` from `ceil(m/2)` to `m`, `b = m - a`: a table of the
`x`-coordinates of all signed sums of `a` classes, built once,
`T_a = f^a 2^(a-1) / a!` additions and as many inserts (`a = 1`: no
additions); per attempt `S_b = (2f)^b / b!` additions and `max(S_b, 1)`
lookups (`b = 0`: one lookup). Study A ran `a = 2` for `m = 2` (`b = 0`) and
`m = 3` (`b = 1`). Memory (`16 T_a` bytes) is reported, not charged.

**(b) Groebner/Macaulay solving** with degree `D(m, k)` on the direct chain
presentation in `mk` chain variables; columns
`|M_<=D| = sum_{d <= D} C(mk, d)` (squarefree monomials). Charged in native
`F_p` multiply-adds per attempt, plus the attempt's walk addition. Three
scenarios:

| scenario | degree `D(m, k)` | solver cost per attempt | role |
|---|---|---|---|
| **LB** (proved floor) | `m = 2`: `max(4, D_growth(k))` on attempts that fail (unit ideal; Theorem `thm:three` needs `3t^2 + a != 0`, Theorem `thm:growth`: `lfd > D` for `k >= k_D = D + Lambda(D - 1)`, which implies `floor((k-10)/14) + 1`); successful attempts (fraction `1 - exp(-lambda)`) are charged 1, since no lower bound is proved for proper ideals. `m = 3`: 7 (`thm:carry`, `k >= 3`, all attempts). `m = 4`: 13 for `k >= 4` (kernel-checked carry statement; its hypotheses `deg g = 12` and `h != c g` are established by computation at `k = 4, 5` only, **extrapolated** in `k`). `m >= 5`: 2 (no proved constant; the chain equations have degree 2). | `C(mk, D)` (the top layer only, linear in the columns) | lower bound on the attack cost that the proved results allow |
| **LAW** (measured law, extrapolated) | `m = 2`: `floor((k+14)/4)` (exact for `k <= 13`, 7 at 14). `m >= 3`: `mk + 1` for `k < m - 1` (forced, measured for `m <= 8`); else `m(m-1) + 1 + floor((k+14)/4) - 4`, i.e. the onset value `deg g + 1` and the two-summand slope `1/4`. This reproduces `m = 3`: 7 at `k = 2..5`, 8 at `k = 6` (Study B, exact); `m = 4`: 9, 13, 13 at `k = 2, 3, 4` (exact). Beyond `k = 6` (`m = 3`), `k = 4` (`m = 4`) and for all `m >= 5` at `k >= m - 1` it is an **assumption**. Capped at the ceiling. | `|M_<=D|^2` (attacker-optimistic: sparse solver, no `omega > 2`, constant 1) | central estimate |
| **CEIL** (proved ceiling) | `m = 2`: `k + 2` (`thm:sd`, `sd <= k + 2`); `m >= 3`: `mk + 1` (`cor:allm`, `thm:sd`). | `7.5e-3 |M_<=D|^3` (Study B's measured constant for a dense closure run to the end) | an achievable attack: upper bound on its cost (modulo the yield model) |

Sensitivity rows (reported, not headline): LAW with Study B's dense
`4.5e-3 |M_<=D|^3`; LB with route (b) removed.

### Total and optimization

`IC(m, k, route) = 94 n + A (1 + per-attempt GA) + [route a] T_a + m R_m +
pi * (3A + per-attempt native + [route a] T_a + LA)`, minimized over route
(`a` included), `k` and `m`, separately per scenario and price.

### Rho baseline

`sqrt(pi N / 4)` additions (r-adding walk with the negation map), plus 2.25
native operations per addition (Study A's measured walk-index,
distinguished-point and negation bookkeeping), so `rho(pi) = sqrt(pi N/4)
(1 + 2.25 pi)`. Setup and the distinguished-point tail are dropped (they
favour rho at large `N` by a negligible amount).

## What counts as "IC beats rho"

At a given `n` and scenario, the minimum IC cost is below `rho` at **all
three** prices `0, 0.1, 1` (the order holds across the band). A win at a
single price is reported but is not a win. At price 0 route (b) charges
nothing for the solver, so price-0 rows are expected to favour IC
artificially; this is why the band is required.

Each reported row states the ratio `log2(IC) - log2(rho)` (the distance
remaining) and the minimizing `(m, k, route)`.

## Validation and falsification (before any extrapolation is reported)

The model, with the constants above (taken from Study A's pooled medians,
so this is an in-sample consistency check, not an out-of-sample test), is
evaluated on every end-to-end cell of Study A (`m = 2, 3`, `log2 p` 20 to 40,
33 cells) with Study A's own route (`a = 2`) and nominal `f = 2^(k-1)`,
and compared with the measured mean total at each price.

* **Pass** if every cell with `log2 p in [28, 40]` (19 cells, 57
  cell-prices) has `predicted / measured` in `[1/2, 2]`, and the median of
  `|log2(predicted/measured)|` over those cell-prices is at most `0.5`; and
  the rho model is within `[0.75, 1.33]` of the measured rho mean at
  `log2 p >= 32`.
* Cells at `log2 p = 20, 24` are reported, not scored (fixed costs).
* **Fail**: the extrapolation is still computed but every table is headed
  "falsified model", and the failure is reported as the result. A fix is
  allowed only if it is a stated accounting error (a term Study A charged and
  the model omitted, or vice versa), reported with the before/after numbers.

Secondary check: the predicted best `k` per size and the predicted
`IC/rho` of the best cell against Study A's table.

## Predictions (written before computing)

* P1: the validation passes.
* P2: under LAW and CEIL, at prices 0.1 and 1, route (a) is the minimum
  for every `n`, and IC does not beat rho at any `n` (the collision argument:
  route (a) costs at least `2 sqrt(R_m N)`).
* P3: at price 0, route (b) costs only `A` additions, so IC is below rho for
  `m >= 3` at every `n` in every scenario (an artifact of the band's lower end).
* P4: under LB, IC beats rho (all three prices) at large `n` only through
  `m >= 5`, where no degree lower bound is proved, i.e. the proved results do
  not exclude IC beating rho.

## Separation of evidence (to be kept in the report)

* Proved: the `m = 2` per-failed-attempt lower bounds (`thm:three`,
  `thm:growth`, written proof), `m = 3` `lfd >= 7`, the ceilings
  `k + 2` and `mk + 1`.
* Measured: the yield law (`m = 2, 3`), relations needed, the LA constants,
  `m = 3` degrees to `k = 6`, `m = 4` to `k = 4`, the closure-cost constant.
* Extrapolated: every degree beyond the measured `k`, every `m >= 4` yield,
  relation and LA constant, the LAW slope `1/4` for `m >= 3` beyond `k = 6`.

## Outputs

Rust crate `experiments/ic-costmodel/` (no new Python): `validate`,
`sweep`, `figure`. `results/validation.md`, `results/sweep.jsonl`,
`results/headline.md`, `results/fig-cost.svg`, `README.md`.
