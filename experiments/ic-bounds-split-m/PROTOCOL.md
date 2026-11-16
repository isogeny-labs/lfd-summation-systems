# Protocol: exact last fall degree of the split presentation with m summands

Frozen 2026-11-16, before any cell with m >= 4 or (m, k) = (3, 5) was run with this engine.
The only runs before the freeze are those of `control.json` (m = 3 at k = 2, 3 against recorded
cases, and five planted negative controls at m = 3, 4, cap 7).

## Question

Does the engine of this repository reproduce, with its own certificates, the stage B readings of
the second instrument (`experiments/crossrepo-pkm-stage-b/`) for the split presentation with
m summands and m - 2 free auxiliaries?

## System and instrument

`csd/ic_bounds_split_m.py`: inputs `S_3(x_1, x_2, z_1)`, `S_3(z_{j-1}, x_{j+1}, z_j)`,
`S_3(z_{m-2}, x_m, t)` in normal form modulo the chains of the subgroup of order `2^k`; closure
`V_D` by `csd.ic_bounds_m3_batched.BatchedClosure` on the monomial space `SplitSpace`.
Exact value of a unit instance: least `D` with `1 in V_D`, every smaller closure verified closed
(all variable shifts of rows of degree `< D` reduce to zero; every input in the span).

## Instances

Per cell `(m, k)`: primes 7681, 12289, 40961 (each `= 1 mod 2^k` for `k <= 9`), seeds 0 and 1, a
seeded random nonsingular curve and a random rational target, redrawn until `S_{m+1}` has no
zero on the grid `H^m` (enumeration by `csd.ic_bounds_msum.prime_grid`, independent of the
closure). Six instances per cell where the cell is feasible; fewer is reported as fewer.
These are subgroup bases; stage B used coset bases at primes 786433 and 2013265921. The
instances differ from stage B's; the comparison is of values per cell, not of systems.

## Predictions (the stage B readings)

| m | k | predicted exact value |
|---|---|---|
| 4 | 2 | 6 |
| 4 | 3 | 7 |
| 4 | 4 | 7 |
| 5 | 2 | 6 |
| 5 | 3 | 8 |
| 6 | 2 | 7 |
| 3 | 5 | 7 |

Stage B's `m = 5, 6` at `k = 2` used a base of four points; here the base is the subgroup of
order 4 as well, so those cells are comparable.

## Reading rule

A cell **confirms** if every completed unit instance has the predicted exact value, **refutes**
if any completed instance has a different exact value, and is **undecided** if no instance
completes (resource cap). A lower bound from a capped ladder is reported as a lower bound.
No cell is dropped: every launched instance is in the record file, including capped ones.

## Controls (`control.json`)

- Positive: 12 recorded unit cases of `experiments/ic-bounds-m3` (six at `k = 2`, six at `k = 3`,
  chosen by a seeded shuffle): generator lists equal to `csd.ic_bounds_m3.split_generators`, exact
  value equal, chain-quotient rank equal at every degree.
- Negative: 5 planted targets (proper ideals): no closure through degree 7 contains 1, above the
  unit values 5 and 6 of those cells, so the engine distinguishes the two cases.
- Agreement of the two seeds and three primes within a cell is the repeat measure.

## Not claimed

No solving cost, no statement about cells that do not complete, no statement about the
2-isogeny arm of stage B.
