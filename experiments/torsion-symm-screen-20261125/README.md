# Matched screen: two-torsion symmetrized summation systems

**Status:** finite-field screen complete; no end-to-end relation-finding win established.

This is an instance of the published torsion-equivariant and invariant-coordinate approach of [Faugère, Huot, Joux, Renault and Vitse](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf). The study uses the paper repository's `csd/ic_dreg.py` group law and `S_3` polynomial, with a fresh coordinate derivation and exact finite checks. It does not claim a new summation polynomial family.

## Frozen matched fixture

For `p = 97, 193, 257`, use `E: y² = x³ + x` and the multiplicative subgroup `H = μ_(2^k)` with `k = 3`; additionally use `k = 4` for `p = 193, 257`. The factor-base abscissas are **exactly** `H \ {1}` in every arm. The omitted `x = 1` is the pole of the affine coordinate below. For each cell, the targets are the first 12 rational points with nonzero ordinate in ascending `x`, using the first square root returned by the frozen group-law code. There are 60 target attempts, including 50 with no scalar `S_3` pair and 10 with at least one pair. The same target, curve, and factor-base membership are used in all comparisons.

The curve has rational `T = (0,0)` of order two. Translation by `T` sends `x` to `1/x`, which preserves `H \ {1}`. Set `u = (x+1)/(x-1)`, so this translation sends `u` to `-u`. For two summands let `e = u+v`, `r = uv`, `q = e²`, and `t = x(R)`. Clearing the two Möbius denominators in `S_3` gives, up to the nonzero scalar 4,

`P(q,r,t) = (1+t²)q - 4t²r - 2t(r²-1)`.

The original `S_3` has degree four in the two free abscissas; `P` has degree two in `(e,r)` or `(q,r)`. The factor-base equations retain a higher degree and **are included** in every closure and Gröbner computation. Specifically `F(u) = (u+1)^N - (u-1)^N` vanishes precisely on the transformed factor-base abscissas, where `N=2^k`. The permutation-symmetric arm imposes `F(u)+F(v)=0` and `F(u)F(v)=0` in `(e,r)`. Since `F` is odd, the first equation is `e H(e²,r)=0`; the torsion quotient imposes `q H(q,r)=0` and the second equation in `(q,r)`. These equations are equivalent over the algebraic closure. The script exhaustively verifies the quotient factor-base image for every frozen field and compares the relation equation on each target. Every retained point decomposition is checked with the elliptic-curve group law.

## Results

| Field, subgroup | Targets without a pair | Exact unit-ideal LFD: direct / permutation / torsion quotient | Unordered abscissa pairs admitting a point decomposition across 12 targets |
|---|---:|---|---:|
| `p=97`, `N=8` | 11 | `9 / 8 / 7` (all 11) | 2 |
| `p=193`, `N=8` | 11 | `9 / 8 / 7` (all 11) | 2 |
| `p=257`, `N=8` | 10 | `9 / 8 / 7` (all 10) | 4 |
| `p=193`, `N=16` | 8 | `17 / 16 / 15` (seven); `17 / 15 / 15` (one) | 7 |
| `p=257`, `N=16` | 10 | `17 / 16 / 15` (all 10) | 4 |

The exact unit-ideal LFD is the first degree whose iterated closure contains `1`; every lower degree was closed without `1`. `screen.py` computes all 50 thresholds. `verify_sage.py` independently constructs dense finite-field spaces and checks both sides of the threshold for representative `N=8` and `N=16` cases; its output is `sage-verification.txt`. The 10 targets with decompositions are positive controls for relation preservation; this record makes no exact LFD claim for those systems. The transformed systems have different presentations, so their LFD values can differ even though their solution sets correspond.

The `N=16` same-target solver timing uses SymPy Gröbner bases over `GF(p)` with grevlex, five repetitions per target and rotated arm order. Each timing includes polynomial construction and the Gröbner computation. The medians below pool all 12 targets per field, including failed attempts; wall and CPU agree within the displayed precision.

| Field | Direct GB | Permutation GB | Torsion quotient GB | Exact pair enumeration and group-law check |
|---|---:|---:|---:|---:|
| `p=193` | 39.2 ms | 4.5 ms | 5.1 ms | 0.121 ms |
| `p=257` | 39.6 ms | 4.5 ms | 5.2 ms | 0.118 ms |

Raw wall and CPU records are in `groebner-timing.json` and `enumeration-timing.json`. The timing records use different repetition counts (five for Gröbner, fifty for enumeration). Their gap is much larger than observed run-to-run variation, but these small fields are a screening fixture, not a scaling estimate.

**Decision:** the torsion quotient lowers the exact degree needed to refute failed attempts, and the invariant presentations make the measured Gröbner step faster on these fixtures. The coordinate map is bijective on the fixed factor base, so it does not increase the number of useful point decompositions. Direct enumeration is faster at these sizes. The torsion quotient's lower LFD also does not translate into the fastest Gröbner timing: permutation symmetry alone is faster here. The proposed combined gain in *relation yield and total relation-finding cost* is therefore unproved; these results support a larger matched crossover test rather than a paper speedup claim.

The timed Gröbner step excludes root extraction and downstream relation rank and linear algebra. We did not measure a same-target ECDLP comparison, and no result here should be extrapolated to ECC2K-130. The published [Galbraith--Gebregiyorgis characteristic-two study](https://www.math.auckland.ac.nz/~sgal018/shishay.pdf) is a separate binary construction with its own factor-base and solver tradeoffs.

## Replay

From the repository root:

```sh
python3 experiments/torsion-symm-screen-20261125/screen.py
sage -python -u experiments/torsion-symm-screen-20261125/verify_sage.py
python3 experiments/torsion-symm-screen-20261125/bench_groebner.py
python3 experiments/torsion-symm-screen-20261125/bench_enumeration.py
```

`results.json` binds the screen and group-law source hashes. The timing scripts regenerate their JSON records; timings naturally vary by host and load.

## Equal-size binary Koblitz follow-up

The [frozen equal-67 Koblitz panel](PROTOCOL_EQUAL67_20261125.md) compares
distinct `x` and torsion-invariant `u` bases of equal rational size on the
same 24 targets. Its [checked result](EQUAL67_RESULT_20261125.md) finds
12 admitted targets for the direct base and 6 for the symmetrized base.
On the 60 factor-base columns, the candidate's checked equation rank is
12 versus 9 for the control, so it passes the frozen rank-yield criterion
despite lower target admission. Exact Boolean LFD is 3 in both arms on each
of the first four targets, so it does not pass the joint yield-and-degree
gate. The result includes complete raw rows, exact
closure profiles, a separately retained environment-error attempt, and an
independent group-law/subgroup-log verifier. It is separate from the
fixed-base ECC2K-130 recipe and the prime-field subgroup systems.
