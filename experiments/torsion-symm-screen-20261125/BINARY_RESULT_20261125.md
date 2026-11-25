# Binary torsion-equivariant matched screen — 2026-11-25

**Decision:** the tested characteristic-two invariant presentation lowers
exact scalar-coordinate unit LFD on the same failed targets. The published
symmetry-breaking factor-base idea does **not** increase useful relation yield
when its total factor-base set is held equal to the union control. The two
requested gains do not coincide in this matched design.

## Published construction and controlled comparison

[Galbraith and Gebregiyorgis](https://www.math.auckland.ac.nz/~sgal018/shishay.pdf)
use binary Edwards `E_(d,d)` and `t(P)=x(P)+y(P)`: translation by the rational
order-four point sends `t` to `t+1`. Their invariant variables for two
summands include `e=t1+t2` and
`s=(t1²+t1)(t2²+t2)`. For `d=1`, their `f3` becomes
`e²+s(t(R)²+t(R))+t(R)²`. Their symmetry-breaking construction assigns the
summands to distinct translates of a common subspace. The paper's
approximately `m!` relation-probability comparison holds each *individual*
base size fixed, so the total base becomes approximately `m` times larger.

The [frozen protocol](PROTOCOL_BINARY_20261125.md) instead holds the **total**
base constant. The baseline allows both summands from
`W=V_(l+1)`; the colored arm permits one from `A=V_l` and one from
`B=theta^l+V_l`, with `W=A disjoint-union B`. Both use the same canonical
point representatives, signs and order-four torsion corrections. The
independent point enumerator checks all target relations and records actual
matrix rows modulo the selected odd prime subgroup order. It verifies the
Edwards group law identities and `f3` on every pair of rational points in
each cell.

The colored relation set is a subset of the union baseline relation set for
**every** target, independent of any solver. Consequently its admitted-target
count and independent row rank cannot exceed the baseline's under this
equal-total-base comparison. This is a construction-level falsifier, not just
a small-sample observation.

## Measured relation yield

All three groups have order `4r`; the target stream consists of the first
`min(24,r-1)` nonzero multiples of a deterministic prime-order generator.
The small `n=5,7` cells have an empty first colored part after removing the
torsion-only and prime-projection-trivial columns; they are degenerate
controls. The `n=9` cell has two columns in `A` and four in `B`.

| `n`, `r` | Total columns, `A/B` | Targets | Admitted baseline / colored | Distinct checked rows baseline / colored | Final rank baseline / colored |
|:--|:--|--:|:--|:--|:--|
| 5, 11 | 1, 0/1 | 10 | 2 / 0 | 2 / 0 | 1 / 0 |
| 7, 29 | 1, 0/1 | 24 | 2 / 0 | 2 / 0 | 1 / 0 |
| 9, 127 | 6, 2/4 | 24 | **12 / 7** | **13 / 7** | **6 / 5** |

Sage independently recomputed the final ranks from the saved relation
vectors. The point enumeration is an exact finite check, not a solver-success
estimate.

## Exact degree on the same failed targets

For the first four baseline-failed `n=9` targets (scalars 2, 3, 4, 7), the
direct scalar ideal uses the full published `f3` and both factor-base
polynomials. The invariant scalar ideal uses `e,s` and the exact finite image
ideal of allowed pairs, derived once per base; field equations are included
in both presentations. The original degree-12 cap censored baseline-direct;
the pre-recorded extension to degree 18 gave its exact threshold. All four
targets gave the same values:

| Factor base | Direct unit LFD | Invariant unit LFD |
|:--|--:|--:|
| Union baseline | 14 | 8 |
| Colored cross-partition | 8 | 5 |

The degree is the first `D` whose **iterated** degree-`D` closure contains 1;
all lower degrees were closed without 1. An independent Sage dense-space
calculation checked the thresholds on target scalar 2 (including the
immediately lower degrees), confirmed the four ideals are unit ideals by
Gröbner basis, and checked the relation ranks. The `n=5,7` colored ideals
contain 1 already at degree zero because one color is empty; those values
are excluded from the degree comparison above.

These are exact LFD values for the **scalar-coordinate finite-image ideals
defined here**. Building that image ideal enumerates the finite factor-base
pairs and is charged as setup; it is not the Boolean Weil-descent system or
the published implementation's total solving degree. In particular, this
does not establish a lower exact LFD for the repository's frozen binary
`ic-curve-sweep` systems.

## Cost and frozen-system compatibility

The saved `n=9` point enumeration, including both arms, all 24 targets,
failed attempts, group checks and rank updates, took 19.11 seconds under
shared host load. The finite-image setup took 0.420 seconds for the baseline
and 0.074 seconds for the colored arm in the saved degree-18 replay.
Closure wall times for the four failed targets were approximately 0.047/0.008
seconds (baseline direct/invariant) and 0.030/0.007 seconds (colored
direct/invariant). These are component timings under shared host load, with
different relation yield; they do not form a favorable total-cost result.
Root extraction, a complete relation-collection matrix solve, and matched
same-target rho were not measured in this screen.

The separate frozen-system compatibility scan covers 488 distinct binary
curve/base configurations from `experiments/ic-curve-sweep/cases.jsonl`.
For a Weierstrass binary curve, translation by its rational order-two point
maps nonzero `x` to `sqrt(b)/x`, as verified against the repository group law
on sampled points of every distinct curve. Only 24 of the 488 frozen
nonzero-`x` base sets are closed under that map. The scan excludes `x=0`,
which maps to the identity; closure of the *whole* existing factor base is a
separate condition. Most frozen bases therefore cannot be used unchanged as
torsion-equivariant bases. This particular scan did not inspect the fixed
ECC2K-130 base; a [separate exact recipe screen](ECC2K130_T2_RESULT_20261125.md)
now does so. Neither screen establishes a solver-degree benefit on that
target family.

## Raw evidence and replay

- `binary-edwards-yield.json`: all target rows, relation vectors and ranks;
  SHA-256 `d8152ec376bb0a39b6d2e63a59761ef23863007e5bcb46dad5c80ceb17b79381`.
- `binary-edwards-degree.json`: original degree-12 panel, including the
  censored baseline-direct arms; SHA-256
  `3627fbf2fdd275b52787648cd63bf85072c01929609e85bf8f002c609c911465`.
- `binary-edwards-degree-extend.json`: separate degree-18 follow-up;
  SHA-256 `c7d651ce2d625bac676ff14a0ed9113f48ae77e063ecb31841189f6db818c2c8`.
- `binary-sage-verification.json`: independent ranks and dense closure;
  SHA-256 `240d0fc77881186081cdf2a0275ca7e3d2833dae3d5dc37850a15b8be2b7bfa8`.
- `binary-base-compatibility.json`: frozen-base screen and source hash;
  SHA-256 `c1ec45c8050ee7c57bac218adb5d9faa86a9dd5f3f6e5123223dfeafaf51085c`.

```sh
python3 experiments/torsion-symm-screen-20261125/binary_edwards_yield.py
sage -python -u experiments/torsion-symm-screen-20261125/binary_edwards_degree.py
sage -python -u experiments/torsion-symm-screen-20261125/binary_edwards_degree.py --cap 18 --degrees 9 --out experiments/torsion-symm-screen-20261125/binary-edwards-degree-extend.json
sage -python -u experiments/torsion-symm-screen-20261125/verify_binary_sage.py
python3 experiments/torsion-symm-screen-20261125/binary_base_compatibility.py
```

This result rejects the two-part lead for the controlled colored-vs-union
binary Edwards design. The remaining scientifically distinct question is
whether a different equal-size factor-base construction can improve useful
relation rank **and** reduce the Boolean solver's exact degree on the
repository's frozen systems, after charging setup and failures.

The [fixed ECC2K-130 recipe screen](ECC2K130_T2_RESULT_20261125.md)
found zero reciprocal weight-at-most-two coordinates among all 8,646
weight-at-most-two type-II normal-basis candidates. It also identifies the
pre-existing native `koblitz-symmetrised` implementation and its archived
complete-cost comparison, which do not pass the joint controlled-yield and
exact-degree gate. These findings delimit the direct transplant to the
configured ECC2K-130 base without claiming a full-size solve or a general
torsion no-go theorem.

## Frozen Boolean follow-up

The [separate matched frozen-system replay](FROZEN_BOOLEAN_RESULT_20261125.md)
tested a two-torsion invariant quotient on the repository's degree-nine
Koblitz Boolean systems. A boundary branch was necessary to retain the
frozen base's zero-coordinate decompositions. Its direct exact LFD values
were 2, 2, and 3; the equivalent boundary-patched presentation gave 6 on
all three targets, with identical checked relation rows. This is a negative
degree-and-yield result for that frozen-base presentation, not a complete
cost or ECC2K-130 result.
