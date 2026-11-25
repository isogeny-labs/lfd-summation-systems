# Frozen binary Boolean torsion follow-up — 2026-11-25

**Decision:** the tested two-torsion invariant quotient does not lower exact
LFD on the frozen `ic-curve-sweep` systems. At the same factor base and
targets, it cannot change useful relation yield. This is a negative result
for this presentation, not a statement about every torsion-symmetrized base.

## Frozen workload and boundary

The [protocol](PROTOCOL_FROZEN_BOOLEAN_20261125.md) fixes the existing
degree-nine Koblitz curve `y²+xy=x³+x²+1` over `F_(2^9)` with modulus
`0x203`, basis `[302,252,1]`, eight `x` values, and three existing targets.
The rational point `(0,1)` has order two. It sends nonzero `x` to `1/x`;
the seven nonzero values in the frozen base are closed under this action.
The eighth value, `x=0`, translates to the point at infinity. The first
pre-degree truth-table check found two valid planted-target decompositions
using `x=0`, so the quotient equation must retain an explicit boundary
branch. The amendment was committed before measuring degree.

For `u=x1+x2`, `v=x1*x2` and target abscissa `z`, the direct equation is
`S3=(v+zu)²+zv+1`. On `v!=0`, the simultaneous-inversion invariants
`s=v+1/v` and `w=u²/v` give `S3/v=s+z²w+z`. The tested six-variable
presentation uses that quotient equation on `v!=0` and the original `S3`
on `v=0`. It is boundary-patched, not a globally invariant coordinate
system on this frozen affine base. All 64 Boolean assignments were checked:
the direct and patched systems have identical zero sets for every target.

## Exact degree and independent relations

The repository's exact iterated Boolean closure engine supplied LFD values.
An independently coded dense 64-column closure replay matched every rank
through each certification degree. For the proper ideal, it also checked
evaluation rank and codimension against the complete solution list. The
solutions were independently enumerated from `S3` over the base.

| Existing case | Target class | Solutions | Direct input degree / LFD | Patched input degree / LFD | Checked ordered decompositions | Distinct relation rows |
|:--|:--|--:|:--:|:--:|--:|--:|
| `1bac99e1e73e08b00999` | random rational | 0 | 2 / **2** | 6 / **6** | 0 | 0 |
| `a5eb99758876c7f2fee6` | random rational | 0 | 2 / **2** | 6 / **6** | 0 | 0 |
| `7c65f50a88d2c5e980b2` | planted rational decomposition | 12 | 2 / **3** | 6 / **6** | 12 | 7 |

The curve has 518 rational points, or `14·37`. The seven distinct checked
rows from the planted target have rank seven modulo 37; Sage independently
confirmed that matrix rank. The two polynomial presentations use the same
point set and factor base, so their checked relation rows are literally
identical. Coordinate symmetrization alone cannot increase relation yield.
The 12 checked ordered decompositions agree with the frozen case's recorded
count; the enumerator avoids duplicate sign choices for the order-two point.

## Costs, limits and replay

The JSON records polynomial-construction, closure, and complete eight-base
pair enumeration wall and CPU for every target, including the two failed
targets. One saved run, under shared host load, spent about 0.0026 s wall
on relation enumeration and rank work. The planted-target closure used about
0.0011 s wall for the direct system and 0.0493 s for the patched system;
the latter CPU time was about 0.0065 s, illustrating why these short shared-
host timings are not a qualified speed comparison. The matched result is the
exact degree and relation count, not a total runtime ratio. Root extraction,
complete relation collection and matrix solve, setup amortization, and
matched cold rho are not measured here.

- `frozen-boolean-torsion.json`: equations, source hashes, all solutions,
  boundary pairs, checked rows, ranks and per-case costs; SHA-256
  `e6b72dcdb0edf279708773a9e81ff9f71427e854b45d70e4003c0ade8611ef58`.
- `frozen-boolean-verification.json`: independent dense closure profiles;
  SHA-256 `5bea53862518511d6c07cbba4dc31bb38193b47e8fef14d7d4f18822b2aa77b5`.

```sh
python3 experiments/torsion-symm-screen-20261125/frozen_boolean_torsion.py
python3 experiments/torsion-symm-screen-20261125/verify_frozen_boolean.py
```

The order-four binary Edwards screen remains in
[BINARY_RESULT_20261125.md](BINARY_RESULT_20261125.md). Neither this
two-torsion frozen-base follow-up nor that equal-total-base colored design
passes the requested joint yield-and-degree gate. A different equal-size
torsion-stable base on the *same target stream* would be a separate test;
the present result does not prejudge it.
