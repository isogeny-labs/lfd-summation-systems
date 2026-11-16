# Rational parameters: Lean classification and growth extension

The source is `formal/ICFallDegrees/GrowthRationalParameterClassification.lean`.
The module is imported by the default target and its eight headline statements
are printed by `Audit.lean`. No project axiom, placeholder, or native decision
procedure is used.

## Statements

For a rational parameter a, put `E_a: y²=x³+a*x-a-1`, target abscissa 1.

* `grid_zero_iff`: at every k>=2, a complex `2^k`-root grid zero of S3
  exists if and only if a=-3, a=-1, or a=1. The reverse direction constructs
  `(1,1)`, `(i,-i)`, and `(-1,-1)`, respectively.
* `unit_complex`: for every other rational a and every k>=2, 1 belongs to
  the ideal of the actual normalized `PrimeM2General.inputs` over C.
* `exists_lfd_linear_window`: for rational a outside
  {-3,-3/4,-1,1} and every k>=24, there exists an actual last fall degree d
  with `(k-10)/14 < d <= k+2`, where division is natural-number division.
  `lfd_linear_window` proves the same inequalities for any supplied LFD.
* `new_members` instantiates the full window at a=0,2,-2, all outside the
  old norm-five condition. `excluded_members` checks that a=-1 and a=1
  really do have grid zeros at every k>=2.

The singular parameters are -3 and -3/4. The unit theorem does not need a
nonsingularity assumption; the full LFD window excludes both because its
selector upper bound needs nonsingularity. This is a characteristic-zero
theorem for rational coefficients viewed in C. It does not prove a uniform
prime-field family or settle the non-torsion, genus-one case.

## Proof structure

From S3=0, its square factor gives `(x-1)z=x+2+a`. Taking norm squares gives
`(a+3)(a+1+2 Re(x))=0`. Unless a=-3, this identifies
`x+conj(x)=-a-1`. Both summands are integral over Z because they are roots
of unity. Injectivity of Q -> C reflects integrality to Q, and the integral
root theorem makes this rational number an integer. Its absolute value is
at most 2. Values -1 and 1 would produce cubic roots of unity (using -x for
the latter); their orders cannot also divide a power of two. The remaining
values give exactly the three classified parameters.

The unit proof uses the existing finite-grid theorem and explicitly connects
the raw S3 identity to the normalized chain inputs. The growth proof composes
that unit result with `GrowthTwoTorsionFamily.model`, which needs only a+3!=0,
the existing lower-bound reduction, and the nonsingular selector bound.
Neither CurveModel, SpreadBound, nor unit-ideal membership is an assumption
of the public growth theorem.

## Reproduction and evidence

From the repository root:

```sh
python3 scripts/verify_lean_fall_degrees.py
python3 scripts/lean_check_all.py --out /tmp/rational-parameter-all-modules
python3 scripts/check_fixed_curve_obstructions.py
python3 scripts/check_publication_artifacts.py
```

`module-build.log` is the output of elaborating the final module with pinned
Lean 4.19.0 and Mathlib `c44e0c8ee63ca166450922a373c7409c5d26b00b`.
`default-build.log` records the default-target verifier result. The corresponding
full default build, axiom output and source digests are in `experiments/ic-lean/`.
The local run reused compiled dependencies; it is not a fresh empty-cache build.
The existing pull-request CI installs the pinned dependencies, builds all modules,
and checks every printed axiom list on Linux.

The earlier exact replay and written argument remain in
`experiments/fixed-curve-validation-20261114/`. This module formalizes that
classification and its growth consequence; the resultant threshold and smaller
inverse-denominator follow-ups in that note remain written arguments.
