# Constructed quadratic S3 function field

This continues the algebraic bridge in PR #17. The actual field is now
`AdjoinRoot (T^2 - C(X^3+aX+b))` over `RatFunc K`. Irreducibility follows
from the already checked nonsquare cubic and Mathlib's root criterion for
quadratic polynomials. Its finite extension degree is exactly two.

The new module `GrowthNonTorsionFunctionField.lean` constructs the generic
coordinates x and y, proves y²=x³+ax+b, proves x is nonzero and differs
from every constant t, and constructs the actual secant companion z. For
any rational target (t,Y) on the curve, S3(x,z,t)=0. If t is nonzero, the
quadratic constant coefficient has nonzero X² coefficient t²; injectivity
into the field proves that coefficient function, and hence z, nonzero.

These are algebraic results over an arbitrary coefficient field. The
geometric growth objective still requires nonsingularity and the relevant
characteristic hypotheses; this module does not silently add a CurveModel
or an assumed spread conclusion. Non-two-torsion means Y≠0, not infinite
order of a finite-field point.

## Validation and failed attempts

All eight audited new lemmas pass a local pinned Lean 4.19.0 / Mathlib
c44e0c8 build. Their transitive axiom lists contain only propext,
Classical.choice and Quot.sound. Exact output and source hashes are retained
in `experiments/non-two-torsion-field-20261116/local-v1/`.

The initial elaboration failed to turn a root equation into an equality
and to simplify the map of zero. A later attempt left a polynomial
coefficient and numeral-map normalization unresolved. Explicit coefficient
lemmas and map simplification corrected these without changing assumptions.
Failed elaboration is not kernel evidence; only the final successful audit
is supplied to the checker. No numerical experiment or ECDLP result ran.

Sources read: the existing GrowthSecant, GrowthSecantFactorization,
GrowthNonTorsionBridge, GrowthNonTorsionRationalObstruction and
GrowthPoleArithmetic modules; pinned Mathlib's Polynomial/SpecificDegree,
RingTheory/AdjoinRoot and RatFunc degree/injectivity sources. These supply
established quotient-field and secant constructions; no novelty claim.

## Remaining discriminating obligation

The earlier algebraic-model milestone is discharged for t≠0: actual field,
actual nonzero x/z and S3, with no supplied model. Next construct the
distinguished divisor valuations at O and R and prove the two pole orders
and genus-one height inequality used by
`GrowthPoleArithmetic.spread_from_genus_one_height`. In characteristic p,
the coefficient field and scalar tower must handle p-th powers explicitly.
The full circuit-spread bound and general non-two-torsion growth theorem
remain open. The new field construction alone does not prove either.
