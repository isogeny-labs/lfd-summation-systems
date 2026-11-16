# Two constructed valuations for the non-two-torsion secant model

This continues PR #20 using Mathlib's standard coordinate ring
`A = K[X,Y]/(W(X,Y))` and function field `Frac(A)`. Both distinguished
valuations are now constructed on this actual field. An isomorphism with
PR #20's quadratic presentation over `RatFunc K` has not been formalized;
the standard presentation independently supplies the curve equation,
nonzero coordinates and S3 identity needed here.

## Checked result

For `W: y²=x³+ax+b`, target `R=(t,Y)` with `Y²=t³+at+b`, `2≠0`,
`Y≠0`, and `t≠0`, put `z=x(R-P)` using the existing secant expression.
The two actual, distinct valuations give the following ordinary orders:

| Function | At infinity O | At target R |
| --- | ---: | ---: |
| x | -2 | 0 |
| y | -3 | — |
| z | 0 | -2 |
| x^i z^j / (x^i' z^j') | -2(i-i') | -2(j-j') |

The Lean statements use multiplicative values `ofAdd(n)` with ordinary
order `-n`: poles have positive exponents. Uniformizers are `x/y` at O
and `x-t` at R; each has value `ofAdd(-1)`.

At infinity, the coordinate-ring norm degree is the weighted degree
`max(2 deg p, 2 deg q+3)` for `p+qY`. Its addition inequality and
multiplicativity produce a valuation, which is extended to the fraction
field. For `t≠0`, the secant numerator and denominator both have pole
order four, so z has order zero at O.

At R, the evaluation kernel is the point ideal. Mathlib's product identity
for that ideal and the opposite-point ideal becomes a principal ideal
after localization: the opposite-point ideal contains an element with
nonzero evaluation `Y-negY(t,Y)`. Thus x-t generates the maximal ideal.
The local ring is noetherian and is not a field, giving an actual DVR.
Its normalized adic valuation extends to the original function field.
The secant numerator evaluates to `4Y²≠0`, while `(x-t)²` has order two.

Only local nonsingularity is needed for these identities. Global
nonsingularity is still required for the intended genus-one argument.
"Non-two-torsion" means `Y≠0` in characteristic other than two, not
infinite order of a finite-field point. The monomial-ratio identities here
also require `t≠0`; the excluded case is not silently covered.

## Evidence and limits

`formal/N8NonTwoTorsionPolesAudit.lean` prints the exact statements and
transitive axioms of 42 new theorems. The dedicated workflow retains the
earlier 16 theorem audits and adds this exact 42-name audit. Local evidence,
including tested source commit, hashes, successful kernel output and
clearly labeled failed elaboration logs, is retained under
`experiments/non-two-torsion-poles-20261116/local-v1/`.

The proof uses pinned Lean 4.19.0 and Mathlib c44e0c8. The allowed
foundational axioms are `propext`, `Classical.choice` and `Quot.sound`;
no project axioms, placeholders or native evaluation are accepted by the
audit checker. The targeted producer build is not a claim that every
default library module has been rebuilt. No numerical study or ECDLP
performance result ran.

Sources read: the existing GrowthSecant, GrowthNonTorsionFunctionField,
GrowthNonTorsionBridge and GrowthPoleArithmetic modules; pinned Mathlib's
EllipticCurve/Group, Norm/Basic, Valuation/ExtendToLocalization,
Localization/AtPrime, Localization/LocalizationLocalization,
DiscreteValuationRing/TFAE and DedekindDomain/AdicValuation. These are
established constructions; no mathematical novelty is claimed.

## Remaining obligation

The two raw monomial pole identities are discharged. The full
`GrowthPoleArithmetic.spread_from_genus_one_height` application still
needs the genus-one height/Wronskian estimate for circuit-supported
relations, plus the coefficient field/scalar tower handling p-th powers
in positive characteristic and the link to the circuit-spread interface.
Neither the general growth theorem nor an ECDLP speedup follows yet.

Failed elaborations involved namespace/coercion normalization, explicit
field-parameter hypotheses, and an attempted use of a field-only secant
helper in the coordinate ring. The latter was replaced by the same
polynomial numerator written directly in the ring. No theorem hypotheses
were weakened to make a failed proof pass. Failed logs are excluded from
the successful axiom audit.
