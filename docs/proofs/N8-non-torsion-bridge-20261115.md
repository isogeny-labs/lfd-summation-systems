# Non-two-torsion algebraic bridge

Requested scope: advance the growth proof beyond its order-two witness
targets, retaining the original prime-field and circuit-spread obligations.

The new `GrowthNonTorsionBridge.lean` module defines

    y = ((x-t)^2 z - companion(a,b,x,t))/(2Y),

for a target R=(t,Y) with Y nonzero. The existing square-factorization
identity shows that S3(x,z,t)=0 implies y²=x³+ax+b. Away from x=t,
the recovered y and the existing secant companion coordinate are mutually
inverse. The new module also isolates a nonsquare obstruction: if the cubic
is not a square in the chosen field, no rational z in that field can solve S3.

This identifies the next genuine interface rather than reusing the E6
rational section with inappropriate hypotheses. For the generic x in K(X),
the cubic's odd degree supplies the nonsquare argument. The second new module,
`GrowthNonTorsionRationalObstruction.lean`, proves it using integer degree
of rational functions and composes it into `generic_s3_ne_zero` (locally
kernel-checked). The non-two-torsion S3 component
requires a quadratic extension of K(X), followed by the genus-one
valuation/Wronskian spread proof. Neither inverse formulas nor a nonzero
companion coordinate establish that spread bound.

## Proof map and limits

| Obligation | New declaration or remaining work |
|---|---|
| Recover the curve ordinate | `recoverY_square` |
| Inverse secant maps | `recoverY_secantX`, `secantX_recoverY` |
| Exclude rational sections under explicit nonsquare hypothesis | `no_section_of_nonsquare` |
| Nonzero companion from quadratic constant | `secantX_ne_zero` |
| Nonsquare generic cubic in K(X) | `cubicPolynomial_nonsquare`, local kernel check passed |
| No generic rational S3 section for nonzero target ordinate | `generic_s3_ne_zero`, local kernel check passed |
| Construct the quadratic function field | Open |
| Genus-one circuit spread / CurveModel | Open |
| Prime-field unit instances and full growth theorem for non-torsion targets | Open |

The eight algebraic lemmas are a partial formalization, not the requested
general growth theorem. Exact build and axiom-audit status is recorded in
the task validation addendum. No new numerical closure study was run.

## Validation checkpoint

Local Lean 4.19.0 was installed and the pinned Mathlib source resolved at
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. Its prebuilt-cache download failed
with HTTP 403; the documented alternative cache returned 404 for a sample
object. The dependency source build completed, first with two and then with
four Lean runtime threads; the build was also resumed on the combined final
target. The initial final-module elaboration failed on lifting 2≠0 into the
rational-function field. An explicit injectivity proof corrected the cast
without changing the theorem's assumptions. The final target build passed.

All eight new lemmas are locally kernel-checked, with transitive axioms only
`propext`, `Classical.choice` and `Quot.sound`. Exact statements, fresh axiom
output, dependency pins and source hashes are retained under
`experiments/non-two-torsion-bridge-20261115/local-v1/`. This is a producer-local
kernel check, not independent mathematical review. A dedicated CI job builds
both modules and requires the same eight explicit axiom answers; it was still
queued when the local check completed. Existing full CI is preserved.

Sources: `GrowthSecant.lean`, `GrowthSecantFactorization.lean`,
`GrowthMain.lean` and `N8-growth.md` sections 1–3. These identify the actual
`CurveModel.spread` requirement and the exceptional rational section.

## Exact next interface

Here “non-two-torsion” means Y≠0. It does not mean infinite order: every
finite-field point has finite order. The original finite-field objective
concerns targets outside the order-two locus.

Over a coefficient field containing Y, put f=X³+aX+b in K(X), and form
L=AdjoinRoot(T²-C(f)). The new nonsquare argument is the input for
irreducibility of this quadratic. Embed x=X into L, take y to be its
adjoined root, and define z by the existing secant formula. The next
algebraic milestone is an actual field instance and the two coordinate
functions, not an assumed `CurveModel` or spread conclusion.

For t≠0, the companion constant is a polynomial in X with leading
coefficient t², so it is nonzero. Injectivity into L and `secantX_ne_zero`
then supply nonvanishing of z; the existing secant identity supplies S3=0.
These must be proved as Lean declarations before claiming the algebraic
model complete.

The analytic milestone remains exact: for every minimal circuit prove the
two pole identities and height inequality accepted by
`GrowthPoleArithmetic.spread_from_genus_one_height`. In characteristic zero
use the constant coefficient field; in positive characteristic build the
p-th-power coefficient field and its scalar tower explicitly. The latter
cannot be replaced by independence over constants. None of these geometric
or positive-characteristic obligations is discharged by the new inverse
formulas or rational nonsquare argument.
