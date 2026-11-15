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
`GrowthNonTorsionRationalObstruction.lean`, implements it using integer degree
of rational functions and composes it into `generic_s3_ne_zero` (pending
kernel build). The non-two-torsion S3 component
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
| Nonsquare generic cubic in K(X) | `cubicPolynomial_nonsquare`, candidate pending build |
| No generic rational S3 section for nonzero target ordinate | `generic_s3_ne_zero`, candidate pending build |
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
object. A dependency source build is running, first with two and then with
four Lean runtime threads (the first build was interrupted cleanly to raise
parallelism). No new theorem is yet reported kernel-checked. A dedicated CI
job builds both modules and requires exactly eight transitive axiom answers,
allowing only the three standard foundations. Existing full CI is preserved.

Sources: `GrowthSecant.lean`, `GrowthSecantFactorization.lean`,
`GrowthMain.lean` and `N8-growth.md` sections 1–3. These identify the actual
`CurveModel.spread` requirement and the exceptional rational section.
