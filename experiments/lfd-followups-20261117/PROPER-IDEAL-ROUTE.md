# Candidate route: distinguish h from h squared with a first derivative

Status: a proposed reduction, finite falsifiers passed, two elementary
identities checked in Lean. This is NOT an asymptotic LFD theorem.

For E6 and target (1,0), let h=xy-x-y-8. The summation polynomial is
h^2. Put N=2^k. In the finite chain quotient over F_p, if N divides p-1,
Frobenius gives h^p=h, so h=h^2*h^(p-2). Thus h is an ideal element.
This observation does not require that the ideal be the unit ideal.

The finite probes find h outside V4 at (k,p)=(7,257),(7,769),(8,7681),
(9,12289), and outside V5 at (10,12289). These are explicit low-degree
witnesses for the new proper-ideal lower bounds. The closure checks are
completed, not inferred from partial spaces. The Frobenius normal-form
identity and the truncated-shift commutator are also checked in
`jet-probe.json`.

## Proposed invariant

Take the rational function field with x transcendental and
T(x)=(x+8)/(x-1). Let phi evaluate an unwrapped chain normal form at (x,T).
Let delta multiply each monomial x^a y^b by b (Euler derivative in y),
without changing its support. Consider the intersection

    phi(f)=0 and phi(delta(f))=0.

The input h^2 satisfies both conditions. The candidate witness h satisfies
the first, while phi(delta(h))=T(x)*(x-1)=x+8 is nonzero.

For a truncated shift T_j in y by s=2^j, and projection P_j onto terms
with b >= N-s, the exact commutator is

    delta T_j = T_j delta + s T_j - N T_j P_j.

For an x shift, delta commutes with the shift. These formulas have been
checked on every monomial at k=2..6. This finite check is not the all-depth
Lean specialization.

The existing growth proof gives a wrapping-part splitting property for
low-weight relations in ker(phi), conditional on its spread bound. Since
delta preserves support and binary weight, the same splitting property
should apply to delta(f). If the required shift vanishings follow, the
commutator preserves both equations, giving a stronger invariant than
ker(phi) alone. That would exclude h from the bounded closure.

## Checked elementary lemmas

`formal/ICFallDegrees/GrowthJetPrototype.lean` proves:

1. h^p=h and p>=2 imply h=h^2*h^(p-2), in any monoid;
2. the displayed commutator and three explicit scalar vanishings imply
   that the shifted element and its derivative both evaluate to zero.

The second lemma ASSUMES those vanishings. It does not establish the
analytic splitting property, instantiate the concrete chain operators,
or prove the candidate growth bound. Its name is not evidence of a
completed proper-ideal theorem.

## Remaining proof obligations

- Define delta and the wrapping projection on the precise normal-form
  space used in GrowthReduction, including the two coordinate blocks.
- Prove the commutator and degree/support preservation there for every k.
- Instantiate the three shift vanishings from the existing FEI/spread
  theorem at its exact characteristic and depth thresholds. Check the
  coefficient field and all reduction conditions; do not assume a model
  exists at an exceptional prime merely because the finite ideal is proper.
- Prove h^p=h modulo the actual chain ideal and transfer the cofactor
  identity to the ordinary polynomial presentation.
- Check that h's derivative image is nonzero in the chosen model, and
  translate its nonmembership to an LFD bound using deg(h)<=D.
- Establish which exceptional admissible primes satisfy all model
  hypotheses. A conditional invariant alone does not settle that question.

A failure of any concrete commutator, degree preservation, or required
wrapping-part vanishing would falsify this proposed bridge. The prototype
module and the finite data are kept separate from the submitted paper.
