# Frozen screen: unramified affine degree-two quotient

## Question and scope

The existing [cover-of-E review](../ecc2k130-weil-factor/round2/cover-of-E-review/REPORT.md)
leaves non-Galois degree-131 maps in affine monodromy open. Test only the
subcase in which the geometric Galois closure has group
`C_131 semidirect C_2`, the quotient `E' -> E` has degree two over `F_2`,
and the normal `C_131`-cover `D -> E'` is geometrically connected,
unramified, and defined with constant deck group over `F_2`. The target
curve `C -> E` has genus 131. These are hypotheses, not conclusions about
every affine cover.

## Predeclared gate

1. Check `#E(F_2)=4` for `E: y^2+xy=x^3+1` by explicit point enumeration.
2. Derive the different of `C -> E` from the 131-point permutation action
   of `C_131 semidirect C_2`, with unramified normal cover. Inertia in the
   quotient `C_2` fixes one sheet and swaps 65 pairs. Check that
   `Diff(C/E)=65 Diff(E'/E)`, including residue degrees. Riemann–Hurwitz
   must then force `g(E')=3`.
3. Geometric connectedness of the unramified constant `C_131`-cover
   requires `131 | #J(E')(F_2)` by geometric class field theory. The
   degree-two map `E' -> E` makes `J(E')` isogenous over `F_2` to
   `E x Prym(E'/E)`, so `4 | #J(E')(F_2)`.
4. Compare the resulting minimum `4*131=524` with the Weil upper bound
   `#J(E')(F_2) <= (1+sqrt(2))^6 = 99+70 sqrt(2) < 198`.

If each implication is justified, record a conditional exclusion of this
specific unramified affine subcase. If any implication fails, keep it open
and record the exact countercondition. Do not infer a result for ramified
normal `C_131` covers, other complements `C_d`, non-affine monodromy, or
Jacobians without such a map. No ECDLP transfer or speed claim follows
from this screen.

## Verification

The finite-field point count and the strict integer inequality have
independent hand calculations in the result note. A Sage check may be used
as a secondary control. The group-action and class-field implications
require a written proof with their hypotheses, rather than treating a
small-field computation as their validation.
