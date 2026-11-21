# Prime-field growth beyond depth fourteen and four summands

Status: work in progress, PR #45. The higher-summand growth theorem is not
proved. The transfer statements below have passed the local Lean elaborator;
their CI validation is tracked on the PR. No new exact closure is asserted.

## What is already proved at arbitrary depth

`GrowthRationalPrimeNative.not_one_of_large_prime` excludes 1 from the
degree-D closure of the normalized two-summand input on

    E_a: y^2 = x^3 + a*x - 1 - a, target T = (1,0)

when D >= 1, k >= 14D+10, p >= 2*2^k+1, and a+3 != 0 in F_p.
This is an arbitrary-depth theorem, not an extrapolation from depth fourteen.
To deduce LFD > D from this exclusion one also needs the unit-ideal property.
`GrowthRationalPrimeCutoff.rational_prime_growth` supplies a complete window
for its rational parameter family above its explicit cutoff, subject to
2^k dividing p-1 and the listed rational exclusions.

These results do not prove the observed sharp staircase with steps at
k=2,6,10,14,18. They also do not themselves cover S_{m+1} for m>2.

## Degree-controlled specialization: formal result

Let phi be an affine polynomial substitution from a source input F to a
target input G. Fix a preprocessing cost C. Suppose every source input f
whose degree is at most D maps into V_{G,max(D,C)}. Then

    phi(V_{F,D}) is contained in V_{G,max(D,C)}.

Proof: induct on the source derivation. Inputs use the supplied derivations;
addition and scalar multiplication commute with phi. For a product h*f,
the affine substitution does not increase either degree, so its product
budget remains at most D, hence at most max(D,C). The preprocessing cost
is paid once, rather than added at each multiplication step.

Consequently exclusion of 1 from the target space implies exclusion from
V_{F,D}. If the source ideal contains 1, this rules out an LFD bound D.
Unit status of G does not imply unit status of F: that is a separate gate.

`SummandGrowthTransfer.lean` proves these three statements, and
`prime_growth_transfer` instantiates them with the native prime-field
two-summand theorem. Its explicit depth condition is

    k >= 14*max(D,C)+10,  p >= 2*2^k+1.

The source can have arbitrarily many variables, but the theorem explicitly
requires the specialized-input derivations. It is not an unconditional
theorem about elliptic summation polynomials with m>4.

## Order-two specialization and a degree-drop trap

The point T belongs to the coordinate factor base because its abscissa is
1. Fixing an extra summand to T is the affine substitution setting every
variable in its squaring chain to 1. Every equation in that chain becomes
zero. Repeating this for m-2 blocks leaves two free chains.

The universal Semaev polynomial must be specialized *after* forming its
resultant with the prescribed degrees. One must not recompute a resultant
after a leading coefficient vanishes. For two quadratics A*z^2+B*z+C and
d*z^2+e*z+f the formal-degree-(2,2) Sylvester formula is

    (A*f-C*d)^2 - (A*e-B*d)*(B*f-C*e).

Putting d=e=0 gives A^2*f^2. The ordinary resultant computed with the
second polynomial now regarded as a constant would instead be f^2.

For E_a, S3(1,1,z)=(a+3)^2 and the leading coefficient in z of S3(x,y,z)
is (x-y)^2. Thus specialization of the universal S4 gives

    S4(x,y,1,1) = (a+3)^4*(x-y)^4.

The exact substituted Sylvester-formula identity is proved by `ring` in
`order_two_pair_specialization`. Independently, symbolic expansion of the
universal quadratic resultant before substitution gives the same expression.
Substituting first and asking a CAS to compute a reduced-degree resultant
gives only (a+3)^4, losing the diagonal factor. That shortcut is invalid.

This already shows why fixing extra summands to T is not a uniform
same-target reduction: with one extra T, the remaining two summands must
sum to the identity, giving the diagonal. With two extra T's their total
is the identity and the original target is retained. The parity matters.

For odd m, the full system with target T has the rational solution with
every summand equal to T: m*T=T. Every chain variable is 1. Therefore
this family cannot supply unit-ideal instances for odd m. This does not
exclude growing LFD for its proper ideals, but exclusion of 1 alone cannot
prove it. It also does not rule out other targets or families for odd m.

## Remaining proof obligations

1. Establish the universal specialization identity at arbitrary m, retaining
   formal degrees and all leading-coefficient factors, before normalizing
   into the squaring-chain presentation. The scalar S4 calculation alone
   does not establish that identity for m=5 or higher.
2. Give a preprocessing degree C_m independent of k for an admissible
   specialization, and prove the specialized-input derivations required by
   `prime_growth_transfer`. Ideal divisibility without a degree-controlled
   derivation is insufficient.
3. Prove unit status of the higher-summand source family, or supply an ideal
   element witnessing an LFD lower bound for a proper ideal. For odd m and
   target T, a unit-ideal claim is contradicted by the explicit grid point.
4. Preserve the original prime-field and subgroup conditions. Extending a
   theorem over a function field does not by itself give prime-field instances.
5. The sharp slope 1/4 remains a separate obligation. Even successful
   transfer of the existing theorem initially gives the lower slope 1/14.

The prior depth-fourteen certificate, split-D8, and depth-eighteen computations
are not marked complete by these transport lemmas.
