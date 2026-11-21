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


## Written extension for every fixed m: proof candidate

This section is a written proof candidate, not a Lean theorem about Semaev
polynomials. It uses the existing written two-summand growth theorem, the
standard defining zero property of Semaev polynomials, and Hasse's bound.
The normalization, universal-resultant, and counting arguments below still
need independent review and formalization before a headline claim is changed.
The earlier remaining-obligations list records the formalization gates; the
following argument proposes how to discharge the mathematical ones.

### Two universal resultant identities

Write S_n for the universal, symmetric summation polynomial with its usual
normalization, S_2(x,y)=x-y, and formal degree 2^(n-2) in each variable.
These are the polynomials defined before specializing the arguments.
The recursion with the last two arguments grouped together gives

    S_n(U,u,v) = Res_z(S_{n-1}(U,z), S_3(u,v,z)).

There are n-2 arguments in U and the first polynomial has degree
 d=2^(n-3) in z. As a polynomial in v, the coefficient of v^2 in
S_3(u,v,z) is (z-u)^2. Homogeneity of the Sylvester determinant in its
second input therefore gives

    [v^(2d)] S_n(U,u,v) = Res_z(S_{n-1}(U,z), (z-u)^2)
                             = S_{n-1}(U,u)^2.

This is an identity of universal coefficients, so it remains valid when
specialization makes that coefficient zero. Equivalently, the leading
formal coefficient of S_{n-1}(U,z) is S_{n-2}(U)^2.

For E_a and c=(a+3)^2, S_3(1,1,z)=c. The formal-degree Sylvester determinant
with a second input specialized to a constant has value LC(P)^2*c^d.
Consequently, for n>=4,

    S_n(U,1,1) = c^(2^(n-3)) * S_{n-2}(U)^4.                 (A)

For one order-two argument and t!=1, put u=(t+a+2)/(t-1). Directly,

    S_3(1,t,z) = ((t-1)*z-(t+a+2))^2,
    S_4(x,y,1,t) = (t-1)^4 * S_3(x,y,u)^2.                  (B)

The second equality follows by taking the resultant against that squared
linear polynomial. No new curve is introduced: u is the abscissa of R+T
on E_a when t is the abscissa of R.

Applying (A) repeatedly gives these exact formulas:

* m=2r+2: S_{m+1}(x,y,1,...,1,t)
      = c^(r*4^r) * S_3(x,y,t)^(4^r).
* m=2r+3, t!=1: the same specialization equals
      c^(2*r*4^r) * (t-1)^(4^(r+1)) * S_3(x,y,u)^(2*4^r).

The even scalar recurrence is separately checked in Lean by
`pairIterate_closed`; that lemma does not itself define Semaev polynomials.
Exact arithmetic using the existing formal-degree resultant engine agreed
with both formulas in 252 cases: m=2..8, primes 257 and 65537, parameters
0 and 6, (x,y)=(2,3),(1,2),(7,7), and t=3,5,11. These checks include
vanishing intermediate leading coefficients. Separately 126 checks at
t=1, including a=-3, agreed with the torsion/diagonal formulas. These are
polynomial identity checks, not closure computations or an induction proof.

### An explicit preprocessing cost independent of depth

A deliberately loose total-degree bound for S_{m+1} is

    C_m = (m-1)*2^m,  m>=2.

To check it without using a sharp degree formula: S_3 has total degree at
most 4. If S_{n-1} has total degree at most delta, its resultant with S_3
has total degree at most 2*delta+4*2^(n-3), since each determinant term uses
two coefficients of S_{n-1} and 2^(n-3) coefficients of S_3. Induction gives
S_n total degree at most (n-2)*2^(n-1).

Every binary-chain normalization step replaces a square by its successor
(or by 1) and does not raise ordinary degree. The difference between a raw
polynomial of degree at most C_m and its chain normal form is therefore
derivable from the chains at degree C_m. This holds before and after fixing
extra chains to 1, and the resulting normal form commutes with that fixing.

The right side of either specialization formula is a scalar times a power
of a raw S_3. That power has total degree at most 2^m <= C_m. Its derivation
from the normalized two-summand S_3 and the two chains stays within C_m:
first recover raw S_3 using its degree-four normalization identity, multiply
by its remaining power with the literal product-degree budget, then normalize
again. Scalars involving t are constants of the coefficient field.
Thus the specialized-input condition in the transport theorem holds at
max(D,C_m). This does not assert that the source generator has small degree
at every D: if its degree exceeds D it is not an available input, exactly
as the transport hypothesis specifies.

### Conditional lower bound for the actual direct m-summand input

Let E_a be nonsingular, t!=0,1, and t+a+2!=0. Assume

    k >= 14*max(D,C_m)+10,  D>=1,  p > 2^(k+1).

For even m use the two-summand target t; for odd m use u. Both target
abscissae are nonzero. The existing written two-summand growth theorem
excludes 1 at degree max(D,C_m). The specialization and preprocessing
argument just given then exclude 1 from the degree-D closure of the
normalized direct m-summand input. If that input is a unit ideal, its
last fall degree is greater than D. This implication does not assume or
infer source unit status from a two-dimensional slice.

### Prime-field unit instances exist for every such depth

Take E_0: y^2=x^3-1 over F_p, and H of order N=2^k. In the algebraic
closure, there are at most 2N affine points whose abscissa belongs to H.
Therefore the set of sums of m such points has size at most

    B = (2N)^m.

It is closed under negation. By the defining zero property of summation
polynomials, any rational target R outside this sum set has no zero of
S_{m+1}(H^m,x(R)). Notice that the summands here range over the algebraic
closure, not just E(F_p); restricting this count to rational summands would
not establish unit status of the polynomial system.

Exclude also O and targets of abscissa 0,1,-2, at most seven further points.
Hasse gives #E(F_p)>=p+1-2*sqrt(p)>p/2 for p>=16. Hence a suitable rational
target exists whenever

    p > max(16, 4*(B+7)).

For this target, t and (t+2)/(t-1) satisfy all guards above. The chain
quotient splits as the function algebra on H^m, since p is odd and
N divides p-1. The nonvanishing summation input is invertible in that
algebra, so the actual source input generates the unit ideal.

Combining the arguments: for each fixed m>=2, each D>=1, and each depth
k>=14*max(D,C_m)+10, every prime satisfying

    N divides p-1,
    p > max(2^(k+1), 16, 4*((2N)^m+7))

admits a rational target on E_0 whose direct m-summand system has LFD>D.
The standard chain upper bound also gives LFD<=mk+1. Primes in this
progression are unbounded, so these conditions provide prime-field
instances at arbitrarily large depths. This is an existence statement
with targets chosen after p and k; it is not a statement about every target
or the frozen computational instances.

For k>=14*C_m+10, taking D=floor((k-10)/14) yields the written candidate
window floor((k-10)/14)+1 <= LFD <= mk+1. Its lower asymptotic slope is 1/14
for every fixed m. No slope-1/4 claim follows. For example C_5=128 gives
an initial sufficient depth 1802 and C_6=320 gives 4490; this deliberately
loose preprocessing bound is far beyond the observed small-depth cells.

Source for the standard summation-polynomial properties: Igor Semaev,
[Summation polynomials and the discrete logarithm problem on elliptic curves](https://eprint.iacr.org/2004/031),
Theorem 1. The new written argument above must be reviewed separately from
that source and from the repository's existing two-summand proof.
