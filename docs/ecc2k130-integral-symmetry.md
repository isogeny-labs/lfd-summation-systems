# ECC2K-130: integral symmetry and principal polarizability

2026-11-24. Continues the exact-invariant computation in PR #51 and the
point-count tests in PR #53. Reproduce the arithmetic with
`python3 experiments/ecc2k130-weil-factor/integral_symmetry.py`.
The output is `experiments/ecc2k130-weil-factor/integral_symmetry.json`.
The arguments below use established theorems; the script checks their
arithmetic hypotheses, not the theorems themselves.

There are three distinct conclusions:

1. The actual trace-zero variety, with its integral rational order-131
   automorphism, is not a genus-130 Jacobian over F2.
2. Its F2-isogeny class **does contain a principally polarized variety**.
   This follows from an explicit witness for Howe's Theorem 1.1.
3. Whether this isogeny class contains a genus-130 Jacobian remains
   unresolved. The first conclusion does not extend automatically to every
   member of the class: the order-131 action is not forced by Z[F,V].

Continuation: [the isogeny-class obstruction](ecc2k130-isogeny-class-obstruction.md)
shows that every isogeny from A to a principally polarized member has an odd
valuation at 146505763881528721 in its degree. An
[explicit principally polarized Deligne module](ecc2k130-explicit-ppav.md)
has now been constructed, with endomorphism-order index that prime.
Its Jacobian status remains undetermined.

If the letter A is used only for a Weil polynomial or an unspecified
member of its isogeny class, conclusion 1 must not be read as an
isogeny-class exclusion. It refers to the concrete kernel defined below.

## The actual trace-zero variety is excluded

Let E/F2 have equation y²+xy=x³+1, and let
A = ker(Tr: Res(F2^131/F2) E -> E). This kernel is connected: after base
change it is the kernel of the summation map E^131 -> E, hence E^130.
The cyclic permutation of the factors restricts to a rational automorphism
rho of A of exact order 131. On its rational Tate module its eigenvalues
are the nontrivial 131st roots of unity, each twice.

The invariant computation gives End_F2(A) tensor Q =
K = Q(zeta_131, sqrt(-7)). A positive Rosati involution on a CM field is
complex conjugation. Thus every F2-polarization lambda of A satisfies
rho^dagger = rho^(-1), and rho preserves lambda.

Suppose A were F2-isomorphic to Jac(C) for a genus-130 curve C/F2.
Transport the canonical principal polarization to A. By strong Torelli,
the odd-order polarized automorphism rho comes from an automorphism sigma
of C of order 131. It is rational over F2: its action on the Jacobian is
rational and the map from curve automorphisms to Jacobian automorphisms
is injective. The possible sign in Torelli causes no obstruction for an
odd-order element (in the hyperelliptic case the sign is itself induced
by the hyperelliptic involution).

The quotient by sigma is tame. If its genus is h and it has r geometric
branch points, Riemann--Hurwitz gives

    258 = 131(2h - 2) + 130r.

The only nonnegative integral solution is h=0, r=4. But a rational
order-131 action over F2 cannot have four geometric branch points.
At a fixed point, the tangent multiplier of sigma is a primitive
131st root of unity. Frobenius squares that multiplier and commutes with
sigma. A closed orbit of length m therefore requires 2^m = 1 mod 131.
Since ord_131(2)=130, every branch orbit has length divisible by 130.
Consequently r is divisible by 130, contradicting r=4.

This is a proof of nonexistence for the actual A. It also excludes any
member of its isogeny class carrying the corresponding integral
order-131 action. It does not assume the genus >=394 result from the
other repository or require a genus-3 census.

For the Torelli input see J. S. Milne, *Jacobian Varieties*,
[Section 12](https://jmilne.org/math/xnotes/JVs.pdf).

## A principal polarization exists somewhere in the isogeny class

Write alpha²+alpha+2=0, beta=-1-alpha, and pi=alpha*zeta.
The exact recurrence gives

    alpha^131 = 38531015900842053623*alpha + 8123678690673902378,
    38531015900842053623 = 263 * 146505763881528721.

Both displayed factors are prime. Put ell=146505763881528721. We have
ell=-1 mod 131 and the Legendre symbol (-7/ell)=-1.
In the quadratic field F_ell[alpha], the element beta/alpha has exact
order 131. Choose zeta=(beta/alpha)^66; then zeta has order 131 and
zeta²=beta/alpha. This defines a reduction of the full ring of integers
O_K=Z[alpha,zeta] to F_(ell²). The equality of this ring with O_K follows
from the coprime discriminants of Q(sqrt(-7)) and Q(zeta_131).

Residue Frobenius sends alpha to beta and zeta to zeta^(-1), so its
element in Gal(K/Q) is complex conjugation. The prime of K+ under this
prime is therefore inert in K/K+. Meanwhile

    pi = alpha*zeta = beta*zeta^(-1) = V,

so this inert prime divides F-V. Concretely the reduction of pi is
14975904886999252 in F_ell, whose square is 2 mod ell.

Howe's Theorem 1.1 states that a simple finite-field isogeny class with
an inert prime of K+ dividing F-V contains a principally polarized
variety. Its hypotheses are now verified. Thus principal polarizability
of the **isogeny class** is established, rather than left as an untested
obstruction. This theorem supplies existence, not an explicit curve or
an evaluated isogeny, and a principally polarized variety of dimension
130 need not be a Jacobian.

Reference: E. W. Howe, *Kernels of polarizations of abelian varieties
over finite fields*, J. Algebraic Geometry 5 (1996), 583--608,
[Theorem 1.1, author preprint](https://arxiv.org/abs/2001.05111).

## Why the symmetry does not settle the isogeny question

Every member of the class has integral F and V, hence an action of
R=Z[F,V]. Rational endomorphisms give K, but algebraic integrality in K
does not mean an element preserves every R-lattice.

Here is a direct certificate that zeta is not in R. Modulo 263 the two
roots of x²+x+2 are 139 and 123. Take a nontrivial 131st root z1 such
that z2=(139/123)*z1 is nontrivial. There are two maps O_K -> F_263:

    (alpha,zeta) -> (139,z1),
    (alpha,zeta) -> (123,z2).

They agree on F=alpha*zeta and V=beta*zeta^(-1), but disagree on zeta.
No polynomial in F,V can have these two different images. Thus
zeta is not in R. The script supplies explicit z1,z2 and checks both
images. This does not assert that any arbitrary R-lattice is principally
polarized; it identifies the missing implication in the symmetry argument.

## Point, place and effective-divisor tests in every degree

The finite checks can be completed with elementary tail bounds, rather
than leaving a cutoff at degree 400 or 1599. Write N_k for the virtual
point count. The Weil bound is

    2^k + 1 - 260*2^(k/2) <= N_k <= 2^k + 1 + 260*2^(k/2).

For d>=32, every proper divisor e of d satisfies e<=d/2. Overcounting
their number by d gives

    N_d - sum_{e|d,e<d} N_e
      >= 2^d + 1 - 260*2^(d/2)
         - d*(2^(d/2) + 1 + 260*2^(d/4)) > 0.

For the last inequality divide the subtracted terms by 2^d. At d=32
their sum is 664577/134217728 < 1, and each term decreases thereafter.
The absolute values of the Mobius coefficients are at most one, so
the same inequality proves positivity of the place count in every
degree d>=32. The script checks degrees 1 through 31 exactly.
Integrality follows in every degree from the unique Euler product
expansion of L(T)/((1-T)(1-2T)) in Z[[T]], where L(0)=1.

For effective divisors, exact checks through degree 258 suffice.
With r=P_A(1), the functional equation gives

    B_m = r*(2^(m-129) - 1),  m>=259.

Indeed L(1)=r and L(1/2)=r/2^130, and the coefficient formula for the
rational zeta function gives this identity. At m=259 the possible
last numerator term contributes zero, so the formula already holds
there. This proves nonnegative integral virtual point, place and
effective-divisor counts in **all** degrees. These necessary conditions
remain weaker than the existence of a curve.

## Remaining proof obligation

The cryptanalytic question is whether a Jacobian occurs in this isogeny
class, with an explicit usable transfer, or whether a larger-genus
Jacobian carries the factor inside the proposed cost window. Neither
question is answered by these calculations. The exact trace-zero
representative is excluded; principal polarizability of its class is
proved; the remaining genus-130 candidate must fail to carry the
corresponding integral cyclotomic action.

An isogeny-class Jacobian exclusion now requires control of those smaller
endomorphism orders, or a different obstruction. A CM-type comparison
limited to curves with a rational order-131 symmetry cannot decide that
remainder. Exhaustive enumeration over the fixed finite field is a
decision procedure in principle; no feasible enumeration at genus 130
is supplied here.

Separately, PR #53's operation estimate needs a linear-algebra audit:
N*g models one sparse matrix-vector product, whereas generic Wiedemann
or Lanczos requires order N such products. N²*g is the corresponding
first-order estimate. At N around 1.46 million and g=130 this is around
2^48 field operations, before field-operation bit costs. Relation trials,
field operations and bit operations cannot serve as interchangeable
end-to-end attack costs. This does not change any geometric conclusion.
