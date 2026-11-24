# ECC2K-130: the required large-prime change of representative

2026-11-24. Continuation of PR #57 at the user's request to address
Jacobians elsewhere in the isogeny class. This note proves a necessary
condition on every principally polarized representative, hence on every
potential Jacobian. It does **not** prove that a Jacobian exists or that
none exists.

Put ell=146505763881528721. For the concrete trace-zero kernel A, we obtain:

**Proposition.** If B/F2 admits a principal polarization over F2 and
f:A -> B is an F2-isogeny, then v_ell(deg f) is odd. In particular,
ell divides deg f. No maximal-endomorphism-order member of this
isogeny class has a principal polarization over F2.

This explains why an attempt to recover the integral symmetry from a
principal polarization fails: a principally polarized member must change
the integral structure at ell. Merely changing A at 263 cannot suffice.
The positive principal-polarizability result from Howe remains intact.

The arithmetic inputs are reproduced by
`python experiments/ecc2k130-weil-factor/local_orders.py`, with output
`local_orders.json` beside that script. The proof also uses polarization
theory, Tate modules, and global Artin reciprocity; those arguments have
not been formalized in a proof assistant.

## The Frobenius order is maximal away from two primes

Keep alpha²+alpha+2=0, beta=-1-alpha, pi=alpha*zeta, V=beta*zeta^-1,
K=Q(alpha,zeta), O=Z[alpha,zeta], and R=Z[pi,V]. Write

    alpha^131 = a*alpha+b,
    a = 263*ell, b = 8123678690673902378.

Then alpha=(pi^131-b)/a belongs to R[1/a]. Also

    beta^130 = 1+alpha*t, t in Z[alpha],
    zeta = V^130-t*pi.

The first identity holds because beta=-1-alpha and 130 is even.
Consequently R[1/a]=O[1/a], so O/R is supported only at the primes
dividing a: R tensor Z_p = O tensor Z_p for
every rational prime p other than 263 and ell. This includes p=2.
The script checks the exact integer coefficients of t.

For reference, direct discriminant identities give

    disc(K) = 7^130*131^258,
    [O:Z[pi]] = 2^8385*a^129,
    [O:R] = a^129.

One way to check the power-order formula is to write
Q(X)=(X^131-alpha^131)(X^131-beta^131)=P_E(X)P_A(X).
Here Res(P_E,P_A)=131²*2^130*a². Applying the product-discriminant
identity gives the stated result. Passing from Z[pi] to R changes
only its 2-part, since 2 is a unit at every odd prime and V=2/pi.
The preceding reconstruction shows that R is maximal at 2.

## The relative CM extension is unramified at finite primes

Let K+ be the maximal real subfield. The two factors Q(sqrt(-7))
and Q(zeta_131) have coprime discriminants. Complex conjugation in
their product Galois group is c=(c_M,c_L). At 7 the inertia group
lies in the M-factor, and at 131 it lies in the L-factor; neither
contains c. Since these are the only ramified rational primes,
K/K+ is unramified at every finite prime.

At 131, the field Q(sqrt(-7)) is inert and the cyclotomic factor
is totally ramified. There is a unique prime p0 of K+ above 131,
with residue field F_131, and it is inert in K/K+. The prime of
K above it has residue field F_(131²).

At 263 both factors split completely: 263=1 mod 131 and
(-7/263)=1. Thus every prime of K+ above 263 splits in K.

At ell, Frobenius is c because ell=-1 mod 131 and (-7/ell)=-1.
Thus ell splits into 130 degree-one primes in K+, all of them
inert in K/K+. These facts are exact congruence/Legendre-symbol checks.

Primes above 2 split in K/K+: their decomposition group lies entirely
in the cyclotomic factor, since 2 splits in Q(sqrt(-7)), so it does
not contain c. Thus every inert prime used in the Tate-module length
argument below has residue characteristic different from 2.

## The polarization-kernel parity proof

The product principal polarization on Res(E) restricts to a polarization
lambda0 on A with kernel E[131]. Indeed A is the orthogonal complement
of the diagonal E, and their intersection is the diagonal E[131].
Over the algebraic closure its matrix on E^130 is I+J, whose determinant
is 131. Thus deg lambda0=131².

The base-curve alpha action and the cyclic zeta action are integral on A.
They generate O, so End_F2(A)=O. As an O-module the kernel of lambda0
is supported at the unique prime over 131 and has length one there:
zeta acts as 1 on the diagonal, and alpha acts irreducibly on E[131]
because (-7/131)=-1.

Suppose B and f satisfy the proposition, and let lambdaB be a principal
polarization. The pullback lambda=f^vee lambdaB f on A has the form

    lambda=lambda0*u, with u in K+ totally positive.

This uses the Rosati description of polarizations on the CM field K.
For a prime q of K+ inert in K, the length over the corresponding
O-prime of ker(lambda) is

    v_q(u) + delta_(q,p0).

At every inert q not above ell, the kernel of f at that prime is
O-stable: R is locally O, except over 263 where all primes split and
there are no inert primes to consider. Since lambdaB is principal,
ker(lambda) is an extension of ker(f) by ker(f^vee). Cartier duality
fixes the O-prime at an inert q, so these two kernels have the same
O-length there. Therefore

    v_q(u) + delta_(q,p0) is even, for inert q not above ell.

Artin reciprocity in the finite-unramified quadratic extension K/K+
now applies. The ideal (u) is generated by a totally positive element,
so its Artin symbol is trivial. Split primes contribute trivially;
each inert prime contributes (-1)^(v_q(u)). The valuations away from
ell sum to an odd integer because precisely p0 has the extra delta.
Consequently

    sum_(q|ell) v_q(u) is odd.

Finally, the degree identity is

    (deg f)^2 = 131²*N_(K/Q)(u)
              = 131²*N_(K+/Q)(u)^2,
    deg f = 131*N_(K+/Q)(u).

All 130 primes of K+ over ell have residue degree one, so

    v_ell(deg f) = sum_(q|ell) v_q(u),

which is odd. This proves the divisibility and parity assertion.
The argument permits negative valuations of u; the local length of
the actual kernel is still nonnegative, and reciprocity uses ideal
valuations with either sign.

If End_F2(B)=O, every local kernel of f above ell is an O-module.
Its residue fields all have cardinality ell², so v_ell(deg f) is even,
contradicting the assertion just proved. This proves the last part.

The reciprocity principle used here is the same maximal-order Artin
obstruction appearing in Howe's polarization theory. In contrast to
the non-principal-polarizability criterion of Howe--Maisner--Nart--
Ritzenthaler, we do not assume that no inert prime divides F-V:
the exceptional prime ell is retained throughout the proof.

## Exact local arithmetic checks

For P_A modulo the two conductor primes, the factor patterns are:

| prime | simple linear | doubled linear | simple quadratic | doubled quadratic |
|---|---:|---:|---:|---:|
| 263 | 2 | 129 | 0 | 0 |
| ell | 0 | 1 | 1 | 64 |

The certificate checks these without relying on a heuristic factorizer:
gcd(P_A,P_A') has degree 129 and is squarefree, and quotienting by
its square leaves a coprime degree-two polynomial. Gcds with X^p-X
and X^(p²)-X determine the linear/quadratic degrees. Independent full
factorization with FLINT 0.9.0 agrees, as does its integer discriminant
calculation (the discriminant has 35736 bits).

The doubled linear factor at ell is X-14975904886999252. It is the
place where the two-dimensional residue field can contract to F_ell,
allowing odd-degree changes. The other colliding factors have residue
field F_(ell²). A classification of polarized modules, rather than
just their residue factor patterns, is still required to construct B.

## Other attempts to close the Jacobian question

Every principally polarized member here is geometrically indecomposable.
All geometric endomorphisms are defined over F_(2^131): the 131st
power of Frobenius is central on the geometric isogeny E^130.
If a principal polarization decomposed geometrically, Frobenius would
permute its nontrivial polarized factors with order dividing 131.
There are at most 130 factors, so each would be defined over F2,
contradicting F2-simplicity. Thus geometric decomposability cannot
exclude the principally polarized candidates.

All hyperelliptic candidates are excluded. A hyperelliptic curve over F2
has its degree-two quotient over F2; the genus-zero quotient is P1
because a conic over a finite field has a rational point. Each of the
three rational base points has a point above it over F4: its fiber is
a degree-two zero-cycle, either a rational support point or a quadratic
point. Distinct fibers are disjoint, giving #C(F4)>=3. The target
requires #C(F4)=2. This excludes the entire hyperelliptic locus,
independently of the integral-symmetry argument.

These statements leave geometrically indecomposable, nonhyperelliptic
principally polarized varieties with a required nonmaximal order at ell.
They do not show that any such variety is a Jacobian. In particular,
ell dividing the transfer degree is not a proved lower bound of ell
operations for evaluating that transfer; no cryptanalytic complexity
claim follows from the degree assertion alone.

The original isogeny-class Jacobian problem is still unsolved by this
work. No curve or general Torelli-locus exclusion has been produced.
The successor construction has been executed: an explicit principally
polarized Deligne module, with endomorphism order of index ell, is in
[the construction note](ecc2k130-explicit-ppav.md). Its isogeny from A
has degree 263*ell, consistent with the parity assertion.
