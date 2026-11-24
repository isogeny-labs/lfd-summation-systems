# An explicit principally polarized representative of the ECC2K-130 class

2026-11-24, PR #57. This constructs a polarized ordinary Deligne module
in the required F2-isogeny class. It is a constructive refinement of
Howe's existence result, not a construction of a curve.

Reproduce with:

```sh
python experiments/ecc2k130-weil-factor/construct_ppav.py
```

The saved artifact `experiments/ecc2k130-weil-factor/explicit_ppav.json`
contains the lattice generator, an explicit positive-unit expression,
its 260 integer coefficients, and the data reconstructing the principal
alternating form. The output explicitly says that Jacobian status is
not determined.

## The polarized module

Let alpha²+alpha+2=0, beta=-1-alpha, zeta a primitive 131st root,
and O=Z[alpha,zeta]. The trace-zero representative has the ordinary
Deligne module O with F=alpha*zeta and V=beta*zeta^-1. To see the
integral lattice identification, write the augmentation lattice of
Z^131 using the 130 consecutive differences e_j-e_(j-1), with indices
modulo 131. The cyclic shift has the zeta action, and this is a basis
over Z. Tensoring with the base elliptic module Z[alpha] gives O.
Here the base module is Z[alpha] because Q(sqrt(-7)) has class number one.

In the ordered basis zeta^j, alpha*zeta^j, j=0,...,129, the natural
restricted product polarization has alternating matrix

    Omega0 = [[0,G],[-G,0]],

where G is the Cartan matrix A_130, with diagonal 2 and adjacent
entries -1. Its determinant is 131².

Put a=263*ell=38531015900842053623, ell=146505763881528721, and

    m = (alpha*zeta-beta*zeta^-1)/(zeta-zeta^-1).

At the 130 real embeddings, indexed by zeta -> exp(2*pi*i/131)
with the upper-half-plane embedding of alpha fixed, its values are

    m_i = -1/2 + sqrt(7)/2*cot(2*pi*i/131).

Exactly the indices 26,...,65 and 91,...,130 are negative.
The script proves these signs with rational Taylor enclosures,
including rational bounds on pi from Machin's formula. No floating
point signs are used. The uncorrected m is a failed positivity control.

## An explicit unit corrects every sign

Take theta=exp(2*pi*i/917). The field K embeds in Q(theta) via

    alpha=theta^131+theta^262+theta^524, zeta=theta^7.

The subgroup fixing K has exponents {1,394,263}. For every k prime
to 917 let

    eta_k = (theta^k-theta^-k)/(theta-theta^-1),
    epsilon_k = product_(h in {1,394,263}) eta_k(theta^h).

These are real integral units in K: the eta_k are cyclotomic units,
and the product is the norm to K of a real element. Their 130 sign
vectors are computed by integer residues of khb modulo 917. The
sign matrix has rank 129 over F2. Gaussian elimination gives the
explicit 69-element list recorded in the JSON such that

    u = product epsilon_k

has precisely the same signs as m. Hence u*m is totally positive.
The polynomial product is expanded modulo Phi_917 and converted
to the O-basis. The resulting coefficients have at most 388 bits.
The script verifies membership in Z[alpha,zeta], rather than merely
asserting that the norm expression lands in K.

Set Omega=Omega0*mult_(u*m) on O. It is an integral alternating form,
as verified entry by entry. The code also verifies F^t*Omega*F=2*Omega.
It is positive by the certified signs and the CM polarization criterion.

Its determinant is **exactly** a². Indeed

    N_(K+/Q)(m) = a/131,

by multiplying the displayed expression for m over the nontrivial
131st roots: the numerator product is
(alpha^131-beta^131)/(alpha-beta)=a and the denominator product is 131.
The unit has absolute norm one. Therefore

    det(Omega) = 131²*N_(K/Q)(u*m) = a².

This determinant assertion uses the exact norm identity. It is not
reported as a direct dense determinant computation.

## The quotient lattice and principal form

For p=263 and p=ell, the script constructs an idempotent e_p in O/pO
supported on the part of ker(Omega) where F=V is scalar. It checks
e_p²=e_p and F*e_p=V*e_p=s_p*e_p modulo p, with

    s_263=110, s_ell=14975904886999252, s_p²=2 mod p.

At 263, multiply the spectral projector by (alpha-123)/(139-123)
to select one of the split O-primes. At ell take the rational line
in the quadratic residue field. Both lines are maximal isotropic
in their two-dimensional alternating kernel factors.

Combine their generators into

    w = ell*e_263 + 263*e_ell, considered modulo a.

After multiplying by a unit modulo a, the script obtains a vector
v with first coordinate 1. Define the rank-260 lattice

    I = O + Z*v/a.

A basis is the standard O-basis with its first vector replaced by
v/a. Its covolume relative to O is 1/a; I/O is cyclic of order a.
The script checks F*v=V*v=s*v modulo a, with s the CRT combination
of the two scalars. Thus F and V preserve I. Ordinary behavior at 2
is unchanged because a is odd.

The entries of Omega*v are divisible by a, so the form is integral
on I. In the new basis the form is obtained by replacing its first
column by Omega*v/a, its first row by the negative transpose, and
the first diagonal entry by zero. It remains positive and alternating.
Its determinant is now

    det(Omega_I)=det(Omega)/a²=1.

Thus (I,F,V,Omega_I) is an explicit principally polarized ordinary
Deligne module. Deligne's equivalence with ordinary abelian varieties,
and its polarization description, give a principally polarized B/F2
with the required Weil polynomial. The inclusion O -> I gives an
isogeny A -> B of degree a.

At 263 the new lattice is O-stable. At ell, the new line in F_(ell²)
is only F_ell-linear; its local multiplier order is Z_ell+ell*O_local.
Elsewhere the order is maximal. Consequently End_F2(B) has index ell
in O. This explicitly realizes the loss of integral cyclotomic symmetry.

## What has and has not been verified

The exact checks cover the sign enclosures, unit expansion, integrality,
alternation, Frobenius compatibility, quotient stability and determinant
norm/index identity. The matrix SHA-256 is recorded with the decimal-row
serialization convention, so all 260-by-260 entries can be reproduced.
The construction is an algebraic module description; no equations for
B, no curve equation, and no executable pointwise cryptanalytic transfer
are supplied.

An exploratory direct dense determinant check was stopped after about
84 seconds without a result; a general Hermite-normal-form attempt was
also stopped after about 63 seconds. Neither is claimed as a successful
verification. The final implementation instead uses the exact norm/index
identity and the explicit cyclic lattice basis, so neither computation
is needed.

The companion obstruction proves B is geometrically indecomposable and
that a curve in the class would have to be nonhyperelliptic. These do
not imply that B is a Jacobian in dimension 130.

For the remaining Jacobian test, two alternatives were investigated:

- ThetAV's documented curve/theta conversion concerns genus-2
  hyperelliptic curves. It does not provide a genus-130 recognition
  operation for this module. See its
  [documentation](https://thetav.readthedocs.io/en/latest/doc/index.html).
- The numerical Fay-trisecant approach takes a complex Riemann matrix
  and reports numerical evidence, with examples through genus seven.
  It does not certify the special fiber over F2. In particular, failure
  of a chosen characteristic-zero lift to be a Jacobian would not exclude
  its special fiber from being one. See
  [de Leon--Frauendiener--Klein](https://arxiv.org/abs/2303.15249).

No valid characteristic-two Torelli-locus test has been implemented for
this 130-dimensional input. A direct full-coordinate method using
|2*Theta| would have 2^130 sections, about 1.36e39 coordinates; the
module certificate does not supply a compressed theta representation
that avoids materializing that space. This is a concrete size barrier
for that route, not a proof that every possible recognition method
has that cost. Nor would excluding this one B exclude
every polarized member of the isogeny class. The requested isogeny-class
Jacobian verdict remains unresolved; this construction supplies an
explicit polarized candidate for further work.
