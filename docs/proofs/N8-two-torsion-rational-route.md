# N8: a fixed characteristic-zero unit family with a rational-function curve model

**Status:** new written reduction, *not* kernel-checked. This is an independent shortcut for establishing the existence of unbounded last-fall degree; it does **not** replace the characteristic-p, every-nonsingular-curve theorem in `N8-growth.md` or discharge `GrowthMain.CurveModel` in Lean.

## 1. An explicit fixed curve

Work over `K = Q` with

    E: y^2 = x^3 + 10*x - 11,
    t = 1,     R = (1,0).

The cubic vanishes at 1 and its short Weierstrass discriminant factor is

    4*a^3 + 27*b^2 = 4*10^3 + 27*121 = 7267 != 0.

The exact third summation polynomial factors as

    S3(x,z,1) = (x*z - x - z - 12)^2.

This is a doubled rational component. The `GrowthSecantFactorization` and
`GrowthTwoTorsion` Lean files contain algebraic proof candidates for the
factorization and the following Möbius formula.

## 2. Unit ideal at *every* depth, with the same K and E

Put N=2^k. A common zero of the normalized two-chain equations over
an algebraic closure of Q gives x^N=z^N=1 and S3(x,z,1)=0.
Embed the algebraic numbers in C. Then |x|=|z|=1, while

    S3=0  =>  x*z-x-z=12
    |x*z-x-z| <= |x*z| + |x| + |z| = 3 < 12,

a contradiction. By the weak Nullstellensatz, the corresponding finite
polynomial system generates the unit ideal over Q at every k >= 2.

**Formalization obligation:** bridge the exact normalized chain quotient
to the equations `x^N=z^N=1`; prove the complex norm inequality and use
a finite-dimensional/nilpotent-free quotient criterion or the weak
Nullstellensatz. No `UnitInstance.unit_zmod` is used: that lemma assumes
a target with nonzero ordinate.

## 3. Concrete rational-function model

Take F=Qbar(X), the rational function field over an algebraic closure of Q, and F0=Qbar. The formulas below are defined over Q and can first be verified in Q(X). Set

    xi1 = X,
    xi2 = T(X) = (X+12)/(X-1).

Both are nonzero. `S3(xi1,xi2,1)=0` identically. Moreover

    T(T(X)) = X,
    T(X)-1 = 13/(X-1).

The divisor table on P^1 is:

| place | ord(xi1) | ord(xi2) |
|---|---:|---:|
| X=0 | +1 | 0 |
| X=-12 | 0 | +1 |
| X=1 | 0 | -1 |
| infinity | -1 | 0 |

There are exactly four support places, with genus zero. In particular,
the two distinguished poles at infinity and X=1 have independent
coordinate orders. The nonzero denominator/derivative condition is 13 != 0.

## 4. A stronger spread bound in characteristic zero

For a minimal linear dependence of m monomials
`xi1^a xi2^b` over Qbar, the usual genus-zero Wronskian argument
(precisely Theorem 3.1 of N8, with genus=0 and |S|=4) yields

    H <= ((m-1)(m-2)/2) * (-2+4) = (m-1)(m-2).

The height comparison at the two distinguished poles gives

    |a-a'| + |b-b'| <= 2*H <= 2*(m-1)*(m-2).

Therefore

    |2*(a-a')| + |2*(b-b')| <= 4*(m-1)*(m-2)
                             <= 6*(m-1)*(m-2)  (m >= 2).

This is stronger than `GrowthSpread.SpreadBound Qbar xi 0` required by
`GrowthMain.CurveModel`. No elliptic function-field places or
positive-characteristic p-th-power Wronskian argument are necessary
**for this special characteristic-zero family**.

**Formalization obligation:** set up `RatFunc Qbar` (with algebra maps from Q), prove the four
order/divisor identities on P^1, and instantiate a kernel-proved
rational-function Wronskian bound with constant coefficients. Merely
postulating this genus-zero spread bound would not close the gap.

## 5. Conditional consequence and scope

Once both the rational-function spread bound and the fixed-system unit
ideal are formalized, the existing Lean `GrowthMain.growth_theorem`,
`growth_lfd_window` and upper bound would prove for this *single fixed
curve over Q*:

    floor((k-10)/14)+1 <= LFD(k) <= k+2,   k>=38.

Thus the objective `LFD(k)=Theta(k)` for a fixed characteristic-zero
curve could be proved while the substantially harder uniform
characteristic-p theorem remains open. It says nothing about ECDLP
hardness or prime-field runtime improvements.

The written derivation above was independently checked symbolically,
and the Möbius formulas passed 27,306 exact finite-field checks in
`scripts/check_n8_torsion_mobius.py`. Neither substitutes for a Lean
kernel proof.
