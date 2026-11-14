# N8 rational family: elementary polynomial Wronskian reduction

**Status (2026-11-12): complete written algebraic reduction to standard
polynomial/differential-field facts; NOT a Lean theorem.**  This route avoids
general function-field genus, canonical divisors, and Brownawell–Masser.
The verification script checks a nontrivial exact instance; it is not a
substitute for `GrowthSpread.SpreadBound`.

## Statement

Let `K` be a field of characteristic zero or characteristic `p>0`,
let `a+3 != 0`, and write `c=a+2`. In `F=K(X)` set

    xi1 = X,    xi2 = Z = (X+c)/(X-1).

For characteristic zero take `C=K`. For characteristic p take
`C=K(X^p)`. In both cases the kernel of `d/dX` acting on `F` is `C`.
For a minimal `C`-linear circuit `J` of distinct monomials
`f_e=X^(e_0)*Z^(e_1)`, with `m=|J|>=2`, put
`M=(m-1)(m-2)`. For any `e,h` in the circuit,

    ||e_0-h_0||_p <= M,     ||e_1-h_1||_p <= M,

where `||n||_0=|n|` and `||n||_p=min_{z in Z}|n-pz|`.
Consequently the EXACT interface required by `GrowthSpread.SpreadBound`
holds, with `4M` in place of `6M`:

    exists z1 z2 : Z,
      |2*(e_0-h_0)-p*z1| + |2*(e_1-h_1)-p*z2| <= 4*M <= 6*M.

This is for the genuine rational two-torsion component, NOT every
elliptic curve. In characteristic zero it applies to the specific
nonsingular choice `K=Q, a=1/2, b=-3/2, t=1`, independently of
the separate pending group-algebra unit-ideal bridge.

## Polynomial derivation (no valuation or genus machinery)

Write `f_(u,v)=X^u*(X+c)^v/(X-1)^v` and set

    D = X*(X+c)*(X-1)                          (degree 3)
    Q_(u,v) = u*(X+c)*(X-1) + v*X*(X-1) - v*X*(X+c)

so `deg Q_(u,v) <= 2` and `f'/f=Q/D`. Define polynomials

    P_(u,v),0 = 1
    P_(u,v),s+1 = D*P_(u,v),s' + (Q_(u,v) - s*D')*P_(u,v),s.

Induction gives

    (d/dX)^s f_(u,v) = f_(u,v)*P_(u,v),s / D^s
    deg P_(u,v),s <= 2*s.

For a minimal `C`-relation, choose ALL nonzero coefficients
`phi_e in C` with `sum_e phi_e*f_e=0`; put `g_e=phi_e*f_e`.
Every `phi_e'` is zero, so the same derivative identity applies to
`g_e`. Let `r=m-1`. For each deleted index `j`, take the
`r x r` Wronskian `W_j` of all `g_e` except `g_j`, using
derivative orders `0,...,r-1`. It factors as

    W_j = (product_(i != j) g_i) * A_j / D^(r*(r-1)/2),

where `A_j` is the determinant of the polynomial matrix with entries
`P_(u_i,v_i),s` (`i != j`, `s=0,...,r-1`). By the degree bound,

    deg A_j <= sum_(s=0)^(r-1) 2*s = r*(r-1) = M.

**Nonvanishing:** Each deleted subfamily is independent over `C`.
For any differential field whose constant subfield is `C`, a
`C`-independent tuple has nonzero ordinary Wronskian. One elementary
proof is induction on the tuple length: divide by the first function,
expand the first Wronskian column (now 1), and apply induction to the
derivatives of the remaining quotients. Those derivatives are
`C`-independent, since a vanishing `C`-combination would integrate
to a constant and contradict original independence. This argument
works in characteristic p as well as characteristic zero. Therefore
`W_j != 0` and `A_j != 0`.

The single constant-coefficient relation `sum g_e=0` shows by
alternating multilinearity that all `W_j` agree up to sign. Cancelling
common factors gives, for any `j,h`,

    g_j/g_h = (+/-) * A_j/A_h,
    f_j/f_h = (+/-) * (phi_h/phi_j) * A_j/A_h.    (*)

**Infinity:** the rational-function integer degree of `f_j/f_h`
is `u_j-u_h`. For `A_j/A_h` it is `deg A_j-deg A_h`, whose
absolute value is at most `M`. The integer degree of
`phi_h/phi_j` is 0 in characteristic zero and divisible by p
in characteristic p (`phi_i in K(X^p)`). Apply (*) to obtain the
first residue bound.

**X=1:** because `a+3 !=0`, the exact order of
`f_j/f_h` at `X=1` is `-(v_j-v_h)` (there is no numerator
cancellation). For a nonzero polynomial `A`, its multiplicity
`mult_(X-1)(A)` is between 0 and `deg A <= M`.
Thus the order of `A_j/A_h` is a difference of two integers in
`[0,M]`. Again the order of `phi_h/phi_j` is 0 in
characteristic zero and divisible by p in characteristic p, since
`(X^p-1)=(X-1)^p`. This gives the second residue bound.

Multiplying the two residue bounds by 2 and summing yields the claimed
`SpreadBound` constant `4M`.

## Differential-constants lemma (concrete formal target)

For a reduced rational function `f/g` with coprime polynomials,
`(f/g)'=0` implies `f'g-fg'=0`. Coprimality forces
`f | f'` and `g | g'`; strict degree drop for a nonzero derivative
forces `f'=g'=0`. A polynomial has derivative zero exactly when it
is constant in characteristic zero, or a polynomial in `X^p` in
characteristic p. Thus `ker(d/dX)=K` or `K(X^p)`, respectively.
No perfectness assumption on `K` is needed. The pinned Mathlib
`Polynomial.wronskian` library handles pairs only; the required
multi-column matrix determinant and this kernel characterization
must be developed separately.

## Lean integration order

1. Prove the elementary polynomial recurrence and `natDegree`
   estimate. This is uniform over fields and does not need `RatFunc`
   derivatives yet.
2. Construct `RatFunc` differentiation and prove its constants
   kernel is the existing concrete `GrowthRationalFrobeniusBase.base`
   (or `K` at characteristic zero).
3. Prove the independent-tuple Wronskian lemma and determinant
   factorization; derive (*) in the minimal-circuit case.
4. Use polynomial integer degree and `(X-1)` multiplicity to obtain
   the two residue bounds and instantiate
   `GrowthSpread.SpreadBound` with no analytic hypothesis.
5. Combine with `GrowthRatFuncModel.ratFunc_algebraic_data`,
   `GrowthRationalChainUnit.unit_of_primitiveRoot_frobenius`
   in characteristic p or the separate characteristic-zero determinant
   and quotient bridge, then call `GrowthMain.growth_theorem`.

**Do not advertise an unconditional Lean result until every one of these
steps builds and the `#print axioms` audit passes.**
## Easier Lean bound: no derivative degree drop required

For the actual `SpreadBound` interface, the stronger constant `4M`
is unnecessary. Bounding the derivative by
`natDegree(P') <= natDegree(P)` (without proving the strict drop)
gives the simpler recursion estimate

    natDegree(P_(u,v),s) <= 3*s.

The same determinant argument then gives
`natDegree(A_j) <= 3*r*(r-1)/2`. The two distinguished pole
residue bounds each cost at most this number; after doubling both
coordinate differences the sum is at most

    4*(3*r*(r-1)/2) = 6*(m-1)*(m-2).

This *exactly matches* the existing `GrowthSpread.SpreadBound`
constant and avoids a degree-zero derivative exception in Lean.
`GrowthRationalWronskianRec.numeratorRec_natDegree_le` implements
this polynomial recurrence bound independently of the as-yet unproved
rational-function derivative and matrix-cofactor bridges. Its status
must be determined by the pinned CI run, not by this note.
