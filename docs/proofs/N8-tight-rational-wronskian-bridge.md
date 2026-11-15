# N8: exact characteristic-zero rational Wronskian bridge

**Evidence as of 2026-11-13:** the argument below is a mathematical proof
plan, not an end-to-end Lean theorem. Pinned Lean checks for the
`GrowthWronskianTightDegree` lemmas are tracked in stacked PR #176.
In particular, nothing here establishes `GrowthSpread.SpreadBound` in
the kernel until the explicitly enumerated bridge lemmas are proved.
No new axiom or `CurveModel` assumption is proposed.

## The particular curve, and no function-field geometry needed

Fix `a=6` and `c=a+2=8` over `K=ℂ`. The curve is

```
E: y² = x³ + 6x - 7,    t=1.
```

Its discriminant factor `(a+3)²(4a+3)=9²·27=2187` is nonzero.
The two rational functions in `F=RatFunc ℂ=ℂ(X)` are

```
x = X,    z = (X+8)/(X-1).
```

`GrowthRatFuncModel.ratFunc_algebraic_data` already proves their
nonvanishing and the specialized `S₃(x,z,1)=0`; equivalently
`S₃(x,z,1) = ((x-1)z-x-8)²`. For the model take `F₀=ℂ`, `p=0`.
The complex unit-ideal result on PR #170 applies because `‖6‖>5`.

For `u,v ∈ ℕ`, write `f_(u,v)=X^u z^v`. The coordinate functions
are nonzero. The relevant pole data can be expressed using valuations
of rational functions at infinity and at `X=1`:

```
ord_∞ f_(u,v) = -u,
ord_1 f_(u,v) = -v.
```

The first formula holds because `z` has degree zero (tends to 1
at infinity). The second holds because `X+8` is nonzero at `X=1`.
At infinity, a nonzero polynomial `A` has order `-natDegree A`;
at 1 it has nonnegative integer order at most its degree.

## Logarithmic derivative, as a *genuine* RatFunc identity

Put

```
D = X(X+8)(X-1),
Q_(u,v) = u(X+8)(X-1) - 9vX,
P_(u,v),0 = 1,
P_(u,v),s+1 = D P_(u,v),s' + (Q_(u,v) - sD') P_(u,v),s.
```

For ordinary differentiation `δ=d/dX` on `RatFunc ℂ`, the
identities to formalize are

```
δ(f_(u,v)) = f_(u,v) Q_(u,v) / D,
δ^s(f_(u,v)) = f_(u,v) P_(u,v),s / D^s.
```

The second follows by induction and the quotient/product rules.
These are rational-function identities: merely defining polynomial
`P` does not establish them. The existing
`GrowthRationalLogNumerator.logNumerator_eq` and
`GrowthRationalWronskianRec.numeratorRec` supply their polynomial
components. `D` is nonzero in the rational-function field.

For any cubic `D` and quadratic `Q`, the recurrence admits the
**tighter degree budget**

```
natDegree P_s ≤ 2s.
```

The base `P_1=Q` is degree ≤2. For `s≥1`, the derivative has
degree at most `2s−1`; multiplying by cubic `D` gives degree ≤
`2s+2`. In the other summand, `Q−sD'` has degree ≤2, and
multiplying by `P_s` again gives ≤`2s+2`. **Separating `s=0`
is essential:** the inequality `natDegree(P_0')≤natDegree(P_0)−1`
is true using natural-number truncated subtraction, but it does not
by itself yield the needed derivative improvement.

The new Lean module in PR #176 is designed to formalize this bound,
and obtains an `r×r` polynomial determinant degree ≤

```
B_r = ∑_(0≤s<r) 2s = r(r−1).
```

## Nonzero multi-column Wronskian and the cofactor identity

Let `J` be **any** finite circuit of the functions `f_e` over
`ℂ`, of cardinality `m`. Because the functions are all nonzero,
`m≥2`. In fact `m≥3`: two distinct `f_e` cannot be proportional
over constants, as their orders at infinity and 1 distinguish both
exponent coordinates.

Choose an enumeration `e_0,…,e_(m−1)`, and let `r=m−1`. By
minimality there are **nonzero** constants `λ_i∈ℂ` such that
`∑_i λ_i f_(e_i)=0`; every `r`-element subfamily is independent
over `ℂ`.

For every omitted index `i` form the `r×r` Wronskian
`W_i=det[δ^s f_(e_j)]_(j≠i, 0≤s<r)` using fixed row/column
conventions. The characteristic-zero Wronskian theorem guarantees
`W_i ≠0` for a linearly independent family over the field of
constants of the derivation. Here `ker δ=ℂ` inside `ℂ(X)`.

**Nonvanishing needs a real formal proof.** An economical standard
induction uses the invertible Wronskian of the last `r−1` functions:
if the full determinant were zero, a relation
`f_1 +∑_(j≥2)c_j f_j=0` between all derivative rows would exist.
Differentiating the first `r−1` rows and subtracting the next rows
shows `δ(c_j)=0` for every `j≥2`. Since `ker δ=ℂ`, this
contradicts constant linear independence. To prove `ker δ=ℂ`,
one can reduce a fraction `g/h` using gcd, then use
`g'h-gh'=0`, coprimality and the strict degree drop for derivatives.

Each Wronskian factors as

```
W_i = D^{-r(r−1)/2} (∏_(j≠i) f_(e_j)) A_i,
A_i = det[P_(e_j),s]_(j≠i, 0≤s<r) ∈ ℂ[X].
```

As `W_i≠0`, every `A_i≠0`, and the degree lemma yields
`natDegree A_i≤B_r`. Multilinearity of determinants applied to the
minimal relation gives `W_i/W_j ∈ ℂ^*`. Taking the ratio of
the displayed factorizations cancels all common factors and yields

```
f_(e_i)/f_(e_j) = κ_(i,j) A_i/A_j
    for some κ_(i,j) ∈ ℂ^*.
```

This is the **general** cofactor comparison, not just a regression
for the four-term example.

## Infinity and X=1 give the exact required circuit constant

At infinity and at `X=1`, valuation of the cofactor ratio implies

```
  e_i,1 - e_j,1 = natDegree(A_i) - natDegree(A_j),
-(e_i,2 - e_j,2) = ord_1(A_i) - ord_1(A_j).
```

Both degrees are between `0` and `B_r`; for each nonzero
polynomial, `0 ≤ ord_1(A_i) ≤ natDegree(A_i) ≤ B_r`. Therefore

```
|e_i,1 - e_j,1| ≤ B_r,      |e_i,2 - e_j,2| ≤ B_r,
2|e_i,1-e_j,1| + 2|e_i,2-e_j,2|
  ≤ 4r(r−1) ≤ 6r(r−1) = 6(m−1)(m−2).
```

For `GrowthSpread.SpreadBound F₀ ξ 0`, choose the existential
integers `z₁=z₂=0`. This establishes its needed bound for **every
circuit**, once the formal derivative/nonvanishing/valuation bridges
are in place. It uses a stricter constant than the current N8 target
and requires neither a new axiom nor a generic all-curves theorem.

## Formalization dependencies and final composition

A minimal dependency graph, using declarations already in the
repository where possible:

1. Implement `RatFunc.derivative` identification (and its constant
   kernel in characteristic zero), then prove the polynomial
   numerator recurrence computes each genuine higher derivative.
2. Prove the finite general Wronskian determinant is nonzero for
   independent rational functions, with no assumptions about an
   individual circuit beyond its existing `Circuit` predicate.
3. Derive arbitrary-index cofactor ratios from the linear relation
   and factorization; use `GrowthWronskianTightDegree` to bound each
   cofactor polynomial.
4. Prove the two stated valuation equalities at infinity and
   `X=1` for nonzero rational functions; combine with the degree
   budget to prove `GrowthSpread.SpreadBound ℂ
   (GrowthRatFuncModel.xi 6) 0`.
5. Assemble the three `GrowthMain.CurveModel` fields from the
   existing `ratFunc_algebraic_data` and new spread theorem.
   Apply `GrowthMain.growth_theorem` or
   `GrowthMain.growth_lfd_gt_linear` with `p=0`, and use the
   independently verified `GrowthCharZeroSeparatedFamily.unit_complex_of_norm_gt`
   on PR #170. The upper bound comes from
   `SelectorExtension.lfd_le_depth_add_two_of_two_ne`.

This would yield for the **fixed complex family** `E_6` and
`k≥24`: no last-fall-degree bound at most `(k−10)/14`;
combined with the `k+2` upper bound, the last fall degree grows
linearly. No generic prime-field ECDLP runtime improvement is implied.

**Status boundary:** A proof of the polynomial degree recurrence,
even once kernel-audited, proves neither `SpreadBound` nor
`CurveModel`. The exact missing mathematical/Lean interface is the
genuine higher `RatFunc` derivative, nonvanishing `r`-column
Wronskian, cofactor ratio, and two pole orders. This file records an
explicit derivation to implement, not a disguised hypothesis.
