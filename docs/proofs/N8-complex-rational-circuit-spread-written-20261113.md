# N8 characteristic-zero rational circuit spread: complete mathematical proof

**Status:** written mathematical argument only. This file is **not** a Lean theorem, does not discharge `GrowthSpread.SpreadBound` in the kernel, and must not be imported as an axiom. Its purpose is to isolate explicit Lean proof obligations. Base: pinned and CI-successful PR #178 commit `48664cd5d2d24fc60df0a5fa7741e7ee9532f496`; all parent branches unchanged.

## Theorem (constant-field circuit spread)

Let `F = ℂ(X)`, `z = (X+8)/(X-1)`, `f_(u,v) = X^u z^v` for `(u,v) ∈ ℕ²`. Let `J ⊂ ℕ²` be a finite circuit over `ℂ` of the family `f_(u,v)`, and put `m=|J|`. For every `(u,v),(u',v')∈J`,

```
2 |u-u'| + 2 |v-v'| ≤ 4 (m-1)(m-2).
```

In particular, this implies the exact proposition `GrowthSpread.SpreadBound ℂ (GrowthRatFuncModel.xi 6) 0`, since `4≤6` and the existential multipliers `z₁,z₂∈ℤ` can both be zero.

### 1. Differential field constants (genuine derivation)

Give `F=ℂ(X)` its ordinary formal derivative `∂`. It is the unique `ℂ`-derivation with `∂X=1`. If `h=A/B∈F` with coprime `A,B∈ℂ[X]`, `B≠0`, then

```
∂h = (A' B - A B') / B².
```

If `∂h=0`, coprimality gives `A|A'` and `B|B'`: indeed, `A|A'B` and `B|AB'`. Every nonconstant polynomial in characteristic zero has a derivative of strictly smaller degree, so `A` and `B` are constants (with the zero numerator case immediate). Thus `ker ∂=ℂ`.

### 2. Wronskian nonvanishing, proved by induction

For `h₁,…,h_r∈F`, write `W(h₁,…,h_r)=det(∂^s h_j)_{0≤s<r,1≤j≤r}`. We claim

```
W(h₁,…,h_r) ≠ 0  iff  h₁,…,h_r are linearly independent over ℂ.
```

The reverse direction follows from multilinearity and the fact that derivatives annihilate constants. For the forward implication needed here, prove its contrapositive in the other direction by induction: assume `h₁,…,h_r` independent, hence `h₁≠0`, and put `g_j=h_j/h₁` for `j≥2`. Generalized Leibniz yields the **exact determinant identity**

```
W(h₁,h₂,…,h_r) = h₁^r W(1,g₂,…,g_r)
                 = h₁^r W(∂g₂,…,∂g_r).
```

The second Wronskian has size `r-1`: expansion along the first column of `W(1,g₂,…,g_r)` removes row 0 and column 1. If `∂g₂,…,∂g_r` were linearly dependent over `ℂ`, there would be nonzero constants `c₂,…,c_r` such that `∂(∑ c_j g_j)=0`; by step 1, `∑c_j g_j=c₁∈ℂ`. Multiplying by `h₁` contradicts independence of `h₁,…,h_r`. Therefore the `r-1` derivatives are independent; induction makes their Wronskian nonzero; the displayed identity makes the original Wronskian nonzero. The case `r=1` is immediate.

**Lean obligations:** define a derivation `F →ₗ[ℂ] F` with quotient rule, prove `ker ∂=ℂ`, generalized Leibniz determinant identity, and induction on `r`. None of these are discharged by a polynomial recurrence in isolation.

### 3. Actual higher derivatives and bounded polynomial cofactors

For `e=(u,v)`, let

```
D=X(X+8)(X-1),      Q_e=u(X+8)(X-1)-9vX,
P_(e,0)=1,           P_(e,s+1)=D P'_(e,s)+(Q_e-s D')P_(e,s).
```

The quotient rule gives `∂f_e=f_e Q_e/D`. By induction using the product rule,

```
∂^s f_e = f_e P_(e,s)/D^s.
```

The identity is in `F`, not merely `ℂ[X]`; `D≠0`. The kernel-checked polynomial work in PRs #176 and #178 supplies the recurrence degree estimate `deg P_(e,s) ≤ 2s` and the cleared first-derivative identity, but **the transport to a genuine rational derivation is still open in Lean**.

Let `J={e₀,…,e_(m-1)}` be a circuit, and `r=m-1`. Every subfamily of size `r` is linearly independent. Put `W_i=W(f_j : j≠i)`, so step 2 gives `W_i≠0`. Factor `f_j` from each column and `D^{-s}` from row `s`:

```
W_i = (∏_(j≠i) f_j) D^{-r(r-1)/2} A_i,
A_i = det(P_(e_j,s))_(0≤s<r, j≠i) ∈ ℂ[X].
```

Thus `A_i ≠0`; by the Leibniz expansion, every term of `A_i` has degree at most `∑_(s<r)2s=r(r-1)`, so

```
0 ≤ deg A_i ≤ B:=r(r-1).
```

The upper bound is already kernel-checked as `GrowthWronskianTightDegree.cofactor_natDegree_le_two` when the input hypotheses are matched.

### 4. Cofactor comparison (all circuits, not one example)

Minimal dependence yields a relation `∑_(i<m) λ_i f_i=0` with **all** `λ_i∈ℂ^×`. Applying `∂^s` for each `s<r` gives `r` linear equations in `m=r+1` unknowns. The corresponding `r×m` derivative matrix has rank `r`, because every maximal minor `W_i` is nonzero. Its kernel is one-dimensional, with signed maximal minors as a nonzero kernel vector. Hence for any `i,j`, `W_i/W_j∈ℂ^×` (indeed, the signed ratio equals `λ_i/λ_j`). Taking the quotient of the factorizations of step 3 gives

```
f_i/f_j = κ_(i,j) A_i/A_j,     κ_(i,j)∈ℂ^×.
```

Signs and ordering of omitted columns affect only the nonzero scalar `κ`; they do not affect valuations.

### 5. The two orders and sharp constant

For nonzero rational functions put `ord∞(A/B)=deg B - deg A`. Also `ord₁(A/B)` is the multiplicity of `X-1` in `A` minus that in `B` (using reduced fractions or the additive extension of multiplicity). Since `z=(X+8)/(X-1)` has equal numerator/denominator degrees and `X+8` evaluates to `9≠0` at `X=1`,

```
ord∞ f_(u,v) = -u,         ord₁ f_(u,v) = -v.
```

For each nonzero polynomial `A_i`, `ord∞A_i=-deg A_i` and `0≤ord₁A_i≤deg A_i≤B`. Applying `ord∞` and `ord₁` to the cofactor identity yields

```
u_i-u_j = deg A_i - deg A_j,
-(v_i-v_j) = ord₁ A_i - ord₁ A_j.
```

Each right-hand-side difference has absolute value at most `B`. Therefore

```
2|u_i-u_j| + 2|v_i-v_j| ≤ 4B
   = 4(m-1)(m-2) ≤ 6(m-1)(m-2).
```

There cannot be a one-element circuit since all `f_e≠0`. There cannot be a two-element circuit: proportional monomials would have the same orders at infinity and 1 and hence the same exponent pair, contrary to distinctness of circuit indices. Thus `m≥3` and the displayed arithmetic is valid.

### 6. Exact consequence conditional *only* on Lean integration

For the fixed curve `E:y²=x³+6x−7`, target `t=1`, the existing `GrowthRatFuncModel.ratFunc_algebraic_data` identifies nonzero coordinates `X,(X+8)/(X-1)` and the Semaev relation in `RatFunc ℂ`; the independent PR #170 supplies the complex chain unit ideal. Once steps 1–5 are **translated into Lean without extra hypotheses**, the concrete spread theorem instantiates `GrowthMain.CurveModel ℂ 6 (-7) 1 0 (GrowthRatFuncModel.xi 6)`. Existing `GrowthMain.growth_lfd_gt_linear` then gives, for `k≥24`, no LFD bound at most `⌊(k-10)/14⌋`. The already formalized selector bound gives an LFD bound `k+2`. Together this is `LFD(k)=Θ(k)` for this **fixed** complex family.

### Remaining formalization milestones (not yet kernel-checked)

1. Implement/identify the real `RatFunc ℂ` derivation, `ker ∂=ℂ`, and `∂^s f_e=f_e P_(e,s)/D^s`.
2. Prove the arbitrary finite Wronskian identity and independence criterion in Lean (step 2), for every circuit subset.
3. Prove determinant cofactor-ratio identity and both integer-valued order identities, combining the already checked cofactor degree bound.
4. Produce the **concrete, assumption-free** `GrowthSpread.SpreadBound`, the `CurveModel`, and a final `#print axioms` on the LFD growth statement; require pinned Lean 4.19 / Mathlib v4.19 full CI success.

**Evidence classification:** steps 1–6 above constitute a written mathematical proof, using the standard algebra of rational-function derivatives and determinants. A passing audit for prerequisite degree lemmas is not a passing audit of this theorem. No claimed all-curves result, characteristic-p extension, or ECDLP speedup.
