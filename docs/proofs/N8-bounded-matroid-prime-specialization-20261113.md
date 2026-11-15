# N8 prime-field growth: bounded-box matroid specialization route

**Status (13 November 2026):** The bounded-spread ⇒ FEI kernel bridge in `GrowthBoxSpread.lean` is a candidate pending CI. The statements below about rational/prime matroid transfer and rational unit-ideal certificates are a mathematical strategy, **not yet Lean theorems**. This is an independent alternative to transporting every bounded-degree `InClosure` row-space stage.

## 1. The finite integer coefficient matrix

For the fixed E6 section `ξ₁=X`, `ξ₂=(X+8)/(X−1)`, fix depth `k` and `N=2^k`. For every exponent pair `(u,v)` with `0≤u,v<N`, multiply `ξ₁^u ξ₂^v` by the **same nonzero** denominator `(X−1)^(N−1)`. This gives the integral polynomial

```
g_(u,v)(X) = X^u (X+8)^v (X−1)^(N−1−v),
degree(g_(u,v)) ≤ 2N−2.
```

Let `A_k∈ℤ^((2N−1)×N²)` be the coefficient matrix with these columns. Common-denominator multiplication preserves all linear dependences over any coefficient field. Over every field `K`, the column matroid of `A_k ⊗ K` therefore agrees with the matroid of the bounded family `{ξ₁^u ξ₂^v | 0≤u,v<N}` in `K(X)`.

For every integer square minor `δ` of `A_k`, if `δ=0` over ℚ it remains zero in every `𝔽_p`. If `δ≠0`, it survives whenever `p∤δ`. There are finitely many minors for fixed k. Consequently there exists a finite bad-prime set `S_k` (or a natural bound `B_k` exceeding the absolute value of every nonzero minor) for which **every subset of the bounded monomial family has exactly the same linear-independence status over ℚ and over 𝔽_p** outside `S_k`. In particular all minimal circuits coincide as *sets of exponent pairs*. This does not require a finite-field Wronskian theorem.

The actual characteristic-zero circuit spread bound is already kernel-verified for the concrete rational E6 section in `GrowthActualCircuitSpread.actual_spread_bound` (specialized to ℚ), with no supplied `SpreadBound`. Hence, for each good p, the same absolute exponent spread inequality holds for the bounded circuits in `𝔽_p(X)`; the p-modular inequality then follows by choosing `z₁=z₂=0`.

## 2. Restrict the FEI hypothesis rather than prove global SpreadBound

The original `GrowthSpread.SpreadBound` quantifies over exponent pairs of *arbitrary size*. However, the proofs of `circuit_sigma_lt` and `fei_of_spread` invoke it only on circuits contained in an exponent box `[0,2^k)²`, because the input polynomial's support is bounded there. The new `GrowthBoxSpread.BoxSpreadBound` formalizes precisely this finite restriction and `fei_of_box_spread` proves the same FEI conclusion. If the bridge passes Lean CI, only the finite integer-matroid specialization is needed to feed this `BoxSpreadBound` at each fixed k. No false assertion about all positive-characteristic Wronskians is required.

The remaining direct proof chain is:

```
rank matroid(A_k over Q) = rank matroid(A_k over F_p) for p outside S_k
    + verified rational E6 circuit spread
 -> finite-box circuit spread over F_p
 -> FEI_k(D−1) over F_p (new Lean bridge)
 -> 1 ∉ InClosure(F_{p,k}, D), D=floor((k−10)/14).
```

At every fixed k≥24, one also needs a **separate** rational/integral Bezout certificate for `1∈Ideal.span(F_{p,k})`, valid outside a second finite set. Existence over ℂ follows from the kernel-verified E6 unit-ideal theorem, and finite rational linear-system descent yields a rational certificate; clearing denominators gives an integer multiple of 1 and finite bad primes. The selector upper bound (already kernel-verified in all characteristics other than 2) then gives

```
floor((k−10)/14) < LFD(F_{p,k}) ≤ k+2
```

for every fixed k≥24 and all sufficiently large admissible primes p. A growing prime sequence can be selected using the already checked Fermat-number-factor method of PR #232. **This conclusion is not yet claimed as kernel-checked**: the matrix/matroid, Bezout and final composition interfaces remain to be implemented in Lean.

## 3. Positive-characteristic Wronskians over 𝔽_p(X^p)

For `p≥2N+1`, every `g_(u,v)` has degree less than p. The standard basis `1,X,…,X^(p−1)` of `𝔽_p(X)` over `𝔽_p(X^p)` implies that any bounded collection of these polynomials is linearly independent over `𝔽_p(X^p)` **if and only if** it is independent over `𝔽_p`: its coefficient matrix has entries in `𝔽_p` and field extension does not change matrix rank. Thus the usual ordinary-Wronskian obstruction from `1,X^p` cannot arise inside this bounded family. For independent bounded families, size is at most `2N−1<p` and the ordinary Wronskian nonvanishing criterion over derivation constants `𝔽_p(X^p)` is applicable in principle. A general unbounded-characteristic-`p` Wronskian theorem is not asserted.

## 4. Independent exact negative test

At k=2 (N=4), E6 has a *different* admissible-characteristic exception beyond the known unit-ideal obstruction. For the exact seven columns

```
(u,v) = (0,1),(0,3),(1,0),(1,1),(2,3),(3,2),(3,3),
```

the 7×7 determinant is

```
578207808 = 2^6 × 3^12 × 17.
```

These columns have rank seven over ℚ but rank six over 𝔽₁₇. Since `17 ≡ 1 mod 8`, one must exclude p=17 from any naive `all admissible primes` matroid-transfer claim. The same columns recover rank seven over 𝔽₂₅₇. The independent exact regression checks 400 additional deterministic subsets at k=1,2, and the previously established unit-ideal grid determinant counterexample remains in the inherited suite.

**No measured ECDLP performance consequence** follows from this proof strategy.
