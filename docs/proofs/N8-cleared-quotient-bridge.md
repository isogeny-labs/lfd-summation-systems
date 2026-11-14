# N8 cleared quotient-rule bridge (fixed complex two-torsion family)

## Kernel target

The Lean module `ICFallDegrees.GrowthPolynomialQuotientBridge` is an
algebraic prerequisite for the genuine characteristic-zero rational
Wronskian.  It never imports or assumes `GrowthSpread.SpreadBound` or
`GrowthMain.CurveModel`.

For any field `K`, any `c : K`, and `u,v : ℕ`, set

```
A = X, B = X + c, T = X - 1
N = A^u B^v, H = T^v
D = ABT
Q = u BT + v AT - v AB
```

The polynomial identity is

```
D (N' H - N H') = Q N H.
```

No denominator is inverted here, and the proof handles u=0 or v=0
without a case omission.  The key generic equality, for P' = 1, is

```
P (P^n)' = n P^n,
```

which follows from the formal derivative power rule (including n=0).
Expanding the derivative of N and substituting those three Euler identities
establishes the result by polynomial commutative-ring normalization.

For `c=8`, we have `Q = u(X+8)(X-1)-9vX` by the already audited
`GrowthRationalLogNumerator.logNumerator_eq`.

## Exact follow-on RatFunc identity still to prove

As `H ≠ 0` in `K[X]`, the natural map to `RatFunc K` permits
`f = N/H`.  Define or identify the **genuine derivation** `δ` on
`RatFunc K`, with `δ(X)=1` and `δ(a)=0` for constants.  The quotient
rule gives

```
δ(f) = (N'H - NH') / H² = f Q/D.
```

The next step is an actual Lean proof of the quotient rule (not an
assumed map with the required formula), followed by induction showing

```
δ^s(f) = f P_s/D^s,
P_0=1, P_(s+1)=D P_s' + (Q-sD')P_s.
```

Then one must **separately** establish characteristic-zero
constant-kernel Wronskian nonvanishing for arbitrary independent
columns, the cofactor ratio, and the two finite/infinity pole
bounds.  The existing bound `deg P_s≤2s` gives
`deg cofactor≤r(r-1)`; it is not by itself a circuit-spread theorem.

## Status protocol

Exact symbolic algebra independently checked 64 exponent pairs with
the parameter c left free.  Treat the Lean declarations as candidates
until the pinned Lean 4.19.0 / Mathlib v4.19.0 workflow compiles and
`#print axioms` shows only standard axioms.  No new result here asserts
unconditional N8 last-fall-degree growth, general curves, or any ECDLP
runtime speedup.
