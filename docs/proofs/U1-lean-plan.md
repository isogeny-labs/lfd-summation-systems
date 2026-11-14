# U1 in Lean: plan for the selector bound LFD ≤ k + 2

Date: 2026-11-12. Source: `docs/proofs/U1-upper-bounds.md` (Theorem 1, Section 2) with the
hypotheses fixed by `docs/reviews-U1-upper-bound-redteam-20261110.md` (Section 7: nonsingular
curve, 2^k-th roots of unity in the field, k ≥ 2). Status of each item is given in Section 5.

## 1. Target

    theorem lfd_le_depth_add_two {K : Type*} [Field K] (k : ℕ) (depth : 2 ≤ k) {ζ : K}
        (root : IsPrimitiveRoot ζ (2^k)) (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
        HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (k + 2)

and the prime-field entry point (p prime, 2^k ∣ p − 1, which forces p ≥ 5 for k ≥ 2).
The field is any field with a primitive 2^k-th root of unity; the characteristic is then odd.
No hypothesis p ≥ 5 is needed for the algebra: Lemma 2.2 uses only 2 ≠ 0 (Section 3, step 5).

## 2. Route: evaluation on the grid, not normal forms

The written proof works in the group algebra A = K[(Z/2^k)²] and returns through the bridge
`chain_mem_closure_iff`. The formal route stays in the polynomial ring and replaces every
"modulo the chains" step by one lemma:

**Grid lemma (G).** If K contains a primitive 2^k-th root of unity and p vanishes at every chain
point (u_{i,j} = x_i^{2^j}, x_i^{2^k} = 1), then `InClosure (chainInputs K m k) p.totalDegree p`.
Proof: `p − NF p` is derivable at degree deg p (`normalForm_reduction`); `NF p` is squarefree and
still vanishes on the grid, so its expansion has partial degrees below 2^k and vanishes on a
2^k × … grid, hence is zero (`eq_zero_of_eval_zero_at_prod_finset`, `bitExp_injective`).

With (G), a polynomial of degree ≤ D that vanishes on the grid is in the degree-D closure, so
selector identities (idempotent, eigenvector, partition of unity, `E_α·g ≡ E_α·q_α(W)`) are
proved by evaluating at points. No statement about normal forms of products is needed.

## 3. Chain of lemmas

| # | module | statement | uses |
|---|---|---|---|
| 1 | `SelectorGrid` | (G) above, for every m | `ChainClosureBridge.normalForm_*`, `ChainQuotient.expand_*`, `bitExp_*` |
| 2 | `SelectorAssembly` | `lfd_bound_of_projector`: a linear N on the ring with deg N p ≤ σ, σ+1 ≤ D, N(X_i·N p) = N(X_i·p), N = 0 on the ideal, 1 − N 1 and X_i·N p − N(X_i·N p) in the degree-D closure ⇒ `HasLFDBound F D` | the induction of `LocalReductionBound.global_reduction_of_local` (Section 2.5 of the proof) |
| 3 | `SelectorUnivariate` | Lemma 2.3 in an abstract K-algebra: rows `E·(w_j − a_j − b_j w_0)`, `E·(w_0 w_j − …)` from `E·(w_1 − s w_0 + n)`; the remainder recursion a_{j+1} = a_j² − n b_j², b_{j+1} = 2a_j b_j + s b_j²; and `linear_mem`: every c₀ + c₁X vanishing on the roots of Q = X² − sX + n in μ_N gives `E·(c₀ + c₁ w_0)` in the closure | Mathlib `Polynomial` (splitting of X^N − 1) |
| 4 | `SelectorIdempotent` | `E α = N⁻¹ ∏_j (1 + α^{−2^j} u_{0,j})`: degree ≤ k, value δ_{α,x₀} on the grid | geometric product identity |
| 5 | `SelectorFibre` | value of the normalized generator at a chain point = S_3(α, β, t) as a quadratic in β with coefficients c₂ = (α−t)², c₁, c₀; Lemma 2.2: not all three vanish (needs 2 ≠ 0 and 4a³+27b² ≠ 0 only); the fibre lemma: for α ∈ μ_N and j ≤ k, `E α·(w_j − c₀ − c₁w_0)` and `E α·(w_0w_j − c₀ − c₁w_0)` are in the closure at degree k+2 whenever c₀ + c₁β matches on the fibre | 1, 3, 4 |
| 6 | `SelectorBound` | N p = Σ_α E α·lin(Lagrange interpolant of p on the fibre of α); hypotheses of 2; the theorem | 2, 5, Mathlib `Lagrange` |

Differences from the written proof (same mathematics, shorter formal path):

* The case distinction of Lemma 2.3 on deg q is removed: a linear or constant fibre polynomial
  is first multiplied by w_0 (once or twice, allowed at weight ≤ k+1) to a monic quadratic Q
  with the same roots in μ_N. The cases r = 0, 1, 2 become three cases on the 2×2 matrix of
  multiplication by 1 − X^N on K[X]/(Q): det ≠ 0; det = 0 and b ≠ 0 (then (1−a)/b is a root of Q
  in μ_N, because det = b²·Q((1−a)/b)); (a, b) = (1, 0) (then Q ∣ X^N − 1 has two distinct roots
  in μ_N and nothing is to be proved). No gcd is computed.
* The section is Σ_α E_α·(Lagrange interpolant on the fibre), so the projector is defined by
  evaluation and π(a·π(c)) = π(a·c) is immediate.
* `LocalReductionBound.lfd_bound_of_local` is not called (it needs a quotient algebra and a
  section as separate objects); its induction is repeated for a projector on the ring.

## 4. Expected hardest step

Step 5 (fibre lemma): it joins the abstract recursion to the concrete ring (the three leading
coefficient cases, degree bounds of every row, the dependent indices `var 2 k 1 ⟨j, _⟩`).
Step 6 is long but routine. Steps 1–4 are independent of the curve.

## 5. Status (2026-11-12): the full theorem is kernel-checked

Label **kernel** for every statement below. Lean 4.19.0, pinned Mathlib; each module was built
through the lake lock (`lake build ICFallDegrees.<Module>`, about 10 s each); no `sorry`, no new
axiom; `#print axioms` gives `[propext, Classical.choice, Quot.sound]` for every declaration
listed. The modules are not yet imported by `formal/ICFallDegrees.lean` nor listed in
`Audit.lean` (lines to add: Section 7).

| module | lines | main declarations |
|---|---|---|
| `SelectorGrid.lean` | 133 | `grid_mem_closure`, `grid_mem`, `pt`, `eval_expand` |
| `SelectorAssembly.lean` | 99 | `lfd_bound_of_projector` |
| `SelectorUnivariate.lean` | 307 | `rem`, `rem_root`, `rem_dvd`, `two_roots`, `rows`, `linear_mem`, `fibre_rows` |
| `SelectorFibre.lean` | 414 | `E`, `eval_E`, `fibre_nonzero`, `eval_normalized`, `fibre` |
| `SelectorBound.lean` | 469 | `N`, `lfd_le_depth_add_two`, `exists_primitiveRoot_zmod`, `prime_field_lfd_le_depth_add_two` |
| `SelectorExtension.lean` | 98 | `map_inputs`, `lfd_le_of_extension`, `lfd_le_depth_add_two_of_two_ne`, `prime_field_lfd_le_depth_add_two`, `exact_lfd_four_at_depth_two` |

Final statements:

    theorem SelectorBound.lfd_le_depth_add_two {K : Type*} [Field K] (k : ℕ) (depth : 2 ≤ k)
        {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) (a b t : K)
        (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
        HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (k+2)

    theorem SelectorBound.prime_field_lfd_le_depth_add_two (p : ℕ) [Fact p.Prime] (k : ℕ)
        (depth : 2 ≤ k) (root : 2^k ∣ p - 1) (a b t : ZMod p)
        (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
        HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (k+2)

    theorem SelectorExtension.lfd_le_depth_add_two_of_two_ne {K : Type*} [Field K]
        (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
        (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
        HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (k+2)

    theorem SelectorExtension.exact_lfd_four_at_depth_two {K : Type*} [Field K]
        (two_ne : (2 : K) ≠ 0) (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0)
        (nondegenerate : 3*t^2 + a ≠ 0)
        (unit : (1 : MvPolynomial (Fin (2*2)) K) ∈
          Ideal.span (Set.range (PrimeM2General.inputs 2 le_rfl a b t))) :
        IsLastFallDegree (Set.range (PrimeM2General.inputs 2 le_rfl a b t)) 4

The third statement removes the hypothesis 2^k ∣ p − 1: the bound descends from the 2^k-th
cyclotomic extension by `ScalarExtensionTransport.bound_iff` (the reduction that the review,
finding F1, left unchecked). The fourth is the sanity check against the exact data: at k = 2
the bound is attained (with `PrimeM2Lower.lfd_ge_four_of_unit`), in agreement with the
recorded value LFD = 4 off the locus 3t² + a = 0.

Not formalized (scope of the kernel statement): unequal depths k₁ ≠ k₂, coset bases
x_i = s_i·u_{i,0}, the raw generator (bound k + 4), depth k = 1 (the presentation
`PrimeM2General.inputs` requires k ≥ 2), and Theorem 1 for a general generator of x₂-degree
≤ 2. The singular case stays open, as in the review. For these the status remains
"written proof".

## 6. Corrections and remarks on the written proof

No gap was found. Three points where the formal statement differs from the text:

1. **p ≥ 5 is not used.** Lemma 2.2 needs only 2 ≠ 0: from c₂ = (α−t)² = 0, c₁ = −4f(t) = 0 and
   c₀ = f′(t)² − 8t·f(t) = 0 one gets f(t) = f′(t) = 0 and then
   4a³ + 27b² = (27b − 27at − 27t³)·f(t) + (4a² + 15at² + 9t⁴)·f′(t) = 0, an identity over Z.
   The kernel hypothesis is `4*a^3 + 27*b^2 ≠ 0` in any characteristic other than two (in
   characteristic 3 this is a ≠ 0, the case the review mentions in 4(g)).
2. **2^k ∣ p − 1 is not needed** for the bound on the polynomial system (third statement
   above). The table of Section 0 of `U1-upper-bounds.md` and the paper may say "every field
   of characteristic ≠ 2"; without the roots of unity the system no longer encodes a
   subgroup factor base, so the hypothesis remains natural for the application.
3. **Lemma 2.3 does not need the gcd d = gcd(q, X^N − 1) or the case distinction on deg q.**
   The formal proof multiplies a fibre polynomial of degree < 2 by X (once or twice) to a
   monic quadratic Q with the same roots in μ_N and uses det M = b²·Q((1−a)/b) for the
   matrix M of multiplication by 1 − X^N on K[X]/(Q). The statement proved is the one used
   in 2.5: every c₀ + c₁X vanishing on the roots of q_α in μ_N gives the row e_α·(c₀ + c₁x₂).

The budget accounting of the text is confirmed by the kernel: the selector multiplies only
the generator (k + 2), every later product is one variable times a row of degree ≤ k + 1.

## 7. Integration (not done here: existing files were not edited)

Add to `formal/ICFallDegrees.lean`:

    import ICFallDegrees.SelectorGrid
    import ICFallDegrees.SelectorAssembly
    import ICFallDegrees.SelectorUnivariate
    import ICFallDegrees.SelectorFibre
    import ICFallDegrees.SelectorBound
    import ICFallDegrees.SelectorExtension

Add to `formal/ICFallDegrees/Audit.lean` (it must import `ICFallDegrees.SelectorExtension`,
directly or through the root module):

    #print axioms ICFallDegrees.SelectorGrid.grid_mem_closure
    #print axioms ICFallDegrees.SelectorAssembly.lfd_bound_of_projector
    #print axioms ICFallDegrees.SelectorUnivariate.fibre_rows
    #print axioms ICFallDegrees.SelectorFibre.fibre_nonzero
    #print axioms ICFallDegrees.SelectorFibre.fibre
    #print axioms ICFallDegrees.SelectorBound.lfd_le_depth_add_two
    #print axioms ICFallDegrees.SelectorBound.prime_field_lfd_le_depth_add_two
    #print axioms ICFallDegrees.SelectorExtension.lfd_le_depth_add_two_of_two_ne
    #print axioms ICFallDegrees.SelectorExtension.prime_field_lfd_le_depth_add_two
    #print axioms ICFallDegrees.SelectorExtension.exact_lfd_four_at_depth_two

`SelectorExtension` imports `Mathlib.NumberTheory.Cyclotomic.PrimitiveRoots` (already in the
Mathlib cache; no new Mathlib build was needed).

Next steps, in order of cost: (i) the same proof for unequal depths needs a presentation with
two depths, which the repository does not define; (ii) coset bases need only a change of the
evaluation point (x_i = s_i·u_{i,0}) in `SelectorFibre.eval_normalized` and the fibre
coefficients; (iii) Theorem 1 for a general generator of x₂-degree ≤ 2 is `SelectorFibre.fibre`
with `c₂, c₁, c₀` replaced by the coefficient functions and `fibre_nonzero` as a hypothesis.
