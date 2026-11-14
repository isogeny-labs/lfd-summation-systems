import ICFallDegrees.GrowthReduction
import ICFallDegrees.GrowthCounting
import Mathlib.LinearAlgebra.LinearIndependent.Defs
import Mathlib.Algebra.MvPolynomial.Eval

noncomputable section

/-!
# Growth theorem: the splitting property from the spread of circuits
(Lemma 3.4, Lemma 4.3 and Theorem 5.1 of `docs/proofs/N8-growth.md`)

Setting: fields `K → F₀ → F` (a scalar tower) and two elements `ξ₁, ξ₂ ∈ F`.  In the written
proof `F` is the function field of the curve, `ξ₁ = x(P)`, `ξ₂ = x(R - P)`, and `F₀` is the
constant field in characteristic `0` and `F^p` in characteristic `p`.  **None of this is
constructed here.**

* `Circuit F₀ v J`: the family `v` is linearly dependent over `F₀` on the finite set `J` and
  independent on every proper subset.
* `SpreadBound F₀ ξ p`: **the statement of Theorem 3.3 of the written proof, as a hypothesis.**
  It is not proved in Lean (it needs Wronskians on function fields of curves).
* `split_of_no_crossing` (Lemma 3.4, **kernel**): a vanishing `F₀`-combination splits along any
  partition that no circuit crosses.
* `circuit_sigma_lt` (Lemma 4.3, **kernel**, from `SpreadBound`).
* `fei_of_spread` (Theorem 5.1, **kernel**, from `SpreadBound`): `FEI_k(w)` for
  `k ≥ w + lam + 1`, for the ideal `ker (x_i ↦ ξ_i)`.
-/

namespace ICFallDegrees.GrowthSpread

open MvPolynomial BinaryWeight GrowthCounting GrowthReduction

/-! ## Circuits and Lemma 3.4 -/

section Circuits

variable {F₀ F ι : Type*} [Field F₀] [AddCommGroup F] [Module F₀ F]

/-- A circuit of the family `v` over `F₀`: a finite index set on which `v` is linearly
dependent, every proper subset of which is linearly independent. -/
def Circuit (F₀ : Type*) {F ι : Type*} [Field F₀] [AddCommGroup F] [Module F₀ F] (v : ι → F)
    (J : Finset ι) : Prop :=
  ¬ LinearIndepOn F₀ v (J : Set ι) ∧ ∀ J' ⊂ J, LinearIndepOn F₀ v (J' : Set ι)

/-- A finite dependent index set contains a circuit. -/
theorem exists_circuit_subset (v : ι → F) (S : Finset ι)
    (dep : ¬ LinearIndepOn F₀ v (S : Set ι)) : ∃ J ⊆ S, Circuit F₀ v J := by
  classical
  obtain ⟨J, ⟨hJS, hJdep⟩, hmin⟩ := WellFounded.has_min
    (Finset.isWellFounded_ssubset (α := ι)).wf
    {J : Finset ι | J ⊆ S ∧ ¬ LinearIndepOn F₀ v (J : Set ι)} ⟨S, subset_rfl, dep⟩
  refine ⟨J, hJS, hJdep, ?_⟩
  intro J' hJ'
  by_contra h
  exact hmin J' ⟨hJ'.subset.trans hJS, h⟩ hJ'

/-- **Lemma 3.4 (splitting along circuits).** If no circuit inside the support of a vanishing
`F₀`-linear combination meets both `{T}` and its complement, then the part of the combination
on `{T}` vanishes. -/
theorem split_of_no_crossing (v : ι → F) (T : ι → Prop) [DecidablePred T] :
    ∀ (n : ℕ) (c : ι →₀ F₀), c.support.card ≤ n →
      (∀ J : Finset ι, J ⊆ c.support → Circuit F₀ v J →
        (∀ e ∈ J, T e) ∨ (∀ e ∈ J, ¬ T e)) →
      Finsupp.linearCombination F₀ v c = 0 →
      Finsupp.linearCombination F₀ v (c.filter T) = 0 := by
  classical
  intro n
  induction n with
  | zero =>
    intro c hc _ _
    have : c = 0 := by
      rw [← Finsupp.support_eq_empty, ← Finset.card_eq_zero]
      omega
    rw [this, Finsupp.filter_zero, map_zero]
  | succ n ih =>
    intro c hc hcross h0
    by_cases hzero : c = 0
    · rw [hzero, Finsupp.filter_zero, map_zero]
    have dep : ¬ LinearIndepOn F₀ v (c.support : Set ι) :=
      linearDepOn_iff'.mpr ⟨c, Finsupp.mem_supported_support F₀ c, h0, hzero⟩
    obtain ⟨J, hJS, hJ⟩ := exists_circuit_subset v c.support dep
    obtain ⟨γ, hγsup, hγ0, hγne⟩ := linearDepOn_iff'.mp hJ.1
    rw [Finsupp.mem_supported] at hγsup
    obtain ⟨e₀, he₀⟩ := Finsupp.support_nonempty_iff.mpr hγne
    have hγe₀ : γ e₀ ≠ 0 := Finsupp.mem_support_iff.mp he₀
    have he₀J : e₀ ∈ J := by exact_mod_cast hγsup he₀
    have he₀c : e₀ ∈ c.support := hJS he₀J
    set r := c e₀ / γ e₀ with hr
    set c' := c - r • γ with hc'
    have hc'e₀ : c' e₀ = 0 := by
      rw [hc', Finsupp.sub_apply, Finsupp.smul_apply, smul_eq_mul, hr, div_mul_cancel₀ _ hγe₀,
        sub_self]
    have hsub : c'.support ⊆ c.support.erase e₀ := by
      intro x hx
      rw [Finsupp.mem_support_iff] at hx
      rw [Finset.mem_erase]
      refine ⟨fun h => hx (h ▸ hc'e₀), ?_⟩
      by_contra hxc
      rw [Finsupp.not_mem_support_iff] at hxc
      have hxγ : γ x = 0 := by
        by_contra hne
        have hxJ : x ∈ J := by exact_mod_cast hγsup (Finsupp.mem_support_iff.mpr hne)
        exact (Finsupp.mem_support_iff.mp (hJS hxJ)) hxc
      apply hx
      rw [hc', Finsupp.sub_apply, Finsupp.smul_apply, hxc, hxγ, smul_zero, sub_zero]
    have hcard : c'.support.card ≤ n := by
      have h1 := Finset.card_le_card hsub
      rw [Finset.card_erase_of_mem he₀c] at h1
      omega
    have hsub' : c'.support ⊆ c.support := hsub.trans (Finset.erase_subset _ _)
    have h0' : Finsupp.linearCombination F₀ v c' = 0 := by
      rw [hc', map_sub, map_smul, h0, hγ0, smul_zero, sub_zero]
    have step := ih c' hcard (fun J' hJ' hcirc => hcross J' (hJ'.trans hsub') hcirc) h0'
    have hγT : Finsupp.linearCombination F₀ v (γ.filter T) = 0 := by
      rcases hcross J hJS hJ with hall | hnone
      · have : γ.filter T = γ := by
          rw [Finsupp.filter_eq_self_iff]
          intro x hx
          exact hall x (by exact_mod_cast hγsup (Finsupp.mem_support_iff.mpr hx))
        rw [this, hγ0]
      · have : γ.filter T = 0 := by
          rw [Finsupp.filter_eq_zero_iff]
          intro x hx
          by_contra hne
          exact hnone x (by exact_mod_cast hγsup (Finsupp.mem_support_iff.mpr hne)) hx
        rw [this, map_zero]
    have hdecomp : c.filter T = c'.filter T + r • γ.filter T := by
      ext x
      simp only [Finsupp.filter_apply, Finsupp.add_apply, Finsupp.smul_apply, hc',
        Finsupp.sub_apply]
      split_ifs
      · simp
      · simp
    rw [hdecomp, map_add, map_smul, step, hγT, smul_zero, add_zero]

end Circuits

/-! ## The analytic hypothesis -/

section Spread

variable {K F₀ F : Type*} [Field K] [Field F₀] [Field F] [Algebra K F₀] [Algebra F₀ F]
  [Algebra K F] [IsScalarTower K F₀ F]

/-- The monomial `ξ^e = ξ₁^(e 0) · ξ₂^(e 1)`. -/
def mono (ξ : Fin 2 → F) (e : Fin 2 →₀ ℕ) : F := ξ 0 ^ e 0 * ξ 1 ^ e 1

/-- **Theorem 3.3 of `docs/proofs/N8-growth.md` ("circuits are local"), as a hypothesis.**

For every circuit `J` of the monomials `ξ^e` over `F₀` (`m = |J|`) and all `e, e' ∈ J`:
`‖2(e₁ - e'₁)‖_p + ‖2(e₂ - e'₂)‖_p ≤ 6(m-1)(m-2)`, where `‖n‖_p` is the distance from `n` to
`pℤ` (written here with explicit multiples `p·z`; for `p = 0` it is `|n|`).

This is **not proved in Lean**.  The written proof derives it, for the function field `F` of
the curve `y² = x³ + a x + b` (nonsingular), `ξ₁ = x(P)`, `ξ₂ = x(R - P)` with `x(R) = t ≠ 0`,
`F₀` the constants (characteristic `0`) or `F^p` (characteristic `p`), from a Wronskian bound
for vanishing sums on a curve (Theorem 3.1; Brownawell–Masser in characteristic `0`). -/
def SpreadBound (F₀ : Type*) {F : Type*} [Field F₀] [Field F] [Algebra F₀ F] (ξ : Fin 2 → F)
    (p : ℕ) : Prop :=
  ∀ J : Finset (Fin 2 →₀ ℕ), Circuit F₀ (mono ξ) J → ∀ e ∈ J, ∀ e' ∈ J, ∃ z₁ z₂ : ℤ,
    |2 * ((e 0 : ℤ) - e' 0) - p * z₁| + |2 * ((e 1 : ℤ) - e' 1) - p * z₂| ≤
      6 * ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2)

/-- `S_3(ξ₁, ξ₂, t)` evaluated in `F`. -/
def S3at (a b t : K) (ξ : Fin 2 → F) : F :=
  (ξ 0 - ξ 1) ^ 2 * algebraMap K F (t ^ 2) -
      2 * ((ξ 0 + ξ 1) * (ξ 0 * ξ 1 + algebraMap K F a) + algebraMap K F (2 * b)) *
        algebraMap K F t +
    (ξ 0 * ξ 1 - algebraMap K F a) ^ 2 - algebraMap K F (4 * b) * (ξ 0 + ξ 1)

omit [Field F₀] [Algebra K F₀] [Algebra F₀ F] [IsScalarTower K F₀ F] in
theorem aeval_S3poly (a b t : K) (ξ : Fin 2 → F) :
    aeval ξ (S3poly a b t) = S3at a b t ξ := by
  simp only [S3poly, S3at, map_sub, map_add, map_mul, map_pow, aeval_X, aeval_C, map_ofNat]

/-- The coefficient vector of a polynomial, with coefficients in `F₀`. -/
def cf (F₀ : Type*) [Field F₀] [Algebra K F₀] (f : Q K) : (Fin 2 →₀ ℕ) →₀ F₀ :=
  Finsupp.onFinset f.support (fun d => algebraMap K F₀ (coeff d f)) (by
    intro d hd
    rw [mem_support_iff]
    intro h
    exact hd (show algebraMap K F₀ (coeff d f) = 0 by rw [h, map_zero]))

omit [Field F] [Algebra F₀ F] [Algebra K F] [IsScalarTower K F₀ F] in
theorem cf_apply (f : Q K) (d : Fin 2 →₀ ℕ) : cf F₀ f d = algebraMap K F₀ (coeff d f) := rfl

omit [Field F] [Algebra F₀ F] [Algebra K F] [IsScalarTower K F₀ F] in
theorem support_cf_subset (f : Q K) : (cf F₀ f).support ⊆ f.support :=
  Finsupp.support_onFinset_subset

omit [Field F] [Algebra F₀ F] [Algebra K F] [IsScalarTower K F₀ F] in
theorem cf_topPart (k : ℕ) (i : Fin 2) (j : ℕ) (f : Q K) :
    cf F₀ (topPart k i j f) = (cf F₀ f).filter (fun d => 2 ^ k - 2 ^ j ≤ d i) := by
  ext d
  rw [cf_apply, coeff_topPart, Finsupp.filter_apply, cf_apply]
  split_ifs
  · rfl
  · rw [map_zero]

/-- Evaluation at `ξ` is the `F₀`-linear combination of the monomials `ξ^e`. -/
theorem linearCombination_cf (ξ : Fin 2 → F) (f : Q K) :
    Finsupp.linearCombination F₀ (mono ξ) (cf F₀ f) = aeval ξ f := by
  classical
  rw [Finsupp.linearCombination_apply, cf, Finsupp.onFinset_sum]
  swap
  · intro d
    exact zero_smul _ _
  conv_rhs => rw [f.as_sum]
  rw [map_sum]
  apply Finset.sum_congr rfl
  intro d _
  rw [aeval_monomial, Finsupp.prod_fintype _ _ (fun i => pow_zero _), Fin.prod_univ_two,
    IsScalarTower.algebraMap_apply K F₀ F, Algebra.smul_def]
  rfl

/-! ## Lemma 4.3 -/

omit [Field K] [Algebra K F₀] [Algebra K F] [IsScalarTower K F₀ F] in
/-- **Lemma 4.3 (a circuit of low weight is small).** From `SpreadBound`: a circuit of
exponent pairs in `[0, 2^k)²` of weight `≤ w` has `σ = 3(m-1)(m-2) < 2^(lam-1)`, for every
`lam ≥ 1` bounding the solutions of `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2` (for example `lam = Λ(w)`). -/
theorem circuit_sigma_lt (ξ : Fin 2 → F) (p k w lam : ℕ) (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p)
    (spread : SpreadBound F₀ ξ p) (hlam : LambdaBound w lam) (h1 : 1 ≤ lam)
    (J : Finset (Fin 2 →₀ ℕ)) (circ : Circuit F₀ (mono ξ) J)
    (hJ : ∀ d ∈ J, d 0 < 2 ^ k ∧ d 1 < 2 ^ k ∧ wt (d 0) + wt (d 1) ≤ w) :
    3 * (J.card - 1) * (J.card - 2) < 2 ^ (lam - 1) := by
  classical
  apply sigma_lt_of_lambda w lam J.card hlam h1
  intro ℓ hℓ
  rcases J.eq_empty_or_nonempty with rfl | ⟨e, he⟩
  · simp
  have hm : 1 ≤ J.card := Finset.card_pos.mpr ⟨e, he⟩
  have hcast : (6 * ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2)) =
      2 * ((3 * (J.card - 1) * (J.card - 2) : ℕ) : ℤ) := by
    obtain ⟨n, hn⟩ : ∃ n, J.card = n + 1 := ⟨J.card - 1, by omega⟩
    rw [hn]
    cases n with
    | zero => norm_num
    | succ n =>
      simp only [Nat.add_sub_cancel, show n + 1 + 1 - 2 = n by omega]
      push_cast
      ring
  let φ : (Fin 2 →₀ ℕ) → ℕ × ℕ := fun d => (d 0, d 1)
  have φinj : Function.Injective φ := by
    intro d d' h
    simp only [φ, Prod.mk.injEq] at h
    ext i
    fin_cases i
    · exact h.1
    · exact h.2
  rw [← Finset.card_image_of_injective J φinj]
  apply card_le_of_window (2 ^ k) p (3 * (J.card - 1) * (J.card - 2)) ℓ w hp hℓ (J.image φ)
    (e 0) (e 1) (hJ e he).1 (hJ e he).2.1
  intro q hq
  obtain ⟨d, hd, rfl⟩ := Finset.mem_image.mp hq
  obtain ⟨h0, h1', hw⟩ := hJ d hd
  obtain ⟨z₁, z₂, hz⟩ := spread J circ d hd e he
  rw [hcast] at hz
  refine ⟨h0, h1', hw, ⟨z₁, ?_⟩, ⟨z₂, ?_⟩⟩
  · have := abs_nonneg (2 * ((d 1 : ℤ) - e 1) - p * z₂)
    simp only [φ]
    push_cast at hz ⊢
    linarith
  · have := abs_nonneg (2 * ((d 0 : ℤ) - e 0) - p * z₁)
    simp only [φ]
    push_cast at hz ⊢
    linarith

/-! ## Theorem 5.1 -/

/-- **Theorem 5.1, conditional on Theorem 3.3 (`SpreadBound`).** For `p = 0` or
`p ≥ 2·2^k + 1`, `lam ≥ 1` with `LambdaBound w lam`, and `k ≥ w + lam + 1`, the ideal
`ker (K[x₁, x₂] → F, x_i ↦ ξ_i)` has the splitting property `FEI_k(w)`. -/
theorem fei_of_spread (ξ : Fin 2 → F) (p k w lam : ℕ) (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p)
    (spread : SpreadBound F₀ ξ p) (hlam : LambdaBound w lam) (h1 : 1 ≤ lam)
    (hk : w + lam + 1 ≤ k) :
    FEI (RingHom.ker (aeval ξ : Q K →ₐ[K] F)) k w := by
  classical
  intro f hf boxed i j hj
  rw [RingHom.mem_ker] at hf ⊢
  have coord : ∀ d ∈ f.support, ∀ i' : Fin 2, d i' < 2 ^ k ∧ wt (d i') ≤ w := by
    intro d hd i'
    obtain ⟨b0, b1, bw⟩ := boxed d hd
    fin_cases i'
    · exact ⟨b0, by simp only [Fin.zero_eta]; omega⟩
    · exact ⟨b1, by simp only [Fin.mk_one]; omega⟩
  have key := split_of_no_crossing (F₀ := F₀) (mono ξ) (fun d => 2 ^ k - 2 ^ j ≤ d i)
    (cf F₀ f).support.card (cf F₀ f) le_rfl ?_ (by rw [linearCombination_cf]; exact hf)
  · rw [← cf_topPart, linearCombination_cf] at key
    exact key
  · intro J hJS circ
    by_contra hcon
    rw [not_or] at hcon
    obtain ⟨hnall, hnnone⟩ := hcon
    push_neg at hnall hnnone
    obtain ⟨e', he', hnot⟩ := hnall
    obtain ⟨e, he, htop⟩ := hnnone
    have hsupp : ∀ d ∈ J, d ∈ f.support := fun d hd => support_cf_subset f (hJS hd)
    have hσ := circuit_sigma_lt (F₀ := F₀) ξ p k w lam hp spread hlam h1 J circ
      (fun d hd => boxed d (hsupp d hd))
    obtain ⟨z₁, z₂, hz⟩ := spread J circ e he e' he'
    have hm : 1 ≤ J.card := Finset.card_pos.mpr ⟨e, he⟩
    have hcast : (6 * ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2)) =
        2 * ((3 * (J.card - 1) * (J.card - 2) : ℕ) : ℤ) := by
      obtain ⟨n, hn⟩ : ∃ n, J.card = n + 1 := ⟨J.card - 1, by omega⟩
      rw [hn]
      cases n with
      | zero => norm_num
      | succ n =>
        simp only [Nat.add_sub_cancel, show n + 1 + 1 - 2 = n by omega]
        push_cast
        ring
    rw [hcast] at hz
    have hzi : ∃ z : ℤ, |2 * ((e i : ℤ) - e' i) - p * z| ≤
        2 * ((3 * (J.card - 1) * (J.card - 2) : ℕ) : ℤ) := by
      have n1 := abs_nonneg (2 * ((e 0 : ℤ) - e' 0) - p * z₁)
      have n2 := abs_nonneg (2 * ((e 1 : ℤ) - e' 1) - p * z₂)
      fin_cases i
      · exact ⟨z₁, by simp only [Fin.zero_eta]; linarith⟩
      · exact ⟨z₂, by simp only [Fin.mk_one]; linarith⟩
    exact straddle_contra k lam j _ p w (e i) (e' i) hp h1 hk hj hσ
      (coord e (hsupp e he) i).1 htop (not_le.mp hnot)
      (coord e (hsupp e he) i).2 (coord e' (hsupp e' he') i).2 hzi

end Spread

end ICFallDegrees.GrowthSpread
