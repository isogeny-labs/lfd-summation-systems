import ICFallDegrees.GrowthSpread
import Mathlib.Tactic

/-!
# Restricted-box spread is sufficient for the actual FEI splitting theorem

The global SpreadBound quantified over all exponent pairs is stronger than
necessary for the degree-k growth theorem: its proof only invokes circuit
spread on exponents 0 <= e_i < 2^k from a fixed bounded support.
This standalone kernel bridge isolates that exact finite hypothesis.
It will permit transporting finite rational-function coefficient
matroids to all but finitely many prime fields at each fixed k.
No rational/prime matroid transfer or prime LFD growth is assumed here.
-/

noncomputable section
namespace ICFallDegrees.GrowthBoxSpread

open MvPolynomial BinaryWeight GrowthCounting GrowthReduction GrowthSpread

/-- Only circuits contained in the finite exponent box are relevant to
FEI at depth k. Unlike global SpreadBound, this excludes all other circuits. -/
def BoxSpreadBound (F₀ : Type*) {F : Type*}
    [Field F₀] [Field F] [Algebra F₀ F]
    (ξ : Fin 2 → F) (p k : ℕ) : Prop :=
  ∀ J : Finset (Fin 2 →₀ ℕ),
    (∀ d ∈ J, d 0 < 2 ^ k ∧ d 1 < 2 ^ k) →
    Circuit F₀ (mono ξ) J →
    ∀ e ∈ J, ∀ e' ∈ J, ∃ z₁ z₂ : ℤ,
      |2 * ((e 0 : ℤ) - e' 0) - p * z₁| +
        |2 * ((e 1 : ℤ) - e' 1) - p * z₂| ≤
          6 * ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2)

section Spread
variable {K F₀ F : Type*} [Field K] [Field F₀] [Field F]
  [Algebra K F₀] [Algebra F₀ F] [Algebra K F]
  [IsScalarTower K F₀ F]

/-! ## Lemma 4.3 -/

omit [Field K] [Algebra K F₀] [Algebra K F] [IsScalarTower K F₀ F] in
/-- **Lemma 4.3 (a circuit of low weight is small).** From `SpreadBound`: a circuit of
exponent pairs in `[0, 2^k)²` of weight `≤ w` has `σ = 3(m-1)(m-2) < 2^(lam-1)`, for every
`lam ≥ 1` bounding the solutions of `2^(ℓ-1) < 7776·C(2ℓ, ≤w)^2` (for example `lam = Λ(w)`). -/
theorem circuit_sigma_lt_box (ξ : Fin 2 → F) (p k w lam : ℕ) (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p)
    (spread : BoxSpreadBound F₀ ξ p k) (hlam : LambdaBound w lam) (h1 : 1 ≤ lam)
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
  obtain ⟨z₁, z₂, hz⟩ := spread J (fun q hq =>
    ⟨(hJ q hq).1, (hJ q hq).2.1⟩) circ d hd e he
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
theorem fei_of_box_spread (ξ : Fin 2 → F) (p k w lam : ℕ) (hp : p = 0 ∨ 2 * 2 ^ k + 1 ≤ p)
    (spread : BoxSpreadBound F₀ ξ p k) (hlam : LambdaBound w lam) (h1 : 1 ≤ lam)
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
    have hσ := circuit_sigma_lt_box (F₀ := F₀) ξ p k w lam hp spread hlam h1 J circ
      (fun d hd => boxed d (hsupp d hd))
    obtain ⟨z₁, z₂, hz⟩ := spread J (fun q hq =>
      ⟨(boxed q (hsupp q hq)).1, (boxed q (hsupp q hq)).2.1⟩) circ e he e' he'
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

#print axioms circuit_sigma_lt_box
#print axioms fei_of_box_spread

end Spread
end ICFallDegrees.GrowthBoxSpread
