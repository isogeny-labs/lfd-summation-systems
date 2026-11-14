import ICFallDegrees.GrowthReduction

noncomputable section

/-!
# Growth theorem: Proposition 2.1 as an equivalence, and the low-part criterion
(Propositions 2.1 and 7.1 of `docs/proofs/N8-growth.md`)

* `closed_iff_fei` (Proposition 2.1): for a prime `𝔭` of `K[x₁, x₂]` containing neither `x_i` nor
  `x_i^(2^k) - 1`, the no-wrap ideal space `S_D(k)` is closed under the degree-`D` rule if and
  only if `FEI_k(D-1)` holds.
* `not_one_of_crit` (Proposition 7.1): the low-part criterion `Crit_D(k)` implies that `1` is
  not in the degree-`D` closure of the audited system; `closure_normalForm_mem_sharpSpace` is
  `V_D(k) ⊂ S^♯`; `crit_of_fei`: `FEI_k(D-1)` implies `Crit_D(k)`.

These do not use the analytic hypothesis of `GrowthSpread`.
-/

namespace ICFallDegrees.GrowthCriterion

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight ChainTranslation ChainModelMap
open GrowthReduction
open PrimeM2General (inputs normalized)

variable {K : Type*} [Field K]

/-- `S_D(k)` is closed under the degree-`D` rule of the chain quotient. -/
def Closed (𝔭 : Ideal (Q K)) (k D : ℕ) : Prop :=
  ∀ (v : Fin (2 * k)) (w : P (2 * k) K), w ∈ idealSpace 𝔭 k D → w.totalDegree < D →
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X v * w) ∈
      idealSpace 𝔭 k D

/-- The reduced chain polynomial with a given lift. -/
def lower (k : ℕ) (f : Q K) : P (2 * k) K :=
  ∑ d ∈ f.support, monomial (digits k (fun i => d i)) (coeff d f)

theorem lower_reduced (k : ℕ) (f : Q K) : lower k f ∈ reducedSpace (chainInputs K 2 k) :=
  (reducedSpace _).sum_mem (fun _ _ => digits_reduced _ _)

theorem expand_lower {k n : ℕ} (f : Q K) (boxed : Boxed k n f) : expand K 2 k (lower k f) = f := by
  classical
  unfold lower
  rw [map_sum]
  conv_rhs => rw [f.as_sum]
  apply Finset.sum_congr rfl
  intro d hd
  obtain ⟨b0, b1, _⟩ := boxed d hd
  have lt : ∀ i : Fin 2, (fun i => d i) i < 2 ^ k := by
    intro i
    fin_cases i
    · exact b0
    · exact b1
  rw [expand_digits _ lt]
  have e : fs (fun i => d i) = d := Finsupp.ext (fun _ => rfl)
  rw [e]

theorem lower_totalDegree_le {k n : ℕ} (f : Q K) (boxed : Boxed k n f) :
    (lower k f).totalDegree ≤ n := by
  classical
  unfold lower
  apply (totalDegree_finset_sum _ _).trans
  apply Finset.sup_le
  intro d hd
  obtain ⟨b0, b1, bw⟩ := boxed d hd
  have lt : ∀ i : Fin 2, (fun i => d i) i < 2 ^ k := by
    intro i
    fin_cases i
    · exact b0
    · exact b1
  have h1 := totalDegree_monomial_le (digits k (fun i => d i)) (coeff d f)
  have h2 : (digits k (fun i => d i)).degree = wt (d 0) + wt (d 1) := by
    have h := digits_degree (k := k) (fun i => d i) lt
    rw [Fin.sum_univ_two] at h
    exact h
  have h3 : (digits k (fun i => d i)).degree = (digits k (fun i => d i)).sum (fun _ => id) := rfl
  omega

/-- **Proposition 2.1.** `S_D(k)` is closed under the degree-`D` rule if and only if
`FEI_k(D-1)` holds (`𝔭` prime, `x_i ∉ 𝔭`, `x_i^(2^k) - 1 ∉ 𝔭`, `D ≥ 1`). -/
theorem closed_iff_fei (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) [𝔭.IsPrime]
    (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (k : ℕ)
    (hXN : ∀ i : Fin 2, (X i : Q K) ^ 2 ^ k - 1 ∉ 𝔭) (D : ℕ) (hD : 1 ≤ D) :
    Closed 𝔭 k D ↔ FEI 𝔭 k (D - 1) := by
  constructor
  · intro closed f hf boxed i j hj
    have hw : lower k f ∈ idealSpace 𝔭 k D :=
      ⟨lower_reduced k f, (lower_totalDegree_le f boxed).trans (by omega), by
        rw [expand_lower f boxed]; exact hf⟩
    have hT := (closed (var 2 k i ⟨j, hj⟩) (lower k f) hw
      (lt_of_le_of_lt (lower_totalDegree_le f boxed) (by omega))).2.2
    have id := lift_identity two_ne (lower k f) (lower_reduced k f) i ⟨j, hj⟩
    rw [expand_lower f boxed] at id
    simp only [Fin.val_mk] at id
    have key : (X i ^ 2 ^ k - 1) * topPart k i j f ∈ 𝔭 := by
      have e : (X i ^ 2 ^ k - 1) * topPart k i j f =
          X i ^ 2 ^ k * f - (X i ^ 2 ^ k * (f - topPart k i j f) + topPart k i j f) := by ring
      rw [e, ← id]
      exact 𝔭.sub_mem (𝔭.mul_mem_left _ hf) (𝔭.mul_mem_left _ hT)
    rcases Ideal.IsPrime.mem_or_mem inferInstance key with h | h
    · exact absurd h (hXN i)
    · exact h
  · intro fei v w hw hdeg
    exact idealSpace_closed two_ne 𝔭 hX k D fei v w hw hdeg

/-! ## Proposition 7.1 -/

/-- `S^♯ = S' + Σ_v u_v · S'` with `S' = S_{D-1}(k)`, products in the chain quotient. -/
def sharpSpace (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) (k D : ℕ) : Submodule K (P (2 * k) K) :=
  idealSpace 𝔭 k (D - 1) ⊔ ⨆ v : Fin (2 * k), (idealSpace 𝔭 k (D - 1)).map
    ((ChainClosureBridge.nfMap (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (reduced_chain_ideal_trivial 2 k two_ne)).comp (LinearMap.mulLeft K (X v)))

/-- **The low-part criterion `Crit_D(k)`**: `S^♯ ∩ M_{≤D-1} = S_{D-1}(k)`. -/
def Crit (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) (k D : ℕ) : Prop :=
  ∀ w ∈ sharpSpace two_ne 𝔭 k D, w.totalDegree ≤ D - 1 → w ∈ idealSpace 𝔭 k (D - 1)

theorem sharp_products (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) (k D : ℕ)
    (crit : Crit two_ne 𝔭 k D) (v : Fin (2 * k)) (w : P (2 * k) K)
    (hw : w ∈ sharpSpace two_ne 𝔭 k D) (hdeg : w.totalDegree < D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X v * w) ∈
      sharpSpace two_ne 𝔭 k D := by
  have hlow := crit w hw (by omega)
  exact Submodule.mem_sup_right (Submodule.mem_iSup_of_mem v (Submodule.mem_map_of_mem hlow))

theorem sharp_gen (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K) (𝔭 : Ideal (Q K))
    (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (hD : 3 ≤ D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (normalized k depth a b t) ∈ sharpSpace two_ne 𝔭 k D := by
  have hdeg : (normalized k depth a b t).totalDegree ≤ D - 1 := by
    have h := (PrimeM2General.degree_and_top k depth a b t (Fin.last (2 * k))).1
    simp only [inputs, Fin.lastCases_last] at h
    omega
  exact Submodule.mem_sup_left (normalized_mem_idealSpace two_ne k depth a b t 𝔭 gen (D - 1) hdeg)

theorem one_not_mem_sharp (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) (h1 : (1 : Q K) ∉ 𝔭)
    (k D : ℕ) (crit : Crit two_ne 𝔭 k D) : (1 : P (2 * k) K) ∉ sharpSpace two_ne 𝔭 k D := by
  intro h
  have hlow := crit 1 h (by simp)
  have := hlow.2.2
  rw [map_one] at this
  exact h1 this

/-- **Proposition 7.1 (low-part criterion).** If `S_3 ∈ 𝔭`, `1 ∉ 𝔭`, `D ≥ 3` and `Crit_D(k)`
holds, then `1` is not in the degree-`D` closure of the two-summand chain presentation. -/
theorem not_one_of_crit (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) (h1 : (1 : Q K) ∉ 𝔭) (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (hD : 3 ≤ D)
    (crit : Crit two_ne 𝔭 k D) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : P (2 * k) K) :=
  ChainClosureBridge.prime_m2_not_one k depth a b t two_ne D (sharpSpace two_ne 𝔭 k D)
    (fun _ => sharp_gen two_ne k depth a b t 𝔭 gen D hD)
    (sharp_products two_ne 𝔭 k D crit)
    (one_not_mem_sharp two_ne 𝔭 h1 k D crit)

/-- **Proposition 7.1, `V_D(k) ⊂ S^♯`.** Under `Crit_D(k)` the chain normal form of every
element of the degree-`D` closure lies in `S^♯`; if it has degree `≤ D - 1` it lies in
`S_{D-1}(k)`. -/
theorem closure_normalForm_mem_sharpSpace (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k)
    (a b t : K) (𝔭 : Ideal (Q K)) (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (hD : 3 ≤ D)
    (crit : Crit two_ne 𝔭 k D) {p : P (2 * k) K}
    (h : InClosure (Set.range (inputs k depth a b t)) D p) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) p ∈
      sharpSpace two_ne 𝔭 k D := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 k two_ne
  have hq := ChainClosureBridge.closure_normalForm_mem (nonzero := chainInputs_nonzero 2 k) zero h
  refine ChainClosureBridge.quot_closure_minimal _ D (sharpSpace two_ne 𝔭 k D) ?_
    (sharp_products two_ne 𝔭 k D crit) hq
  rintro f ⟨i, rfl⟩ hd
  revert hd
  refine Fin.lastCases ?_ (fun j => ?_) i <;> intro hd
  · simp only [inputs, Fin.lastCases_last] at hd ⊢
    exact sharp_gen two_ne k depth a b t 𝔭 gen D hD
  · rw [ChainClosureBridge.normalForm_chain_eq_zero 2 k
      (ChainClosureBridge.inputs_castSucc_mem_chainInputs k depth a b t j) two_ne]
    exact (sharpSpace two_ne 𝔭 k D).zero_mem

/-- `FEI_k(D-1)` implies `Crit_D(k)` (`D ≥ 1`). -/
theorem crit_of_fei (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) [𝔭.IsPrime]
    (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (k D : ℕ) (hD : 1 ≤ D) (fei : FEI 𝔭 k (D - 1)) :
    Crit two_ne 𝔭 k D := by
  have mono : idealSpace 𝔭 k (D - 1) ≤ idealSpace 𝔭 k D := by
    rintro w ⟨h1, h2, h3⟩
    exact ⟨h1, h2.trans (Nat.sub_le D 1), h3⟩
  have le : sharpSpace two_ne 𝔭 k D ≤ idealSpace 𝔭 k D := by
    apply sup_le mono
    apply iSup_le
    intro v
    rw [Submodule.map_le_iff_le_comap]
    intro w hw
    exact idealSpace_closed two_ne 𝔭 hX k D fei v w (mono hw)
      (lt_of_le_of_lt hw.2.1 (by omega))
  intro w hw hdeg
  have := le hw
  exact ⟨this.1, hdeg, this.2.2⟩

end ICFallDegrees.GrowthCriterion
