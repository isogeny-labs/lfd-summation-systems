import ICFallDegrees.ShiftWitness
import ICFallDegrees.N5CertSound

noncomputable section

/-!
# N7: the growth schema with the explicit generator vector

`genVec a b t` is the exponent-model vector of the normalized two-summand generator
`PrimeM2General.normalized k hk2 a b t` (nine monomials `x₀^i x₁^j`, `i, j ≤ 2`), the same at
every depth.  `Φ_genVec` identifies its image under the model map with the generator at every
depth `k ≥ 2` and every base depth `k0`.  With it the hypotheses of
`ShiftWitness.lfd_gt_of_lowIdeal` become depth-independent, which gives the

**growth schema** (`growth_schema`): for every degree `D` and base depth `k0 ≥ 2`, a subspace
`S` of the depth-`k0` exponent model such that

* the free extension `W_{k0+1}(S)` is closed (`Closed k0 (k0+1) D (W k0 D S (k0+1))`),
* the generator vector lies in `S`,
* some low ideal element `v` (`LowIdeal k0 (genVec a b t) v`) of weight `≤ D` is not in `S`,

forces `LFD > D` at every depth `k ≥ k0+1`: no `e ≤ D` is an LFD bound.  This is the formal
statement behind "each degree `D` is eventually exceeded".  It is conditional on a base for
each `D`: the theorem constructs none, and only the base `D = 4`, `k0 = 7` is certified (in
`N5CertWitnessInstance*`).
-/

namespace ICFallDegrees.GrowthSchema

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight DepthReduction ChainTranslation
open ChainModelMap ChainClosureBridge PrimeM2General GeneralQuadraticTop ShiftWitness

variable {K : Type*} [Field K]

/-- The basis vector of a low exponent. -/
def δ (l : Fin 2 → ℕ) : Lo 2 K := Finsupp.single l 1

/-- The exponent-model vector of the normalized two-summand generator. -/
def genVec (a b t : K) : Lo 2 K :=
  (-2 * (t ^ 2 + a)) • δ ![1, 1] + (-2 * t) • δ ![1, 2] + (-2 * t) • δ ![2, 1] + δ ![2, 2] +
    (t ^ 2) • δ ![2, 0] + (t ^ 2) • δ ![0, 2] +
    (-2 * a * t - 4 * b) • δ ![1, 0] + (-2 * a * t - 4 * b) • δ ![0, 1] +
    (a ^ 2 - 4 * b * t) • δ ![0, 0]

theorem Φ_embed_δ (k0 k : ℕ) (l : Fin 2 → ℕ) :
    Φ (K := K) k0 k (embed 0 (δ l)) = monomial (digits k l) 1 := by
  have h1 : embed (K := K) 0 (δ l) = Finsupp.single (l, (0 : Fin 2 → ℕ)) 1 := by
    unfold embed δ
    exact Finsupp.mapDomain_single
  rw [h1, Φ_single]
  have : encode k0 (l, (0 : Fin 2 → ℕ)) = l := by
    funext i
    simp [encode]
  rw [this]

theorem var00 {k : ℕ} (hk : 2 ≤ k) : var 2 k 0 ⟨0, by omega⟩ = active k hk 0 :=
  Fin.ext (by simp [var_val, active])
theorem var01 {k : ℕ} (hk : 2 ≤ k) : var 2 k 0 ⟨1, by omega⟩ = active k hk 1 :=
  Fin.ext (by simp [var_val, active])
theorem var10 {k : ℕ} (hk : 2 ≤ k) : var 2 k 1 ⟨0, by omega⟩ = active k hk 2 :=
  Fin.ext (by simp [var_val, active])
theorem var11 {k : ℕ} (hk : 2 ≤ k) : var 2 k 1 ⟨1, by omega⟩ = active k hk 3 :=
  Fin.ext (by simp [var_val, active]; omega)

/-- The generator vector maps to the normalized generator at every depth. -/
theorem Φ_genVec (k0 k : ℕ) (hk2 : 2 ≤ k) (a b t : K) :
    Φ (K := K) k0 k (embed 0 (genVec a b t)) = normalized k hk2 a b t := by
  simp only [genVec, map_add, map_smul, Φ_embed_δ]
  rw [N5CertSound.digits_two_bits hk2 ![1, 1] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![1, 2] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![2, 1] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![2, 2] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![2, 0] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![0, 2] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![1, 0] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![0, 1] (by decide) (by decide),
    N5CertSound.digits_two_bits hk2 ![0, 0] (by decide) (by decide)]
  simp only [var00 hk2, var01 hk2, var10 hk2, var11 hk2]
  norm_num [bit]
  simp only [normalized, top, tail, smul_eq_C_mul, monomial_single_add,
    ← X_pow_eq_monomial, pow_one, mul_one]
  simp only [map_add, map_sub, map_mul, map_pow, map_neg, map_ofNat]
  ring

/-- The property shared by the nine monomials of the generator. -/
def Small (l : Fin 2 → ℕ) : Prop := l 0 < 4 ∧ l 1 < 4 ∧ lowWt l ≤ 2

theorem small_mk {i j : ℕ} (hi : i < 3) (hj : j < 3) : Small ![i, j] := by
  refine ⟨by simpa using by omega, by simpa using by omega, ?_⟩
  unfold lowWt
  rw [Fin.sum_univ_two]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
  have hw : ∀ n : ℕ, n < 3 → wt n ≤ 1 := by
    intro n hn
    interval_cases n
    · simp
    · rw [wt_one]
    · rw [wt_eq, wt_one]
  have := hw i hi
  have := hw j hj
  omega

theorem small_δ {l0 : Fin 2 → ℕ} (h : Small l0) : ∀ l ∈ (δ (K := K) l0).support, Small l := by
  intro l hl
  have := Finsupp.support_single_subset hl
  rw [Finset.mem_singleton] at this
  rwa [this]

theorem small_smul (c : K) {v : Lo 2 K} (h : ∀ l ∈ v.support, Small l) :
    ∀ l ∈ (c • v).support, Small l :=
  fun l hl => h l (Finsupp.support_smul hl)

theorem small_add {v v' : Lo 2 K} (h : ∀ l ∈ v.support, Small l)
    (h' : ∀ l ∈ v'.support, Small l) : ∀ l ∈ (v + v').support, Small l := by
  classical
  intro l hl
  rcases Finset.mem_union.mp (Finsupp.support_add hl) with h1 | h1
  · exact h l h1
  · exact h' l h1

theorem genVec_support (a b t : K) (l : Fin 2 → ℕ) (hl : l ∈ (genVec a b t).support) :
    l 0 < 4 ∧ l 1 < 4 ∧ lowWt l ≤ 2 := by
  have key : ∀ l ∈ (genVec a b t).support, Small l := by
    unfold genVec
    exact small_add (small_add (small_add (small_add (small_add (small_add (small_add (small_add
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_δ (small_mk (by norm_num) (by norm_num))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num)))))
      (small_smul _ (small_δ (small_mk (by norm_num) (by norm_num))))
  exact key l hl

/-- **Growth schema.**  For every degree `D` and base depth `k0 ≥ 2`: a base `S` whose free
extension is closed at depth `k0+1`, which contains the generator vector, and which misses one
low ideal element of weight `≤ D`, gives `LFD > D` at every depth `k ≥ k0+1`, for the
two-summand presentation with parameters `a b t` over any field of characteristic `≠ 2`.
Conditional on the base: none is constructed here. -/
theorem growth_schema (two_ne : (2 : K) ≠ 0) (a b t : K) (k0 D : ℕ) (hk0 : 2 ≤ k0)
    (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (gen : genVec a b t ∈ S)
    (v : Lo 2 K) (ideal : LowIdeal k0 (genVec a b t) v)
    (weight : ∀ l ∈ v.support, lowWt l ≤ D) (outside : v ∉ S) :
    ∀ k (hk2 : 2 ≤ k), k0 + 1 ≤ k →
      ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  intro k hk2 hk e he bound
  have h4 : (4 : ℕ) ≤ 2 ^ k0 := by
    calc (4 : ℕ) = 2 ^ 2 := by norm_num
      _ ≤ 2 ^ k0 := Nat.pow_le_pow_right (by norm_num) hk0
  have hgbox : ∀ l ∈ (genVec a b t).support, LowBox k0 l := by
    intro l hl i
    obtain ⟨h0, h1, _⟩ := genVec_support a b t l hl
    fin_cases i
    · exact lt_of_lt_of_le h0 h4
    · exact lt_of_lt_of_le h1 h4
  by_cases hD : 2 ≤ D
  · exact lfd_gt_of_lowIdeal k0 D S base k hk hk2 two_ne a b t (genVec a b t) gen hgbox
      (fun l hl => (genVec_support a b t l hl).2.2.trans hD) (Φ_genVec k0 k hk2 a b t)
      v ideal weight outside e he bound
  · -- `D ≤ 1`: the generator is not an admissible input, so the base hypothesis on it is void
    have gen' : GeneratorInBase' k0 D S k hk2 a b t := by
      intro hd
      have := (degree_and_top k hk2 a b t (Fin.last _)).1
      simp only [inputs, Fin.lastCases_last] at this
      omega
    have hvbox := LowIdeal.lowBox hgbox ideal
    refine lfd_gt_of_base_witness k0 D S base k hk hk2 two_ne a b t gen' hvbox weight outside ?_
      e he bound
    refine lowIdeal_mem_ideal two_ne (by omega) _ ?_ hgbox ?_ ideal
    · exact (chainInputs_two_subset k hk2 a b t).trans Ideal.subset_span
    · rw [Φ_genVec k0 k hk2 a b t]
      apply Ideal.subset_span
      exact ⟨Fin.last _, by simp [inputs]⟩

end ICFallDegrees.GrowthSchema
