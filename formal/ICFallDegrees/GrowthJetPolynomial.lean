import ICFallDegrees.GrowthReduction
import Mathlib.Algebra.MvPolynomial.PDeriv
import ICFallDegrees.GrowthRationalFrobenius
import ICFallDegrees.UnitInstance
import ICFallDegrees.GrowthPrimeE6NativeNonmembership
import Mathlib.Tactic

noncomputable section
namespace ICFallDegrees.GrowthJetPolynomial

open MvPolynomial FiniteQuotient ChainQuotient ChainTranslation ChainModelMap
open GrowthReduction
variable {K : Type*} [Field K]

/-- The actual polynomial Euler derivative, not an abstract coefficient operator. -/
def euler (i : Fin 2) : Q K →ₗ[K] Q K :=
  (LinearMap.mulLeft K (X i)).comp (pderiv i).toLinearMap

theorem euler_apply (i : Fin 2) (f : Q K) :
    euler i f = X i * pderiv i f := rfl

theorem euler_monomial (i : Fin 2) (d : Fin 2 →₀ ℕ) (c : K) :
    euler i (monomial d c) = d i • monomial d c :=
  X_mul_pderiv_monomial

theorem coeff_euler (i : Fin 2) (f : Q K) (d : Fin 2 →₀ ℕ) :
    coeff d (euler i f) = (d i : K) * coeff d f := by
  classical
  conv_lhs => rw [f.as_sum]
  rw [map_sum, coeff_sum]
  simp only [euler_monomial, coeff_smul, smul_eq_mul]
  by_cases hd : d ∈ f.support
  · rw [Finset.sum_eq_single d]
    · simp
    · intro b _ hbd
      rw [coeff_monomial, if_neg hbd, smul_zero]
    · exact fun h => False.elim (h hd)
  · rw [not_mem_support_iff.mp hd, mul_zero]
    apply Finset.sum_eq_zero
    intro b hb
    rw [coeff_monomial, if_neg (by intro h; exact hd (h ▸ hb)), smul_zero]

theorem boxed_euler (i : Fin 2) (k D : ℕ) (f : Q K) (hf : Boxed k D f) :
    Boxed k D (euler i f) := by
  intro d hd
  apply hf d
  rw [mem_support_iff] at hd ⊢
  rw [coeff_euler] at hd
  exact fun h => hd (by rw [h, mul_zero])

theorem euler_topPart (i l : Fin 2) (k j : ℕ) (f : Q K) :
    euler i (topPart k l j f) = topPart k l j (euler i f) := by
  ext d
  rw [coeff_euler, coeff_topPart, coeff_topPart, coeff_euler]
  split_ifs <;> simp

theorem euler_mul (i : Fin 2) (f g : Q K) :
    euler i (f * g) = euler i f * g + f * euler i g := by
  simp only [euler_apply, Derivation.leibniz, smul_eq_mul]
  ring

/-- Differentiated lift identity: FEI suffices to preserve vanishing of the
Euler derivative after an admissible chain shift. No separate wrapping
vanishing assumption is introduced here. FEI itself remains a hypothesis. -/
theorem derivative_shift_mem (two_ne : (2 : K) ≠ 0)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (k D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (i l : Fin 2) (j : Fin k) (w : P (2 * k) K)
    (hw : w ∈ idealSpace 𝔭 k D) (hdeg : w.totalDegree < D)
    (hderiv : euler l (expand K 2 k w) ∈ 𝔭) :
    euler l (expand K 2 k
      (ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
        (X (var 2 k i j) * w))) ∈ 𝔭 := by
  let f := expand K 2 k w
  let t := topPart k i j f
  let z := expand K 2 k
    (ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (X (var 2 k i j) * w))
  have hf : f ∈ 𝔭 := hw.2.2
  have boxed : Boxed k (D - 1) f := boxed_expand w hw.1 (D - 1) (by omega)
  have ht : t ∈ 𝔭 := fei f hf boxed i j j.isLt
  have hdt : euler l t ∈ 𝔭 := by
    rw [show t = topPart k i j f from rfl, euler_topPart]
    exact fei _ hderiv (boxed_euler l k (D - 1) f boxed) i j j.isLt
  have hz : z ∈ 𝔭 := by
    have h := idealSpace_closed two_ne 𝔭 hX k D fei
      (var 2 k i j) w hw hdeg
    exact h.2.2
  have ident := lift_identity two_ne w hw.1 i j
  change X i ^ (2 ^ k - 2 ^ (j : ℕ)) * z = X i ^ 2 ^ k * (f - t) + t at ident
  have rhs : euler l (X i ^ 2 ^ k * (f - t) + t) ∈ 𝔭 := by
    rw [map_add, euler_mul, map_sub]
    exact 𝔭.add_mem
      (𝔭.add_mem (𝔭.mul_mem_left _ (𝔭.sub_mem hf ht))
        (𝔭.mul_mem_left _ (𝔭.sub_mem hderiv hdt))) hdt
  rw [← ident, euler_mul] at rhs
  have prod : X i ^ (2 ^ k - 2 ^ (j : ℕ)) * euler l z ∈ 𝔭 := by
    have hh := 𝔭.sub_mem rhs (𝔭.mul_mem_left (euler l (X i ^ (2 ^ k - 2 ^ (j : ℕ)))) hz)
    simpa only [add_sub_cancel_left] using hh
  rcases Ideal.IsPrime.mem_or_mem inferInstance prod with h | h
  · exact False.elim (hX i (Ideal.IsPrime.mem_of_pow_mem inferInstance _ h))
  · exact h

/-- The first-jet invariant on the actual reduced chain polynomial space. -/
def jetSpace (𝔭 : Ideal (Q K)) (k D : ℕ) (l : Fin 2) :
    Submodule K (P (2 * k) K) :=
  idealSpace 𝔭 k D ⊓
    (𝔭.restrictScalars K).comap ((euler l).comp (expand K 2 k).toLinearMap)

theorem jetSpace_closed (two_ne : (2 : K) ≠ 0)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (k D : ℕ) (fei : FEI 𝔭 k (D - 1)) (l : Fin 2)
    (v : Fin (2 * k)) (w : P (2 * k) K)
    (hw : w ∈ jetSpace 𝔭 k D l) (hdeg : w.totalDegree < D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (X v * w) ∈ jetSpace 𝔭 k D l := by
  constructor
  · exact idealSpace_closed two_ne 𝔭 hX k D fei v w hw.1 hdeg
  · have id := derivative_shift_mem two_ne 𝔭 hX k D fei
      (row v) l (depth v) w hw.1 hdeg hw.2
    rw [var_row_depth] at id
    exact id

theorem normalized_mem_jetSpace (two_ne : (2 : K) ≠ 0)
    (k : ℕ) (depth : 2 ≤ k) (a b t : K) (𝔭 : Ideal (Q K))
    (gen : S3poly a b t ∈ 𝔭) (l : Fin 2)
    (dgen : euler l (S3poly a b t) ∈ 𝔭) (D : ℕ)
    (hD : (PrimeM2General.normalized k depth a b t).totalDegree ≤ D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (PrimeM2General.normalized k depth a b t) ∈ jetSpace 𝔭 k D l := by
  constructor
  · exact normalized_mem_idealSpace two_ne k depth a b t 𝔭 gen D hD
  · change euler l (expand K 2 k
      (ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
        (PrimeM2General.normalized k depth a b t))) ∈ 𝔭
    rw [normalForm_normalized two_ne k depth a b t, expand_normalized]
    exact dgen

/-- Actual InClosure consequence. No unit-ideal hypothesis is used. -/
theorem closure_normalForm_mem_jetSpace (two_ne : (2 : K) ≠ 0)
    (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (gen : S3poly a b t ∈ 𝔭) (l : Fin 2)
    (dgen : euler l (S3poly a b t) ∈ 𝔭)
    (D : ℕ) (fei : FEI 𝔭 k (D - 1)) {f : P (2 * k) K}
    (hf : InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D f) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) f ∈
      jetSpace 𝔭 k D l := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 k two_ne
  have hq := ChainClosureBridge.closure_normalForm_mem
    (nonzero := chainInputs_nonzero 2 k) zero hf
  refine ChainClosureBridge.quot_closure_minimal _ D (jetSpace 𝔭 k D l) ?_
    (fun v w hw hd => jetSpace_closed two_ne 𝔭 hX k D fei l v w hw hd) hq
  rintro f ⟨i, rfl⟩ hd
  revert hd
  refine Fin.lastCases ?_ (fun j => ?_) i <;> intro hd
  · simp only [PrimeM2General.inputs, Fin.lastCases_last] at hd ⊢
    exact normalized_mem_jetSpace two_ne k depth a b t 𝔭 gen l dgen D hd
  · rw [ChainClosureBridge.normalForm_chain_eq_zero 2 k
      (ChainClosureBridge.inputs_castSucc_mem_chainInputs k depth a b t j) two_ne]
    exact (jetSpace 𝔭 k D l).zero_mem

/-- A derivative witness refutes every LFD bound at most D without
requiring that the input ideal be the unit ideal. The ideal witness and
FEI must both be established by the application. -/
theorem lfd_gt_of_derivative_witness (two_ne : (2 : K) ≠ 0)
    (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (gen : S3poly a b t ∈ 𝔭) (l : Fin 2)
    (dgen : euler l (S3poly a b t) ∈ 𝔭)
    (D : ℕ) (fei : FEI 𝔭 k (D - 1)) (w : P (2 * k) K)
    (hw : w ∈ Ideal.span (Set.range (PrimeM2General.inputs k depth a b t)))
    (hdeg : w.totalDegree ≤ D)
    (outside : euler l (expand K 2 k
      (ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) w)) ∉ 𝔭) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) e := by
  intro e he bound
  have derived := (bound w hw).mono (max_le he hdeg)
  have h := closure_normalForm_mem_jetSpace two_ne k depth a b t 𝔭 hX gen l dgen D fei derived
  exact outside h.2

/-- The E6 square root in the two unwrapped coordinates. -/
def e6Factor : Q K := X 0 * X 1 - X 0 - X 1 - 8

theorem e6_square : S3poly (6 : K) (-7) 1 = (e6Factor : Q K) ^ 2 := by
  norm_num [S3poly, e6Factor]
  simp only [map_ofNat, map_natCast, map_intCast]
  ring

/-- Both input jet conditions follow from membership of the factor. -/
theorem e6_input_jet (𝔭 : Ideal (Q K)) (h : (e6Factor : Q K) ∈ 𝔭)
    (l : Fin 2) : S3poly (6 : K) (-7) 1 ∈ 𝔭 ∧
      euler l (S3poly (6 : K) (-7) 1) ∈ 𝔭 := by
  rw [e6_square, pow_two, euler_mul]
  exact ⟨𝔭.mul_mem_left _ h,
    𝔭.add_mem (𝔭.mul_mem_left _ h) (𝔭.mul_mem_right _ h)⟩

/-- The factor in the actual chain variables. -/
def e6Witness (k : ℕ) (depth : 2 ≤ k) : P (2 * k) K :=
  X (PrimeM2General.active k depth 0) * X (PrimeM2General.active k depth 2) -
    X (PrimeM2General.active k depth 0) - X (PrimeM2General.active k depth 2) - 8

/-- Frobenius gives actual ideal membership of the low-degree factor,
even at proper-ideal exceptional primes. It gives no degree budget. -/
theorem e6Witness_mem_ideal (p k : ℕ) [Fact p.Prime] [CharP K p]
    (depth : 2 ≤ k) (hdvd : 2 ^ k ∣ p - 1) :
    (e6Witness k depth : P (2 * k) K) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k depth (6 : K) (-7) 1)) := by
  let I := Ideal.span (Set.range (PrimeM2General.inputs k depth (6 : K) (-7) 1))
  let π := Ideal.Quotient.mk I
  let x : P (2 * k) K := X (PrimeM2General.active k depth 0)
  let y : P (2 * k) K := X (PrimeM2General.active k depth 2)
  have root (v : Fin (2 * k)) : π (X v) ^ (2 ^ k) = 1 := by
    have hh := UnitInstance.var_pow_mem k depth (6 : K) (-7) 1 v
    have hz := (Ideal.Quotient.eq_zero_iff_mem).mpr hh
    simpa only [map_sub, map_pow, map_one, sub_eq_zero] using hz
  have xp := GrowthRationalFrobenius.pow_fixed_of_root (2 ^ k) p
    (Fact.out : p.Prime).one_le hdvd (π x) (root _)
  have yp := GrowthRationalFrobenius.pow_fixed_of_root (2 ^ k) p
    (Fact.out : p.Prime).one_le hdvd (π y) (root _)
  have eight : (8 : P (2 * k) K) ^ p = 8 := by
    exact (map_ofNat (frobenius (P (2 * k) K) p) 8)
  have fixed : π (e6Witness k depth : P (2 * k) K) ^ p = π (e6Witness k depth) := by
    rw [← map_pow]
    change π ((x * y - x - y - 8) ^ p) = π (x * y - x - y - 8)
    rw [sub_pow_char, sub_pow_char, sub_pow_char, mul_pow, eight]
    simp only [map_sub, map_mul, map_pow, xp, yp]
  have step0 := UnitInstance.chain_step_mem k depth (6 : K) (-7) 1
    (PrimeM2General.active k depth 0) (by simp [PrimeM2General.active]; omega)
  have step2 := UnitInstance.chain_step_mem k depth (6 : K) (-7) 1
    (PrimeM2General.active k depth 2) (by simp [PrimeM2General.active]; omega)
  have active1 : (⟨(PrimeM2General.active k depth 0).val + 1, by
      simp [PrimeM2General.active]; omega⟩ : Fin (2*k)) = PrimeM2General.active k depth 1 := by
    apply Fin.ext; simp [PrimeM2General.active]
  have active3 : (⟨(PrimeM2General.active k depth 2).val + 1, by
      simp [PrimeM2General.active]; omega⟩ : Fin (2*k)) = PrimeM2General.active k depth 3 := by
    apply Fin.ext; simp [PrimeM2General.active]
  rw [active1] at step0
  rw [active3] at step2
  have rawmem : PrimeM2General.rawS3 k depth (6 : K) (-7) 1 ∈ I := by
    have id := PrimeM2General.raw_normalization_identity k depth (6 : K) (-7) 1
    have hh : PrimeM2General.rawS3 k depth (6 : K) (-7) 1 -
        PrimeM2General.normalized k depth (6 : K) (-7) 1 ∈ I := by
      rw [id]
      exact I.add_mem (I.mul_mem_left _ step0) (I.mul_mem_left _ step2)
    have ng : PrimeM2General.normalized k depth (6 : K) (-7) 1 ∈ I :=
      Ideal.subset_span ⟨Fin.last _, by simp [PrimeM2General.inputs]⟩
    simpa only [sub_add_cancel] using I.add_mem hh ng
  have sq : (e6Witness k depth : P (2*k) K) ^ 2 ∈ I := by
    convert rawmem using 1
    norm_num [e6Witness, PrimeM2General.rawS3]
    simp only [map_ofNat, map_natCast, map_intCast]
    ring
  have zero : π (e6Witness k depth : P (2*k) K) ^ 2 = 0 := by
    rw [← map_pow]
    exact Ideal.Quotient.eq_zero_iff_mem.mpr sq
  have factored := congrArg π (pow_add (e6Witness k depth : P (2*k) K) 2 (p-2))
  have hp : 2 + (p - 2) = p := by have := (Fact.out : p.Prime).two_le; omega
  rw [hp, map_pow, fixed, map_mul, map_pow, zero, zero_mul] at factored
  exact Ideal.Quotient.eq_zero_iff_mem.mp factored

theorem e6Witness_reduced (k : ℕ) (depth : 2 ≤ k) :
    (e6Witness k depth : P (2*k) K) ∈ reducedSpace (chainInputs K 2 k) := by
  let S := reducedSpace (chainInputs K 2 k)
  have hxy : (X (PrimeM2General.active k depth 0) * X (PrimeM2General.active k depth 2) : P (2*k) K) ∈ S :=
    X_mul_X_reduced ((PrimeM2General.active_injective k depth).ne (by decide))
  have hx : (X (PrimeM2General.active k depth 0) : P (2*k) K) ∈ S := X_reduced _
  have hy : (X (PrimeM2General.active k depth 2) : P (2*k) K) ∈ S := X_reduced _
  have h8 : (8 : P (2*k) K) ∈ S := by
    simpa only [map_ofNat] using (C_reduced (m := 2) (k := k) (8 : K))
  exact S.sub_mem (S.sub_mem (S.sub_mem hxy hx) hy) h8

theorem expand_e6Witness (k : ℕ) (depth : 2 ≤ k) :
    expand K 2 k (e6Witness k depth) = (e6Factor : Q K) := by
  obtain ⟨hx, _, hy, _⟩ := expand_active (K := K) k depth
  simp only [e6Witness, e6Factor, map_sub, map_mul, map_ofNat, hx, hy]

theorem e6Witness_degree (k : ℕ) (depth : 2 ≤ k) :
    (e6Witness k depth : P (2*k) K).totalDegree ≤ 2 := by
  unfold e6Witness
  apply (totalDegree_sub_le _ _).trans
  apply max_le
  · apply (totalDegree_sub_le _ _).trans
    apply max_le
    · apply (totalDegree_sub_le _ _).trans
      apply max_le
      · exact (totalDegree_mul _ _).trans (by simp)
      · simp
    · simp
  · rw [show (8 : P (2*k) K) = C (8 : K) by simp only [map_ofNat]]
    simp

theorem euler_e6Factor : euler 1 (e6Factor : Q K) = X 1 * (X 0 - 1) := by
  have h8 : (pderiv (1 : Fin 2)) (8 : Q K) = 0 := by
    rw [show (8 : Q K) = C (8 : K) by simp only [map_ofNat], pderiv_C]
  simp [euler_apply, e6Factor, Derivation.leibniz, smul_eq_mul, h8]

/-- The witness derivative is nonzero on the actual rational-function model. -/
theorem e6_derivative_outside :
    euler 1 (e6Factor : Q K) ∉ RingHom.ker
      (aeval (GrowthRatFuncModel.xi (6 : K)) : Q K →ₐ[K] RatFunc K) := by
  intro h
  change aeval (GrowthRatFuncModel.xi (6 : K)) (euler 1 (e6Factor : Q K)) = 0 at h
  rw [euler_e6Factor] at h
  simp only [map_mul, map_sub, map_one, aeval_X, GrowthRatFuncModel.xi,
    show (1 : Fin 2) ≠ 0 by decide, if_false, if_pos rfl] at h
  simp only [ite_true] at h
  unfold GrowthRationalCoordinates.rationalZ at h
  have hd : (RatFunc.X : RatFunc K) - 1 ≠ 0 := sub_ne_zero.mpr GrowthRatFuncModel.X_ne_one
  rw [div_mul_cancel₀ _ hd] at h
  exact GrowthRatFuncModel.numerator_ne_zero (6 : K) h

/-- Explicit E6 lower bound over every sufficiently large admissible prime,
including proper-ideal exceptional primes. No unit-ideal hypothesis and
no supplied CurveModel, FEI or spread hypothesis appears in this statement. -/
theorem e6_prime_lfd_gt (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (D : ℕ) (hD : 2 ≤ D) (hk : 14 * D + 10 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p) (hdvd : 2 ^ k ∣ p - 1) :
    ∀ e ≤ D, ¬ HasLFDBound
      (Set.range (PrimeM2General.inputs k depth (6 : ZMod p) (-7) 1)) e := by
  let ξ := GrowthRatFuncModel.xi (6 : ZMod p)
  let 𝔭 := RingHom.ker (aeval ξ : Q (ZMod p) →ₐ[ZMod p] RatFunc (ZMod p))
  have data := GrowthRatFuncModel.ratFunc_algebraic_data (6 : ZMod p)
  have curve : GrowthSpread.S3at (6 : ZMod p) (-7) 1 ξ = 0 := by
    simpa only [ξ, show (-1 - (6 : ZMod p)) = (-7 : ZMod p) by ring] using data.2
  have gen : S3poly (6 : ZMod p) (-7) 1 ∈ 𝔭 := GrowthMain.S3poly_mem_ker curve
  have factor : (e6Factor : Q (ZMod p)) ∈ 𝔭 := by
    rw [e6_square] at gen
    exact Ideal.IsPrime.mem_of_pow_mem inferInstance _ gen
  have dgen := (e6_input_jet 𝔭 factor 1).2
  have hlarge : 5 ≤ p := by
    have hpow : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
    omega
  have htwo : (2 : ZMod p) ≠ 0 := PrimeFieldConclusions.two_ne_zero p hlarge
  have spread := GrowthPrimeE6NativeNonmembership.e6_native_box_spread_large_prime
    p k (by omega) hp
  have fei : FEI 𝔭 k (D - 1) := by
    apply GrowthBoxSpread.fei_of_box_spread
      (K := ZMod p) (F₀ := ZMod p) (F := RatFunc (ZMod p))
      ξ p k (D - 1) (13 * (D - 1) + 23) (Or.inr hp)
      spread (GrowthCounting.lambdaBound_linear (D - 1)) (by omega) (by omega)
  apply lfd_gt_of_derivative_witness htwo k depth (6 : ZMod p) (-7) 1
    𝔭 (GrowthMain.X_not_mem_ker data.1) gen 1 dgen D fei (e6Witness k depth)
    (e6Witness_mem_ideal p k depth hdvd) ((e6Witness_degree k depth).trans hD)
  rw [ChainClosureBridge.normalForm_of_reduced _ _
    (reduced_chain_ideal_trivial (K := ZMod p) 2 k htwo) (e6Witness_reduced k depth),
    expand_e6Witness]
  exact e6_derivative_outside

/-- The explicit linear lower bound at every admissible prime, including
proper-ideal exceptional primes: LFD > floor((k-10)/14), for k>=38. -/
theorem e6_prime_lfd_linear (p k : ℕ) [Fact p.Prime] (hk : 38 ≤ k)
    (admissible : 2 ^ (k + 1) ∣ p - 1) :
    ∀ e ≤ (k - 10) / 14, ¬ HasLFDBound
      (Set.range (PrimeM2General.inputs k (by omega) (6 : ZMod p) (-7) 1)) e := by
  have hprime : 2 ≤ p := (Fact.out : p.Prime).two_le
  have hdvd : 2 ^ k ∣ p - 1 := by
    apply dvd_trans _ admissible
    rw [pow_succ]
    exact dvd_mul_right _ _
  have hbound : 2 ^ (k + 1) ≤ p - 1 := Nat.le_of_dvd (by omega) admissible
  rw [pow_succ] at hbound
  have hp : 2 * 2 ^ k + 1 ≤ p := by omega
  exact e6_prime_lfd_gt p k (by omega) ((k - 10) / 14) (by omega) (by omega) hp hdvd

#print axioms coeff_euler
#print axioms boxed_euler
#print axioms euler_topPart
#print axioms euler_mul
#print axioms derivative_shift_mem
#print axioms jetSpace_closed
#print axioms normalized_mem_jetSpace
#print axioms closure_normalForm_mem_jetSpace
#print axioms lfd_gt_of_derivative_witness
#print axioms e6_square
#print axioms e6_input_jet
#print axioms e6Witness_mem_ideal
#print axioms e6Witness_reduced
#print axioms expand_e6Witness
#print axioms e6Witness_degree
#print axioms euler_e6Factor
#print axioms e6_derivative_outside
#print axioms e6_prime_lfd_gt
#print axioms e6_prime_lfd_linear
end ICFallDegrees.GrowthJetPolynomial
