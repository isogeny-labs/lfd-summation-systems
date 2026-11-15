import ICFallDegrees.GrowthTwoTorsionFamily
import Mathlib.RingTheory.Polynomial.RationalRoot

/-! # Rational parameters in the two-torsion growth family -/

noncomputable section
namespace ICFallDegrees.GrowthRationalParameterClassification

open MvPolynomial
open ICFallDegrees.GrowthSecant ICFallDegrees.GrowthRationalUnit
open ICFallDegrees.GrowthCharZeroComplexUnit

theorem root_integral (k : ℕ) (x : ℂ) (hx : x ^ (2 ^ k) = 1) :
    IsIntegral ℤ x := by
  apply IsIntegral.of_pow (n := 2 ^ k) (by positivity)
  rw [hx]
  exact isIntegral_one

theorem cubic_root_of_two_power (k : ℕ) (x : ℂ)
    (hx : x ^ (2 ^ k) = 1) (h3 : x ^ 3 = 1) : x = 1 := by
  have hd := Nat.dvd_gcd (orderOf_dvd_of_pow_eq_one hx)
    (orderOf_dvd_of_pow_eq_one h3)
  have hc : Nat.Coprime (2 ^ k) 3 := (by decide : Nat.Coprime 2 3).pow_left k
  rw [hc.gcd_eq_one] at hd
  exact orderOf_eq_one_iff.mp (Nat.dvd_one.mp hd)

theorem trace_of_s3_zero (k : ℕ) (a : ℚ) (x z : ℂ)
    (ha : a ≠ -3) (_hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1)
    (h : s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 = 0) :
    x + (starRingEnd ℂ) x = ((-a - 1 : ℚ) : ℂ) := by
  have hp := parameter_of_s3_zero (a : ℂ) x z h
  have hlin : (x - 1) * z = x + 2 + (a : ℂ) := by linear_combination -hp
  have hzn := norm_eq_one_of_root k z hz
  have hs : Complex.normSq z = 1 := by
    rw [Complex.normSq_eq_norm_sq, hzn]; norm_num
  have hn := congrArg Complex.normSq hlin
  rw [Complex.normSq_mul, hs, mul_one] at hn
  simp [Complex.normSq_apply] at hn
  have hprod : ((a : ℝ) + 3) * ((a : ℝ) + 1 + 2 * x.re) = 0 := by
    nlinarith [hn]
  have hne : (a : ℝ) + 3 ≠ 0 := by
    intro he
    apply ha
    exact_mod_cast (show (a : ℝ) = -3 by linarith)
  have ht := (mul_eq_zero.mp hprod).resolve_left hne
  apply Complex.ext <;> simp
  · linarith

theorem trace_is_integer (k : ℕ) (a : ℚ) (x z : ℂ)
    (ha : a ≠ -3) (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1)
    (h : s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 = 0) :
    ∃ n : ℤ, (n : ℚ) = -a - 1 := by
  have hc : ((starRingEnd ℂ) x) ^ (2 ^ k) = 1 := by
    rw [← map_pow, hx, map_one]
  have hi := (root_integral k x hx).add (root_integral k _ hc)
  rw [trace_of_s3_zero k a x z ha hx hz h] at hi
  have hq : IsIntegral ℤ (-a - 1 : ℚ) :=
    (isIntegral_algebraMap_iff (R := ℤ) (A := ℚ) (B := ℂ)
      (algebraMap ℚ ℂ).injective).mp hi
  exact UniqueFactorizationMonoid.integer_of_integral hq

theorem parameter_of_grid_zero (k : ℕ) (hk : 2 ≤ k) (a : ℚ) (x z : ℂ)
    (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1)
    (h : s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 = 0) :
    a = -3 ∨ a = -1 ∨ a = 1 := by
  by_cases ha : a = -3
  · exact Or.inl ha
  obtain ⟨n, hn⟩ := trace_is_integer k a x z ha hx hz h
  have ht : x + (starRingEnd ℂ) x = (n : ℂ) := by
    rw [trace_of_s3_zero k a x z ha hx hz h]
    exact_mod_cast hn.symm
  have hxn := norm_eq_one_of_root k x hx
  have hr := congrArg Complex.re ht
  simp at hr
  have habs : |x.re| ≤ 1 := by simpa [hxn] using Complex.abs_re_le_norm x
  have hb := abs_le.mp habs
  have hlo : (-2 : ℤ) ≤ n := by
    exact_mod_cast (show (-2 : ℝ) ≤ (n : ℝ) by linarith)
  have hhi : n ≤ (2 : ℤ) := by
    exact_mod_cast (show (n : ℝ) ≤ 2 by linarith)
  have hne : x ≠ 0 := by intro he; rw [he, zero_pow (by positivity)] at hx; norm_num at hx
  have hquad : x ^ 2 - (n : ℂ) * x + 1 = 0 := by
    rw [← Complex.inv_eq_conj hxn] at ht
    have ht' := congrArg (fun y : ℂ => y * x) ht
    simp only [add_mul, inv_mul_cancel₀ hne] at ht'
    linear_combination ht'
  interval_cases n
  · norm_num at hn
    exact Or.inr (Or.inr (by linarith))
  · have h3 : x ^ 3 = 1 := by
      norm_num at hquad
      linear_combination (x - 1) * hquad
    have he := cubic_root_of_two_power k x hx h3
    rw [he] at hquad
    norm_num at hquad
  · norm_num at hn
    exact Or.inr (Or.inl (by linarith))
  · have heven : Even (2 ^ k) := by
      exact even_two.pow_of_ne_zero (by omega)
    have hneg : (-x) ^ (2 ^ k) = 1 := by rw [heven.neg_pow, hx]
    have h3 : (-x) ^ 3 = 1 := by
      norm_num at hquad
      linear_combination -(x + 1) * hquad
    have he := cubic_root_of_two_power k (-x) hneg h3
    have hxm : x = -1 := by linear_combination -he
    rw [hxm] at hquad
    norm_num at hquad
  · norm_num at hn
    exact False.elim (ha (by linarith))

theorem fourth_root_pow_two (k : ℕ) (hk : 2 ≤ k) (x : ℂ) (h4 : x ^ 4 = 1) :
    x ^ (2 ^ k) = 1 := by
  obtain ⟨m, hm⟩ := Nat.exists_eq_add_of_le hk
  rw [hm, pow_add, show (2 : ℕ) ^ 2 = 4 by norm_num, pow_mul, h4, one_pow]

/-- Exact classification, including explicit witnesses for all exceptions. -/
theorem grid_zero_iff (k : ℕ) (hk : 2 ≤ k) (a : ℚ) :
    (∃ x z : ℂ, x ^ (2 ^ k) = 1 ∧ z ^ (2 ^ k) = 1 ∧
      s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 = 0) ↔
      a = -3 ∨ a = -1 ∨ a = 1 := by
  constructor
  · rintro ⟨x, z, hx, hz, h⟩
    exact parameter_of_grid_zero k hk a x z hx hz h
  · rintro (rfl | rfl | rfl)
    · exact ⟨1, 1, by simp, by simp, by norm_num [s3]⟩
    · refine ⟨Complex.I, -Complex.I, fourth_root_pow_two k hk _ ?_,
        fourth_root_pow_two k hk _ ?_, ?_⟩
      · norm_num [neg_pow, Complex.I_pow_four]
      · norm_num [neg_pow, Complex.I_pow_four]
      · rw [s3_square_family]
        ring_nf
        norm_num [Complex.I_sq, Complex.I_pow_four]
    · exact ⟨-1, -1, fourth_root_pow_two k hk _ (by norm_num),
        fourth_root_pow_two k hk _ (by norm_num), by norm_num [s3]⟩

theorem s3_ne_zero (k : ℕ) (hk : 2 ≤ k) (a : ℚ)
    (h3 : a ≠ -3) (hm : a ≠ -1) (hp : a ≠ 1) (x z : ℂ)
    (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1) :
    s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 ≠ 0 := by
  intro h
  rcases parameter_of_grid_zero k hk a x z hx hz h with h | h | h
  · exact h3 h
  · exact hm h
  · exact hp h

/-- The actual normalized chain ideal is the unit ideal for each allowed rational a. -/
theorem unit_complex (k : ℕ) (depth : 2 ≤ k) (a : ℚ)
    (h3 : a ≠ -3) (hm : a ≠ -1) (hp : a ≠ 1) :
    (1 : PrimeM2General.P k ℂ) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k depth
        (a : ℂ) (-1 - (a : ℂ)) 1)) := by
  classical
  have pos : 0 < 2 ^ k := by positivity
  have root := Complex.isPrimitiveRoot_exp (2 ^ k) (by positivity : 2 ^ k ≠ 0)
  refine UnitInstance.one_mem_of_no_grid_zero
    (fun _ => Polynomial.nthRootsFinset (2 ^ k) (1 : ℂ)) _ ?_ ?_
  · intro i
    have identity := congrArg (Polynomial.aeval (X i : PrimeM2General.P k ℂ))
      (Polynomial.X_pow_sub_one_eq_prod pos root)
    simp only [map_sub, map_pow, Polynomial.aeval_X, map_one, map_prod,
      Polynomial.aeval_C, MvPolynomial.algebraMap_eq] at identity
    rw [← identity]
    exact UnitInstance.var_pow_mem k depth (a : ℂ) (-1 - (a : ℂ)) 1 i
  · intro x hx
    have roots : ∀ i, x i ^ 2 ^ k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (PrimeM2General.normalized k depth
        (a : ℂ) (-1 - (a : ℂ)) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [PrimeM2General.inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth (a : ℂ) (-1 - (a : ℂ)) 1
      (PrimeM2General.active k depth 0) (by simp [PrimeM2General.active]; omega)
    have second := UnitInstance.chain_step_mem k depth (a : ℂ) (-1 - (a : ℂ)) 1
      (PrimeM2General.active k depth 2) (by simp [PrimeM2General.active]; omega)
    have first_eq :
        (⟨(PrimeM2General.active k depth 0).val+1,
          by simp [PrimeM2General.active]; omega⟩ : Fin (2*k)) =
        PrimeM2General.active k depth 1 := by
      apply Fin.ext; simp [PrimeM2General.active]
    have second_eq :
        (⟨(PrimeM2General.active k depth 2).val+1,
          by simp [PrimeM2General.active]; omega⟩ : Fin (2*k)) =
        PrimeM2General.active k depth 3 := by
      apply Fin.ext; simp [PrimeM2General.active]
    rw [first_eq] at first
    rw [second_eq] at second
    by_cases c₁ : eval x
        (X (PrimeM2General.active k depth 0)^2 -
          X (PrimeM2General.active k depth 1) : PrimeM2General.P k ℂ) = 0
    swap
    · exact ⟨_, first, c₁⟩
    by_cases c₂ : eval x
        (X (PrimeM2General.active k depth 2)^2 -
          X (PrimeM2General.active k depth 3) : PrimeM2General.P k ℂ) = 0
    swap
    · exact ⟨_, second, c₂⟩
    exfalso
    have identity := congrArg (eval x)
      (PrimeM2General.raw_normalization_identity k depth
        (a : ℂ) (-1 - (a : ℂ)) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, GrowthRationalChainUnit.eval_rawS3_rational] at identity
    exact s3_ne_zero k depth a h3 hm hp _ _
      (roots (PrimeM2General.active k depth 0))
      (roots (PrimeM2General.active k depth 2)) identity

theorem add_three_ne_zero (a : ℚ) (h3 : a ≠ -3) : (a : ℂ) + 3 ≠ 0 := by
  have hc : (a : ℂ) ≠ -3 := by exact_mod_cast h3
  intro h
  apply hc
  linear_combination h

theorem nonsingular (a : ℚ) (h3 : a ≠ -3) (hs : a ≠ -(3 / 4 : ℚ)) :
    4 * (a : ℂ) ^ 3 + 27 * (-1 - (a : ℂ)) ^ 2 ≠ 0 := by
  apply nonsingular_of_parameter_conditions _ (add_three_ne_zero a h3)
  intro h
  have hq : 4 * a + 3 = 0 := by exact_mod_cast h
  apply hs
  linarith

/-- Linear lower bound for rational members outside the exact grid exceptions. -/
theorem lfd_gt_linear (a : ℚ) (h3 : a ≠ -3) (hm : a ≠ -1) (hp : a ≠ 1)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∀ e ≤ (k - 10) / 14,
      ¬ HasLFDBound (Set.range (PrimeM2General.inputs k depth
        (a : ℂ) (-1 - (a : ℂ)) 1)) e := by
  exact GrowthMain.growth_lfd_gt_linear
    (K := ℂ) (F₀ := ℂ) (F := RatFunc ℂ) (by norm_num)
    k depth (a : ℂ) (-1 - (a : ℂ)) 1 (GrowthRatFuncModel.xi (a : ℂ)) 0
    (GrowthTwoTorsionFamily.model (a : ℂ) (add_three_ne_zero a h3))
    (Or.inl rfl) hk (unit_complex k depth a h3 hm hp)

/-- The full growth window, with no norm-five or supplied model/unit hypothesis. -/
theorem lfd_linear_window (a : ℚ)
    (h3 : a ≠ -3) (hs : a ≠ -(3 / 4 : ℚ)) (hm : a ≠ -1) (hp : a ≠ 1)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) (d : ℕ)
    (hd : IsLastFallDegree (Set.range (PrimeM2General.inputs k depth
      (a : ℂ) (-1 - (a : ℂ)) 1)) d) :
    (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  constructor
  · by_contra h
    exact lfd_gt_linear a h3 hm hp k depth hk d (by omega) hd.1
  · have hu := SelectorExtension.lfd_le_depth_add_two_of_two_ne
      (K := ℂ) (by norm_num) k depth (a : ℂ) (-1 - (a : ℂ)) 1
      (nonsingular a h3 hs)
    by_contra h
    exact hd.2 (k + 2) (by omega) hu

theorem exists_lfd_linear_window (a : ℚ)
    (h3 : a ≠ -3) (hs : a ≠ -(3 / 4 : ℚ)) (hm : a ≠ -1) (hp : a ≠ 1)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∃ d : ℕ, IsLastFallDegree (Set.range (PrimeM2General.inputs k depth
        (a : ℂ) (-1 - (a : ℂ)) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  classical
  let F := Set.range (PrimeM2General.inputs k depth (a : ℂ) (-1 - (a : ℂ)) 1)
  have hex : ∃ n, HasLFDBound F n :=
    ⟨k + 2, SelectorExtension.lfd_le_depth_add_two_of_two_ne
      (K := ℂ) (by norm_num) k depth (a : ℂ) (-1 - (a : ℂ)) 1
      (nonsingular a h3 hs)⟩
  let d := Nat.find hex
  have hd : IsLastFallDegree F d :=
    ⟨Nat.find_spec hex, fun e he => Nat.find_min hex he⟩
  exact ⟨d, hd, lfd_linear_window a h3 hs hm hp k depth hk d hd⟩

/-- Three small parameters genuinely outside the old norm-five hypothesis. -/
theorem new_members (a : ℚ) (ha : a = 0 ∨ a = 2 ∨ a = -2)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∃ d : ℕ, IsLastFallDegree (Set.range (PrimeM2General.inputs k depth
        (a : ℂ) (-1 - (a : ℂ)) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  rcases ha with rfl | rfl | rfl <;>
    exact exists_lfd_linear_window _ (by norm_num) (by norm_num)
      (by norm_num) (by norm_num) k depth hk

/-- Both nonsingular excluded parameters really have grid zeros at every depth. -/
theorem excluded_members (k : ℕ) (hk : 2 ≤ k) (a : ℚ)
    (ha : a = -1 ∨ a = 1) :
    ∃ x z : ℂ, x ^ (2 ^ k) = 1 ∧ z ^ (2 ^ k) = 1 ∧
      s3 (a : ℂ) (-1 - (a : ℂ)) x z 1 = 0 := by
  apply (grid_zero_iff k hk a).mpr
  exact Or.inr ha

#print axioms parameter_of_grid_zero
#print axioms grid_zero_iff
#print axioms unit_complex
#print axioms lfd_gt_linear
#print axioms lfd_linear_window
#print axioms exists_lfd_linear_window
#print axioms new_members
#print axioms excluded_members

end ICFallDegrees.GrowthRationalParameterClassification
