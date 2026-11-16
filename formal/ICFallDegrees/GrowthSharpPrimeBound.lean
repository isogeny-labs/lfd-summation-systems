import ICFallDegrees.GrowthExplicitPrimeBound

noncomputable section
namespace ICFallDegrees.GrowthSharpPrimeBound
open Matrix Finset GrowthExplicitPrimeBound

section Columns
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem abs_pow_le (U V : Matrix ι ι ℤ)
    (h : ∀ i j, |V i j| ≤ U i j) (n : ℕ) :
    ∀ i j, |(V ^ n) i j| ≤ (U ^ n) i j := by
  induction n with
  | zero => intro i j; simp [Matrix.one_apply]; split_ifs <;> norm_num
  | succ n ih =>
      intro i j
      rw [pow_succ, pow_succ, Matrix.mul_apply, Matrix.mul_apply]
      apply (Finset.abs_sum_le_sum_abs _ _).trans
      apply Finset.sum_le_sum
      intro l _
      rw [abs_mul]
      exact mul_le_mul (ih i l) (h l j) (abs_nonneg _) ((abs_nonneg _).trans (ih i l))

omit [DecidableEq ι] in
theorem column_sum_mul (M N : Matrix ι ι ℤ) (c : ℤ)
    (h : ∀ j, ∑ i, M i j = c) (j : ι) :
    ∑ i, (M * N) i j = c * ∑ i, N i j := by
  simp only [Matrix.mul_apply]
  rw [Finset.sum_comm]
  simp_rw [← Finset.sum_mul, h]
  exact (Finset.mul_sum _ _ _).symm

theorem column_sum_pow (M : Matrix ι ι ℤ) (c : ℤ)
    (h : ∀ j, ∑ i, M i j = c) (n : ℕ) (j : ι) :
    ∑ i, (M ^ n) i j = c ^ n := by
  induction n generalizing j with
  | zero => simp [Matrix.one_apply]
  | succ n ih =>
      rw [pow_succ, column_sum_mul _ _ _ (fun j => ih j), h, pow_succ]

end Columns

theorem eliminant_nonnegative (k : ℕ) :
    ∀ i j, 0 ≤ eliminantMatrix k ℤ i j := by
  have h : ∀ i j : Fin (2^k),
      |(shift (2^k) ℤ - 1) i j| ≤ (shift (2^k) ℤ + 8) i j := by
    intro i j
    simp only [Matrix.sub_apply, Matrix.add_apply, shift, Matrix.one_apply]
    change |(if finRotate (2^k) i = j then (1:ℤ) else 0) -
      (if i = j then 1 else 0)| ≤
      (if finRotate (2^k) i = j then (1:ℤ) else 0) +
      (if i = j then 8 else 0)
    split_ifs <;> norm_num
  intro i j
  exact sub_nonneg.mpr ((le_abs_self _).trans (abs_pow_le _ _ h (2^k) i j))

theorem shift_column_sum (n : ℕ) (j : Fin n) : ∑ i, shift n ℤ i j = 1 := by
  simp only [shift, Equiv.apply_eq_iff_eq_symm_apply]
  simp

theorem eliminant_column_sum (k : ℕ) (j : Fin (2^k)) :
    ∑ i, eliminantMatrix k ℤ i j = 9 ^ (2^k) := by
  have hU : ∀ j : Fin (2^k), ∑ i, (shift (2^k) ℤ + 8) i j = 9 := by
    intro j
    rw [show (∑ i, (shift (2^k) ℤ + 8) i j) =
      (∑ i, shift (2^k) ℤ i j) + ∑ i, (8 : Matrix (Fin (2^k)) (Fin (2^k)) ℤ) i j
      by simp [Finset.sum_add_distrib]]
    rw [shift_column_sum]
    change 1 + (∑ i : Fin (2^k), if i = j then (8:ℤ) else 0) = 9
    simp
  have hV : ∀ j : Fin (2^k), ∑ i, (shift (2^k) ℤ - 1) i j = 0 := by
    intro j
    simp only [Matrix.sub_apply, Finset.sum_sub_distrib, shift_column_sum]
    simp [Matrix.one_apply]
  simp only [eliminantMatrix, Matrix.sub_apply, Finset.sum_sub_distrib,
    column_sum_pow _ 9 hU, column_sum_pow _ 0 hV]
  simp

def bound (k : ℕ) : ℕ := 9 ^ (4 ^ k)

theorem eliminant_det_bound (k : ℕ) :
    |(eliminantMatrix k ℤ).det| ≤ (bound k : ℤ) := by
  have hc : ColBound (eliminantMatrix k ℤ) (9 ^ (2^k)) := by
    intro j
    simp_rw [abs_of_nonneg (eliminant_nonnegative k _ _)]
    exact (eliminant_column_sum k j).le
  have h := abs_det_le_of_colBound _ _ hc
  have he : (2 ^ k) * (2 ^ k) = 4 ^ k := by rw [← mul_pow]; norm_num
  simpa only [Fintype.card_fin, ← pow_mul, he, bound, Nat.cast_pow, Nat.cast_ofNat] using h


theorem eliminant_det_mod_prime (k p : ℕ) [Fact p.Prime] (hp : bound k < p) :
    (eliminantMatrix k (ZMod p)).det ≠ 0 := by
  have hnz := GrowthExplicitPrimeBound.eliminant_det_ne_zero k
  have hb : (eliminantMatrix k ℤ).det.natAbs ≤ bound k := by
    have h := eliminant_det_bound k
    rw [Int.abs_eq_natAbs] at h
    exact_mod_cast h
  intro hz
  have hcast : ((eliminantMatrix k ℤ).det : ZMod p) = 0 := by
    change (Int.castRingHom (ZMod p)) (eliminantMatrix k ℤ).det = 0
    rw [(Int.castRingHom (ZMod p)).map_det, map_eliminant]
    exact hz
  have hd : (p : ℤ) ∣ (eliminantMatrix k ℤ).det :=
    (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp hcast
  have hd' : p ∣ (eliminantMatrix k ℤ).det.natAbs :=
    Int.natCast_dvd_natCast.mp (Int.dvd_natAbs.mpr hd)
  have hp' := Nat.le_of_dvd (Int.natAbs_pos.mpr hnz) hd'
  omega
theorem no_grid_zero (k p : ℕ) [Fact p.Prime] (hdet : (eliminantMatrix k (ZMod p)).det ≠ 0)
    (x z : ZMod p) (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1) :
    GrowthSecant.s3 (6 : ZMod p) (-1 - 6) x z 1 ≠ 0 := by
  intro h
  have hlin : (x - 1) * z = x + 8 := by
    have hh := GrowthRationalUnit.parameter_of_s3_zero (6 : ZMod p) x z h
    linear_combination -hh
  have he : (x + 8) ^ (2 ^ k) = (x - 1) ^ (2 ^ k) := by
    rw [← hlin, mul_pow, hz, mul_one]
  have hv := eliminant_mulVec k x hx
  rw [he, sub_self, zero_smul] at hv
  have hzero := Matrix.eq_zero_of_mulVec_eq_zero hdet hv
  have hc := congrFun hzero (⟨0, by positivity⟩ : Fin (2 ^ k))
  simp at hc

open MvPolynomial PrimeM2General in
theorem unit_prime (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (hdet : (eliminantMatrix k (ZMod p)).det ≠ 0) (hadmissible : 2 ^ k ∣ p - 1) :
    (1 : P k (ZMod p)) ∈
      Ideal.span (Set.range (inputs k depth (6 : ZMod p) (-7) 1)) := by
  classical
  have hb : (-1 - (6 : ZMod p)) = -7 := by ring
  rw [← hb]
  obtain ⟨ζ, root⟩ := UnitInstance.exists_primitiveRoot_zmod p (2 ^ k) hadmissible
  have pos : 0 < 2 ^ k := by positivity
  refine UnitInstance.one_mem_of_no_grid_zero
    (fun _ => Polynomial.nthRootsFinset (2 ^ k) (1 : ZMod p)) _ ?_ ?_
  · intro i
    have identity := congrArg (Polynomial.aeval (X i : P k (ZMod p)))
      (Polynomial.X_pow_sub_one_eq_prod pos root)
    simp only [map_sub, map_pow, Polynomial.aeval_X, map_one, map_prod,
      Polynomial.aeval_C, MvPolynomial.algebraMap_eq] at identity
    rw [← identity]
    exact UnitInstance.var_pow_mem k depth (6 : ZMod p) (-1 - 6) 1 i
  · intro x hx
    have roots : ∀ i, x i ^ 2 ^ k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth (6 : ZMod p) (-1 - 6) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth (6 : ZMod p) (-1 - 6) 1
      (active k depth 0) (by simp [active]; omega)
    have second := UnitInstance.chain_step_mem k depth (6 : ZMod p) (-1 - 6) 1
      (active k depth 2) (by simp [active]; omega)
    have first_eq :
        (⟨(active k depth 0).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
        active k depth 1 := by apply Fin.ext; simp [active]
    have second_eq :
        (⟨(active k depth 2).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
        active k depth 3 := by apply Fin.ext; simp [active]
    rw [first_eq] at first
    rw [second_eq] at second
    by_cases c₁ : eval x
        (X (active k depth 0)^2 - X (active k depth 1) : P k (ZMod p)) = 0
    swap
    · exact ⟨_, first, c₁⟩
    by_cases c₂ : eval x
        (X (active k depth 2)^2 - X (active k depth 3) : P k (ZMod p)) = 0
    swap
    · exact ⟨_, second, c₂⟩
    exfalso
    have identity := congrArg (eval x)
      (raw_normalization_identity k depth (6 : ZMod p) (-1 - 6) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, GrowthRationalChainUnit.eval_rawS3_rational] at identity
    exact no_grid_zero k p hdet _ _ (roots (active k depth 0)) (roots (active k depth 2)) identity


/-- Divisibility gives a useful exact sufficient criterion below the size cutoff. -/
theorem unit_of_not_dvd (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (hadmissible : 2 ^ k ∣ p - 1)
    (hgood : ¬ (p : ℤ) ∣ (eliminantMatrix k ℤ).det) :
    (1 : PrimeM2General.P k (ZMod p)) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k depth (6 : ZMod p) (-7) 1)) := by
  apply unit_prime k p depth _ hadmissible
  intro he
  apply hgood
  apply (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp
  change (Int.castRingHom (ZMod p)) (eliminantMatrix k ℤ).det = 0
  rw [(Int.castRingHom (ZMod p)).map_det, map_eliminant]
  exact he

theorem proper_prime_le_bound (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (hadmissible : 2 ^ k ∣ p - 1)
    (hproper : (1 : PrimeM2General.P k (ZMod p)) ∉
      Ideal.span (Set.range (PrimeM2General.inputs k depth (6 : ZMod p) (-7) 1))) :
    p ≤ bound k := by
  by_contra h
  exact hproper (unit_prime k p depth (eliminant_det_mod_prime k p (by omega)) hadmissible)

theorem prime_lfd_window_of_det (k p : ℕ) [Fact p.Prime] (hk : 24 ≤ k)
    (hpDisc : 2187 < p) (hpLarge : 2 * 2 ^ k + 1 ≤ p)
    (hadmissible : 2 ^ k ∣ p - 1)
    (hdet : (eliminantMatrix k (ZMod p)).det ≠ 0) :
    ∃ d : ℕ, IsLastFallDegree
      (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
        (6 : ZMod p) (-7) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  have depth : 2 ≤ k := by omega
  have hu := unit_prime k p depth hdet hadmissible
  have ho := GrowthPrimeE6NativeNonmembership.e6_not_one_of_large_prime
    p k depth ((k - 10) / 14) (by omega) (by omega) hpLarge
  obtain ⟨d, hd, hlo, hhi, _, _⟩ :=
    GrowthPrimeNativeWindow.fixed_e6_prime_theta_of_native_certificates
      p hpDisc k depth hk hu ho
  exact ⟨d, hd, hlo, hhi⟩

theorem characteristic_threshold (k : ℕ) (hk : 2 ≤ k) :
    2 * 2 ^ k ≤ bound k ∧ 2187 ≤ bound k := by
  constructor
  · exact (GrowthExplicitPrimeBound.two_mul_le_nine_pow (2 ^ k)).trans
      (by unfold bound; gcongr <;> norm_num)
  · have hh : (4 : ℕ) ^ 1 ≤ 4 ^ k := Nat.pow_le_pow_right (by norm_num) (by omega)
    have hh' : (9 : ℕ) ^ 4 ≤ 9 ^ (4 ^ k) := Nat.pow_le_pow_right (by norm_num) (by simpa using hh)
    exact (by norm_num : 2187 ≤ 9 ^ 4).trans hh'

theorem prime_lfd_window (k p : ℕ) [Fact p.Prime] (hk : 24 ≤ k)
    (hp : bound k < p) (hadmissible : 2 ^ k ∣ p - 1) :
    ∃ d : ℕ, IsLastFallDegree
      (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
        (6 : ZMod p) (-7) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  have hth := characteristic_threshold k (by omega)
  exact prime_lfd_window_of_det k p hk (by omega) (by omega)
    hadmissible (eliminant_det_mod_prime k p hp)

theorem strict_improvement (k : ℕ) : bound k < GrowthExplicitPrimeBound.bound k := by
  have hp : 0 < bound k := by unfold bound; positivity
  have ht : 1 < (2 : ℕ) ^ (2 ^ k) := one_lt_pow₀ (by norm_num) (by positivity)
  change bound k < 2 ^ (2 ^ k) * bound k
  nlinarith

#print axioms eliminant_nonnegative
#print axioms eliminant_column_sum
#print axioms eliminant_det_bound
#print axioms unit_of_not_dvd
#print axioms proper_prime_le_bound
#print axioms prime_lfd_window_of_det
#print axioms prime_lfd_window
#print axioms strict_improvement
end ICFallDegrees.GrowthSharpPrimeBound

