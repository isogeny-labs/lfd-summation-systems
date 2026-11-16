import ICFallDegrees.GrowthPrimeIntegralUnit
import ICFallDegrees.GrowthFixedE6PrimeTheta
import Mathlib.LinearAlgebra.Matrix.Nondegenerate
import Mathlib.Logic.Equiv.Fin.Rotate

/-! # Explicit exceptional-prime bound for the fixed E6 chain system -/

noncomputable section
namespace ICFallDegrees.GrowthExplicitPrimeBound

open Matrix Finset

section Bounds
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

def ColBound (M : Matrix ι ι ℤ) (c : ℤ) : Prop :=
  ∀ j, ∑ i, |M i j| ≤ c

theorem colBound_one : ColBound (1 : Matrix ι ι ℤ) 1 := by
  intro j
  simp only [Matrix.one_apply, apply_ite abs, abs_one, abs_zero]
  simp

theorem colBound_natCast (n : ℕ) : ColBound (n : Matrix ι ι ℤ) n := by
  intro j
  change (∑ i, |if i = j then (n : ℤ) else 0|) ≤ n
  simp only [apply_ite abs, abs_zero, abs_of_nonneg (Int.natCast_nonneg n)]
  simp

omit [DecidableEq ι] in
theorem colBound_add {M N : Matrix ι ι ℤ} {a b : ℤ}
    (hM : ColBound M a) (hN : ColBound N b) : ColBound (M + N) (a + b) := by
  intro j
  calc
    ∑ i, |(M + N) i j| ≤ ∑ i, (|M i j| + |N i j|) :=
      Finset.sum_le_sum (fun i _ => abs_add_le _ _)
    _ = (∑ i, |M i j|) + ∑ i, |N i j| := Finset.sum_add_distrib
    _ ≤ a + b := add_le_add (hM j) (hN j)

omit [DecidableEq ι] in
theorem colBound_sub {M N : Matrix ι ι ℤ} {a b : ℤ}
    (hM : ColBound M a) (hN : ColBound N b) : ColBound (M - N) (a + b) := by
  intro j
  calc
    ∑ i, |(M - N) i j| ≤ ∑ i, (|M i j| + |N i j|) :=
      Finset.sum_le_sum (fun i _ => abs_sub _ _)
    _ = (∑ i, |M i j|) + ∑ i, |N i j| := Finset.sum_add_distrib
    _ ≤ a + b := add_le_add (hM j) (hN j)

omit [DecidableEq ι] in
theorem colBound_mul {M N : Matrix ι ι ℤ} {a b : ℤ}
    (ha : 0 ≤ a) (hM : ColBound M a) (hN : ColBound N b) :
    ColBound (M * N) (a * b) := by
  intro j
  calc
    ∑ i, |(M * N) i j| ≤ ∑ i, ∑ l, |M i l| * |N l j| := by
      apply Finset.sum_le_sum
      intro i _
      simpa only [Matrix.mul_apply, abs_mul] using
        Finset.abs_sum_le_sum_abs (fun l => M i l * N l j) Finset.univ
    _ = ∑ l, (∑ i, |M i l|) * |N l j| := by
      rw [Finset.sum_comm]
      simp_rw [Finset.sum_mul]
    _ ≤ ∑ l, a * |N l j| := Finset.sum_le_sum
      (fun l _ => mul_le_mul_of_nonneg_right (hM l) (abs_nonneg _))
    _ = a * ∑ l, |N l j| := (Finset.mul_sum _ _ _).symm
    _ ≤ a * b := mul_le_mul_of_nonneg_left (hN j) ha

theorem colBound_pow {M : Matrix ι ι ℤ} {c : ℤ}
    (hc : 0 ≤ c) (hM : ColBound M c) (n : ℕ) : ColBound (M ^ n) (c ^ n) := by
  induction n with
  | zero => simpa using (colBound_one (ι := ι))
  | succ n ih =>
      simpa only [pow_succ] using colBound_mul (pow_nonneg hc n) ih hM

/-- The determinant is bounded by the product of the column absolute sums. -/
theorem abs_det_le_of_colBound (M : Matrix ι ι ℤ) (c : ℤ) (hM : ColBound M c) :
    |M.det| ≤ c ^ Fintype.card ι := by
  let f : (ι → ι) → ℤ := fun p => ∏ i, |M (p i) i|
  have hf : ∀ p, 0 ≤ f p := fun p => Finset.prod_nonneg (fun _ _ => abs_nonneg _)
  have hi : Function.Injective (fun p : Equiv.Perm ι => (p : ι → ι)) :=
    fun _ _ h => Equiv.ext (congrFun h)
  calc
    |M.det| ≤ ∑ p : Equiv.Perm ι, f p := by
      rw [Matrix.det_apply']
      apply (Finset.abs_sum_le_sum_abs _ _).trans
      apply Finset.sum_le_sum
      intro p _
      simp [f, abs_mul, Finset.abs_prod, Equiv.Perm.sign_abs]
    _ = ∑ p ∈ Finset.univ.image (fun p : Equiv.Perm ι => (p : ι → ι)), f p := by
      rw [Finset.sum_image]
      exact fun _ _ _ _ h => hi h
    _ ≤ ∑ p : ι → ι, f p := by
      apply Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      intro p _ _
      exact hf p
    _ = ∏ i, ∑ j, |M j i| := (Fintype.prod_sum (fun i j => |M j i|)).symm
    _ ≤ ∏ _i : ι, c := Finset.prod_le_prod
      (fun i _ => Finset.sum_nonneg (fun _ _ => abs_nonneg _)) (fun i _ => hM i)
    _ = c ^ Fintype.card ι := by simp

end Bounds

def shift (n : ℕ) (R : Type*) [CommRing R] : Matrix (Fin n) (Fin n) R :=
  fun i j => if finRotate n i = j then 1 else 0

def eliminantMatrix (k : ℕ) (R : Type*) [CommRing R] :
    Matrix (Fin (2 ^ k)) (Fin (2 ^ k)) R :=
  (shift (2 ^ k) R + 8) ^ (2 ^ k) - (shift (2 ^ k) R - 1) ^ (2 ^ k)

def bound (k : ℕ) : ℕ := 2 ^ (2 ^ k) * 9 ^ (4 ^ k)

theorem shift_colBound (n : ℕ) : ColBound (shift n ℤ) 1 := by
  intro j
  simp only [shift, apply_ite abs, abs_one, abs_zero, Equiv.apply_eq_iff_eq_symm_apply]
  simp

theorem eliminant_colBound (k : ℕ) :
    ColBound (eliminantMatrix k ℤ) (9 ^ (2 ^ k) + 2 ^ (2 ^ k)) := by
  apply colBound_sub
  · apply colBound_pow (by norm_num)
    simpa using colBound_add (shift_colBound (2 ^ k)) (colBound_natCast (ι := Fin (2 ^ k)) 8)
  · apply colBound_pow (by norm_num)
    simpa using colBound_sub (shift_colBound (2 ^ k)) (colBound_one (ι := Fin (2 ^ k)))

theorem eliminant_det_bound (k : ℕ) :
    |(eliminantMatrix k ℤ).det| ≤ (bound k : ℤ) := by
  have h := abs_det_le_of_colBound _ _ (eliminant_colBound k)
  simp only [Fintype.card_fin] at h
  calc
    |(eliminantMatrix k ℤ).det| ≤ (9 ^ (2 ^ k) + 2 ^ (2 ^ k)) ^ (2 ^ k) := h
    _ ≤ (2 * 9 ^ (2 ^ k)) ^ (2 ^ k) := by
      gcongr
      have hh : (2 : ℤ) ^ (2 ^ k) ≤ 9 ^ (2 ^ k) := pow_le_pow_left₀ (by norm_num) (by norm_num) _
      linarith
    _ = (bound k : ℤ) := by
      have he : (2 ^ k) * (2 ^ k) = 4 ^ k := by rw [← mul_pow]; norm_num
      simp only [bound, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat, mul_pow, ← pow_mul, he]

theorem map_shift {R S : Type*} [CommRing R] [CommRing S] (f : R →+* S) (n : ℕ) :
    (f.mapMatrix (shift n R)) = shift n S := by
  ext i j
  simp [shift, RingHom.mapMatrix_apply, apply_ite f]

theorem map_eliminant {R S : Type*} [CommRing R] [CommRing S] (f : R →+* S) (k : ℕ) :
    f.mapMatrix (eliminantMatrix k R) = eliminantMatrix k S := by
  simp only [eliminantMatrix, map_sub, map_add, map_pow, map_ofNat, map_one, map_shift]

theorem eliminant_mod_two (k : ℕ) : eliminantMatrix k (ZMod 2) = 1 := by
  have h := congrArg (fun f : Polynomial (ZMod 2) =>
    Polynomial.aeval (shift (2 ^ k) (ZMod 2)) f)
    (GrowthPrimeIntegralUnit.eliminant_mod_two_eq_one k)
  simpa only [GrowthPrimeIntegralUnit.eliminantModTwo, eliminantMatrix,
    map_sub, map_pow, map_add, Polynomial.aeval_X, Polynomial.aeval_C,
    map_ofNat, map_one] using h

theorem eliminant_det_ne_zero (k : ℕ) : (eliminantMatrix k ℤ).det ≠ 0 := by
  apply GrowthCharZeroDeterminant.det_ne_zero_of_mod_two_one
  change (Int.castRingHom (ZMod 2)).mapMatrix (eliminantMatrix k ℤ) = 1
  rw [map_eliminant, eliminant_mod_two]

theorem eliminant_det_mod_prime (k p : ℕ) [Fact p.Prime] (hp : bound k < p) :
    (eliminantMatrix k (ZMod p)).det ≠ 0 := by
  have hnz := eliminant_det_ne_zero k
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

section Eigenvector
variable {K : Type*} [Field K]

theorem shift_mulVec (n : ℕ) (v : Fin n → K) :
    (shift n K).mulVec v = fun i => v (finRotate n i) := by
  ext i
  simp [Matrix.mulVec, dotProduct, shift, ite_mul]

theorem root_rotate (n : ℕ) (hn : 0 < n) (x : K) (hx : x ^ n = 1) (i : Fin n) :
    x ^ (finRotate n i).val = x * x ^ i.val := by
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hn.ne'
  by_cases hi : i = Fin.last m
  · subst i
    simpa only [finRotate_last, Fin.val_zero, pow_zero, Fin.val_last, pow_succ, mul_comm] using hx.symm
  · rw [coe_finRotate_of_ne_last hi, pow_succ, mul_comm]

theorem power_mulVec {n : ℕ} (M : Matrix (Fin n) (Fin n) K)
    (v : Fin n → K) (t : K) (h : M.mulVec v = t • v) (r : ℕ) :
    (M ^ r).mulVec v = t ^ r • v := by
  induction r with
  | zero => simp
  | succ r ih =>
      rw [pow_succ, ← Matrix.mulVec_mulVec, h, Matrix.mulVec_smul, ih]
      simp [smul_smul, pow_succ, mul_comm]

theorem eliminant_mulVec (k : ℕ) (x : K) (hx : x ^ (2 ^ k) = 1) :
    (eliminantMatrix k K).mulVec (fun i => x ^ i.val) =
      ((x + 8) ^ (2 ^ k) - (x - 1) ^ (2 ^ k)) • (fun i => x ^ i.val) := by
  let v : Fin (2 ^ k) → K := fun i => x ^ i.val
  have he : (shift (2 ^ k) K).mulVec v = x • v := by
    rw [shift_mulVec]
    ext i
    exact root_rotate _ (by positivity) x hx i
  have h8 : (shift (2 ^ k) K + 8).mulVec v = (x + 8) • v := by
    rw [Matrix.add_mulVec, he]
    simp only [Matrix.ofNat_mulVec, ← add_smul]
  have h1 : (shift (2 ^ k) K - 1).mulVec v = (x - 1) • v := by
    rw [Matrix.sub_mulVec, he]
    simp [sub_smul]
  change ((shift (2 ^ k) K + 8) ^ (2 ^ k) -
    (shift (2 ^ k) K - 1) ^ (2 ^ k)).mulVec v = _
  rw [Matrix.sub_mulVec, power_mulVec _ _ _ h8, power_mulVec _ _ _ h1, sub_smul]

end Eigenvector

theorem no_grid_zero (k p : ℕ) [Fact p.Prime] (hp : bound k < p)
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
  have hzero := Matrix.eq_zero_of_mulVec_eq_zero (eliminant_det_mod_prime k p hp) hv
  have hc := congrFun hzero (⟨0, by positivity⟩ : Fin (2 ^ k))
  simp at hc

open MvPolynomial PrimeM2General in
theorem unit_prime (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (hp : bound k < p) (hadmissible : 2 ^ k ∣ p - 1) :
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
    exact no_grid_zero k p hp _ _ (roots (active k depth 0)) (roots (active k depth 2)) identity

theorem proper_prime_le_bound (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (hadmissible : 2 ^ k ∣ p - 1)
    (hproper : (1 : PrimeM2General.P k (ZMod p)) ∉
      Ideal.span (Set.range (PrimeM2General.inputs k depth (6 : ZMod p) (-7) 1))) :
    p ≤ bound k := by
  by_contra h
  exact hproper (unit_prime k p depth (by omega) hadmissible)

theorem two_mul_le_nine_pow (n : ℕ) : 2 * n ≤ 9 ^ n := by
  induction n with
  | zero => norm_num
  | succ n ih =>
      rw [pow_succ]
      have hp : 0 < 9 ^ n := by positivity
      omega

theorem nine_pow_le_bound (k : ℕ) : 9 ^ (4 ^ k) ≤ bound k := by
  unfold bound
  have h : (1 : ℕ) ≤ 2 ^ (2 ^ k) := by
    have ht : 0 < (2 : ℕ) ^ (2 ^ k) := by positivity
    omega
  simpa only [one_mul] using Nat.mul_le_mul_right (9 ^ (4 ^ k)) h

theorem characteristic_threshold (k : ℕ) (hk : 2 ≤ k) :
    2 * 2 ^ k ≤ bound k ∧ 2187 ≤ bound k := by
  constructor
  · have h1 := two_mul_le_nine_pow (2 ^ k)
    have h2 : 9 ^ (2 ^ k) ≤ 9 ^ (4 ^ k) := by gcongr <;> norm_num
    exact h1.trans (h2.trans (nine_pow_le_bound k))
  · have h4 : 4 ≤ 4 ^ k := by
      have hh : (4 : ℕ) ^ 1 ≤ 4 ^ k := Nat.pow_le_pow_right (by norm_num) (by omega)
      simpa using hh
    have hh : (9 : ℕ) ^ 4 ≤ 9 ^ (4 ^ k) := Nat.pow_le_pow_right (by norm_num) h4
    exact (by norm_num : 2187 ≤ 9 ^ 4).trans (hh.trans (nine_pow_le_bound k))

/-- Every admissible prime above the displayed computable threshold has the window. -/
theorem prime_lfd_window (k p : ℕ) [Fact p.Prime] (hk : 24 ≤ k)
    (hp : bound k < p) (hadmissible : 2 ^ k ∣ p - 1) :
    ∃ d : ℕ, IsLastFallDegree
      (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
        (6 : ZMod p) (-7) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  have depth : 2 ≤ k := by omega
  have hth := characteristic_threshold k depth
  have hpDisc : 2187 < p := lt_of_le_of_lt hth.2 hp
  have hpLarge : 2 * 2 ^ k + 1 ≤ p := by omega
  have hu := unit_prime k p depth hp hadmissible
  have ho := GrowthPrimeE6NativeNonmembership.e6_not_one_of_large_prime
    p k depth ((k - 10) / 14) (by omega) (by omega) hpLarge
  obtain ⟨d, hd, hlo, hhi, _, _⟩ :=
    GrowthPrimeNativeWindow.fixed_e6_prime_theta_of_native_certificates
      p hpDisc k depth hk hu ho
  exact ⟨d, hd, hlo, hhi⟩

/-- The formerly existential unit threshold is replaced by max(B, bound k). -/
theorem exists_prime_window_above (k B : ℕ) (hk : 24 ≤ k) :
    ∃ p : ℕ, Nat.Prime p ∧ max B (bound k) < p ∧
      p ≡ 1 [MOD 2 ^ (k + 1)] ∧ ∃ d : ℕ,
        (∀ [Fact p.Prime], IsLastFallDegree
          (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
            (6 : ZMod p) (-7) 1)) d) ∧
        (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  obtain ⟨p, hprime, hgt, hmod⟩ :=
    GrowthPrimeAdmissible.exists_admissible_prime_above k (max B (bound k))
  haveI : Fact p.Prime := ⟨hprime⟩
  have hd : 2 ^ (k + 1) ∣ p - 1 :=
    (Nat.modEq_iff_dvd' hprime.one_le).mp hmod.symm
  have hd' : 2 ^ k ∣ 2 ^ (k + 1) := by rw [pow_succ]; exact dvd_mul_right _ _
  obtain ⟨d, hlfd, hlo, hhi⟩ := prime_lfd_window k p hk
    (lt_of_le_of_lt (le_max_right _ _) hgt) (hd'.trans hd)
  refine ⟨p, hprime, hgt, hmod, d, ?_, hlo, hhi⟩
  intro _
  simpa using hlfd

#print axioms eliminant_det_bound
#print axioms eliminant_mod_two
#print axioms eliminant_det_ne_zero
#print axioms eliminant_det_mod_prime
#print axioms no_grid_zero
#print axioms unit_prime
#print axioms proper_prime_le_bound
#print axioms prime_lfd_window
#print axioms exists_prime_window_above

end ICFallDegrees.GrowthExplicitPrimeBound
