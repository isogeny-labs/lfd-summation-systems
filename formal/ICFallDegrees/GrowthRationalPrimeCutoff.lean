import ICFallDegrees.GrowthSharpPrimeBound
import ICFallDegrees.GrowthRationalPrimeNative
import ICFallDegrees.GrowthRationalParameterClassification
import Mathlib.LinearAlgebra.Vandermonde

/-!
# Rational-parameter prime-field growth with a closed cutoff

For a=m/q with q>0, use the integer cyclic matrix for
(q*X+m+2*q)^N-(q*(X-1))^N, N=2^k. Complex root-grid classification
and a genuine Vandermonde diagonalization prove its determinant nonzero.
Column bounds give |det| <= 2^N*(3*q+|m+2*q|)^(N*N).

The public cutoff also dominates 2*N, q, |m+3*q| and |4*m+3*q|.
These terms prove the required characteristic, denominator and discriminant
conditions after reduction. Rational reduction is explicitly num/den in Fp;
there is no assumed field homomorphism from Q to Fp.

The headline theorem is for the actual normalized PrimeM2General.inputs,
not a supplied unit certificate or CurveModel. The rational exclusions
are -3,-3/4,-1,1; the prime hypothesis is universal above the cutoff,
subject to the usual factor-base admissibility 2^k dividing p-1.
-/
noncomputable section
namespace ICFallDegrees.GrowthRationalPrimeCutoff
open Matrix Finset GrowthExplicitPrimeBound

def matrix (k : ℕ) (m : ℤ) (q : ℕ) (R : Type*) [CommRing R] :
    Matrix (Fin (2^k)) (Fin (2^k)) R :=
  ((q : Matrix (Fin (2^k)) (Fin (2^k)) R) * shift (2^k) R +
    ((m + 2 * (q:ℤ) : ℤ) : Matrix (Fin (2^k)) (Fin (2^k)) R)) ^ (2^k) -
  ((q : Matrix (Fin (2^k)) (Fin (2^k)) R) * (shift (2^k) R - 1)) ^ (2^k)

def value {R : Type*} [CommRing R] (k : ℕ) (m : ℤ) (q : ℕ) (x : R) : R :=
  ((q:R) * x + ((m + 2*(q:ℤ) : ℤ):R)) ^ (2^k) -
    ((q:R) * (x-1)) ^ (2^k)

def height (m : ℤ) (q : ℕ) : ℕ := 3*q + (m+2*(q:ℤ)).natAbs
def detBound (k : ℕ) (m : ℤ) (q : ℕ) : ℕ :=
  2 ^ (2^k) * height m q ^ (4^k)
def cutoff (k : ℕ) (m : ℤ) (q : ℕ) : ℕ :=
  detBound k m q + 2 * 2^k + q +
    (m+3*(q:ℤ)).natAbs + (4*m+3*(q:ℤ)).natAbs

theorem map_matrix {R S : Type*} [CommRing R] [CommRing S]
    (f : R →+* S) (k : ℕ) (m : ℤ) (q : ℕ) :
    f.mapMatrix (matrix k m q R) = matrix k m q S := by
  simp only [matrix, map_sub, map_pow, map_add, map_mul, map_natCast,
    map_intCast, map_one, map_shift]

theorem colBound_intCast {ι : Type*} [Fintype ι] [DecidableEq ι] (z : ℤ) :
    ColBound (z : Matrix ι ι ℤ) |z| := by
  intro j
  change (∑ i, |if i = j then z else 0|) ≤ |z|
  simp only [apply_ite abs, abs_zero]
  simp

theorem matrix_det_bound (k : ℕ) (m : ℤ) (q : ℕ) :
    |(matrix k m q ℤ).det| ≤ (detBound k m q : ℤ) := by
  have hH : 0 ≤ (height m q : ℤ) := by positivity
  have ha : (q:ℤ) + |m+2*(q:ℤ)| ≤ (height m q : ℤ) := by
    simp only [height, Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat, ← Int.abs_eq_natAbs]
    omega
  have hb : (q:ℤ) * 2 ≤ (height m q : ℤ) := by
    simp only [height, Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat, ← Int.abs_eq_natAbs]
    have := abs_nonneg (m+2*(q:ℤ))
    omega
  have hu := colBound_add
    (colBound_mul (by positivity : 0 ≤ (q:ℤ))
      (colBound_natCast (ι := Fin (2^k)) q) (shift_colBound (2^k)))
    (colBound_intCast (ι := Fin (2^k)) (m+2*(q:ℤ)))
  have hv := colBound_mul (by positivity : 0 ≤ (q:ℤ))
    (colBound_natCast (ι := Fin (2^k)) q)
    (colBound_sub (shift_colBound (2^k)) (colBound_one (ι := Fin (2^k))))
  have hu' : ColBound
      ((q : Matrix (Fin (2^k)) (Fin (2^k)) ℤ) * shift (2^k) ℤ +
        ((m+2*(q:ℤ):ℤ) : Matrix (Fin (2^k)) (Fin (2^k)) ℤ)) (height m q) :=
    fun j => (hu j).trans (by simpa using ha)
  have hv' : ColBound
      ((q : Matrix (Fin (2^k)) (Fin (2^k)) ℤ) * (shift (2^k) ℤ - 1)) (height m q) :=
    fun j => (hv j).trans (by simpa using hb)
  have hc := colBound_sub (colBound_pow hH hu' (2^k)) (colBound_pow hH hv' (2^k))
  have hd := abs_det_le_of_colBound _ _ hc
  have he : (2^k)*(2^k) = 4^k := by rw [← mul_pow]; norm_num
  simpa only [matrix, Fintype.card_fin, ← two_mul, mul_pow, ← pow_mul,
    he, detBound, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat] using hd

theorem matrix_mulVec {K : Type*} [Field K] (k : ℕ) (m : ℤ) (q : ℕ)
    (x : K) (hx : x ^ (2^k) = 1) :
    (matrix k m q K).mulVec (fun i => x ^ i.val) =
      value k m q x • (fun i => x ^ i.val) := by
  let v : Fin (2^k) → K := fun i => x ^ i.val
  have he : (shift (2^k) K).mulVec v = x • v := by
    rw [shift_mulVec]
    ext i
    exact root_rotate _ (by positivity) x hx i
  have hU :
      ((q : Matrix (Fin (2^k)) (Fin (2^k)) K) * shift (2^k) K +
        ((m+2*(q:ℤ):ℤ) : Matrix (Fin (2^k)) (Fin (2^k)) K)).mulVec v =
      ((q:K)*x + ((m+2*(q:ℤ):ℤ):K)) • v := by
    rw [Matrix.add_mulVec, ← Matrix.mulVec_mulVec, he]
    simp only [Matrix.natCast_mulVec, Matrix.intCast_mulVec, smul_smul, ← add_smul]
  have hV :
      ((q : Matrix (Fin (2^k)) (Fin (2^k)) K) * (shift (2^k) K - 1)).mulVec v =
      ((q:K)*(x-1)) • v := by
    rw [← Matrix.mulVec_mulVec, Matrix.sub_mulVec, he]
    rw [Matrix.one_mulVec, Matrix.natCast_mulVec]
    ext i
    simp only [Pi.smul_apply, smul_eq_mul, Pi.sub_apply]
    ring
  change (matrix k m q K).mulVec v = _
  unfold matrix
  rw [Matrix.sub_mulVec, power_mulVec _ _ _ hU, power_mulVec _ _ _ hV, ← sub_smul]
  rfl

theorem value_ne_zero_complex (k : ℕ) (hk : 2 ≤ k) (a : ℚ)
    (m : ℤ) (q : ℕ) (hq : 0 < q) (hrep : (m:ℚ) = (q:ℚ)*a)
    (h3 : a ≠ -3) (hm : a ≠ -1) (hp : a ≠ 1)
    (x : ℂ) (hx : x ^ (2^k) = 1) : value k m q x ≠ 0 := by
  have hqC : (q:ℂ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hrepC : (m:ℂ) = (q:ℂ)*(a:ℂ) := by exact_mod_cast hrep
  have htop : (q:ℂ)*x + ((m+2*(q:ℤ):ℤ):ℂ) =
      (q:ℂ)*(x+2+(a:ℂ)) := by push_cast; rw [hrepC]; ring
  intro hzero
  have he : (x+2+(a:ℂ)) ^ (2^k) = (x-1) ^ (2^k) := by
    unfold value at hzero
    rw [htop, mul_pow, mul_pow, ← mul_sub] at hzero
    exact sub_eq_zero.mp ((mul_eq_zero.mp hzero).resolve_left (pow_ne_zero _ hqC))
  have hx1 : x-1 ≠ 0 := by
    intro hh
    have hx' : x=1 := sub_eq_zero.mp hh
    rw [hx'] at he
    have hz : (3+(a:ℂ)) ^ (2^k) = 0 := by norm_num at he ⊢; exact he
    have hz' : (a:ℂ) = -3 := by
      have := (pow_eq_zero hz)
      linear_combination this
    apply h3
    exact_mod_cast hz'
  let z : ℂ := (x+2+(a:ℂ))/(x-1)
  have hz : z ^ (2^k) = 1 := by
    dsimp [z]
    rw [div_pow, he, div_self (pow_ne_zero _ hx1)]
  apply GrowthRationalParameterClassification.s3_ne_zero k hk a h3 hm hp x z hx hz
  rw [GrowthRationalUnit.s3_square_family]
  have hh : (x-1)*z-(x+2+(a:ℂ)) = 0 := by
    dsimp [z]
    field_simp
  rw [hh]
  simp

theorem matrix_det_ne_zero (k : ℕ) (hk : 2 ≤ k) (a : ℚ)
    (m : ℤ) (q : ℕ) (hq : 0 < q) (hrep : (m:ℚ) = (q:ℚ)*a)
    (h3 : a ≠ -3) (hm : a ≠ -1) (hp : a ≠ 1) :
    (matrix k m q ℤ).det ≠ 0 := by
  let ζ : ℂ := Complex.exp (2*↑Real.pi*Complex.I / ↑(2^k : ℕ))
  have root : IsPrimitiveRoot ζ (2^k) :=
    Complex.isPrimitiveRoot_exp (2^k) (by positivity : 2^k ≠ 0)
  let r : Fin (2^k) → ℂ := fun j => ζ ^ j.val
  have hr : ∀ j, r j ^ (2^k) = 1 := by
    intro j
    change (ζ ^ j.val) ^ (2^k) = 1
    rw [← pow_mul, Nat.mul_comm, pow_mul, root.pow_eq_one, one_pow]
  have hinj : Function.Injective r := by
    intro i j hij
    exact Fin.ext (root.pow_inj i.isLt j.isLt (by simpa only [r] using hij))
  let V : Matrix (Fin (2^k)) (Fin (2^k)) ℂ := (Matrix.vandermonde r).transpose
  let g : Fin (2^k) → ℂ := fun j => value k m q (r j)
  have hv : V.det ≠ 0 := by
    rw [Matrix.det_transpose]
    exact Matrix.det_vandermonde_ne_zero_iff.mpr hinj
  have hg : (Matrix.diagonal g).det ≠ 0 := by
    rw [Matrix.det_diagonal]
    exact Finset.prod_ne_zero_iff.mpr (fun j _ =>
      value_ne_zero_complex k hk a m q hq hrep h3 hm hp (r j) (hr j))
  have he : matrix k m q ℂ * V = V * Matrix.diagonal g := by
    ext i j
    have hh := congrFun (matrix_mulVec k m q (r j) (hr j)) i
    simpa [Matrix.mul_apply, Matrix.mulVec, dotProduct, V,
      Matrix.vandermonde, g, Matrix.diagonal_apply, mul_comm] using hh
  have hnz : (matrix k m q ℂ).det ≠ 0 := by
    intro hz
    have hh := congrArg Matrix.det he
    rw [Matrix.det_mul, Matrix.det_mul, hz, zero_mul] at hh
    exact (mul_ne_zero hv hg) hh.symm
  intro hz
  apply hnz
  rw [← map_matrix (Int.castRingHom ℂ), ← RingHom.map_det, hz, map_zero]

theorem intCast_ne_zero (z : ℤ) (hz : z ≠ 0) (p : ℕ) [Fact p.Prime]
    (hp : z.natAbs < p) : (z : ZMod p) ≠ 0 := by
  intro h
  have hd : (p:ℤ) ∣ z := (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp h
  have hd' : p ∣ z.natAbs := Int.natCast_dvd_natCast.mp (Int.dvd_natAbs.mpr hd)
  have hh := Nat.le_of_dvd (Int.natAbs_pos.mpr hz) hd'
  omega

theorem matrix_det_mod_prime (k : ℕ) (hk : 2 ≤ k) (a : ℚ)
    (m : ℤ) (q : ℕ) (hq : 0 < q) (hrep : (m:ℚ) = (q:ℚ)*a)
    (h3 : a ≠ -3) (hm : a ≠ -1) (ha1 : a ≠ 1)
    (p : ℕ) [Fact p.Prime] (hp : detBound k m q < p) :
    (matrix k m q (ZMod p)).det ≠ 0 := by
  have hb : (matrix k m q ℤ).det.natAbs ≤ detBound k m q := by
    have h := matrix_det_bound k m q
    rw [Int.abs_eq_natAbs] at h
    exact_mod_cast h
  have hh := intCast_ne_zero _ (matrix_det_ne_zero k hk a m q hq hrep h3 hm ha1) p
    (lt_of_le_of_lt hb hp)
  intro hz
  apply hh
  change (Int.castRingHom (ZMod p)) (matrix k m q ℤ).det = 0
  rw [RingHom.map_det, map_matrix]
  exact hz

theorem no_grid_zero (k p : ℕ) [Fact p.Prime] (m : ℤ) (q : ℕ) (a : ZMod p)
    (hrep : (m:ZMod p) = (q:ZMod p)*a)
    (hdet : (matrix k m q (ZMod p)).det ≠ 0)
    (x z : ZMod p) (hx : x ^ (2^k) = 1) (hz : z ^ (2^k) = 1) :
    GrowthSecant.s3 a (-1-a) x z 1 ≠ 0 := by
  intro h
  have hl := GrowthRationalUnit.parameter_of_s3_zero a x z h
  have ht : (q:ZMod p)*x + ((m+2*(q:ℤ):ℤ):ZMod p) =
      ((q:ZMod p)*(x-1))*z := by
    push_cast
    rw [hrep]
    linear_combination (q:ZMod p)*hl
  have he : value k m q x = 0 := by
    unfold value
    rw [ht, mul_pow, hz, mul_one, sub_self]
  have hv := matrix_mulVec k m q x hx
  rw [he, zero_smul] at hv
  have hzero := Matrix.eq_zero_of_mulVec_eq_zero hdet hv
  have hc := congrFun hzero (⟨0, by positivity⟩ : Fin (2^k))
  simp at hc

open MvPolynomial PrimeM2General in
theorem unit_prime (k p : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (m : ℤ) (q : ℕ) (a : ZMod p)
    (hrep : (m:ZMod p) = (q:ZMod p)*a)
    (hdet : (matrix k m q (ZMod p)).det ≠ 0) (hadmissible : 2 ^ k ∣ p - 1) :
    (1 : P k (ZMod p)) ∈
      Ideal.span (Set.range (inputs k depth a (-1-a) 1)) := by
  classical
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
    exact UnitInstance.var_pow_mem k depth a (-1-a) 1 i
  · intro x hx
    have roots : ∀ i, x i ^ 2 ^ k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth a (-1-a) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth a (-1-a) 1
      (active k depth 0) (by simp [active]; omega)
    have second := UnitInstance.chain_step_mem k depth a (-1-a) 1
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
      (raw_normalization_identity k depth a (-1-a) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, GrowthRationalChainUnit.eval_rawS3_rational] at identity
    exact no_grid_zero k p m q a hrep hdet _ _ (roots (active k depth 0)) (roots (active k depth 2)) identity


/-- Reduction is defined by numerator/denominator; no map Q -> Fp is assumed. -/
def reduction (a : ℚ) (p : ℕ) [Fact p.Prime] : ZMod p :=
  (a.num : ZMod p) / (a.den : ZMod p)

theorem denominator_and_parameter_guards (k : ℕ) (a : ℚ) (m : ℤ) (q : ℕ)
    (hq : 0 < q) (hrep : (m:ℚ) = (q:ℚ)*a)
    (h3 : a ≠ -3) (hs : a ≠ -(3/4:ℚ))
    (p : ℕ) [Fact p.Prime] (hp : cutoff k m q < p) :
    (q:ZMod p) ≠ 0 ∧
      ((m:ZMod p)/(q:ZMod p))+3 ≠ 0 ∧
      4*((m:ZMod p)/(q:ZMod p))+3 ≠ 0 := by
  have hqQ : (q:ℚ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hqZ : (q:ℤ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hqFp : (q:ZMod p) ≠ 0 := by
    have hqp : q < p := by unfold cutoff at hp; omega
    exact_mod_cast intCast_ne_zero (q:ℤ) hqZ p (by simpa using hqp)
  have hm3 : m+3*(q:ℤ) ≠ 0 := by
    intro hz
    have hzQ : (m:ℚ)+3*(q:ℚ) = 0 := by exact_mod_cast hz
    have hh : (q:ℚ)*(a+3) = 0 := by linear_combination hzQ - hrep
    have hh' := (mul_eq_zero.mp hh).resolve_left hqQ
    exact h3 (by linarith)
  have hm4 : 4*m+3*(q:ℤ) ≠ 0 := by
    intro hz
    have hzQ : 4*(m:ℚ)+3*(q:ℚ) = 0 := by exact_mod_cast hz
    have hh : (q:ℚ)*(4*a+3) = 0 := by linear_combination hzQ - 4*hrep
    have hh' := (mul_eq_zero.mp hh).resolve_left hqQ
    exact hs (by linarith)
  have h3Fp := intCast_ne_zero _ hm3 p
    (show (m+3*(q:ℤ)).natAbs < p by unfold cutoff at hp; omega)
  have h4Fp := intCast_ne_zero _ hm4 p
    (show (4*m+3*(q:ℤ)).natAbs < p by unfold cutoff at hp; omega)
  refine ⟨hqFp, ?_, ?_⟩
  · intro hz
    apply h3Fp
    push_cast
    field_simp [hqFp] at hz
    linear_combination hz
  · intro hz
    apply h4Fp
    push_cast
    field_simp [hqFp] at hz
    linear_combination hz

/-- Explicit all-rational prime-field growth for a=m/q, with q>0. -/
theorem prime_lfd_window (k : ℕ) (hk : 24 ≤ k) (a : ℚ)
    (m : ℤ) (q : ℕ) (hq : 0 < q) (hrep : (m:ℚ) = (q:ℚ)*a)
    (h3 : a ≠ -3) (hs : a ≠ -(3/4:ℚ)) (hm : a ≠ -1) (ha1 : a ≠ 1)
    (p : ℕ) [Fact p.Prime] (hp : cutoff k m q < p)
    (hadmissible : 2 ^ k ∣ p - 1) :
    (q:ZMod p) ≠ 0 ∧ ∃ d : ℕ, IsLastFallDegree
      (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
        ((m:ZMod p)/(q:ZMod p)) (-1-(m:ZMod p)/(q:ZMod p)) 1)) d ∧
      (k-10)/14 < d ∧ d ≤ k+2 := by
  have depth : 2 ≤ k := by omega
  obtain ⟨hqFp,h3Fp,h4Fp⟩ :=
    denominator_and_parameter_guards k a m q hq hrep h3 hs p hp
  have hpDet : detBound k m q < p := by unfold cutoff at hp; omega
  have hpLarge : 2*2^k+1 ≤ p := by unfold cutoff at hp; omega
  have hdet := matrix_det_mod_prime k depth a m q hq hrep h3 hm ha1 p hpDet
  have hrepFp : (m:ZMod p) = (q:ZMod p)*((m:ZMod p)/(q:ZMod p)) := by
    field_simp [hqFp]
  have hu := unit_prime k p depth m q _ hrepFp hdet hadmissible
  have ho := GrowthRationalPrimeNative.not_one_of_large_prime p k _ depth
    ((k-10)/14) (by omega) (by omega) hpLarge h3Fp
  have htwo : (2:ZMod p) ≠ 0 := by
    apply PrimeFieldConclusions.two_ne_zero p
    have hh : (2:ℕ)^2 ≤ 2^k := Nat.pow_le_pow_right (by norm_num) depth
    omega
  have hns := GrowthRationalUnit.nonsingular_of_parameter_conditions _ h3Fp h4Fp
  exact ⟨hqFp, GrowthPrimeNativeWindow.exists_lfd_window_of_native_certificates
    k depth hk _ _ _ htwo hns hu ho⟩

/-- Canonical numerator and positive denominator give an explicit cutoff for every a. -/
theorem rational_prime_growth (a : ℚ)
    (h3 : a ≠ -3) (hs : a ≠ -(3/4:ℚ)) (hm : a ≠ -1) (ha1 : a ≠ 1)
    (k : ℕ) (hk : 24 ≤ k) (p : ℕ) [Fact p.Prime]
    (hp : cutoff k a.num a.den < p) (hadmissible : 2 ^ k ∣ p - 1) :
    (a.den:ZMod p) ≠ 0 ∧ ∃ d : ℕ, IsLastFallDegree
      (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
        (reduction a p) (-1-reduction a p) 1)) d ∧
      (k-10)/14 < d ∧ d ≤ k+2 := by
  apply prime_lfd_window k hk a a.num a.den a.den_pos _ h3 hs hm ha1 p hp hadmissible
  have hh := a.num_div_den
  have hq : (a.den:ℚ) ≠ 0 := by exact_mod_cast a.den_nz
  apply (div_eq_iff hq).mp at hh
  linear_combination hh

#print axioms matrix_det_bound
#print axioms value_ne_zero_complex
#print axioms matrix_det_ne_zero
#print axioms matrix_det_mod_prime
#print axioms unit_prime
#print axioms denominator_and_parameter_guards
#print axioms prime_lfd_window
#print axioms rational_prime_growth

end ICFallDegrees.GrowthRationalPrimeCutoff
