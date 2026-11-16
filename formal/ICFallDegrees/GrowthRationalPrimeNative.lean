import ICFallDegrees.GrowthPrimeE6NativeNonmembership
import ICFallDegrees.GrowthPrimeE6BoxCofactorRatio
import ICFallDegrees.GrowthPrimeE6BoxReplacementGap
import Mathlib.Tactic

/-!
# Native bounded-box nonmembership for every two-torsion parameter

Here c=a+2, so the rational section is (X,(X+c)/(X-1)).
The finite-box denominator, degree and cofactor arguments work for arbitrary c;
the pole-gap bound needs c+1 != 0. This generalizes the fixed-E6 proof path
without changing that path or assuming specialization of circuits from Q.
The final theorem excludes 1 from the actual native InClosure for every
prime characteristic greater than 2^(k+1), when a+3 != 0.
-/
noncomputable section
namespace ICFallDegrees.GrowthRationalPrimeNative
variable {K : Type*} [Field K] (c : K)

open Polynomial
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncDenominatorIndependence
open ICFallDegrees.GrowthPrimeBoundedRationalWronskian



def boxDenominator (N : ℕ) : Polynomial K :=
  (Polynomial.X - 1) ^ (N - 1)

def boxNumerator (N : ℕ) (e : ℕ × ℕ) : Polynomial K :=
  Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C c) ^ e.2 *
    (Polynomial.X - 1) ^ (N - 1 - e.2)

theorem boxDenominator_ne_zero (N : ℕ) :
    boxDenominator (K := K) N ≠ 0 := by
  have ht : (Polynomial.X - 1 : Polynomial K) ≠ 0 := by
    simpa only [Polynomial.C_1] using
      (Polynomial.X_sub_C_ne_zero (1 : K))
  exact pow_ne_zero _ ht

theorem box_monomial_cleared (N : ℕ) (e : ℕ × ℕ)
    (he : e.2 < N) :
    algebraMap (Polynomial K) (RatFunc K) (boxDenominator (K := K) N) *
        monomial c e =
      algebraMap (Polynomial K) (RatFunc K) (boxNumerator c (K := K) N e) := by
  let T : Polynomial K := Polynomial.X - 1
  let U : Polynomial K :=
    Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C c) ^ e.2
  have hle : e.2 ≤ N - 1 := by omega
  have ht : (algebraMap (Polynomial K) (RatFunc K) T) ≠ 0 := by
    have hT : T ≠ 0 := by
      change (Polynomial.X - 1 : Polynomial K) ≠ 0
      simpa only [Polynomial.C_1] using
        (Polynomial.X_sub_C_ne_zero (1 : K))
    exact RatFunc.algebraMap_ne_zero hT
  change algebraMap (Polynomial K) (RatFunc K) (T ^ (N - 1)) *
      RatFunc.mk U (T ^ e.2) =
    algebraMap (Polynomial K) (RatFunc K) (U * T ^ (N - 1 - e.2))
  calc
    algebraMap (Polynomial K) (RatFunc K) (T ^ (N - 1)) *
        RatFunc.mk U (T ^ e.2) =
      (algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1) *
        (algebraMap (Polynomial K) (RatFunc K) U /
          (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) := by
          rw [RatFunc.mk_eq_div, map_pow, map_pow]
    _ = ((algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1 - e.2) *
        (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) *
        (algebraMap (Polynomial K) (RatFunc K) U /
          (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) := by
            rw [← pow_add, Nat.sub_add_cancel hle]
    _ = algebraMap (Polynomial K) (RatFunc K) U *
        (algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1 - e.2) := by
          field_simp [pow_ne_zero _ ht]
          ring
    _ = algebraMap (Polynomial K) (RatFunc K)
        (U * T ^ (N - 1 - e.2)) := by
          simp only [map_mul, map_pow, U]

theorem boxNumerator_linearIndependent (r N : ℕ)
    (e : Fin r → ℕ × ℕ)
    (he : ∀ j, (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j))) :
    LinearIndependent K (fun j => boxNumerator c (K := K) N (e j)) := by
  exact polynomial_linearIndependent_of_clear r
    (fun j => monomial c (e j))
    (boxDenominator (K := K) N)
    (fun j => boxNumerator c (K := K) N (e j))
    (boxDenominator_ne_zero N)
    (fun j => box_monomial_cleared c N (e j) (he j)) hind

theorem box_wronskian_ne_zero_of_degree_bound
    (r N k : ℕ) (e : Fin r → ℕ × ℕ)
    (he : ∀ j, (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hdeg : ∀ j, (boxNumerator c (K := K) N (e j)).natDegree < k)
    (hchar : ringChar K = 0 ∨ k ≤ ringChar K) :
    (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0 := by
  exact rational_wronskian_ne_zero_of_bounded_polynomial_representatives
    r k (fun j => monomial c (e j))
    (boxDenominator (K := K) N)
    (fun j => boxNumerator c (K := K) N (e j))
    (boxDenominator_ne_zero N)
    (fun j => box_monomial_cleared c N (e j) (he j))
    (boxNumerator_linearIndependent c r N e he hind)
    hdeg hchar

theorem boxNumerator_natDegree_lt_two_mul (N : ℕ) (e : ℕ × ℕ)
    (hu : e.1 < N) (hv : e.2 < N) :
    (boxNumerator c (K := K) N e).natDegree < 2 * N := by
  have hx : (Polynomial.X : Polynomial K).natDegree ≤ 1 :=
    Polynomial.natDegree_X_le
  have ha : (Polynomial.X + Polynomial.C c).natDegree ≤ 1 :=
    Polynomial.natDegree_add_le_of_degree_le hx (by simp)
  have ht : (Polynomial.X - 1 : Polynomial K).natDegree ≤ 1 := by
    compute_degree
  have h0 : (Polynomial.X ^ e.1 : Polynomial K).natDegree ≤ e.1 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le e.1 hx)
  have h1 : ((Polynomial.X + Polynomial.C c) ^ e.2).natDegree ≤ e.2 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le e.2 ha)
  have h2 : ((Polynomial.X - 1 : Polynomial K) ^ (N - 1 - e.2)).natDegree ≤
      N - 1 - e.2 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le (N - 1 - e.2) ht)
  have hleft :
      (Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C c) ^ e.2).natDegree ≤
        e.1 + e.2 :=
    Polynomial.natDegree_mul_le_of_le h0 h1
  have htot :
      (boxNumerator c (K := K) N e).natDegree ≤
        (e.1 + e.2) + (N - 1 - e.2) := by
    change (Polynomial.X ^ e.1 *
      (Polynomial.X + Polynomial.C c) ^ e.2 *
      (Polynomial.X - 1) ^ (N - 1 - e.2)).natDegree ≤
        (e.1 + e.2) + (N - 1 - e.2)
    exact Polynomial.natDegree_mul_le_of_le hleft h2
  omega

theorem box_wronskian_ne_zero_of_box_char_bound
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0 := by
  exact box_wronskian_ne_zero_of_degree_bound c r N (2 * N) e
    (fun j => (hbox j).2) hind
    (fun j => boxNumerator_natDegree_lt_two_mul c N (e j)
      (hbox j).1 (hbox j).2) hchar



open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree



theorem parameter_cofactor_ne_zero_of_box
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (polynomialMatrix r c e).det ≠ 0 := by
  have hW : (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0 :=
    box_wronskian_ne_zero_of_box_char_bound c r N e hbox hind hchar
  have hmap := embedded_cofactor_ne_zero_of_wronskian_ne_zero r c e hW
  let φ : Polynomial K →+* RatFunc K :=
    algebraMap (Polynomial K) (RatFunc K)
  have hmatrix :
      (Matrix.of fun i j : Fin r => φ (polynomialMatrix r c e i j)) =
        φ.mapMatrix (polynomialMatrix r c e) := by
    ext i j
    rfl
  rw [hmatrix, ← RingHom.map_det] at hmap
  intro hzero
  exact hmap (by rw [hzero, map_zero])

theorem parameter_cofactor_nonzero_and_degree_of_box
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (polynomialMatrix r c e).det ≠ 0 ∧
    (polynomialMatrix r c e).det.natDegree ≤
      ∑ j : Fin r, 2 * j.val := by
  exact ⟨parameter_cofactor_ne_zero_of_box c r N e hbox hind hchar,
    concrete_cofactor_natDegree_le_two r c e⟩



open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthActualMinorReplacement
open ICFallDegrees.GrowthActualCofactorCross
open ICFallDegrees.GrowthActualCofactorQuotient
open ICFallDegrees.GrowthRatFuncInfinityDegree



theorem parameter_cofactor_ratio_of_box_replacement
    (r N : ℕ) (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
    (polynomialMatrix r c er).det ≠ 0 ∧
    monomial c enew / monomial c (e i) =
      RatFunc.C (coeff i) *
        RatFunc.mk ((polynomialMatrix r c e).det)
          ((polynomialMatrix r c er).det) := by
  dsimp only
  let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
  let A : Polynomial K := (polynomialMatrix r c e).det
  let B : Polynomial K := (polynomialMatrix r c er).det
  have hA : A ≠ 0 :=
    parameter_cofactor_ne_zero_of_box c r N e hbox hind hchar
  have hα : algebraMap (Polynomial K) (RatFunc K) A ≠ 0 :=
    RatFunc.algebraMap_ne_zero hA
  have hγ : RatFunc.C (coeff i) ≠ 0 := by
    rw [← RatFunc.algebraMap_C]
    exact RatFunc.algebraMap_ne_zero
      (show Polynomial.C (coeff i) ≠ 0 by simpa using hcoeff)
  have hf : monomial c (e i) ≠ 0 :=
    monomial_ne_zero c (e i).1 (e i).2
  have hcross := cofactor_cross_of_const_relation
    r c e enew coeff i hrel
  dsimp only at hcross
  simp only [cofactorValue_eq_map_det, RatFunc.algebraMap_eq_C] at hcross
  have hcross' :
      algebraMap (Polynomial K) (RatFunc K) B *
          monomial c enew =
        RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A *
            monomial c (e i) := by
    simpa only [A, B, er] using hcross
  have hβ : algebraMap (Polynomial K) (RatFunc K) B ≠ 0 := by
    intro hz
    rw [hz, zero_mul] at hcross'
    exact (mul_ne_zero (mul_ne_zero hγ hα) hf) hcross'.symm
  have hB : B ≠ 0 := by
    intro hz
    apply hβ
    rw [hz, map_zero]
  constructor
  · exact hB
  · calc
      monomial c enew / monomial c (e i) =
        (algebraMap (Polynomial K) (RatFunc K) B *
          monomial c enew) /
          (algebraMap (Polynomial K) (RatFunc K) B *
            monomial c (e i)) := by
              field_simp [hβ, hf]; ring
      _ = (RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A *
            monomial c (e i)) /
          (algebraMap (Polynomial K) (RatFunc K) B *
            monomial c (e i)) := by
              rw [hcross']
      _ = RatFunc.C (coeff i) * RatFunc.mk A B := by
              rw [RatFunc.mk_eq_div]
              field_simp [hβ, hf]; ring



open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthFiniteCofactorGap



theorem parameter_box_replacement_l1_le_four_budget
    (r N : ℕ) (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K)
    (h9 : c + 1 ≠ 0) :
    |2 * ((enew.1 : ℤ) - (e i).1)| +
      |2 * ((enew.2 : ℤ) - (e i).2)| ≤
        4 * ((∑ j : Fin r, 2 * j.val : ℕ) : ℤ) := by
  let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
  let A : Polynomial K := (polynomialMatrix r c e).det
  let B : Polynomial K := (polynomialMatrix r c er).det
  have hA : A ≠ 0 :=
    parameter_cofactor_ne_zero_of_box c r N e hbox hind hchar
  have hq := parameter_cofactor_ratio_of_box_replacement c r N e enew
    coeff i hbox hind hcoeff hrel hchar
  dsimp only at hq
  obtain ⟨hB, hratio⟩ := hq
  have hdA : A.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r c e
  have hdB : B.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r c er
  exact exponent_pair_gap_l1_le_four_budget
    c (coeff i) h9 enew.1 enew.2 (e i).1 (e i).2
    A B (∑ j : Fin r, 2 * j.val)
    hA hB hcoeff hdA hdB hratio



open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthCircuitInterface
open ICFallDegrees.GrowthCircuitErasedBasis
open ICFallDegrees.GrowthCircuitReplacementFin



theorem parameter_bounded_circuit_pair_le_four_budget
    (N : ℕ)
    (J : Finset (Fin 2 →₀ ℕ))
    (hbox : ∀ d ∈ J, d 0 < N ∧ d 1 < N)
    (hcirc : Circuit K (mono (xi (c - 2))) J)
    (x y : Fin 2 →₀ ℕ)
    (hx : x ∈ J) (hy : y ∈ J) (hxy : y ≠ x)
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K)
    (h9 : c + 1 ≠ 0) :
    |2 * ((x 0 : ℤ) - y 0)| +
      |2 * ((x 1 : ℤ) - y 1)| ≤
        4 * ((∑ j : Fin (J.erase x).card, 2 * j.val : ℕ) : ℤ) := by
  classical
  let er : Fin (J.erase x).card → (Fin 2 →₀ ℕ) :=
    fun i => (((Finset.equivFin (J.erase x)).symm i).1)
  let e : Fin (J.erase x).card → ℕ × ℕ :=
    fun i => (er i 0, er i 1)
  let enew : ℕ × ℕ := (x 0, x 1)
  have h8 : (c - 2) + 2 = c := by ring
  have hind0 :=
    erase_reindex_linearIndependent (mono (xi (c - 2))) J hcirc x hx
  have hind : LinearIndependent K
      (fun j : Fin (J.erase x).card => monomial c (e j)) := by
    simpa only [e, er, mono_xi_eq_monomial, h8] using hind0
  have hboxe : ∀ j : Fin (J.erase x).card, (e j).1 < N ∧ (e j).2 < N := by
    intro j
    have hj : er j ∈ J.erase x :=
      (((Finset.equivFin (J.erase x)).symm j).2)
    exact hbox (er j) (Finset.mem_of_mem_erase hj)
  obtain ⟨coeff, hcoeff, hrel⟩ :=
    exists_replacement_fin (mono (xi (c - 2))) J hcirc x hx
  have hrel' : monomial c enew =
      ∑ j : Fin (J.erase x).card, coeff j • monomial c (e j) := by
    simpa only [enew, e, er, mono_xi_eq_monomial, h8] using hrel
  obtain ⟨i, hi⟩ := erase_reindex_contains J x y hx hy hxy
  have hbound := parameter_box_replacement_l1_le_four_budget c
    (J.erase x).card N e enew coeff i hboxe hind
    (hcoeff i) hrel' hchar h9
  simpa only [enew, e, er, hi] using hbound



open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthBoxSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthActualCircuitSpread



theorem parameter_native_box_spread
    (p k : ℕ)
    (hchar : ringChar K = 0 ∨ 2 * 2 ^ k ≤ ringChar K)
    (h9 : c + 1 ≠ 0) :
    BoxSpreadBound K (xi (c - 2)) p k := by
  intro J hbox hcirc x hx y hy
  refine ⟨0, 0, ?_⟩
  simp only [Nat.cast_zero, zero_mul, mul_zero, sub_zero]
  have hcard : 1 ≤ J.card := Finset.card_pos.mpr ⟨x, hx⟩
  by_cases heq : y = x
  · subst y
    simp only [sub_self, mul_zero, abs_zero, zero_add]
    by_cases hone : J.card = 1
    · simp [hone]
    · have htwo : 2 ≤ J.card := by omega
      have htwoz : (2 : ℤ) ≤ (J.card : ℤ) := by exact_mod_cast htwo
      have hp : 0 ≤ ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2) :=
        mul_nonneg (by omega) (by omega)
      nlinarith
  · have hxy : y ≠ x := heq
    have hyerase : y ∈ J.erase x := Finset.mem_erase.mpr ⟨hxy, hy⟩
    have hrpos : 0 < (J.erase x).card := Finset.card_pos.mpr ⟨y, hyerase⟩
    have hr1 : (J.erase x).card = J.card - 1 :=
      Finset.card_erase_of_mem hx
    have hr2 : (J.erase x).card - 1 = J.card - 2 := by omega
    have htwo : 2 ≤ J.card := by omega
    have htwoz : (2 : ℤ) ≤ (J.card : ℤ) := by exact_mod_cast htwo
    have hbudget :=
      parameter_bounded_circuit_pair_le_four_budget c
        (2 ^ k) J hbox hcirc x y hx hy hxy hchar h9
    rw [sum_twice_fin_eq_budget, Nat.cast_mul, hr2, hr1] at hbudget
    have hcast1 : ((J.card - 1 : ℕ) : ℤ) = (J.card : ℤ) - 1 := by omega
    have hcast2 : ((J.card - 2 : ℕ) : ℤ) = (J.card : ℤ) - 2 := by omega
    rw [hcast1, hcast2] at hbudget
    have hp : 0 ≤ ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2) :=
      mul_nonneg (by omega) (by omega)
    nlinarith



open MvPolynomial GrowthBoxSpread GrowthReduction
theorem not_one_of_box_spread
    (p k : ℕ) [Fact p.Prime] (a : ZMod p) (depth : 2 ≤ k)
    (D : ℕ) (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p)
    (hspread : BoxSpreadBound (ZMod p)
      (GrowthRatFuncModel.xi a) p k) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth
        a (-1-a) 1)) D
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) := by
  let ξ : Fin 2 → RatFunc (ZMod p) := GrowthRatFuncModel.xi a
  have hdata := GrowthRatFuncModel.ratFunc_algebraic_data a
  have hcurve : GrowthSpread.S3at
      a (-1-a) 1 ξ = 0 := by
    simpa only [ξ, show (-1 - a) = (-1-a) by ring]
      using hdata.2
  have hlarge : 5 ≤ p := by
    have hpow : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
    omega
  have htwo : (2 : ZMod p) ≠ 0 :=
    PrimeFieldConclusions.two_ne_zero p hlarge
  have hfei : FEI
      (RingHom.ker (aeval ξ :
        MvPolynomial (Fin 2) (ZMod p) →ₐ[ZMod p] RatFunc (ZMod p)))
      k (D - 1) := by
    apply GrowthBoxSpread.fei_of_box_spread
      (K := ZMod p) (F₀ := ZMod p) (F := RatFunc (ZMod p))
      ξ p k (D - 1) (13 * (D - 1) + 23) (Or.inr hp)
      hspread (GrowthCounting.lambdaBound_linear (D - 1))
      (by omega) (by omega)
  exact GrowthReduction.not_one_of_fei htwo k depth
    a (-1-a) 1
    (RingHom.ker (aeval ξ :
       MvPolynomial (Fin 2) (ZMod p) →ₐ[ZMod p] RatFunc (ZMod p)))
    (GrowthMain.X_not_mem_ker hdata.1)
    (GrowthMain.S3poly_mem_ker hcurve)
    D hfei

theorem not_one_of_large_prime (p k : ℕ) [Fact p.Prime] (a : ZMod p)
    (depth : 2 ≤ k) (D : ℕ) (hD : 1 ≤ D) (hk : 14*D+10 ≤ k)
    (hp : 2*2^k+1 ≤ p) (ha : a+3 ≠ 0) :
    ¬ InClosure (Set.range (PrimeM2General.inputs k depth a (-1-a) 1)) D
      (1 : MvPolynomial (Fin (2*k)) (ZMod p)) := by
  apply not_one_of_box_spread p k a depth D hD hk hp
  have hh : ringChar (ZMod p) = 0 ∨ 2*2^k ≤ ringChar (ZMod p) := by
    right
    rw [ZMod.ringChar_zmod_n]
    omega
  have hc : a+2+1 ≠ 0 := by convert ha using 1; ring
  have h := parameter_native_box_spread (a+2) p k hh hc
  simpa only [add_sub_cancel_right] using h

#print axioms parameter_native_box_spread
#print axioms not_one_of_large_prime

end ICFallDegrees.GrowthRationalPrimeNative
