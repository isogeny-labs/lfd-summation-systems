import ICFallDegrees.GrowthActualCofactorCross
import ICFallDegrees.GrowthActualCofactorNonzero
import Mathlib.Tactic

/-!
# N8: true polynomial cofactor quotient from a replaced rational Wronskian

This bridge combines the actual all-rank cofactor cross identity with
the characteristic-zero independent-family nonvanishing theorem.
In particular the denominator cofactor is PROVED nonzero, not postulated.

The input here is an actual constant-linear replacement relation with
a nonzero replacement coefficient. Converting arbitrary minimal circuits
to such a replacement and then proving SpreadBound remain separate tasks.
-/

noncomputable section

namespace ICFallDegrees.GrowthActualCofactorQuotient

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthActualMinorReplacement
open ICFallDegrees.GrowthActualCofactorCross
open ICFallDegrees.GrowthActualCofactorNonzero
open ICFallDegrees.GrowthRatFuncInfinityDegree

variable {K : Type*} [Field K] [CharZero K]

omit [CharZero K] in
/-- The genuine rational cofactor is the image of the true polynomial
determinant under the injective polynomial-to-rational-function map. -/
theorem cofactorValue_eq_map_det (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) :
    cofactorValue r c e =
      algebraMap (Polynomial K) (RatFunc K) ((polynomialMatrix r c e).det) := by
  let φ : Polynomial K →+* RatFunc K :=
    algebraMap (Polynomial K) (RatFunc K)
  have hmatrix :
      (Matrix.of (fun i j : Fin r =>
        φ (polynomialMatrix r c e i j))) =
        φ.mapMatrix (polynomialMatrix r c e) := by
    ext i j
    rfl
  unfold cofactorValue
  rw [hmatrix, ← RingHom.map_det]

/-- A constant-linear replacement of an actually independent family
with nonzero replacement coefficient produces a nonzero denominator
cofactor, and forces the exact fraction equality of actual rational
monomials and actual polynomial Wronskian cofactors.

The only relation input is the genuine linear combination, not an
assumed Wronskian identity, cofactor ratio, or spread hypothesis. -/
theorem cofactor_ratio_of_independent_replacement (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j)) :
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
    polynomial_cofactor_ne_zero_of_linearIndependent r c e hind
  have hα : algebraMap (Polynomial K) (RatFunc K) A ≠ 0 :=
    RatFunc.algebraMap_ne_zero hA
  have hγ : RatFunc.C (coeff i) ≠ 0 := by
    rw [← RatFunc.algebraMap_C]
    exact RatFunc.algebraMap_ne_zero
      (show Polynomial.C (coeff i) ≠ 0 by simpa using hcoeff)
  have hf : monomial c (e i) ≠ 0 :=
    monomial_ne_zero c (e i).1 (e i).2
  have hcross := cofactor_cross_of_const_relation r c e enew coeff i hrel
  dsimp only at hcross
  simp only [cofactorValue_eq_map_det, RatFunc.algebraMap_eq_C] at hcross
  have hcross' :
      algebraMap (Polynomial K) (RatFunc K) B * monomial c enew =
        RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A * monomial c (e i) := by
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
          (algebraMap (Polynomial K) (RatFunc K) B * monomial c enew) /
            (algebraMap (Polynomial K) (RatFunc K) B * monomial c (e i)) := by
              field_simp [hβ, hf] <;> ring
      _ = (RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A * monomial c (e i)) /
            (algebraMap (Polynomial K) (RatFunc K) B * monomial c (e i)) := by
              rw [hcross']
      _ = RatFunc.C (coeff i) * RatFunc.mk A B := by
              rw [RatFunc.mk_eq_div]
              field_simp [hβ, hf] <;> ring

end ICFallDegrees.GrowthActualCofactorQuotient
