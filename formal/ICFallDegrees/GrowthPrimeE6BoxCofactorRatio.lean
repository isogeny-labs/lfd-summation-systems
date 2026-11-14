import ICFallDegrees.GrowthPrimeE6BoxCofactor
import ICFallDegrees.GrowthActualCofactorQuotient
import Mathlib.Tactic

/-!
# N8: characteristic-bounded genuine E6 cofactor ratio

An independent finite E6 exponent-box family with nonzero replacement
coefficient has a proved nonzero polynomial cofactor and its exact
replacement ratio in the genuine rational function field, even in
positive characteristic above 2N. This reuses the characteristic-free
Wronskian replacement and polynomial normalization lemmas.

Unlike the previous global characteristic-zero version, no [CharZero K]
instance is required. In particular this theorem does not postulate
nonzero Wronskians, cofactor fractions, or finite matroid transfer.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6BoxCofactorRatio

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthPrimeE6BoxCofactor
open ICFallDegrees.GrowthActualMinorReplacement
open ICFallDegrees.GrowthActualCofactorCross
open ICFallDegrees.GrowthActualCofactorQuotient
open ICFallDegrees.GrowthRatFuncInfinityDegree

variable {K : Type*} [Field K]

/-- Bounded positive-characteristic E6 replacement implies the exact
true polynomial-cofactor fraction identity, with both its numerator
and denominator proved nonzero. -/
theorem e6_cofactor_ratio_of_box_replacement
    (r N : ℕ) (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial (8 : K) enew =
      ∑ j : Fin r, coeff j • monomial (8 : K) (e j))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
    (polynomialMatrix r (8 : K) er).det ≠ 0 ∧
    monomial (8 : K) enew / monomial (8 : K) (e i) =
      RatFunc.C (coeff i) *
        RatFunc.mk ((polynomialMatrix r (8 : K) e).det)
          ((polynomialMatrix r (8 : K) er).det) := by
  dsimp only
  let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
  let A : Polynomial K := (polynomialMatrix r (8 : K) e).det
  let B : Polynomial K := (polynomialMatrix r (8 : K) er).det
  have hA : A ≠ 0 :=
    e6_cofactor_ne_zero_of_box r N e hbox hind hchar
  have hα : algebraMap (Polynomial K) (RatFunc K) A ≠ 0 :=
    RatFunc.algebraMap_ne_zero hA
  have hγ : RatFunc.C (coeff i) ≠ 0 := by
    rw [← RatFunc.algebraMap_C]
    exact RatFunc.algebraMap_ne_zero
      (show Polynomial.C (coeff i) ≠ 0 by simpa using hcoeff)
  have hf : monomial (8 : K) (e i) ≠ 0 :=
    monomial_ne_zero (8 : K) (e i).1 (e i).2
  have hcross := cofactor_cross_of_const_relation
    r (8 : K) e enew coeff i hrel
  dsimp only at hcross
  simp only [cofactorValue_eq_map_det, RatFunc.algebraMap_eq_C] at hcross
  have hcross' :
      algebraMap (Polynomial K) (RatFunc K) B *
          monomial (8 : K) enew =
        RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A *
            monomial (8 : K) (e i) := by
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
      monomial (8 : K) enew / monomial (8 : K) (e i) =
        (algebraMap (Polynomial K) (RatFunc K) B *
          monomial (8 : K) enew) /
          (algebraMap (Polynomial K) (RatFunc K) B *
            monomial (8 : K) (e i)) := by
              field_simp [hβ, hf] <;> ring
      _ = (RatFunc.C (coeff i) *
          algebraMap (Polynomial K) (RatFunc K) A *
            monomial (8 : K) (e i)) /
          (algebraMap (Polynomial K) (RatFunc K) B *
            monomial (8 : K) (e i)) := by
              rw [hcross']
      _ = RatFunc.C (coeff i) * RatFunc.mk A B := by
              rw [RatFunc.mk_eq_div]
              field_simp [hβ, hf] <;> ring

#print axioms e6_cofactor_ratio_of_box_replacement

end ICFallDegrees.GrowthPrimeE6BoxCofactorRatio
