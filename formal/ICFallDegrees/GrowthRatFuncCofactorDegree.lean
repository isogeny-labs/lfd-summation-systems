import ICFallDegrees.GrowthRatFuncWronskianNormalization
import ICFallDegrees.GrowthWronskianTightDegree
import Mathlib.Tactic

/-!
# N8: tight polynomial cofactor degree for the actual rational Wronskian

The previous determinant normalization used an actual rational derivative,
but its cofactor matrix was not yet identified with the independently
degree-bounded transposed polynomial matrix. This file makes the bridge
explicit for *every* number of columns and every exponent family.

Nonvanishing of the rational Wronskian remains a separate theorem; no
assumption that the circuit spread holds appears here.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncCofactorDegree

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRationalLogNumerator
open ICFallDegrees.GrowthWronskianCofactor

variable {K : Type*} [Field K]

/-- The exact polynomial matrix obtained from the genuine higher
rational derivatives, transposed into the existing cofactor convention. -/
theorem polynomialMatrix_transpose_eq (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) :
    (polynomialMatrix r c e).transpose =
      matrix r (denominator c)
        (fun j => logNumerator c (e j).1 (e j).2) := by
  ext i j
  rfl

/-- A degree bound for the polynomial determinant in the *actual*
all-size rational Wronskian normalization: at most Σ 2s, with no
hypothesis on coefficients, independence, or the shape of a circuit. -/
theorem concrete_cofactor_natDegree_le_two (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) :
    (polynomialMatrix r c e).det.natDegree ≤
      ∑ j : Fin r, 2 * j.val := by
  have h := ICFallDegrees.GrowthWronskianTightDegree.cofactor_natDegree_le_two
    r (denominator c)
      (fun j => logNumerator c (e j).1 (e j).2)
    (denominator_natDegree_le c)
    (fun j => logNumerator_natDegree_le c (e j).1 (e j).2)
  rw [← polynomialMatrix_transpose_eq r c e] at h
  simpa only [Matrix.det_transpose] using h

/-- The common cubic polynomial denominator is genuinely nonzero,
uniformly over every coefficient field and rational parameter. -/
theorem denominator_ne_zero (c : K) : denominator c ≠ 0 := by
  unfold denominator
  have hx : (X : Polynomial K) ≠ 0 := X_ne_zero
  have hc : (X + C c : Polynomial K) ≠ 0 := by
    simpa using (X_sub_C_ne_zero (-c))
  have h1 : (X - 1 : Polynomial K) ≠ 0 := by
    simpa using (X_sub_C_ne_zero (1 : K))
  exact mul_ne_zero (mul_ne_zero hx hc) h1

/-- In the rational function field the common denominator remains
nonzero: this is a consequence of polynomial-to-fraction injectivity. -/
theorem rational_denominator_ne_zero (c : K) :
    algebraMap (Polynomial K) (RatFunc K) (denominator c) ≠ 0 :=
  by
    simpa only [map_zero] using
      ((RatFunc.algebraMap_injective K).ne (denominator_ne_zero c))

/-- Any nonzero true rational Wronskian forces its actual embedded
polynomial cofactor determinant to be nonzero, without assuming an
independent-family Wronskian theorem. This direction will be used
once such nonvanishing is formally established. -/
theorem embedded_cofactor_ne_zero_of_wronskian_ne_zero
    (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ)
    (hW : (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0) :
    (Matrix.of (fun i j =>
      algebraMap (Polynomial K) (RatFunc K)
        (polynomialMatrix r c e i j))).det ≠ 0 := by
  let d : RatFunc K := algebraMap (Polynomial K) (RatFunc K) (denominator c)
  have hd : d ≠ 0 := rational_denominator_ne_zero c
  have hprod : (∏ i : Fin r, d ^ i.val) ≠ 0 := by
    apply Finset.prod_ne_zero_iff.mpr
    intro i hi
    exact pow_ne_zero _ hd
  have hl : (∏ i : Fin r, d ^ i.val) *
      (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0 :=
    mul_ne_zero hprod hW
  have hfactor := concrete_determinant_factorization r c e
  intro hzero
  rw [hzero, zero_mul] at hfactor
  exact hl hfactor

end ICFallDegrees.GrowthRatFuncCofactorDegree
