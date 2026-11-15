import ICFallDegrees.GrowthGeneralRationalWronskian
import ICFallDegrees.GrowthRatFuncCofactorDegree
import Mathlib.Tactic

/-!
# N8: unconditional nonvanishing of actual polynomial Wronskian cofactors

This composes the now kernel-checked arbitrary-column characteristic-zero
rational Wronskian nonvanishing criterion with the genuine rational-to-
polynomial cofactor determinant factorization and tight degree budget.

There is NO separate assumption that the Wronskian or cofactor is nonzero.
The only hypothesis is actual K-linear independence of the listed rational
monomial functions, exactly what minimality of an N8 circuit supplies when
any one circuit index is omitted. All omitted-index enumeration and two
pole comparisons still remain to be formalized.
-/

noncomputable section

namespace ICFallDegrees.GrowthActualCofactorNonzero

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthGeneralRationalWronskian

variable {K : Type*} [Field K] [CharZero K]

/-- A K-independent family of true rational two-torsion monomials has
a nonzero determinant of the concrete polynomial cofactor matrix.
This uses genuine quotient-rule Wronskians and does not assume their
nonvanishing as a separate input. -/
theorem polynomial_cofactor_ne_zero_of_linearIndependent
    (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ)
    (hind : LinearIndependent K (fun j => monomial c (e j))) :
    (polynomialMatrix r c e).det ≠ 0 := by
  have hW : (derivativeMatrix r (fun j => monomial c (e j))).det ≠ 0 :=
    rational_wronskian_ne_zero_of_linearIndependent r
      (fun j => monomial c (e j)) hind
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

/-- For the *same actual* polynomial cofactor, independence supplies
both nonvanishing and the sharp polynomial degree bound r(r-1).
There is no cofactor or spread assumption hidden in these hypotheses. -/
theorem polynomial_cofactor_nonzero_and_degree
    (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ)
    (hind : LinearIndependent K (fun j => monomial c (e j))) :
    (polynomialMatrix r c e).det ≠ 0 ∧
      (polynomialMatrix r c e).det.natDegree ≤
        ∑ j : Fin r, 2 * j.val := by
  exact ⟨polynomial_cofactor_ne_zero_of_linearIndependent r c e hind,
    concrete_cofactor_natDegree_le_two r c e⟩

end ICFallDegrees.GrowthActualCofactorNonzero
