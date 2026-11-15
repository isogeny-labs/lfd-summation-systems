import ICFallDegrees.GrowthRatFuncWronskianNormalization
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

/-!
# N8: constant basis changes for the *genuine* rational Wronskian

These statements apply to the quotient-rule derivation on RatFunc K.
A constant matrix U acts on columns, and the Wronskian determinant
transforms by its (embedded) determinant. This is a necessary bridge
for using a degree-echelon polynomial basis after common-denominator
clearing. We do *not* claim arbitrary-column Wronskian nonvanishing,
any valuation identity, or a circuit-spread theorem.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncWronskianLinearChange

open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K]

/-- The actual iterated derivative matrix respects any finite
K-constant change of columns. -/
theorem derivativeMatrix_linear_change
    (r : ℕ) (f : Fin r → RatFunc K)
    (U : Matrix (Fin r) (Fin r) K) :
    derivativeMatrix r (fun j => ∑ i : Fin r, U i j • f i) =
      derivativeMatrix r f * U.map (algebraMap K (RatFunc K)) := by
  classical
  ext i j
  simp only [derivativeMatrix, Matrix.of_apply, Matrix.mul_apply, Matrix.map_apply]
  rw [iterDerivative_const_sum]
  apply Finset.sum_congr rfl
  intro k hk
  simp only [Algebra.smul_def]
  ring

/-- The rational Wronskian determinant transforms by the
determinant of the *constant* coefficient matrix. Unlike an
abstract linear-independence assumption, this is a proved
identity for the actual iterated rational derivatives. -/
theorem determinant_linear_change
    (r : ℕ) (f : Fin r → RatFunc K)
    (U : Matrix (Fin r) (Fin r) K) :
    (derivativeMatrix r (fun j => ∑ i : Fin r, U i j • f i)).det =
      (derivativeMatrix r f).det *
        algebraMap K (RatFunc K) U.det := by
  have hmap : U.map (algebraMap K (RatFunc K)) =
      (algebraMap K (RatFunc K)).mapMatrix U := by
    ext i j
    rfl
  rw [derivativeMatrix_linear_change, Matrix.det_mul, hmap, ← RingHom.map_det]

/-- A singular constant matrix necessarily produces a zero
genuine rational Wronskian, in every characteristic. -/
theorem determinant_linear_change_zero_of_singular
    (r : ℕ) (f : Fin r → RatFunc K)
    (U : Matrix (Fin r) (Fin r) K) (hU : U.det = 0) :
    (derivativeMatrix r (fun j => ∑ i : Fin r, U i j • f i)).det = 0 := by
  rw [determinant_linear_change, hU, map_zero, mul_zero]

end ICFallDegrees.GrowthRatFuncWronskianLinearChange
