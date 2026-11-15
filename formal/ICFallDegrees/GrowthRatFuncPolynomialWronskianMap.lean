import ICFallDegrees.GrowthRatFuncWronskianCommonScalar
import ICFallDegrees.GrowthComplexRatFuncDerivative
import Mathlib.Tactic

/-!
# N8: exact polynomial-to-rational Wronskian comparison

Iterating genuine quotient-rule differentiation of the embedded polynomial
agrees with iterating the polynomial derivative. Therefore the all-size
genuine rational Wronskian is the image of the classical polynomial
Wronskian under the *injective* algebra map K[X] -> K(X).

Together with the proven common-scalar Wronskian law, this reduces the
remaining arbitrary-column rational Wronskian nonvanishing lemma to a
polynomial Wronskian independence theorem. It does NOT assert that
independence theorem or assume any circuit spread bound.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncPolynomialWronskianMap

open Polynomial
open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K]

/-- Genuine repeated differentiation commutes with embedding a polynomial
in K(X), for *every* order, not merely the specialized N8 monomials. -/
theorem iterDerivative_algebraMap (n : ℕ) (p : Polynomial K) :
    iterDerivative n (algebraMap (Polynomial K) (RatFunc K) p) =
      algebraMap (Polynomial K) (RatFunc K) ((Polynomial.derivative^[n]) p) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      calc
        iterDerivative (n + 1) (algebraMap (Polynomial K) (RatFunc K) p) =
          rationalDerivative
            (iterDerivative n (algebraMap (Polynomial K) (RatFunc K) p)) := rfl
        _ = algebraMap (Polynomial K) (RatFunc K)
            (Polynomial.derivative ((Polynomial.derivative^[n]) p)) := by
              rw [ih, rationalDerivative_algebraMap]
        _ = algebraMap (Polynomial K) (RatFunc K)
            ((Polynomial.derivative^[n + 1]) p) := by
              rw [Function.iterate_succ_apply']

/-- The ordinary polynomial Wronskian matrix, with derivative order
specified by its row. -/
def polynomialDerivativeMatrix (r : ℕ) (p : Fin r → Polynomial K) :
    Matrix (Fin r) (Fin r) (Polynomial K) :=
  Matrix.of (fun i j => (Polynomial.derivative^[i.val]) (p j))

/-- The *entire* rational Wronskian matrix on embedded polynomials
is the entrywise injection of the polynomial Wronskian matrix. -/
theorem derivativeMatrix_polynomialMap (r : ℕ) (p : Fin r → Polynomial K) :
    derivativeMatrix r (fun j => algebraMap (Polynomial K) (RatFunc K) (p j)) =
      (polynomialDerivativeMatrix r p).map
        (algebraMap (Polynomial K) (RatFunc K)) := by
  ext i j
  simp only [derivativeMatrix, polynomialDerivativeMatrix, Matrix.of_apply, Matrix.map_apply]
  exact iterDerivative_algebraMap i.val (p j)

/-- The rational determinant of a polynomial family equals
the injection of its classical polynomial Wronskian determinant. -/
theorem determinant_polynomialMap (r : ℕ) (p : Fin r → Polynomial K) :
    (derivativeMatrix r (fun j =>
      algebraMap (Polynomial K) (RatFunc K) (p j))).det =
    algebraMap (Polynomial K) (RatFunc K) (polynomialDerivativeMatrix r p).det := by
  have hmap :
      (polynomialDerivativeMatrix r p).map
          (algebraMap (Polynomial K) (RatFunc K)) =
        (algebraMap (Polynomial K) (RatFunc K)).mapMatrix
          (polynomialDerivativeMatrix r p) := by
    ext i j
    rfl
  rw [derivativeMatrix_polynomialMap, hmap, ← RingHom.map_det]

/-- Exact nonvanishing equivalence: there is no loss when carrying
a polynomial Wronskian into the actual rational-function field. -/
theorem determinant_polynomialMap_ne_zero_iff (r : ℕ) (p : Fin r → Polynomial K) :
    (derivativeMatrix r (fun j =>
      algebraMap (Polynomial K) (RatFunc K) (p j))).det ≠ 0 ↔
        (polynomialDerivativeMatrix r p).det ≠ 0 := by
  rw [determinant_polynomialMap]
  constructor
  · intro h hz
    exact h (by rw [hz, map_zero])
  · intro h hz
    apply h
    apply RatFunc.algebraMap_injective K
    simpa only [map_zero] using hz

end ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
