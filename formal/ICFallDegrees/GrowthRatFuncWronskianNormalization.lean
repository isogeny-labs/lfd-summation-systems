import ICFallDegrees.GrowthRatFuncPolynomialRecurrence
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

/-!
# N8 genuine rational Wronskian determinant normalization

For any finite set of rational monomials, every derivative row of the genuine
Wronskian has a common denominator power, and every column has a common
monomial factor. Multiplying by the diagonal row/column matrices and taking
determinants gives an exact **all-size** rational-vs-polynomial Wronskian
identity, with no asserted Wronskian nonvanishing.

The arbitrary-size *nonvanishing* theorem is still required before the
cofactor comparison and two-pole spread estimates can be concluded.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncWronskianNormalization

open Polynomial
open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncPolynomialRecurrence
open ICFallDegrees.GrowthRationalLogNumerator
open ICFallDegrees.GrowthRationalWronskianRec

variable {K : Type*} [Field K]

/-- The Wronskian matrix uses genuine iterated quotient-rule derivatives
as rows, and the given rational functions as columns. -/
def derivativeMatrix (r : ℕ) (f : Fin r → RatFunc K) :
    Matrix (Fin r) (Fin r) (RatFunc K) :=
  Matrix.of (fun i j => iterDerivative i.val (f j))

/-- Factor determinant of any Wronskian satisfying an exact
row-and-column scaled numerator relation. All row and column indices
are arbitrary Fin r; no characteristic-zero or independence assumption
is involved in the factorization. -/
theorem det_scaled (r : ℕ) (f : Fin r → RatFunc K) (d : RatFunc K)
    (P : Matrix (Fin r) (Fin r) (RatFunc K))
    (h : ∀ i j, iterDerivative i.val (f j) * d ^ i.val = f j * P i j) :
    (∏ i : Fin r, d ^ i.val) * (derivativeMatrix r f).det =
      P.det * (∏ j : Fin r, f j) := by
  classical
  have hM :
      Matrix.diagonal (fun i : Fin r => d ^ i.val) * derivativeMatrix r f =
        P * Matrix.diagonal f := by
    ext i j
    simp only [Matrix.diagonal_mul, Matrix.mul_diagonal]
    change d ^ i.val * iterDerivative i.val (f j) = P i j * f j
    have hij := h i j
    calc
      d ^ i.val * iterDerivative i.val (f j) =
        iterDerivative i.val (f j) * d ^ i.val := mul_comm _ _
      _ = f j * P i j := hij
      _ = P i j * f j := mul_comm _ _
  have hd := congrArg Matrix.det hM
  simpa only [Matrix.det_mul, Matrix.det_diagonal] using hd

/-- The exact polynomial-recursion cofactor matrix for a family
of rational two-torsion monomials. -/
def polynomialMatrix (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ) :
    Matrix (Fin r) (Fin r) (Polynomial K) :=
  Matrix.of (fun i j =>
    numeratorRec (denominator c) (logNumerator c (e j).1 (e j).2) i.val)

/-- The concrete rational monomial associated to exponent pair e. -/
def monomial (c : K) (e : ℕ × ℕ) : RatFunc K :=
  RatFunc.mk
    (X ^ e.1 * (X + C c) ^ e.2)
    ((X - 1 : Polynomial K) ^ e.2)

/-- **All sizes:** the determinant of the genuine higher-derivative
Wronskian, multiplied by its common row-denominator, equals the
embedded polynomial cofactor determinant times the product of all
rational monomial columns. No new algebraic assumption is introduced.
-/
theorem concrete_determinant_factorization (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) :
    (∏ i : Fin r,
      (algebraMap (Polynomial K) (RatFunc K) (denominator c)) ^ i.val) *
        (derivativeMatrix r (fun j => monomial c (e j))).det =
    (Matrix.of (fun i j =>
      algebraMap (Polynomial K) (RatFunc K)
        (polynomialMatrix r c e i j))).det *
      (∏ j : Fin r, monomial c (e j)) := by
  apply det_scaled r (fun j => monomial c (e j))
    (algebraMap (Polynomial K) (RatFunc K) (denominator c))
    (Matrix.of (fun i j =>
      algebraMap (Polynomial K) (RatFunc K)
        (polynomialMatrix r c e i j)))
  intro i j
  exact concrete_higher_derivative c (e j).1 (e j).2 i.val

end ICFallDegrees.GrowthRatFuncWronskianNormalization
