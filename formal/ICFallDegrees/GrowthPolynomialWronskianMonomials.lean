import ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Data.Nat.Factorial.NatCast
import Mathlib.Tactic

/-!
# N8: characteristic-zero polynomial monomial Wronskians

Proves, in the pinned kernel, that the genuine polynomial Wronskian
of X^0,...,X^(r-1) is nonzero for every r in characteristic zero.
The derivative matrix is upper triangular with diagonal j!, giving an
exact determinant, including r=0. This is the indispensable base
family for a degree-echelon proof of arbitrary polynomial Wronskian
nonvanishing. We do NOT assume a general Wronskian theorem, cofactor
identities, valuations, or SpreadBound here.
-/

noncomputable section
namespace ICFallDegrees.GrowthPolynomialWronskianMonomials

open Polynomial
open ICFallDegrees.GrowthRatFuncPolynomialWronskianMap

variable {K : Type*} [Field K]

/-- Every entry below the diagonal is zero: the i-th derivative of X^j
vanishes whenever j<i. -/
theorem derivativeMatrix_monomials_upper (r : ℕ) :
    (polynomialDerivativeMatrix r
      (fun j : Fin r => (X : Polynomial K) ^ j.val)).BlockTriangular id := by
  intro i j hij
  change ((Polynomial.derivative^[i.val]) ((X : Polynomial K) ^ j.val)) = 0
  rw [Polynomial.iterate_derivative_X_pow_eq_natCast_mul]
  have hz : (j.val.descFactorial i.val : ℕ) = 0 :=
    Nat.descFactorial_eq_zero_iff_lt.mpr hij
  simp [hz]

/-- Diagonal entry j is precisely the constant polynomial j!. -/
theorem derivativeMatrix_monomials_diag (r : ℕ) (j : Fin r) :
    polynomialDerivativeMatrix r
      (fun k : Fin r => (X : Polynomial K) ^ k.val) j j =
      C (Nat.factorial j.val : K) := by
  change ((Polynomial.derivative^[j.val]) ((X : Polynomial K) ^ j.val)) =
    C (Nat.factorial j.val : K)
  rw [Polynomial.iterate_derivative_X_pow_eq_natCast_mul]
  simp [Nat.descFactorial_self]

/-- Exact determinant of the full r-column polynomial Wronskian. -/
theorem monomial_wronskian_det (r : ℕ) :
    (polynomialDerivativeMatrix r
      (fun j : Fin r => (X : Polynomial K) ^ j.val)).det =
      ∏ j : Fin r, C (Nat.factorial j.val : K) := by
  classical
  rw [Matrix.det_of_upperTriangular (derivativeMatrix_monomials_upper r)]
  exact Finset.prod_congr rfl (fun j _ => derivativeMatrix_monomials_diag r j)

/-- Characteristic zero makes every factorial diagonal coefficient a unit,
so the polynomial Wronskian of the standard monomial basis never vanishes. -/
theorem monomial_wronskian_ne_zero [CharZero K] (r : ℕ) :
    (polynomialDerivativeMatrix r
      (fun j : Fin r => (X : Polynomial K) ^ j.val)).det ≠ 0 := by
  rw [monomial_wronskian_det]
  apply Finset.prod_ne_zero_iff.mpr
  intro j _
  have hj : (Nat.factorial j.val : K) ≠ 0 := by
    exact_mod_cast Nat.factorial_ne_zero j.val
  exact (Polynomial.C_ne_zero).mpr hj

end ICFallDegrees.GrowthPolynomialWronskianMonomials
