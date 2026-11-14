import ICFallDegrees.GrowthPolynomialWronskianBasis
import ICFallDegrees.GrowthRatFuncDenominatorIndependence
import ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
import ICFallDegrees.GrowthRatFuncWronskianCommonScalar
import Mathlib.Tactic

/-!
# N8 bounded rational Wronskian nonvanishing at large characteristic

This theorem uses the real characteristic-bounded polynomial degree-echelon
Vandermonde argument and genuine RatFunc differentiation. It avoids
asserting that ordinary derivatives detect independence over F_p
unconditionally: the polynomial representative degree < ringChar(K)
is an explicit essential premise.

For N8 bounded exponent-box monomials, common clearing by (X-1)^(N-1)
produces polynomial representatives of degree at most 2N-2. The
application to each circuit and resulting spread theorem are not
established in this module.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeBoundedRationalWronskian

open Polynomial
open ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncWronskianCommonScalar

variable {K : Type*} [Field K]

/-- A linearly independent polynomial family with degrees strictly
below the positive characteristic, or in characteristic zero, has a
genuine nonzero arbitrary-column polynomial Wronskian. -/
theorem polynomial_wronskian_ne_zero_of_independent_degree_bound
    (r k : ℕ) (p : Fin r → Polynomial K)
    (hind : LinearIndependent K p)
    (hdeg : ∀ j, (p j).natDegree < k)
    (hchar : ringChar K = 0 ∨ k ≤ ringChar K) :
    (polynomialDerivativeMatrix r p).det ≠ 0 := by
  classical
  let B : Submodule K (Polynomial K) := Submodule.span K (Set.range p)
  let bas : Basis (Fin r) K B := Basis.span hind
  have hbas (j : Fin r) : ((bas j : B) : Polynomial K) = p j := by
    exact Basis.span_apply hind j
  have hdeg' : ∀ j : Fin r, ((bas j : B) : Polynomial K).natDegree < k := by
    intro j
    rw [hbas j]
    exact hdeg j
  have hwr :=
    ICFallDegrees.GrowthPolynomialWronskianBasis.polynomial_wronskian_ne_zero_of_basis_char_bound
      bas hdeg' hchar
  change (polynomialDerivativeMatrix r
    (fun j => ((bas j : B) : Polynomial K))).det ≠ 0 at hwr
  have heq : (fun j : Fin r => ((bas j : B) : Polynomial K)) = p := by
    funext j
    exact hbas j
  rw [heq] at hwr
  exact hwr

/-- Genuine positive-characteristic rational Wronskian nonvanishing
from a *proved* polynomial representation with degree below the
characteristic and constant-linear independence. Neither Wronskian
nonvanishing nor SpreadBound is assumed. -/
theorem rational_wronskian_ne_zero_of_bounded_polynomial_representatives
    (r k : ℕ) (f : Fin r → RatFunc K)
    (H : Polynomial K) (p : Fin r → Polynomial K)
    (hH : H ≠ 0)
    (hrep : ∀ j,
      algebraMap (Polynomial K) (RatFunc K) H * f j =
      algebraMap (Polynomial K) (RatFunc K) (p j))
    (hind : LinearIndependent K p)
    (hdeg : ∀ j, (p j).natDegree < k)
    (hchar : ringChar K = 0 ∨ k ≤ ringChar K) :
    (derivativeMatrix r f).det ≠ 0 := by
  have hpoly : (polynomialDerivativeMatrix r p).det ≠ 0 :=
    polynomial_wronskian_ne_zero_of_independent_degree_bound r k p hind hdeg hchar
  have hrat :
      (derivativeMatrix r
        (fun j => algebraMap (Polynomial K) (RatFunc K) (p j))).det ≠ 0 :=
    (determinant_polynomialMap_ne_zero_iff r p).2 hpoly
  let h : RatFunc K := algebraMap (Polynomial K) (RatFunc K) H
  have hh : h ≠ 0 := RatFunc.algebraMap_ne_zero hH
  have hscaled :
      (derivativeMatrix r (fun j => h * f j)).det ≠ 0 := by
    have heq : (fun j : Fin r => h * f j) =
        (fun j => algebraMap (Polynomial K) (RatFunc K) (p j)) := by
      funext j
      exact hrep j
    rw [heq]
    exact hrat
  exact (determinant_mul_scalar_ne_zero_iff r h f hh).mp hscaled

#print axioms polynomial_wronskian_ne_zero_of_independent_degree_bound
#print axioms rational_wronskian_ne_zero_of_bounded_polynomial_representatives

end ICFallDegrees.GrowthPrimeBoundedRationalWronskian
