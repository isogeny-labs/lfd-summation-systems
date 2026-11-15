import ICFallDegrees.GrowthComplexRatFuncDerivative
import ICFallDegrees.GrowthPolynomialQuotientBridge
import Mathlib.Tactic

/-!
# N8: genuine RatFunc logarithmic derivative of the two-torsion coordinates

The polynomial identity of GrowthPolynomialQuotientBridge is transported
into the true quotient-rule derivation on RatFunc. No SpreadBound,
Wronskian nonvanishing or assumed differential identity is used.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncLogDerivative

open Polynomial ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthPolynomialQuotientBridge
open ICFallDegrees.GrowthRationalLogNumerator

variable {K : Type*} [Field K]

/-- Transport a denominator-cleared polynomial derivative equality through
the genuine RatFunc quotient derivation. The only hypothesis besides the
polynomial identity is that its rational denominator H is nonzero. -/
theorem derivative_mk_mul_of_cleared
    (N H D Q : Polynomial K)
    (hH : H ≠ 0)
    (hpoly :
      D * (N.derivative * H - N * H.derivative) = Q * N * H) :
    rationalDerivative (RatFunc.mk N H) *
        algebraMap (Polynomial K) (RatFunc K) D =
      RatFunc.mk N H * algebraMap (Polynomial K) (RatFunc K) Q := by
  let φ : Polynomial K →+* RatFunc K := algebraMap (Polynomial K) (RatFunc K)
  have hφH : φ H ≠ 0 := RatFunc.algebraMap_ne_zero hH
  have hmap : φ D * (φ N.derivative * φ H - φ N * φ H.derivative) =
      φ Q * φ N * φ H := by
    have h := congrArg φ hpoly
    simpa only [map_mul, map_sub] using h
  have hmap' :
      (φ N.derivative * φ H - φ N * φ H.derivative) * φ D =
        ((φ N / φ H) * φ Q) * (φ H) ^ 2 := by
    calc
      _ = φ D * (φ N.derivative * φ H - φ N * φ H.derivative) := by ring
      _ = φ Q * φ N * φ H := hmap
      _ = ((φ N / φ H) * φ Q) * (φ H) ^ 2 := by
        field_simp [hφH]
        ring
  rw [rationalDerivative_mk]
  simp only [RatFunc.mk_eq_div, map_mul, map_sub, map_pow]
  calc
    ((φ N.derivative * φ H - φ N * φ H.derivative) / (φ H) ^ 2) * φ D =
        ((φ N.derivative * φ H - φ N * φ H.derivative) * φ D) / (φ H) ^ 2 := by ring
    _ = (((φ N / φ H) * φ Q) * (φ H) ^ 2) / (φ H) ^ 2 := by rw [hmap']
    _ = (φ N / φ H) * φ Q := by
      field_simp [hφH]
      ring

/-- The true rational derivative identity for
X^u * ((X+c)/(X-1))^v, written without dividing by the monomial
(which may simplify a future extension to positive characteristic). -/
theorem concrete_rational_derivative (c : K) (u v : ℕ) :
    rationalDerivative
        (RatFunc.mk
          (X ^ u * (X + C c) ^ v)
          ((X - 1 : Polynomial K) ^ v)) *
      algebraMap (Polynomial K) (RatFunc K) (denominator c) =
    RatFunc.mk
        (X ^ u * (X + C c) ^ v)
        ((X - 1 : Polynomial K) ^ v) *
      algebraMap (Polynomial K) (RatFunc K) (logNumerator c u v) := by
  have hT : (X - 1 : Polynomial K) ≠ 0 := by
    simpa only [Polynomial.C_1] using (Polynomial.X_sub_C_ne_zero (1 : K))
  exact derivative_mk_mul_of_cleared
    (X ^ u * (X + C c) ^ v) ((X - 1 : Polynomial K) ^ v)
    (denominator c) (logNumerator c u v)
    (pow_ne_zero _ hT)
    (concrete_cleared_log_derivative c u v)

end ICFallDegrees.GrowthRatFuncLogDerivative
