import ICFallDegrees.GrowthRationalLogNumerator
import Mathlib.Tactic

/-!
# N8: a cleared rational logarithmic-derivative identity

This is an **unconditional polynomial identity** needed before defining
and identifying the genuine derivative on RatFunc K.  It proves the
quotient-rule numerator for the two-torsion rational coordinates:
  f = X^u * (X+c)^v / (X-1)^v,
  D = X*(X+c)*(X-1),
  Q = u*(X+c)*(X-1) + v*X*(X-1) - v*X*(X+c).

The proof does not appeal to any proposed Wronskian nonvanishing,
SpreadBound, or CurveModel.  The equality D*(N'H-NH')=Q*N*H
holds even when some exponent is zero; it can be transported through
the injection K[X] -> K(X) and used with an actual quotient derivative.
-/

noncomputable section

namespace ICFallDegrees.GrowthPolynomialQuotientBridge

open Polynomial
open ICFallDegrees.GrowthRationalLogNumerator

variable {K : Type*} [Field K]

/-- Euler's identity for the formal derivative of a polynomial power.
The n=0 case is essential (there is no informal division by n). -/
theorem mul_derivative_pow (P : Polynomial K) (n : ℕ) :
    P * (P ^ n).derivative = C (n : K) * P ^ n * P.derivative := by
  cases n with
  | zero =>
      simp
  | succ n =>
      rw [Polynomial.derivative_pow]
      simp only [Nat.succ_eq_add_one, Nat.add_sub_cancel]
      rw [pow_succ]
      ring

/-- The denominator-cleared quotient rule for three polynomials
with derivative equal to one.  The polynomial factors need not be
nonzero for this identity (nonzero denominators are needed only when
passing to RatFunc). -/
theorem cleared_quotient_rule
    (A B T : Polynomial K) (u v : ℕ)
    (hA : A.derivative = 1) (hB : B.derivative = 1)
    (hT : T.derivative = 1) :
    (A * B * T) *
      ((A ^ u * B ^ v).derivative * T ^ v -
        (A ^ u * B ^ v) * (T ^ v).derivative) =
      (C (u : K) * B * T + C (v : K) * A * T -
        C (v : K) * A * B) * (A ^ u * B ^ v) * T ^ v := by
  have hAu : A * (A ^ u).derivative = C (u : K) * A ^ u := by
    rw [mul_derivative_pow, hA, mul_one]
  have hBv : B * (B ^ v).derivative = C (v : K) * B ^ v := by
    rw [mul_derivative_pow, hB, mul_one]
  have hTv : T * (T ^ v).derivative = C (v : K) * T ^ v := by
    rw [mul_derivative_pow, hT, mul_one]
  calc
    (A * B * T) *
        ((A ^ u * B ^ v).derivative * T ^ v -
          (A ^ u * B ^ v) * (T ^ v).derivative) =
      (B * T * B ^ v * T ^ v) * (A * (A ^ u).derivative) +
        (A * T * A ^ u * T ^ v) * (B * (B ^ v).derivative) -
        (A * B * A ^ u * B ^ v) * (T * (T ^ v).derivative) := by
          rw [Polynomial.derivative_mul]
          ring
    _ = (C (u : K) * B * T + C (v : K) * A * T -
          C (v : K) * A * B) * (A ^ u * B ^ v) * T ^ v := by
          rw [hAu, hBv, hTv]
          ring

/-- Exact polynomial identity for X^u ((X+c)/(X-1))^v.
Unlike numerical regression, it is uniform in the coefficient c,
both exponents, and the coefficient field. -/
theorem concrete_cleared_log_derivative (c : K) (u v : ℕ) :
    denominator c *
      ((X ^ u * (X + C c) ^ v).derivative * (X - 1) ^ v -
        (X ^ u * (X + C c) ^ v) * (((X - 1 : Polynomial K) ^ v).derivative)) =
      logNumerator c u v * (X ^ u * (X + C c) ^ v) * (X - 1) ^ v := by
  have hX : (X : Polynomial K).derivative = 1 := by simp
  have hplus : (X + C c : Polynomial K).derivative = 1 := by simp
  have hminus : (X - 1 : Polynomial K).derivative = 1 := by simp
  simpa only [denominator, logNumerator] using
    (cleared_quotient_rule (X : Polynomial K) (X + C c) (X - 1)
      u v hX hplus hminus)

end ICFallDegrees.GrowthPolynomialQuotientBridge
