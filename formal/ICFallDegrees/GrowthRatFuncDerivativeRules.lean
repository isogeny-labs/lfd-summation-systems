import ICFallDegrees.GrowthRatFuncLogDerivative
import Mathlib.Tactic

/-!
# N8: actual rational-field inverse and quotient differentiation

Two general quotient-rule consequences of the genuine derivation from
GrowthComplexRatFuncDerivative. In contrast to the polynomial-level
cleared identity, these statements identify true rational derivatives.

This file does not assume Wronskian nonvanishing, cofactor ratios,
pole valuations, or GrowthSpread.SpreadBound.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncDerivativeRules

open ICFallDegrees.GrowthComplexRatFuncDerivative

variable {K : Type*} [Field K]

/-- The inverse rule for the actual derivative of K(X).
The nonzero denominator hypothesis is essential in the field. -/
theorem rationalDerivative_inv (x : RatFunc K) (hx : x ≠ 0) :
    rationalDerivative x⁻¹ = -(rationalDerivative x) / x ^ 2 := by
  have hmul : x * rationalDerivative x⁻¹ +
        x⁻¹ * rationalDerivative x = 0 := by
    calc
      x * rationalDerivative x⁻¹ + x⁻¹ * rationalDerivative x =
          rationalDerivative (x * x⁻¹) := (rationalDerivative_mul x x⁻¹).symm
      _ = 0 := by rw [mul_inv_cancel₀ hx, rationalDerivative_one]
  apply (eq_div_iff (pow_ne_zero 2 hx)).2
  calc
    rationalDerivative x⁻¹ * x ^ 2 =
        x * (x * rationalDerivative x⁻¹ +
          x⁻¹ * rationalDerivative x) - rationalDerivative x := by
            field_simp [hx]
            ring
    _ = -(rationalDerivative x) := by rw [hmul]; ring

/-- True quotient rule for arbitrary elements of the rational function
field. This does not rely on identifying the pair with polynomial
representatives, only on the already constructed derivation. -/
theorem rationalDerivative_div (x y : RatFunc K) (hy : y ≠ 0) :
    rationalDerivative (x / y) =
      (rationalDerivative x * y - x * rationalDerivative y) / y ^ 2 := by
  rw [div_eq_mul_inv, rationalDerivative_mul, rationalDerivative_inv y hy]
  field_simp [hy]
  ring

/-- Formal derivative of every natural power in the genuine rational
function field. This specializes Mathlib's general derivation power
Leibniz rule and supplies the needed denominator-power differentiation
step for the future higher-derivative numerator recurrence. -/
theorem rationalDerivative_pow (x : RatFunc K) (n : ℕ) :
    rationalDerivative (x ^ n) =
      (n : RatFunc K) * x ^ (n - 1) * rationalDerivative x := by
  simpa only [rationalDerivativeDerivation_apply, nsmul_eq_mul,
    smul_eq_mul, mul_assoc] using
    ((rationalDerivativeDerivation (K := K)).leibniz_pow x n)

end ICFallDegrees.GrowthRatFuncDerivativeRules
