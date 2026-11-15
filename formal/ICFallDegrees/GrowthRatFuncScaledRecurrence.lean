import ICFallDegrees.GrowthRatFuncIteratedLinear
import ICFallDegrees.GrowthRationalWronskianRec
import Mathlib.Tactic

/-!
# N8: genuine iterated rational derivatives satisfy a scaled numerator recurrence

This proves the recurrence for the *actual* derivative on K(X), not a
polynomial-only or formal symbolic operation. No Wronskian nonvanishing,
cofactor identity or circuit-spread hypothesis is assumed.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncScaledRecurrence

open ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthRatFuncDerivativeRules
open ICFallDegrees.GrowthRatFuncIteratedLinear

variable {K : Type*} [Field K]

/-- For elements of the rational-function field, the scaled numerator
recurrence associated to a logarithmic derivative equation f' d = f q. -/
def rationalNumeratorRec (d q : RatFunc K) : ℕ → RatFunc K
  | 0 => 1
  | n + 1 =>
      d * rationalDerivative (rationalNumeratorRec d q n) +
        (q - (n : RatFunc K) * rationalDerivative d) *
          rationalNumeratorRec d q n

@[simp]
theorem rationalNumeratorRec_zero (d q : RatFunc K) :
    rationalNumeratorRec d q 0 = 1 := rfl

@[simp]
theorem rationalNumeratorRec_succ (d q : RatFunc K) (n : ℕ) :
    rationalNumeratorRec d q (n + 1) =
      d * rationalDerivative (rationalNumeratorRec d q n) +
        (q - (n : RatFunc K) * rationalDerivative d) *
          rationalNumeratorRec d q n := rfl

/-- One exact step of the scaled higher derivative recurrence, in the
genuine rational function field. The equation f' d = f q is a *real
differential equation hypothesis*, not a circuit-spread hypothesis. -/
theorem scaled_derivative_step (f d q p dp : RatFunc K) (n : ℕ)
    (hlog : rationalDerivative f * d = f * q)
    (hn : iterDerivative n f * d ^ n = f * p)
    (hpd : rationalDerivative p = dp) :
    iterDerivative (n + 1) f * d ^ (n + 1) =
      f * (d * dp + (q - (n : RatFunc K) * rationalDerivative d) * p) := by
  have hpow : (n : RatFunc K) * d ^ (n - 1) * d =
      (n : RatFunc K) * d ^ n := by
    cases n with
    | zero => norm_num
    | succ k =>
        simp only [Nat.succ_eq_add_one, Nat.add_sub_cancel, pow_succ]
        ring
  have hr :
      rationalDerivative (iterDerivative n f * d ^ n) =
        iterDerivative (n + 1) f * d ^ n +
          iterDerivative n f *
            ((n : RatFunc K) * d ^ (n - 1) * rationalDerivative d) := by
    rw [rationalDerivative_mul, rationalDerivative_pow]
    simp only [iterDerivative_succ]
    ring
  have hscale :
      iterDerivative (n + 1) f * d ^ (n + 1) +
        (iterDerivative n f * d ^ n) * (n : RatFunc K) * rationalDerivative d =
        rationalDerivative (iterDerivative n f * d ^ n) * d := by
    rw [hr, pow_succ]
    calc
      iterDerivative (n + 1) f * (d ^ n * d) +
          (iterDerivative n f * d ^ n) * (n : RatFunc K) * rationalDerivative d =
        iterDerivative (n + 1) f * d ^ n * d +
          iterDerivative n f *
            ((n : RatFunc K) * d ^ (n - 1) * d) * rationalDerivative d := by
              rw [hpow]
              ring
      _ = (iterDerivative (n + 1) f * d ^ n +
             iterDerivative n f *
               ((n : RatFunc K) * d ^ (n - 1) * rationalDerivative d)) * d := by ring
  have hright : rationalDerivative (iterDerivative n f * d ^ n) * d =
      (rationalDerivative f * d) * p + f * d * dp := by
    rw [hn, rationalDerivative_mul, hpd]
    ring
  calc
    iterDerivative (n + 1) f * d ^ (n + 1) =
        ((rationalDerivative f * d) * p + f * d * dp) -
          (iterDerivative n f * d ^ n) * (n : RatFunc K) * rationalDerivative d := by
            rw [← hright, ← hscale]
            ring
    _ = f * (d * dp + (q - (n : RatFunc K) * rationalDerivative d) * p) := by
      rw [hn, hlog]
      ring

/-- Authentic arbitrary-order formula in K(X), valid without assuming
Wronskian nonvanishing, pole orders, or a curve model. -/
theorem iterDerivative_scaled (f d q : RatFunc K)
    (hlog : rationalDerivative f * d = f * q) (n : ℕ) :
    iterDerivative n f * d ^ n = f * rationalNumeratorRec d q n := by
  induction n with
  | zero => simp
  | succ n ih =>
      exact scaled_derivative_step f d q (rationalNumeratorRec d q n)
        (rationalDerivative (rationalNumeratorRec d q n)) n hlog ih rfl

end ICFallDegrees.GrowthRatFuncScaledRecurrence
