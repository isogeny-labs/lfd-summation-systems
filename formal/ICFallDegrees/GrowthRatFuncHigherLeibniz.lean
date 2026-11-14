import ICFallDegrees.GrowthRatFuncWronskianLinearChange
import Mathlib.Tactic

/-!
# N8: full higher Leibniz rule for genuine rational derivatives

This is the missing all-order product formula for the actual
quotient-rule derivation on K(X), rather than a formal numerator
recurrence. It will supply the lower-triangular row operations in the
common-factor Wronskian identity W(hf_1,...,hf_r)=h^r W(f_1,...,f_r).

The inductive Pascal identity follows the technique of Mathlib's
`Polynomial.iterate_derivative_mul` (Mathlib, Apache 2.0).

No Wronskian nonvanishing, pole valuation, or spread hypothesis is
introduced.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncHigherLeibniz

open Finset
open ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthRatFuncIteratedLinear

variable {K : Type*} [Field K]

/-- The genuine quotient-rule derivation commutes with finite sums. -/
theorem rationalDerivative_sum {ι : Type*} (s : Finset ι)
    (f : ι → RatFunc K) :
    rationalDerivative (∑ i ∈ s, f i) = ∑ i ∈ s, rationalDerivative (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert i s hi ih =>
      simp [Finset.sum_insert hi, rationalDerivative_add, ih]

/-- The genuine derivative commutes with natural-number scaling. -/
theorem rationalDerivative_nsmul (n : ℕ) (f : RatFunc K) :
    rationalDerivative (n • f) = n • rationalDerivative f := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [succ_nsmul, rationalDerivative_add, ih, succ_nsmul]

/-- All-order Leibniz identity for the *actual* rational derivation.
Unlike the special logarithmic-numerator recurrence, this holds for
arbitrary rational functions f and g over every coefficient field. -/
theorem iterDerivative_mul (n : ℕ) (f g : RatFunc K) :
    iterDerivative n (f * g) =
      ∑ k ∈ Finset.range (n + 1),
        n.choose k • (iterDerivative (n - k) f * iterDerivative k g) := by
  induction n with
  | zero => simp [Finset.range]
  | succ n ih =>
      calc
        iterDerivative (n + 1) (f * g) =
            rationalDerivative
              (∑ k ∈ Finset.range (n + 1),
                n.choose k • (iterDerivative (n - k) f * iterDerivative k g)) := by
                  rw [iterDerivative_succ, ih]
        _ = (∑ k ∈ Finset.range (n + 1),
              n.choose k • (iterDerivative (n - k + 1) f * iterDerivative k g)) +
            ∑ k ∈ Finset.range (n + 1),
              n.choose k • (iterDerivative (n - k) f * iterDerivative (k + 1) g) := by
                simp_rw [rationalDerivative_sum, rationalDerivative_nsmul,
                  rationalDerivative_mul, smul_add, Finset.sum_add_distrib,
                  ← iterDerivative_succ]
                conv_lhs => rw [add_comm]
                congr 1 <;> apply Finset.sum_congr rfl <;> intro k hk
                all_goals ring
        _ = (∑ k ∈ Finset.range (n + 1),
                n.choose (k + 1) •
                  (iterDerivative (n - k) f * iterDerivative (k + 1) g)) +
              1 • (iterDerivative (n + 1) f * iterDerivative 0 g) +
            ∑ k ∈ Finset.range (n + 1),
              n.choose k • (iterDerivative (n - k) f * iterDerivative (k + 1) g) := by
                congr 1
                refine (Finset.sum_range_succ' _ _).trans
                  (congrArg₂ (· + ·) ?_ ?_)
                · rw [Finset.sum_range_succ, Nat.choose_succ_self,
                    zero_smul, add_zero]
                  refine Finset.sum_congr rfl fun k hk => ?_
                  rw [Finset.mem_range] at hk
                  congr
                  omega
                · rw [Nat.choose_zero_right, tsub_zero]
        _ = ((∑ k ∈ Finset.range (n + 1),
                n.choose k • (iterDerivative (n - k) f * iterDerivative (k + 1) g)) +
              ∑ k ∈ Finset.range (n + 1),
                n.choose (k + 1) •
                  (iterDerivative (n - k) f * iterDerivative (k + 1) g)) +
              1 • (iterDerivative (n + 1) f * iterDerivative 0 g) := by
                rw [add_comm, add_assoc]
        _ = (∑ i ∈ Finset.range (n + 1),
              (n + 1).choose (i + 1) •
                (iterDerivative (n + 1 - (i + 1)) f * iterDerivative (i + 1) g)) +
              1 • (iterDerivative (n + 1) f * iterDerivative 0 g) := by
                simp_rw [Nat.choose_succ_succ, Nat.succ_sub_succ,
                  add_smul, Finset.sum_add_distrib]
        _ = ∑ k ∈ Finset.range ((n + 1) + 1),
              (n + 1).choose k •
                (iterDerivative (n + 1 - k) f * iterDerivative k g) := by
                  rw [Finset.sum_range_succ' _ (n + 1),
                      Nat.choose_zero_right, tsub_zero]

end ICFallDegrees.GrowthRatFuncHigherLeibniz
