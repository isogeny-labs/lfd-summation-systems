import ICFallDegrees.GrowthCharZeroDeterminant
import Mathlib.Tactic

/-!
# A verified characteristic-two identity for E6's integral unit resultant

This is the Frobenius residue of the explicit eliminant
F_N(X)=(X+8)^N-(X-1)^N, N=2^k. It is an ingredient for the
resultant oddness / finite-bad-prime unit-ideal argument.

This module does not yet prove resultant naturality, an integral Bezout
identity, or the prime-field N8 last-fall-degree theorem.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeIntegralUnit

open Polynomial

def eliminantModTwo (k : ℕ) : Polynomial (ZMod 2) :=
  (X + C (8 : ZMod 2)) ^ (2 ^ k) -
    (X - C (1 : ZMod 2)) ^ (2 ^ k)

/-- In any characteristic-two polynomial ring, a²-b²=(a-b)². -/
theorem square_diff_char_two (a b : Polynomial (ZMod 2)) :
    a ^ 2 - b ^ 2 = (a - b) ^ 2 := by
  have h2scalar : (2 : ZMod 2) = 0 := by decide
  have h2 : (2 : Polynomial (ZMod 2)) = 0 := by
    calc
      (2 : Polynomial (ZMod 2)) = C (2 : ZMod 2) := by
        exact (map_ofNat (C : ZMod 2 →+* Polynomial (ZMod 2)) 2).symm
      _ = 0 := by simp [h2scalar]
  linear_combination (a * b - b ^ 2) * h2

/-- The E6 eliminant is 1 modulo two, at every power-of-two depth.
This is stronger than merely proving one sample resultant is nonzero. -/
theorem eliminant_mod_two_eq_one (k : ℕ) :
    eliminantModTwo k = 1 := by
  induction k with
  | zero =>
      have h8 : (8 : ZMod 2) = 0 := by decide
      simp [eliminantModTwo, h8]
  | succ k ih =>
      change (X + C (8 : ZMod 2)) ^ (2 ^ (k + 1)) -
        (X - C (1 : ZMod 2)) ^ (2 ^ (k + 1)) = 1
      rw [show 2 ^ (k + 1) = 2 ^ k * 2 by rw [pow_succ]]
      rw [pow_mul, pow_mul, square_diff_char_two]
      have hsquare : ((X + C (8 : ZMod 2)) ^ (2 ^ k) -
          (X - C (1 : ZMod 2)) ^ (2 ^ k)) ^ 2 = 1 ^ 2 := by
        exact congrArg (fun f : Polynomial (ZMod 2) => f ^ 2) ih
      exact hsquare.trans (by norm_num)

#print axioms square_diff_char_two
#print axioms eliminant_mod_two_eq_one

end ICFallDegrees.GrowthPrimeIntegralUnit
