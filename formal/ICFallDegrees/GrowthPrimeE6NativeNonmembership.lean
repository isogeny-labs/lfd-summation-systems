import ICFallDegrees.GrowthPrimeE6NativeBoxSpread
import ICFallDegrees.GrowthPrimeBoxNonmembership
import Mathlib.Tactic

/-!
# Actual E6 native prime-field closure nonmembership

The finite-box spread in characteristic larger than the bounded polynomial
degree is proved over the native field itself, without a circuit-matroid
transfer premise. This file composes that result with the existing
BoxSpreadBound -> FEI -> InClosure exclusion chain. The large-prime
condition is essential; this does not claim all admissible primes.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6NativeNonmembership

open MvPolynomial
open ICFallDegrees.GrowthBoxSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GrowthReduction

/-- At sufficiently large prime characteristic, the *actual* fixed E6
bounded exponent box has native spread, without any rational-circuit
specialization hypothesis. -/
theorem e6_native_box_spread_large_prime
    (p k : ℕ) [Fact p.Prime] (hk : 3 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p) :
    BoxSpreadBound (ZMod p) (xi (6 : ZMod p)) p k := by
  have hchar : ringChar (ZMod p) = 0 ∨ 2 * 2 ^ k ≤ ringChar (ZMod p) := by
    right
    rw [ZMod.ringChar_zmod_n]
    omega
  have hpos : 0 < (9 : ℕ) := by norm_num
  have hgt : 9 < p := by
    have hpow : 2 ^ 3 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) hk
    omega
  have hnine : (9 : ZMod p) ≠ 0 := by
    intro heq
    have hdvd : p ∣ 9 := (ZMod.natCast_zmod_eq_zero_iff_dvd 9 p).mp heq
    have hle : p ≤ 9 := Nat.le_of_dvd hpos hdvd
    omega
  have h9 : (8 : ZMod p) + 1 ≠ 0 := by
    convert hnine using 1 <;> ring
  exact GrowthPrimeE6NativeBoxSpread.e6_native_box_spread
    (K := ZMod p) p k hchar h9

/-- No finite-matroid circuit-transfer hypothesis: true native E6
closure nonmembership for every prime above the depth-dependent bound. -/
theorem e6_not_one_of_large_prime
    (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (D : ℕ) (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth
        (6 : ZMod p) (-7 : ZMod p) 1)) D
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) := by
  exact GrowthPrimeBoxNonmembership.e6_not_one_of_box_spread
    p k depth D hD hk hp (e6_native_box_spread_large_prime p k (by omega) hp)

#print axioms e6_native_box_spread_large_prime
#print axioms e6_not_one_of_large_prime

end ICFallDegrees.GrowthPrimeE6NativeNonmembership
