import ICFallDegrees.GrowthPrimeE6BoxCircuitPair
import ICFallDegrees.GrowthBoxSpread
import ICFallDegrees.GrowthActualCircuitSpread
import Mathlib.Tactic

/-!
# N8: native E6 bounded-box spread, without rational circuit transfer

This proves the actual finite-box spread bound for E6 over the native
coefficient field by using *its own* characteristic-bounded rational
Wronskians, polynomial cofactors, and true pole-gap estimates.

Unlike #287, there is no rational-vs-mod-p circuit transfer hypothesis.
The only extra arithmetic guards are characteristic zero or >=2^(k+1)
and (8+1)!=0; native finite primes above 2^(k+1) with p !=3 satisfy
them. This is not a statement about every admissible prime.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6NativeBoxSpread

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthBoxSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthPrimeE6BoxCircuitPair
open ICFallDegrees.GrowthActualCircuitSpread

variable {K : Type*} [Field K]

/-- Genuine finite-box E6 spread is established by native
characteristic-bounded polynomial/rational Wronskians and a cofactor
pair bound, without finite-matroid circuit specialization. -/
theorem e6_native_box_spread
    (p k : ℕ)
    (hchar : ringChar K = 0 ∨ 2 * 2 ^ k ≤ ringChar K)
    (h9 : (8 : K) + 1 ≠ 0) :
    BoxSpreadBound K (xi (6 : K)) p k := by
  intro J hbox hcirc x hx y hy
  refine ⟨0, 0, ?_⟩
  simp only [Nat.cast_zero, zero_mul, mul_zero, sub_zero]
  have hcard : 1 ≤ J.card := Finset.card_pos.mpr ⟨x, hx⟩
  by_cases heq : y = x
  · subst y
    simp only [sub_self, mul_zero, abs_zero, zero_add]
    by_cases hone : J.card = 1
    · simp [hone]
    · have htwo : 2 ≤ J.card := by omega
      have htwoz : (2 : ℤ) ≤ (J.card : ℤ) := by exact_mod_cast htwo
      have hp : 0 ≤ ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2) :=
        mul_nonneg (by omega) (by omega)
      nlinarith
  · have hxy : y ≠ x := heq
    have hyerase : y ∈ J.erase x := Finset.mem_erase.mpr ⟨hxy, hy⟩
    have hrpos : 0 < (J.erase x).card := Finset.card_pos.mpr ⟨y, hyerase⟩
    have hr1 : (J.erase x).card = J.card - 1 :=
      Finset.card_erase_of_mem hx
    have hr2 : (J.erase x).card - 1 = J.card - 2 := by omega
    have htwo : 2 ≤ J.card := by omega
    have htwoz : (2 : ℤ) ≤ (J.card : ℤ) := by exact_mod_cast htwo
    have hbudget :=
      e6_bounded_circuit_pair_le_four_budget
        (2 ^ k) J hbox hcirc x y hx hy hxy hchar h9
    rw [sum_twice_fin_eq_budget, Nat.cast_mul, hr2, hr1] at hbudget
    have hcast1 : ((J.card - 1 : ℕ) : ℤ) = (J.card : ℤ) - 1 := by omega
    have hcast2 : ((J.card - 2 : ℕ) : ℤ) = (J.card : ℤ) - 2 := by omega
    rw [hcast1, hcast2] at hbudget
    have hp : 0 ≤ ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2) :=
      mul_nonneg (by omega) (by omega)
    nlinarith

#print axioms e6_native_box_spread

end ICFallDegrees.GrowthPrimeE6NativeBoxSpread
