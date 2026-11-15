import ICFallDegrees.GrowthCircuitPairBound
import Mathlib.Tactic

/-!
# N8: concrete unconditional circuit spread for rational two-torsion section

This module attempts the final arithmetic transport of the genuine
Wronskian/cofactor estimate on every minimal circuit to the exact
GrowthSpread.SpreadBound predicate at characteristic zero.

No replacement relation, cofactor ratio, analytic estimate or
SpreadBound is assumed. The result is specialized to the rational section
xi a, with a+3 nonzero (for E6, a=6).
-/

noncomputable section
namespace ICFallDegrees.GrowthActualCircuitSpread

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthCircuitPairBound

variable {K : Type*} [Field K] [CharZero K]

theorem sum_twice_fin_eq_budget (r : ℕ) :
    (∑ j : Fin r, 2 * j.val) = r * (r - 1) := by
  calc
    (∑ j : Fin r, 2 * j.val) =
        ∑ j ∈ Finset.range r, 2 * j :=
      Fin.sum_univ_eq_sum_range (fun j : ℕ => 2 * j) r
    _ = (∑ j ∈ Finset.range r, j) * 2 := by
      rw [Finset.sum_mul]
      exact Finset.sum_congr rfl (by intro j hj; omega)
    _ = r * (r - 1) := Finset.sum_range_id_mul_two r

theorem actual_spread_bound
    (a : K) (ha : a + 3 ≠ 0) :
    SpreadBound K (xi a) 0 := by
  intro J hcirc x hx y hy
  refine ⟨0, 0, ?_⟩
  simp only [Nat.cast_zero, zero_mul, sub_zero]
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
      distinct_circuit_pair_le_four_budget a ha J hcirc x y hx hy hxy
    rw [sum_twice_fin_eq_budget, Nat.cast_mul, hr2, hr1] at hbudget
    have hcast1 : ((J.card - 1 : ℕ) : ℤ) = (J.card : ℤ) - 1 := by omega
    have hcast2 : ((J.card - 2 : ℕ) : ℤ) = (J.card : ℤ) - 2 := by omega
    rw [hcast1, hcast2] at hbudget
    have hp : 0 ≤ ((J.card : ℤ) - 1) * ((J.card : ℤ) - 2) :=
      mul_nonneg (by omega) (by omega)
    nlinarith

theorem fixed_E6_spread_bound :
    SpreadBound ℂ (xi (6 : ℂ)) 0 :=
  actual_spread_bound (6 : ℂ) (by norm_num)

end ICFallDegrees.GrowthActualCircuitSpread
