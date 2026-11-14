import ICFallDegrees.GrowthFixedE6
import Mathlib.Tactic

/-!
# Fixed E6: existence of the last fall degree and unconditional Theta bounds

The previous fixed E6 theorem bounds every actual last fall degree but quantifies
over a supplied witness `IsLastFallDegree`. This module closes that final
logical interface: the proven selector upper bound gives a nonempty set of
degrees satisfying `HasLFDBound`, and the least such degree exists by
well-ordering of the natural numbers. Its window and genuine asymptotic
constants then follow from `GrowthFixedE6.lfd_linear_window`.

No additional spread, unit-ideal, model, or existence hypothesis.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6

open ICFallDegrees.PrimeM2General
open ICFallDegrees.GrowthMain

/-- For every depth k ≥ 24 the fixed E6 system possesses an actual least
last-fall-degree, satisfying the already kernel-checked two-sided window. -/
theorem exists_lfd_linear_window (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∃ d : ℕ,
      IsLastFallDegree
        (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  classical
  let F : Set (MvPolynomial (Fin (2 * k)) ℂ) :=
    Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)
  have hex : ∃ n : ℕ, HasLFDBound F n := by
    refine ⟨k + 2, ?_⟩
    exact SelectorExtension.lfd_le_depth_add_two_of_two_ne
      (K := ℂ) (by norm_num) k depth
      (6 : ℂ) (-7 : ℂ) (1 : ℂ) (by norm_num)
  let d : ℕ := Nat.find hex
  have hd : IsLastFallDegree F d := by
    constructor
    · exact Nat.find_spec hex
    · intro e he hbound
      exact (Nat.find_min hex he) hbound
  exact ⟨d, hd, lfd_linear_window k depth hk d hd⟩

/-- Explicit unconditional Theta(k) constants, with existence, for the
fixed curve y²=x³+6x-7 and target t=1: for all k≥24,
there is an actual last-fall-degree d with k ≤ 28d and d ≤ 2k.
The sharper floor-window is available in `exists_lfd_linear_window`. -/
theorem exists_lfd_theta (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∃ d : ℕ,
      IsLastFallDegree
        (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)) d ∧
      k ≤ 28 * d ∧ d ≤ 2 * k := by
  obtain ⟨d, hd, hlo, hhi⟩ := exists_lfd_linear_window k depth hk
  refine ⟨d, hd, ?_, ?_⟩ <;> omega

end ICFallDegrees.GrowthFixedE6
