import ICFallDegrees.GrowthFixedE6

/-!
# N8: the growth window on a one-parameter family of complex curves

`GrowthFixedE6` proves the window `⌊(k-10)/14⌋ < LFD ≤ k+2` for the single
curve `y² = x³ + 6x - 7` with target abscissa `t = 1`. Its three ingredients
are already parametric in the coefficient `a` of the curve

    E_a : y² = x³ + a x - (1 + a),        (1, 0) ∈ E_a of order two,

namely the rational-function section (`GrowthRatFuncModel.xi a`), the spread
bound (`GrowthActualCircuitSpread.actual_spread_bound`, for `a + 3 ≠ 0`) and
the unit ideal over `ℂ` (`GrowthCharZeroSeparatedFamily`, for `5 < ‖a‖`).
This file composes them for every complex `a` with `5 < ‖a‖`.

No curve, spread bound or unit ideal is assumed. Only characteristic zero is
covered here: the prime-field statements of the development remain specific
to `a = 6`.
-/

noncomputable section
namespace ICFallDegrees.GrowthTwoTorsionFamily

open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthActualCircuitSpread
open ICFallDegrees.GrowthMain
open ICFallDegrees.GrowthCharZeroSeparatedFamily
open ICFallDegrees.GrowthReduction
open ICFallDegrees.PrimeM2General

/-- `a + 3 ≠ 0` on the family: it is the value of `3t² + a` at `t = 1`. -/
theorem add_three_ne_zero (a : ℂ) (ha : 5 < ‖a‖) : a + 3 ≠ 0 := by
  intro h
  have heq : a = -(3 : ℂ) := by linear_combination h
  rw [heq] at ha
  norm_num at ha

/-- The rational-function section is a curve model for every `a` with
`a + 3 ≠ 0`, over the complex numbers. -/
theorem model (a : ℂ) (h3 : a + 3 ≠ 0) :
    GrowthMain.CurveModel ℂ a (-1 - a) (1 : ℂ) 0 (xi a) := by
  obtain ⟨hne, hcurve⟩ := ratFunc_algebraic_data a
  exact ⟨hne, hcurve, actual_spread_bound a h3⟩

/-- No degree at most `⌊(k-10)/14⌋` bounds the last fall degree. -/
theorem lfd_gt_linear (a : ℂ) (ha : 5 < ‖a‖)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∀ e ≤ (k - 10) / 14,
      ¬ HasLFDBound
        (Set.range (PrimeM2General.inputs k depth a (-1 - a) 1)) e := by
  exact GrowthMain.growth_lfd_gt_linear
    (K := ℂ) (F₀ := ℂ) (F := RatFunc ℂ)
    (by norm_num) k depth a (-1 - a) (1 : ℂ)
    (xi a) 0 (model a (add_three_ne_zero a ha)) (Or.inl rfl) hk
    (unit_complex_of_norm_gt k depth a ha)

/-- The window `⌊(k-10)/14⌋ < d ≤ k+2` for every last fall degree `d` of the
system of `E_a` with target abscissa `1`, for every complex `a` with
`5 < ‖a‖` and every depth `k ≥ 24`. -/
theorem lfd_linear_window (a : ℂ) (ha : 5 < ‖a‖)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (d : ℕ)
    (hd : IsLastFallDegree
      (Set.range (PrimeM2General.inputs k depth a (-1 - a) 1)) d) :
    (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  constructor
  · by_contra h
    have hle : d ≤ (k - 10) / 14 := by omega
    exact lfd_gt_linear a ha k depth hk d hle hd.1
  · have hupper :
        HasLFDBound
          (Set.range (PrimeM2General.inputs k depth a (-1 - a) 1))
          (k + 2) :=
      SelectorExtension.lfd_le_depth_add_two_of_two_ne
        (K := ℂ) (by norm_num) k depth
        a (-1 - a) (1 : ℂ) (nonsingular_of_norm_gt a ha)
    by_contra h
    exact hd.2 (k + 2) (by omega) hupper

/-- Existence included: the system of `E_a` has a last fall degree, and it
lies in the window. -/
theorem exists_lfd_linear_window (a : ℂ) (ha : 5 < ‖a‖)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∃ d : ℕ,
      IsLastFallDegree
        (Set.range (PrimeM2General.inputs k depth a (-1 - a) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  classical
  let F : Set (MvPolynomial (Fin (2 * k)) ℂ) :=
    Set.range (PrimeM2General.inputs k depth a (-1 - a) 1)
  have hex : ∃ n : ℕ, HasLFDBound F n := by
    refine ⟨k + 2, ?_⟩
    exact SelectorExtension.lfd_le_depth_add_two_of_two_ne
      (K := ℂ) (by norm_num) k depth
      a (-1 - a) (1 : ℂ) (nonsingular_of_norm_gt a ha)
  let d : ℕ := Nat.find hex
  have hd : IsLastFallDegree F d := by
    constructor
    · exact Nat.find_spec hex
    · intro e he hbound
      exact (Nat.find_min hex he) hbound
  exact ⟨d, hd, lfd_linear_window a ha k depth hk d hd⟩

/-- The fixed curve of `GrowthFixedE6` is the member `a = 6` of the family. -/
theorem fixed_E6_is_member :
    (-1 - (6 : ℂ)) = (-7 : ℂ) ∧ 5 < ‖(6 : ℂ)‖ := by
  constructor <;> norm_num

#print axioms model
#print axioms lfd_gt_linear
#print axioms lfd_linear_window
#print axioms exists_lfd_linear_window

end ICFallDegrees.GrowthTwoTorsionFamily
