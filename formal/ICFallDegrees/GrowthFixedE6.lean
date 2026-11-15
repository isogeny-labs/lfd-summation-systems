import ICFallDegrees.GrowthActualCircuitSpread
import ICFallDegrees.GrowthCharZeroSeparatedFamily
import ICFallDegrees.GrowthMain
import Mathlib.Tactic

/-!
# N8: proposed unconditional fixed complex two-torsion growth theorem

Target fixed nonsingular E: y²=x³+6x−7, t=1. Compose:
- actual RatFunc algebraic model and a genuinely proved circuit SpreadBound;
- the independently verified norm-separated complex unit ideal from PR 170;
- kernel GrowthMain lower bound and SelectorExtension upper bound.

There is no assumed CurveModel or SpreadBound in the theorem signatures.
All claims conditional on final pinned Lean 4.19.0/Mathlib v4.19.0 CI passing.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6

open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthActualCircuitSpread
open ICFallDegrees.GrowthMain
open ICFallDegrees.GrowthCharZeroSeparatedFamily
open ICFallDegrees.GrowthReduction
open ICFallDegrees.PrimeM2General

/-- The rational-function section is an actual CurveModel for E_6,
    without a supplied analytic spread hypothesis. -/
theorem model :
    GrowthMain.CurveModel ℂ (6 : ℂ) (-7 : ℂ) (1 : ℂ) 0 (xi (6 : ℂ)) := by
  obtain ⟨hne, hcurve⟩ := ratFunc_algebraic_data (6 : ℂ)
  refine ⟨hne, ?_, fixed_E6_spread_bound⟩
  simpa only [show (-1 - (6 : ℂ)) = (-7 : ℂ) by norm_num] using hcurve

/-- Actual unit ideal from the independent complex proof, not a premise. -/
theorem unit (k : ℕ) (depth : 2 ≤ k) :
    (1 : MvPolynomial (Fin (2 * k)) ℂ) ∈
      Ideal.span (Set.range
        (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)) := by
  have hnorm : 5 < ‖(6 : ℂ)‖ := by norm_num
  have hun := unit_complex_of_norm_gt k depth (6 : ℂ) hnorm
  simpa only [show (-1 - (6 : ℂ)) = (-7 : ℂ) by norm_num] using hun

/-- No degree at most floor((k-10)/14) can bound the last fall degree. -/
theorem lfd_gt_linear
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ∀ e ≤ (k - 10) / 14,
      ¬ HasLFDBound
        (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)) e := by
  exact GrowthMain.growth_lfd_gt_linear
    (K := ℂ) (F₀ := ℂ) (F := RatFunc ℂ)
    (by norm_num) k depth (6 : ℂ) (-7 : ℂ) (1 : ℂ)
    (xi (6 : ℂ)) 0 model (Or.inl rfl) hk (unit k depth)

/-- Explicit unconditional linear window for any last fall degree:
    floor((k-10)/14) < d ≤ k+2. -/
theorem lfd_linear_window
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (d : ℕ)
    (hd : IsLastFallDegree
      (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1)) d) :
    (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  constructor
  · by_contra h
    have hle : d ≤ (k - 10) / 14 := by omega
    exact lfd_gt_linear k depth hk d hle hd.1
  · have hupper :
        HasLFDBound
          (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1))
          (k + 2) :=
      SelectorExtension.lfd_le_depth_add_two_of_two_ne
        (K := ℂ) (by norm_num) k depth
        (6 : ℂ) (-7 : ℂ) (1 : ℂ) (by norm_num)
    by_contra h
    exact hd.2 (k + 2) (by omega) hupper

end ICFallDegrees.GrowthFixedE6
