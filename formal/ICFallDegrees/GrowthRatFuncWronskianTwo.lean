import ICFallDegrees.GrowthRatFuncIteratedLinear
import Mathlib.Tactic

/-!
# N8: genuine two-column rational-function Wronskian criterion

Over any characteristic-zero field K, the determinant of the first two
genuine derivative rows of f and g in K(X) vanishes exactly when g
is a constant multiple of nonzero f.

This is a proof for the actual differential field, *not* an assumption
that independent arbitrary collections have nonzero Wronskians.  The
general multi-column induction and cofactor/pole bridges remain open.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncWronskianTwo

open ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthRatFuncDerivativeRules

variable {K : Type*} [Field K] [CharZero K]

/-- Actual 2-by-2 Wronskian of two rational functions in K(X). -/
def wronskianTwo (f g : RatFunc K) : RatFunc K :=
  f * rationalDerivative g - g * rationalDerivative f

/-- Characteristic-zero two-column Wronskian criterion, derived
from the *kernel-proved* rational quotient rule and derivative
constant-field theorem. -/
theorem wronskianTwo_eq_zero_iff (f g : RatFunc K) (hf : f ≠ 0) :
    wronskianTwo f g = 0 ↔
      ∃ c : K, g = RatFunc.C c * f := by
  constructor
  · intro h
    have hn : rationalDerivative g * f - g * rationalDerivative f = 0 := by
      calc
        _ = wronskianTwo f g := by
          unfold wronskianTwo
          ring
        _ = 0 := h
    have hr : rationalDerivative (g / f) = 0 := by
      rw [rationalDerivative_div g f hf, hn]
      simp
    obtain ⟨c, hc⟩ := (rationalDerivative_eq_zero_iff (g / f)).mp hr
    exact ⟨c, (div_eq_iff hf).mp hc⟩
  · rintro ⟨c, rfl⟩
    have hc : rationalDerivative (RatFunc.C c) = 0 := by
      rw [← RatFunc.algebraMap_eq_C]
      exact rationalDerivative_algebraMap_base c
    unfold wronskianTwo
    rw [rationalDerivative_mul, hc]
    ring

/-- For any nonzero first column, the two-column Wronskian
is nonzero precisely for constant-linearly independent columns. -/
theorem wronskianTwo_ne_zero_iff (f g : RatFunc K) (hf : f ≠ 0) :
    wronskianTwo f g ≠ 0 ↔
      ∀ c : K, g ≠ RatFunc.C c * f := by
  constructor
  · intro h c hc
    exact h ((wronskianTwo_eq_zero_iff f g hf).mpr ⟨c, hc⟩)
  · intro h hz
    obtain ⟨c, hc⟩ := (wronskianTwo_eq_zero_iff f g hf).mp hz
    exact h c hc

end ICFallDegrees.GrowthRatFuncWronskianTwo
