import ICFallDegrees.Prime257System
import ICFallDegrees.GeneratedPrimeUpper
import ICFallDegrees.GeneratedPrimeFirstFall
import ICFallDegrees.Prime257Homogeneous

namespace ICFallDegrees.Prime257System

/-- Exact full ordinary-ring LFD for the actual normalized prime curve
system, from the separately checked sparse upper DAG and lower space. -/
theorem exact_lfd_four : IsLastFallDegree normalizedInputs 4 := by
  rw [←input_set_semantics]
  exact GeneratedPrimeUpper.exact_lfd_four

/-- The affine strict-fall convention and the full ordinary-ring last fall
are separated on the same curve-derived, explicitly fixed presentation. -/
theorem affine_first_three_lfd_four :
    IsFirstStrictFall normalizedInputs 3 ∧ IsLastFallDegree normalizedInputs 4 :=
  ⟨GeneratedPrimeFirstFall.actual_curve_first_strict_fall_three,exact_lfd_four⟩

/-- The ordinary homogeneous-syzygy convention is verified independently of
its affine counterpart, and separated from the exact full LFD. -/
theorem homogeneous_ffd_three_lfd_four :
    HomogeneousFirstFall.IsKoszulFirstFall indexedInputs 3 ∧
      IsLastFallDegree normalizedInputs 4 :=
  ⟨homogeneous_first_fall_three,exact_lfd_four⟩

end ICFallDegrees.Prime257System
