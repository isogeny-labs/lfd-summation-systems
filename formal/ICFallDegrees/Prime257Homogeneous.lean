import ICFallDegrees.Prime257System
import ICFallDegrees.PrimeM2Homogeneous

noncomputable section
namespace ICFallDegrees.Prime257System
open MvPolynomial HomogeneousFirstFall

/-- Index the exact previously verified curve-derived input list. -/
def indexedInputs (i : Fin 7) : P := (chains++[g]).get i

theorem parametric_input_identity : PrimeM2Homogeneous.inputs (105 : K) 196 42=indexedInputs := by
  funext i
  fin_cases i <;> rfl

theorem indexed_inputs_range : Set.range indexedInputs=normalizedInputs := by
  ext p
  change (∃ i : Fin 7, (chains++[g]).get i=p) ↔ p∈chains++[g]
  constructor
  · rintro ⟨i,rfl⟩
    exact List.get_mem _ _
  · intro member
    obtain ⟨i,equal⟩ := List.mem_iff_get.mp member
    exact ⟨i,equal⟩

/-- Exact ordinary-polynomial homogeneous-syzygy FFD, modulo the full Koszul
submodule, for the same actual normalized curve system as the LFD proof. -/
theorem homogeneous_first_fall_three : IsKoszulFirstFall indexedInputs 3 := by
  rw [←parametric_input_identity]
  exact PrimeM2Homogeneous.derivative_zero_ffd_three (105 : K) 196 42 derivative_zero

end ICFallDegrees.Prime257System
