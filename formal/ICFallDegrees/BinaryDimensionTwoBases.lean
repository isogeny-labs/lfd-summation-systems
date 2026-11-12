import ICFallDegrees.BinaryDimensionTwoCurves
import Mathlib.LinearAlgebra.Dimension.Constructions

noncomputable section
namespace ICFallDegrees.BinaryDimensionTwoCurves
open BinaryExtensionFields

theorem ordinary_base_independent : LinearIndependent k ![(1 : K8),alpha8] := by
  let embed : Fin 2 → Fin 3 := fun i => ⟨i.val,by omega⟩
  have injective : Function.Injective embed := by
    intro i j equal
    exact Fin.ext (congrArg (fun a : Fin 3 => a.val) equal)
  have linear := basis8.linearIndependent.comp embed injective
  have same : (fun i => basis8 (embed i))= ![(1 : K8),alpha8] := by
    funext i
    fin_cases i <;> simp [embed,basis8_eq_pow]
  change LinearIndependent k (fun i => basis8 (embed i)) at linear
  rw [same] at linear
  exact linear

theorem koblitz_base_independent : LinearIndependent k ![(1 : K16),alpha16] := by
  let embed : Fin 2 → Fin 4 := fun i => ⟨i.val,by omega⟩
  have injective : Function.Injective embed := by
    intro i j equal
    exact Fin.ext (congrArg (fun a : Fin 4 => a.val) equal)
  have linear := basis16.linearIndependent.comp embed injective
  have same : (fun i => basis16 (embed i))= ![(1 : K16),alpha16] := by
    funext i
    fin_cases i <;> simp [embed,basis16_eq_pow]
  change LinearIndependent k (fun i => basis16 (embed i)) at linear
  rw [same] at linear
  exact linear

theorem ordinary_factor_base_dimension : Module.finrank k factorBase8=2 := by
  have same : Set.range (![(1 : K8),alpha8] : Fin 2 → K8)={1,alpha8} := by
    ext a
    simp [Fin.exists_fin_succ,or_comm]
  rw [factorBase8,←same]
  simpa using finrank_span_eq_card ordinary_base_independent

theorem koblitz_factor_base_dimension : Module.finrank k factorBase16=2 := by
  have same : Set.range (![(1 : K16),alpha16] : Fin 2 → K16)={1,alpha16} := by
    ext a
    simp [Fin.exists_fin_succ,or_comm]
  rw [factorBase16,←same]
  simpa using finrank_span_eq_card koblitz_base_independent

end ICFallDegrees.BinaryDimensionTwoCurves
