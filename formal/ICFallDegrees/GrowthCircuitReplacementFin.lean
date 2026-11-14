import ICFallDegrees.GrowthCircuitReplacementSum
import Mathlib.Tactic

/-!
# Reindex a minimal circuit replacement to an actual Fin family

The deleted circuit is enumerated by the finite equivalence, preserving
both the replacement equality and the nonzero coefficients. This is
precisely the statement type accepted by the actual Wronskian
replacement theorem. No circuit spread axiom is assumed.
-/

noncomputable section

namespace ICFallDegrees.GrowthCircuitReplacementFin

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthCircuitReplacementSum

variable {K : Type*} [Field K]

theorem erase_sum_equivFin
    {ι F : Type*} [DecidableEq ι] [AddCommMonoid F]
    (J : Finset ι) (x : ι) (g : ι → F) :
    (J.erase x).sum g =
      ∑ i : Fin (J.erase x).card,
        g (((Finset.equivFin (J.erase x)).symm i).1) := by
  classical
  let e := (Finset.equivFin (J.erase x)).symm
  have hsum :
      (∑ i : Fin (J.erase x).card, g ((e i).1)) =
        (∑ y : (J.erase x : Finset ι), g y.1) :=
    Equiv.sum_comp e (fun y : (J.erase x : Finset ι) => g y.1)
  calc
    (J.erase x).sum g =
        ∑ y : (J.erase x : Finset ι), g y.1 := by
      simpa only using (Finset.sum_attach (J.erase x) g).symm
    _ = ∑ i : Fin (J.erase x).card,
          g (((Finset.equivFin (J.erase x)).symm i).1) := by
      simpa only [e] using hsum.symm

theorem exists_replacement_fin
    {F ι : Type*} [AddCommGroup F] [Module K F] [DecidableEq ι]
    (v : ι → F) (J : Finset ι) (hcirc : Circuit K v J)
    (x : ι) (hx : x ∈ J) :
    ∃ coeff : Fin (J.erase x).card → K,
      (∀ i, coeff i ≠ 0) ∧
      v x = ∑ i : Fin (J.erase x).card,
        coeff i • v (((Finset.equivFin (J.erase x)).symm i).1) := by
  classical
  obtain ⟨coeff, hall, hsum⟩ :=
    exists_replacement_coeff v J hcirc x hx
  refine ⟨fun i => coeff (((Finset.equivFin (J.erase x)).symm i).1), ?_, ?_⟩
  · intro i
    exact hall _ (((Finset.equivFin (J.erase x)).symm i).2)
  · calc
      v x = (J.erase x).sum (fun y => coeff y • v y) := hsum
      _ = ∑ i : Fin (J.erase x).card,
            coeff (((Finset.equivFin (J.erase x)).symm i).1) •
              v (((Finset.equivFin (J.erase x)).symm i).1) :=
        erase_sum_equivFin J x (fun y => coeff y • v y)

end ICFallDegrees.GrowthCircuitReplacementFin
