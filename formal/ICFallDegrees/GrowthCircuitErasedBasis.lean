import ICFallDegrees.GrowthCircuitInterface
import ICFallDegrees.GrowthCircuitFullSupport
import Mathlib.Tactic

/-!
# N8: enumerating the independent complement of a minimal circuit

Every circuit becomes an independent family after deleting one index.
The finite erased set can be enumerated by `Fin (card - 1)`, and
every other circuit index appears in that enumeration.

This is a genuine general finite-circuit interface, with no analytic
assumption, cofactor ratio, or SpreadBound assumption.
-/

noncomputable section

namespace ICFallDegrees.GrowthCircuitErasedBasis

open ICFallDegrees.GrowthSpread

variable {K : Type*} [Field K]

/-- A minimal circuit minus any one index is genuinely independent
when the remaining finite set is reindexed by Fin of its cardinality. -/
theorem erase_reindex_linearIndependent
    {F ι : Type*} [AddCommGroup F] [Module K F]
    [DecidableEq ι]
    (v : ι → F) (J : Finset ι)
    (hcirc : Circuit K v J) (x : ι) (hx : x ∈ J) :
    LinearIndependent K
      (fun i : Fin (J.erase x).card =>
        v (((Finset.equivFin (J.erase x)).symm i).1)) := by
  classical
  have hind : LinearIndepOn K v ((J.erase x : Finset ι) : Set ι) :=
    hcirc.2 (J.erase x) (Finset.erase_ssubset hx)
  have hind' : LinearIndependent K
      (fun y : (J.erase x : Finset ι) => v y.1) := hind
  have hfin := hind'.comp (Finset.equivFin (J.erase x)).symm
    (Finset.equivFin (J.erase x)).symm.injective
  simpa only [Function.comp_def] using hfin

/-- Every other circuit index has a preimage in the erased-family Fin
enumeration. This also identifies the chosen replacement column index. -/
theorem erase_reindex_contains
    {ι : Type*} [DecidableEq ι] (J : Finset ι) (x y : ι)
    (hx : x ∈ J) (hy : y ∈ J) (hxy : y ≠ x) :
    ∃ i : Fin (J.erase x).card,
      (((Finset.equivFin (J.erase x)).symm i).1) = y := by
  classical
  have hmem : y ∈ J.erase x := Finset.mem_erase.mpr ⟨hxy, hy⟩
  let ysub : (J.erase x : Finset ι) := ⟨y, hmem⟩
  refine ⟨(Finset.equivFin (J.erase x)) ysub, ?_⟩
  simp [ysub]

/-- The independent replacement family has exactly m-1 columns. -/
theorem erase_reindex_card {ι : Type*}
    [DecidableEq ι]
    (J : Finset ι) (x : ι) (hx : x ∈ J) :
    (J.erase x).card = J.card - 1 :=
  Finset.card_erase_of_mem hx

end ICFallDegrees.GrowthCircuitErasedBasis
