import ICFallDegrees.GrowthCircuitErasedBasis
import ICFallDegrees.GrowthCircuitFullSupport
import Mathlib.Tactic

/-!
# Genuine circuit relation, solved for an omitted column

A minimal circuit has a full-support relation. Dividing by any chosen
nonzero coefficient provides a genuine replacement equation for that
monomial, with every retained coefficient nonzero. This is valid for
any vector space over a field; no Wronskian or SpreadBound is assumed.
-/

noncomputable section

namespace ICFallDegrees.GrowthCircuitReplacementSum

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthCircuitFullSupport

variable {K : Type*} [Field K]

theorem exists_replacement_coeff
    {F ι : Type*} [AddCommGroup F] [Module K F] [DecidableEq ι]
    (v : ι → F) (J : Finset ι) (hcirc : Circuit K v J)
    (x : ι) (hx : x ∈ J) :
    ∃ coeff : ι → K,
      (∀ y ∈ J.erase x, coeff y ≠ 0) ∧
      v x = (J.erase x).sum (fun y => coeff y • v y) := by
  classical
  obtain ⟨c, hfull, hrel, hne⟩ :=
    exists_full_support_relation v J hcirc
  have hcx : c x ≠ 0 := hne x hx
  let coeff : ι → K := fun y => (-(c x)⁻¹) * c y
  refine ⟨coeff, ?_, ?_⟩
  · intro y hy
    exact mul_ne_zero (neg_ne_zero.mpr (inv_ne_zero hcx))
      (hne y ((Finset.erase_subset _ _) hy))
  · rw [Finsupp.linearCombination_apply] at hrel
    change c.support.sum (fun y => c y • v y) = 0 at hrel
    rw [hfull] at hrel
    have hsum :
        (J.erase x).sum (fun y => c y • v y) + c x • v x = 0 := by
      rw [Finset.sum_erase_add J (fun y => c y • v y) hx]
      exact hrel
    have hscaled := congrArg (fun z : F => (c x)⁻¹ • z) hsum
    simp only [smul_add, smul_zero, smul_smul,
      inv_mul_cancel₀ hcx, one_smul] at hscaled
    have hsolve :
        v x = -((c x)⁻¹ • (J.erase x).sum (fun y => c y • v y)) :=
      eq_neg_of_add_eq_zero_right hscaled
    calc
      v x = -((c x)⁻¹ • (J.erase x).sum (fun y => c y • v y)) := hsolve
      _ = (-(c x)⁻¹) • (J.erase x).sum (fun y => c y • v y) := by
        rw [neg_smul]
      _ = (J.erase x).sum (fun y => coeff y • v y) := by
        rw [Finset.smul_sum]
        apply Finset.sum_congr rfl
        intro y hy
        simp only [coeff, mul_smul]

end ICFallDegrees.GrowthCircuitReplacementSum
