import ICFallDegrees.GrowthSpread
import Mathlib.Tactic

/-!
# N8: a minimal circuit has a relation with no zero coefficients

A GrowthSpread.Circuit F₀ v J is dependent, but every proper subset
is independent. We prove the exact finite-support consequence: every
nonzero constant-coefficient relation supported on J necessarily uses
all indices of J, and such a relation exists. In particular, all
cofactor ratios can use nonzero constant coefficients.

No Wronskian, cofactor nonvanishing or pole valuation is assumed.
-/

noncomputable section

namespace ICFallDegrees.GrowthCircuitFullSupport

variable {F₀ F ι : Type*} [Field F₀] [AddCommGroup F] [Module F₀ F]

/-- An actual nonzero relation on a circuit cannot omit even one circuit
column; otherwise its proper support would be linearly dependent. -/
theorem relation_support_eq_circuit
    (v : ι → F) (J : Finset ι)
    (hcirc : ICFallDegrees.GrowthSpread.Circuit F₀ v J)
    (c : ι →₀ F₀)
    (hsubset : c.support ⊆ J)
    (hrel : Finsupp.linearCombination F₀ v c = 0)
    (hne : c ≠ 0) :
    c.support = J := by
  classical
  by_contra hneq
  have hproper : c.support ⊂ J :=
    Finset.ssubset_iff_subset_ne.mpr ⟨hsubset, hneq⟩
  have hlin := hcirc.2 c.support hproper
  have hdep : ¬ LinearIndepOn F₀ v (c.support : Set ι) :=
    linearDepOn_iff'.mpr
      ⟨c, Finsupp.mem_supported_support F₀ c, hrel, hne⟩
  exact hdep hlin

/-- Every coefficient of a nontrivial circuit relation is nonzero. -/
theorem relation_coeff_ne_zero
    (v : ι → F) (J : Finset ι)
    (hcirc : ICFallDegrees.GrowthSpread.Circuit F₀ v J)
    (c : ι →₀ F₀)
    (hsubset : c.support ⊆ J)
    (hrel : Finsupp.linearCombination F₀ v c = 0)
    (hne : c ≠ 0)
    (e : ι) (he : e ∈ J) :
    c e ≠ 0 := by
  have hs := relation_support_eq_circuit v J hcirc c hsubset hrel hne
  have hmem : e ∈ c.support := by simpa only [hs] using he
  exact Finsupp.mem_support_iff.mp hmem

/-- A circuit admits a genuine constant-coefficient relation using every
one of its columns. The existence follows from dependent-on-J, and
the all-nonzero assertion follows from minimality of J. -/
theorem exists_full_support_relation
    (v : ι → F) (J : Finset ι)
    (hcirc : ICFallDegrees.GrowthSpread.Circuit F₀ v J) :
    ∃ c : ι →₀ F₀, c.support = J ∧
      Finsupp.linearCombination F₀ v c = 0 ∧
      ∀ e ∈ J, c e ≠ 0 := by
  classical
  obtain ⟨c, hcJ, hrel, hne⟩ := linearDepOn_iff'.mp hcirc.1
  rw [Finsupp.mem_supported] at hcJ
  have hfull := relation_support_eq_circuit v J hcirc c hcJ hrel hne
  refine ⟨c, hfull, hrel, ?_⟩
  intro e he
  exact relation_coeff_ne_zero v J hcirc c hcJ hrel hne e he

end ICFallDegrees.GrowthCircuitFullSupport
