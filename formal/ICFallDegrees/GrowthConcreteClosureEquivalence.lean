import ICFallDegrees.GrowthConcreteClosureStep
import Mathlib.Algebra.Module.Submodule.Map
import Mathlib.Tactic

/-!
# N8: completeness of the actual finite-stage bounded-degree closure

This closes the major soundness/completeness interface: the least fixed point
of the *concrete* operation adjoining all degree-budget-permitted products of
current rows equals the original inductive `InClosure`, over every field.

The proof uses `Closure.closure_minimal` on the ambient polynomial module,
not an illicit "multiply only the initial rows" shortcut. Since the bounded
degree window is finite-dimensional, this gives a finite row-space algorithm
that exactly computes `InClosure` after finitely many stages.

It does NOT yet show this finite algorithm commutes with prime specialization,
or establish any unconditional prime-field LFD lower bound.
-/

noncomputable section
namespace ICFallDegrees.GrowthConcreteClosureEquivalence

open MvPolynomial ICFallDegrees.GrowthConcreteClosureStep
open ICFallDegrees.GrowthBoundedDegreeLeastClosure

variable (K : Type*) [Field K] (k D : ℕ)

/-- Embed a subspace of the degree window back in the full polynomial ring. -/
def embedded (U : Submodule K (Window K k D)) :
    Submodule K (Poly K k) :=
  U.map ((MvPolynomial.restrictTotalDegree (Fin (2 * k)) K D).subtype)

/-- Membership in the embedded submodule is witnessed by an actual bounded
polynomial, including its membership in the original stage. -/
theorem mem_embedded_iff (U : Submodule K (Window K k D))
    (f : Poly K k) :
    f ∈ embedded K k D U ↔
      ∃ q : Window K k D, q ∈ U ∧ (q : Poly K k) = f := by
  change f ∈ U.map ((MvPolynomial.restrictTotalDegree (Fin (2 * k)) K D).subtype) ↔ _
  rw [Submodule.mem_map]
  rfl

/-- Any concrete step-closed stage containing the input rows also contains
all the inductive closure. The essential point is that multiplication can
use an arbitrary linear combination already in U, subject to the real
total-degree budget, not merely a raw input polynomial. -/
theorem actual_le_of_step_closed (F : Set (Poly K k))
    (U : Submodule K (Window K k D))
    (hseed : seed K k D F ≤ U)
    (hclosed : step K k D U ≤ U) :
    actual K k D F ≤ U := by
  have hmin : closureSpace F D ≤ embedded K k D U := by
    apply closure_minimal F D (embedded K k D U)
    · intro f hf hd
      apply (mem_embedded_iff K k D U f).2
      refine ⟨⟨f, (MvPolynomial.mem_restrictTotalDegree
        (Fin (2 * k)) D f).2 hd⟩, ?_, rfl⟩
      apply hseed
      apply Submodule.subset_span
      exact hf
    · intro a f hf hbudget
      obtain ⟨q, hq, hqf⟩ := (mem_embedded_iff K k D U f).1 hf
      subst f
      have hprodDegree : (a * (q : Poly K k)).totalDegree ≤ D :=
        (MvPolynomial.totalDegree_mul a (q : Poly K k)).trans hbudget
      let z : Window K k D :=
        ⟨a * (q : Poly K k), (MvPolynomial.mem_restrictTotalDegree
          (Fin (2 * k)) D (a * (q : Poly K k))).2 hprodDegree⟩
      have hzProducts : z ∈ permittedProducts K k D U := by
        exact ⟨a, q, hq, hbudget, rfl⟩
      have hzStep : z ∈ step K k D U := by
        change z ∈ U ⊔ Submodule.span K (permittedProducts K k D U)
        exact Submodule.mem_sup_right (Submodule.subset_span hzProducts)
      exact (mem_embedded_iff K k D U (a * (q : Poly K k))).2
        ⟨z, hclosed hzStep, rfl⟩
  intro q hq
  have hfull : (q : Poly K k) ∈ closureSpace F D := hq
  obtain ⟨z, hz, hzq⟩ := (mem_embedded_iff K k D U (q : Poly K k)).1 (hmin hfull)
  have heq : z = q := Subtype.ext hzq
  simpa [heq] using hz

/-- The actual textbook bounded closure is *the* least stable concrete
row-space; no hidden closure assumption or field-characteristic assumption. -/
theorem actual_least_stable (F : Set (Poly K k)) :
    ∃ C : Submodule K (Window K k D),
      seed K k D F ≤ C ∧ step K k D C = C ∧
      C = actual K k D F := by
  obtain ⟨C, hs, hf, hmin⟩ := exists_concrete_least_stable K k D F
  refine ⟨C, hs, hf, ?_⟩
  apply le_antisymm
  · exact hmin (actual K k D F) (seed_le_actual K k D F)
      (le_of_eq (step_actual K k D F))
  · exact actual_le_of_step_closed K k D F C hs (le_of_eq hf)

/-- In particular the actual inductive N8 degree-D closure admits a
finite-stage computation by linear spaces over any field K. -/
theorem exists_finite_exact_closure (F : Set (Poly K k)) :
    ∃ C : Submodule K (Window K k D),
      step K k D C = C ∧
      ∀ q : Window K k D, q ∈ C ↔ InClosure F D (q : Poly K k) := by
  obtain ⟨C, _hs, hc, heq⟩ := actual_least_stable K k D F
  refine ⟨C, hc, ?_⟩
  intro q
  rw [heq]
  rfl

#print axioms mem_embedded_iff
#print axioms actual_le_of_step_closed
#print axioms actual_least_stable
#print axioms exists_finite_exact_closure

end ICFallDegrees.GrowthConcreteClosureEquivalence
