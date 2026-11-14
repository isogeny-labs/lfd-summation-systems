import ICFallDegrees.GrowthConcreteClosureEquivalence
import Mathlib.Tactic

/-!
# N8: a genuine finite number of stages computes the actual InClosure

A depth-k, degree-D N8 system has at most
`Module.finrank K (DegreeWindow K k D)` strict growth steps.
The stage is *explicitly* iterated from the original input row-span using
the real degree-budget multiplication operator, and at some stage within
this bound it is exactly the inductively generated `InClosure`.

This is an exact finite-dimensional termination theorem for every field K,
not a p-specialization or LFD-growth theorem. The finite-stage closure
still requires a separate integer-matrix/rank-specialization proof.
-/

noncomputable section
namespace ICFallDegrees.GrowthConcreteClosureStages

open ICFallDegrees.GrowthConcreteClosureStep
open ICFallDegrees.GrowthConcreteClosureEquivalence
open ICFallDegrees.GrowthFiniteStabilization
open ICFallDegrees.GrowthBoundedDegreeLeastClosure

variable (K : Type*) [Field K] (k D : ℕ)

/-- Deterministic row-space iteration: first input span, then all
degree-admissible products of every row in the previous stage. -/
def stages (F : Set (Poly K k)) : ℕ → Submodule K (Window K k D) :=
  Nat.rec (seed K k D F) (fun _ W => step K k D W)

theorem stages_zero (F : Set (Poly K k)) :
    stages K k D F 0 = seed K k D F := rfl

theorem stages_succ (F : Set (Poly K k)) (n : ℕ) :
    stages K k D F (n + 1) = step K k D (stages K k D F n) := rfl

theorem stages_increasing (F : Set (Poly K k)) (n : ℕ) :
    stages K k D F n ≤ stages K k D F (n + 1) := by
  rw [stages_succ]
  exact step_inflationary K k D _

theorem seed_le_stages (F : Set (Poly K k)) :
    ∀ n, seed K k D F ≤ stages K k D F n := by
  intro n
  induction n with
  | zero => exact le_refl _
  | succ n ih =>
      exact ih.trans (stages_increasing K k D F n)

theorem stages_sound (F : Set (Poly K k)) :
    ∀ n, stages K k D F n ≤ actual K k D F := by
  intro n
  induction n with
  | zero => exact seed_le_actual K k D F
  | succ n ih =>
      rw [stages_succ]
      exact (step_monotone K k D ih).trans
        (le_of_eq (step_actual K k D F))

/-- **Exact degree-D closure by stage ≤ finite dimension**, uniformly
over arbitrary characteristic. In particular no recursion beyond this
finite bound is needed to represent the textbook InClosure. -/
theorem stages_exact_by_finrank (F : Set (Poly K k)) :
    ∃ n ≤ Module.finrank K (Window K k D),
      stages K k D F n = actual K k D F := by
  letI : FiniteDimensional K (Window K k D) :=
    degreeWindow_finite (K := K) k D
  obtain ⟨n, hn, hfixed, _⟩ :=
    eventual_stabilization (step K k D) (stages K k D F)
      (stages_increasing K k D F) (stages_succ K k D F)
  refine ⟨n, hn, ?_⟩
  apply le_antisymm
  · exact stages_sound K k D F n
  · exact actual_le_of_step_closed K k D F (stages K k D F n)
      (seed_le_stages K k D F n) (le_of_eq hfixed)

/-- Finite-stage characterization directly in the original predicate. -/
theorem inClosure_iff_finite_stage (F : Set (Poly K k)) :
    ∃ n ≤ Module.finrank K (Window K k D),
      ∀ f : Window K k D,
        f ∈ stages K k D F n ↔ InClosure F D (f : Poly K k) := by
  obtain ⟨n, hn, hstage⟩ := stages_exact_by_finrank K k D F
  refine ⟨n, hn, ?_⟩
  intro f
  rw [hstage]
  rfl

#print axioms stages_increasing
#print axioms stages_sound
#print axioms stages_exact_by_finrank
#print axioms inClosure_iff_finite_stage

end ICFallDegrees.GrowthConcreteClosureStages
