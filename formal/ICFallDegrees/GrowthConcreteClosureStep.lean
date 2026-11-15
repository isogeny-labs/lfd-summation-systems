import ICFallDegrees.GrowthBoundedDegreeLeastClosure
import Mathlib.RingTheory.MvPolynomial.Basic
import Mathlib.Tactic

/-!
# N8 concrete bounded-degree closure step

For the actual degree-D window in 2k variables, adjoin every permitted
polynomial multiple of every current row and take its K-linear span.
This includes multiplication of rows obtained from cancellation, not just
multiples of the original inputs. In particular it is not the usual
one-shot Macaulay matrix.

A next lemma must establish that its least fixed point is the inductive
`InClosure`. No prime-field theorem is assumed here.
-/

noncomputable section

namespace ICFallDegrees.GrowthConcreteClosureStep

open MvPolynomial ICFallDegrees.GrowthBoundedDegreeLeastClosure

variable (K : Type*) [Field K] (k D : ℕ)

abbrev Poly := MvPolynomial (Fin (2 * k)) K
abbrev Window := DegreeWindow K k D

/-- Each candidate is a permitted product of an *actual row* in U,
including a row obtained by cancellations within U. -/
def permittedProducts (U : Submodule K (Window K k D)) : Set (Window K k D) :=
  {q | ∃ (a : Poly K k) (f : Window K k D),
    f ∈ U ∧ a.totalDegree + (f : Poly K k).totalDegree ≤ D ∧
      (q : Poly K k) = a * (f : Poly K k)}

/-- The exact bounded-degree span step (not a one-shot input-multiple span). -/
def step (U : Submodule K (Window K k D)) : Submodule K (Window K k D) :=
  U ⊔ Submodule.span K (permittedProducts K k D U)

theorem permittedProducts_mono {U V : Submodule K (Window K k D)}
    (h : U ≤ V) :
    permittedProducts K k D U ⊆ permittedProducts K k D V := by
  intro q hq
  rcases hq with ⟨a, f, hf, budget, heq⟩
  exact ⟨a, f, h hf, budget, heq⟩

theorem step_inflationary (U : Submodule K (Window K k D)) :
    U ≤ step K k D U :=
  le_sup_left

theorem step_monotone :
    Monotone (step K k D) := by
  intro U V h
  change U ⊔ Submodule.span K (permittedProducts K k D U) ≤
    V ⊔ Submodule.span K (permittedProducts K k D V)
  refine sup_le (h.trans le_sup_left) ?_
  exact (Submodule.span_mono
    (permittedProducts_mono K k D h)).trans le_sup_right

/-- The actual inductive closure, regarded inside the finite-dimensional
degree-D window. No generators or multiplication rules are lost. -/
def actual (F : Set (Poly K k)) : Submodule K (Window K k D) where
  carrier := {q | InClosure F D (q : Poly K k)}
  zero_mem' := by
    change InClosure F D 0
    exact InClosure.zero
  add_mem' := by
    intro f g hf hg
    simpa using (InClosure.add hf hg)
  smul_mem' := by
    intro c f hf
    simpa using (InClosure.smul c hf)

/-- Permitted multiplication preserves the inductive closure. -/
theorem permittedProducts_actual (F : Set (Poly K k)) :
    permittedProducts K k D (actual K k D F) ⊆
      (actual K k D F : Set (Window K k D)) := by
  intro q hq
  rcases hq with ⟨a, f, hf, budget, heq⟩
  change InClosure F D (q : Poly K k)
  change InClosure F D (f : Poly K k) at hf
  rw [heq]
  exact InClosure.mul a hf budget

/-- The inductively defined textbook closure is a fixed point of the
concrete span step. This direction requires no artificial degree truncation
and includes cancellation-derived rows. -/
theorem step_actual (F : Set (Poly K k)) :
    step K k D (actual K k D F) = actual K k D F := by
  apply le_antisymm
  · change (actual K k D F ⊔
      Submodule.span K (permittedProducts K k D (actual K k D F))) ≤
      actual K k D F
    refine sup_le le_rfl ?_
    exact Submodule.span_le.mpr (permittedProducts_actual K k D F)
  · exact step_inflationary K k D _

/-- Initial rows: exactly the degree-at-most-D members of the original
input family, represented in the bounded window. -/
def inputRows (F : Set (Poly K k)) : Set (Window K k D) :=
  {q | (q : Poly K k) ∈ F}

def seed (F : Set (Poly K k)) : Submodule K (Window K k D) :=
  Submodule.span K (inputRows K k D F)

theorem seed_le_actual (F : Set (Poly K k)) :
    seed K k D F ≤ actual K k D F := by
  apply Submodule.span_le.mpr
  intro q hq
  change InClosure F D (q : Poly K k)
  apply InClosure.input hq
  exact (MvPolynomial.mem_restrictTotalDegree (Fin (2 * k)) D
    (q : Poly K k)).mp q.property

/-- Unconditional existence of a finite least stable *concrete* row space.
The remaining theorem will identify it with `actual`. -/
theorem exists_concrete_least_stable (F : Set (Poly K k)) :
    ∃ C : Submodule K (Window K k D),
      seed K k D F ≤ C ∧ step K k D C = C ∧
      ∀ U : Submodule K (Window K k D),
        seed K k D F ≤ U → step K k D U ≤ U → C ≤ U := by
  exact exists_least_bounded_degree_stable (K := K) k D
    (step K k D) (step_monotone K k D) (step_inflationary K k D)
    (seed K k D F)

/-- The stabilized stage cannot contain a spurious row: it is always
contained in the actual inductive degree-D closure. -/
theorem least_stable_sound (F : Set (Poly K k)) :
    ∃ C : Submodule K (Window K k D),
      seed K k D F ≤ C ∧ step K k D C = C ∧ C ≤ actual K k D F := by
  obtain ⟨C, hseed, hfixed, hmin⟩ :=
    exists_concrete_least_stable K k D F
  refine ⟨C, hseed, hfixed, ?_⟩
  exact hmin (actual K k D F) (seed_le_actual K k D F)
    (le_of_eq (step_actual K k D F))

#print axioms step_monotone
#print axioms step_actual
#print axioms exists_concrete_least_stable
#print axioms least_stable_sound

end ICFallDegrees.GrowthConcreteClosureStep
