import ICFallDegrees.GrowthFiniteStabilization
import Mathlib.RingTheory.MvPolynomial.Basic
import Mathlib.Tactic

/-!
# N8: finite degree truncation and the least stable row space

The ambient space of polynomials in 2k variables with total degree ≤ D
is finite dimensional over any field, including Q and F_p.

For any monotone inflationary row-extension operator T on a
finite-dimensional space, iteration from a seed F stabilizes at the least
T-closed submodule containing F. This is not yet the N8-specific T, nor
does it assert equality with the inductive InClosure.
-/

noncomputable section
namespace ICFallDegrees.GrowthBoundedDegreeLeastClosure

open ICFallDegrees.GrowthFiniteStabilization

variable {K V : Type*} [Field K] [AddCommGroup V] [Module K V]
  [FiniteDimensional K V]

/-- Monotone inflationary finite-dimensional span iterations give a
least fixed point, with no chosen echelon basis or pivot-order hypothesis. -/
theorem exists_least_stable
    (T : Submodule K V → Submodule K V)
    (hmon : Monotone T)
    (hext : ∀ U, U ≤ T U)
    (F : Submodule K V) :
    ∃ C : Submodule K V,
      F ≤ C ∧ T C = C ∧
      ∀ U : Submodule K V, F ≤ U → T U ≤ U → C ≤ U := by
  let W : ℕ → Submodule K V := Nat.rec F (fun _ U => T U)
  have hstep (n : ℕ) : W (n + 1) = T (W n) := rfl
  have hg (n : ℕ) : W n ≤ W (n + 1) := by
    rw [hstep]
    exact hext _
  have hseed : ∀ n, F ≤ W n := by
    intro n
    induction n with
    | zero => exact le_rfl
    | succ n ih => exact ih.trans (hg n)
  have hminimal (U : Submodule K V) (hF : F ≤ U) (hT : T U ≤ U) :
      ∀ n, W n ≤ U := by
    intro n
    induction n with
    | zero => exact hF
    | succ n ih =>
        rw [hstep]
        exact (hmon ih).trans hT
  obtain ⟨n, _hn, hfixed, _hstable⟩ :=
    eventual_stabilization T W hg hstep
  exact ⟨W n, hseed n, hfixed, fun U hF hT => hminimal U hF hT n⟩

/-- The actual N8 ambient degree-window as a finite-dimensional space. -/
abbrev DegreeWindow (K : Type*) [Field K] (k D : ℕ) :=
  ↥(MvPolynomial.restrictTotalDegree (Fin (2 * k)) K D)

/-- Bounded-degree polynomials in finitely many variables form a finite
dimensional coefficient-field vector space at every depth and degree. -/
theorem degreeWindow_finite (k D : ℕ) :
    FiniteDimensional K (DegreeWindow K k D) := by
  infer_instance

/-- Apply the least-closed finite-stage theorem specifically to the
2k-variable N8 polynomial degree window, over an arbitrary field K. -/
theorem exists_least_bounded_degree_stable (k D : ℕ)
    (T : Submodule K (DegreeWindow K k D) →
        Submodule K (DegreeWindow K k D))
    (hmon : Monotone T)
    (hext : ∀ U, U ≤ T U)
    (F : Submodule K (DegreeWindow K k D)) :
    ∃ C : Submodule K (DegreeWindow K k D),
      F ≤ C ∧ T C = C ∧
      ∀ U : Submodule K (DegreeWindow K k D),
        F ≤ U → T U ≤ U → C ≤ U := by
  letI : FiniteDimensional K (DegreeWindow K k D) := degreeWindow_finite k D
  exact exists_least_stable T hmon hext F

#print axioms exists_least_stable
#print axioms degreeWindow_finite
#print axioms exists_least_bounded_degree_stable

end ICFallDegrees.GrowthBoundedDegreeLeastClosure
