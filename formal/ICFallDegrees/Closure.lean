import Mathlib.Algebra.MvPolynomial.Degrees
import Mathlib.Algebra.MvPolynomial.CommRing
import Mathlib.RingTheory.Ideal.Span
import Lean.Elab.Tactic.Omega

/-!
# Degree-bounded ideal closure

This is the textbook `V_(F,D)` closure, using the literal multiplier budget
`totalDegree a + totalDegree f ≤ D`. In particular, multiplication applies
to newly derived rows, not merely original inputs.
-/

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K]

theorem totalDegree_neg_eq (f : MvPolynomial σ K) :
    (-f).totalDegree = f.totalDegree := by
  classical
  simp only [MvPolynomial.totalDegree, MvPolynomial.support_neg]

theorem totalDegree_sub_le (f g : MvPolynomial σ K) :
    (f-g).totalDegree ≤ max f.totalDegree g.totalDegree := by
  simpa only [sub_eq_add_neg, totalDegree_neg_eq] using totalDegree_add f (-g)

inductive InClosure (F : Set (MvPolynomial σ K)) (D : ℕ) :
    MvPolynomial σ K → Prop
  | input {f} (member : f ∈ F) (degree : f.totalDegree ≤ D) : InClosure F D f
  | zero : InClosure F D 0
  | add {f g} : InClosure F D f → InClosure F D g → InClosure F D (f + g)
  | smul (c : K) {f} : InClosure F D f → InClosure F D (c • f)
  | mul (a : MvPolynomial σ K) {f} : InClosure F D f →
      a.totalDegree + f.totalDegree ≤ D → InClosure F D (a * f)

namespace InClosure

theorem degree_le {F : Set (MvPolynomial σ K)} {D f}
    (h : InClosure F D f) : f.totalDegree ≤ D := by
  induction h with
  | input _ hd => exact hd
  | zero => simp
  | add _ _ ih₁ ih₂ => exact (totalDegree_add _ _).trans (max_le ih₁ ih₂)
  | smul c _ ih => exact (totalDegree_smul_le c _).trans ih
  | mul a _ budget _ => exact (totalDegree_mul a _).trans budget

theorem sound {F : Set (MvPolynomial σ K)} {D f}
    (h : InClosure F D f) : f ∈ Ideal.span F := by
  induction h with
  | input member _ => exact Ideal.subset_span member
  | zero => exact (Ideal.span F).zero_mem
  | add _ _ ih₁ ih₂ => exact (Ideal.span F).add_mem ih₁ ih₂
  | smul c _ ih =>
      rw [smul_eq_C_mul]
      exact (Ideal.span F).mul_mem_left _ ih
  | mul a _ _ ih => exact (Ideal.span F).mul_mem_left a ih

theorem mono {F : Set (MvPolynomial σ K)} {D E f}
    (hDE : D ≤ E) (h : InClosure F D f) : InClosure F E f := by
  induction h with
  | input member hd => exact .input member (hd.trans hDE)
  | zero => exact .zero
  | add _ _ ih₁ ih₂ => exact .add ih₁ ih₂
  | smul c _ ih => exact .smul c ih
  | mul a _ budget ih => exact .mul a ih (budget.trans hDE)

/-- Replace available inputs by checked derivations in a second system,
retaining an explicit preprocessing degree E. -/
theorem transport {F G : Set (MvPolynomial σ K)} {D E f}
    (hDE : D ≤ E)
    (inputs : ∀ p ∈ F, p.totalDegree ≤ D → InClosure G E p)
    (h : InClosure F D f) : InClosure G E f := by
  induction h with
  | input member hd => exact inputs _ member hd
  | zero => exact .zero
  | add _ _ ih₁ ih₂ => exact .add ih₁ ih₂
  | smul c _ ih => exact .smul c ih
  | mul a _ budget ih => exact .mul a ih (budget.trans hDE)

theorem neg {F : Set (MvPolynomial σ K)} {D f}
    (h : InClosure F D f) : InClosure F D (-f) := by
  simpa using InClosure.smul (-1 : K) h

theorem sub {F : Set (MvPolynomial σ K)} {D f g}
    (hf : InClosure F D f) (hg : InClosure F D g) : InClosure F D (f - g) := by
  simpa [sub_eq_add_neg] using hf.add hg.neg

end InClosure

def closureSpace (F : Set (MvPolynomial σ K)) (D : ℕ) :
    Submodule K (MvPolynomial σ K) where
  carrier := InClosure F D
  zero_mem' := .zero
  add_mem' := InClosure.add
  smul_mem' := InClosure.smul

/-- The inductive closure is contained in any input-containing submodule
closed under permitted multiplication. This is its least-space property. -/
theorem closure_minimal (F : Set (MvPolynomial σ K)) (D : ℕ)
    (W : Submodule K (MvPolynomial σ K))
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → f ∈ W)
    (products : ∀ a f, f ∈ W → a.totalDegree + f.totalDegree ≤ D → a * f ∈ W) :
    closureSpace F D ≤ W := by
  intro f hf
  induction hf with
  | input member degree => exact inputs _ member degree
  | zero => exact W.zero_mem
  | add _ _ ih₁ ih₂ => exact W.add_mem ih₁ ih₂
  | smul c _ ih => exact W.smul_mem c ih
  | mul a _ budget ih => exact products a _ ih budget

/-- The defining last-fall-degree property, expressed without a choice of a
minimum or an implicit assumption that a finite bound has been established. -/
def HasLFDBound (F : Set (MvPolynomial σ K)) (d : ℕ) : Prop :=
  ∀ f ∈ Ideal.span F, InClosure F (max d f.totalDegree) f

def IsLastFallDegree (F : Set (MvPolynomial σ K)) (d : ℕ) : Prop :=
  HasLFDBound F d ∧ ∀ e < d, ¬ HasLFDBound F e

theorem hasLFDBound_of_one {F : Set (MvPolynomial σ K)} {d}
    (one : InClosure F d 1) : HasLFDBound F d := by
  intro f _
  have hone := one.mono (Nat.le_max_left d f.totalDegree)
  simpa using InClosure.mul f hone (by simp)

theorem one_of_hasLFDBound {F : Set (MvPolynomial σ K)} {d}
    (unit : (1 : MvPolynomial σ K) ∈ Ideal.span F) (bound : HasLFDBound F d) :
    InClosure F d 1 := by
  simpa using bound 1 unit

/-- A refutation at degree d together with exclusion one degree earlier
establishes the exact LFD of a unit ideal. -/
theorem exactLFD_of_first_one {F : Set (MvPolynomial σ K)} {d : ℕ}
    (positive : 0 < d) (upper : InClosure F d 1)
    (lower : ¬ InClosure F (d - 1) 1) : IsLastFallDegree F d := by
  refine ⟨hasLFDBound_of_one upper, ?_⟩
  intro e he bound
  have hone := one_of_hasLFDBound upper.sound bound
  exact lower (hone.mono (by omega))

end ICFallDegrees
