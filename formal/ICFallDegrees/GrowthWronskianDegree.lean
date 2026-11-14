import Mathlib.Algebra.Polynomial.Degree.Definitions
import Mathlib.Tactic

/-!
# N8: degree filtration under finite polynomial products and sums

These elementary lemmata are the input to the Leibniz expansion of the
multi-column Wronskian cofactor. They do not assume Wronskian nonvanishing,
a spread bound, or CurveModel. They are candidates pending pinned CI.
-/

namespace ICFallDegrees.GrowthWronskianDegree

open Polynomial

variable {K ι : Type*} [Field K]

/-- Finite products of polynomials have degree at most the sum of
individual degree budgets, even when some entries vanish. -/
theorem natDegree_prod_le (s : Finset ι) (f : ι → Polynomial K) (d : ι → ℕ)
    (h : ∀ i ∈ s, (f i).natDegree ≤ d i) :
    (∏ i ∈ s, f i).natDegree ≤ ∑ i ∈ s, d i := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert a s ha ih =>
    simp only [Finset.prod_insert ha, Finset.sum_insert ha]
    exact Polynomial.natDegree_mul_le_of_le
      (h a (Finset.mem_insert_self a s))
      (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

/-- Finite sums do not increase the maximum degree bound. -/
theorem natDegree_sum_le (s : Finset ι) (f : ι → Polynomial K) (D : ℕ)
    (h : ∀ i ∈ s, (f i).natDegree ≤ D) :
    (∑ i ∈ s, f i).natDegree ≤ D := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert a s ha ih =>
    rw [Finset.sum_insert ha]
    simpa only [max_self] using
      (Polynomial.natDegree_add_le_of_le
        (h a (Finset.mem_insert_self a s))
        (ih (fun i hi => h i (Finset.mem_insert_of_mem hi))))

end ICFallDegrees.GrowthWronskianDegree
