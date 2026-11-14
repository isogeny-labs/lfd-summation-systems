import ICFallDegrees.GrowthWronskianDegree
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

/-!
# N8: determinant degree bound for multi-column polynomial Wronskians

A cofactor Wronskian is a determinant of an r-by-r polynomial matrix
whose row j has degree at most 3*j. The Leibniz expansion bounds its
degree by the sum of the row budgets, uniformly in characteristic.

This lemma requires no analytic CurveModel, Wronskian nonvanishing,
or S-unit theorem. Proof candidate pending pinned Lean CI.
-/

namespace ICFallDegrees.GrowthWronskianDeterminant

open Polynomial
open ICFallDegrees.GrowthWronskianDegree

variable {K ι : Type*} [Field K] [Fintype ι] [DecidableEq ι]

/-- A determinant has degree no larger than the sum of its *column*
degree bounds. Each Leibniz term selects exactly one entry per column. -/
theorem natDegree_det_le (M : Matrix ι ι (Polynomial K)) (d : ι → ℕ)
    (h : ∀ i j, (M i j).natDegree ≤ d j) :
    M.det.natDegree ≤ ∑ j : ι, d j := by
  classical
  rw [Matrix.det_apply']
  apply ICFallDegrees.GrowthWronskianDegree.natDegree_sum_le (s := Finset.univ)
  intro σ _
  have hp : (∏ j : ι, M (σ j) j).natDegree ≤ ∑ j : ι, d j := by
    exact natDegree_prod_le Finset.univ (fun j => M (σ j) j) d
      (fun j _ => h (σ j) j)
  have hc : (((Equiv.Perm.sign σ : ℤ) : Polynomial K)).natDegree ≤ 0 := by
    simp
  have hm := Polynomial.natDegree_mul_le_of_le hc hp
  simpa using hm

end ICFallDegrees.GrowthWronskianDeterminant
