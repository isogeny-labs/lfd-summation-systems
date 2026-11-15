import ICFallDegrees.GrowthWronskianDeterminant
import ICFallDegrees.GrowthRationalWronskianRec
import Mathlib.Tactic

/-!
# N8 rational Wronskian cofactor matrix and degree budget

A column j of the polynomial cofactor matrix consists of j-th
derivative numerators, each bounded by degree 3*j. Combining the
genuine polynomial recurrence and the determinant degree lemma
establishes the uniform polynomial cofactor degree bound.

This does NOT assert that the Wronskian is nonzero or derive
exponent spread; these are remaining proof obligations. Candidates
await pinned Lean CI and #print axioms.
-/

noncomputable section

namespace ICFallDegrees.GrowthWronskianCofactor

open Polynomial
open ICFallDegrees.GrowthRationalWronskianRec
open ICFallDegrees.GrowthWronskianDeterminant

variable {K : Type*} [Field K]

/-- The polynomial matrix in the factorization of a rational
Wronskian. Every column is a fixed derivative order. -/
def matrix (r : ℕ) (D : Polynomial K) (Q : Fin r → Polynomial K) :
    Matrix (Fin r) (Fin r) (Polynomial K) :=
  fun i j => numeratorRec D (Q i) j.val

/-- Determinant degree bound for arbitrary choices of numerator rows,
provided D is cubic and each Q_i is quadratic. -/
theorem cofactor_natDegree_le (r : ℕ) (D : Polynomial K)
    (Q : Fin r → Polynomial K)
    (hD : D.natDegree ≤ 3)
    (hQ : ∀ i, (Q i).natDegree ≤ 2) :
    (matrix r D Q).det.natDegree ≤
      ∑ j : Fin r, 3 * j.val := by
  apply natDegree_det_le (d := fun j : Fin r => 3 * j.val)
  intro i j
  exact numeratorRec_natDegree_le D (Q i) hD (hQ i) j.val

end ICFallDegrees.GrowthWronskianCofactor
