import ICFallDegrees.GrowthWronskianCofactor
import ICFallDegrees.GrowthRationalLogNumerator
import Mathlib.Tactic

/-!
# N8: a tight polynomial degree budget for rational Wronskians

For D of degree at most three and Q of degree at most two, the
recurrence P₀ = 1, Pₙ₊₁ = D Pₙ' + (Q - nD') Pₙ actually satisfies
deg(Pₙ) ≤ 2n, not just the inherited 3n bound.

The base case n=0 of the derivative induction is treated separately,
because natDegree 0 = 0. This tightens an r-by-r polynomial
Wronskian cofactor to degree ≤ Σ_(j<r) 2j = r(r-1).

This is polynomial algebra; Wronskian nonvanishing, the RatFunc
derivative identification, and GrowthSpread.SpreadBound remain open.
-/

noncomputable section

namespace ICFallDegrees.GrowthWronskianTightDegree

open Polynomial
open ICFallDegrees.GrowthRationalWronskianRec
open ICFallDegrees.GrowthWronskianCofactor

variable {K : Type*} [Field K]

/-- The two-degree-per-derivative recurrence budget. In the successor
step the derivative has one fewer degree (once n ≥ 1), exactly
compensating for the cubic D. -/
theorem numeratorRec_natDegree_le_two (D Q : Polynomial K)
    (hD : D.natDegree ≤ 3) (hQ : Q.natDegree ≤ 2) (s : ℕ) :
    (numeratorRec D Q s).natDegree ≤ 2 * s := by
  induction s with
  | zero =>
      simp [numeratorRec]
  | succ n ih =>
      cases n with
      | zero =>
          simpa [numeratorRec] using hQ
      | succ n =>
          change (D * (numeratorRec D Q (n + 1)).derivative +
            (Q - C (((n + 1 : ℕ) : K)) * D.derivative) *
              numeratorRec D Q (n + 1)).natDegree ≤ 2 * ((n + 1) + 1)
          have hprev : (numeratorRec D Q (n + 1)).natDegree ≤ 2 * (n + 1) := by
            simpa only [Nat.succ_eq_add_one] using ih
          have hder₀ := Polynomial.natDegree_derivative_le (numeratorRec D Q (n + 1))
          have hder : (numeratorRec D Q (n + 1)).derivative.natDegree ≤
              2 * (n + 1) - 1 := by omega
          have hDder₀ := Polynomial.natDegree_derivative_le D
          have hDder : D.derivative.natDegree ≤ 2 := by omega
          have hC : (C (((n + 1 : ℕ) : K))).natDegree ≤ 0 :=
            Mathlib.Tactic.ComputeDegree.natDegree_C_le _
          have hCD : (C (((n + 1 : ℕ) : K)) * D.derivative).natDegree ≤ 2 := by
            simpa only [Nat.zero_add] using
              (Polynomial.natDegree_mul_le_of_le hC hDder)
          have hQsub : (Q - C (((n + 1 : ℕ) : K)) * D.derivative).natDegree ≤ 2 := by
            simpa only [max_self] using
              (Polynomial.natDegree_sub_le_of_le hQ hCD)
          have hleft : (D * (numeratorRec D Q (n + 1)).derivative).natDegree ≤
              2 * ((n + 1) + 1) := by
            have h := Polynomial.natDegree_mul_le_of_le hD hder
            omega
          have hright :
              ((Q - C (((n + 1 : ℕ) : K)) * D.derivative) *
                numeratorRec D Q (n + 1)).natDegree ≤ 2 * ((n + 1) + 1) := by
            have h := Polynomial.natDegree_mul_le_of_le hQsub hprev
            omega
          simpa only [max_self] using
            (Polynomial.natDegree_add_le_of_le hleft hright)

/-- An r-column polynomial cofactor is bounded by r(r-1)
rather than the inherited 3r(r-1)/2. -/
theorem cofactor_natDegree_le_two (r : ℕ) (D : Polynomial K)
    (Q : Fin r → Polynomial K) (hD : D.natDegree ≤ 3)
    (hQ : ∀ i, (Q i).natDegree ≤ 2) :
    (matrix r D Q).det.natDegree ≤ ∑ j : Fin r, 2 * j.val := by
  apply ICFallDegrees.GrowthWronskianDeterminant.natDegree_det_le
    (d := fun j : Fin r => 2 * j.val)
  intro i j
  exact numeratorRec_natDegree_le_two D (Q i) hD (hQ i) j.val

end ICFallDegrees.GrowthWronskianTightDegree
