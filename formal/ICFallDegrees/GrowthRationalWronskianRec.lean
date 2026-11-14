import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.Degree.Definitions
import Mathlib.Tactic

/-!
# N8 rational spread: universal polynomial Wronskian recurrence

This is a *polynomial* calculation in K[X], not an assumed CurveModel.
For the rational monomials X^u*((X+c)/(X-1))^v the actual
denominator D has degree at most 3 and logarithmic numerator Q
has degree at most 2. The recursive numerator for the s-th derivative
then has degree at most 3*s, which is already sufficient for the
6*(m-1)*(m-2) constant in GrowthSpread.SpreadBound.

This file deliberately DOES NOT assert the missing Wronskian
nonvanishing, cofactor comparison or rational residue spread.
-/

namespace ICFallDegrees.GrowthRationalWronskianRec

open Polynomial

variable {K : Type*} [Field K]

/-- Derivative-numerator recurrence, with a polynomial logarithmic
denominator D and numerator Q. -/
noncomputable def numeratorRec (D Q : Polynomial K) : ℕ → Polynomial K
  | 0 => 1
  | s + 1 =>
      D * (numeratorRec D Q s).derivative +
        (Q - Polynomial.C (s : K) * D.derivative) * numeratorRec D Q s

@[simp] theorem numeratorRec_zero (D Q : Polynomial K) :
    numeratorRec D Q 0 = 1 := rfl

@[simp] theorem numeratorRec_one (D Q : Polynomial K) :
    numeratorRec D Q 1 = Q := by
  simp [numeratorRec]

/-- The coarse 3*s bound is enough for the EXACT existing N8 spread
constant 6*(m-1)*(m-2), and avoids a special degree-zero derivative case. -/
theorem numeratorRec_natDegree_le (D Q : Polynomial K)
    (hD : D.natDegree ≤ 3) (hQ : Q.natDegree ≤ 2) (s : ℕ) :
    (numeratorRec D Q s).natDegree ≤ 3 * s := by
  induction s with
  | zero =>
      simp [numeratorRec]
  | succ n ih =>
      change (D * (numeratorRec D Q n).derivative +
        (Q - Polynomial.C (n : K) * D.derivative) * numeratorRec D Q n).natDegree
          ≤ 3 * (n + 1)
      have hdP : (numeratorRec D Q n).derivative.natDegree ≤ 3 * n :=
        (Polynomial.natDegree_derivative_le _).trans
          ((Nat.sub_le _ _).trans ih)
      have hdD : D.derivative.natDegree ≤ 3 :=
        (Polynomial.natDegree_derivative_le _).trans
          ((Nat.sub_le _ _).trans hD)
      have hc : (Polynomial.C (n : K)).natDegree ≤ 0 := by simp
      have hcD : (Polynomial.C (n : K) * D.derivative).natDegree ≤ 3 := by
        have h := Polynomial.natDegree_mul_le_of_le hc hdD
        simpa using h
      have hqD : (Q - Polynomial.C (n : K) * D.derivative).natDegree ≤ 3 := by
        have hq' : Q.natDegree ≤ 3 := by omega
        have h := Polynomial.natDegree_sub_le_of_le hq' hcD
        simpa using h
      have hleft : (D * (numeratorRec D Q n).derivative).natDegree ≤ 3 + 3 * n :=
        Polynomial.natDegree_mul_le_of_le hD hdP
      have hright : ((Q - Polynomial.C (n : K) * D.derivative) *
          numeratorRec D Q n).natDegree ≤ 3 + 3 * n :=
        Polynomial.natDegree_mul_le_of_le hqD ih
      have hsum := Polynomial.natDegree_add_le_of_le hleft hright
      simpa only [max_self] using (hsum.trans (by omega))

end ICFallDegrees.GrowthRationalWronskianRec
