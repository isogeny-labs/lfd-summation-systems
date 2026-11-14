import ICFallDegrees.GrowthRationalWronskianRec
import Mathlib.Tactic

/-!
# N8 rational logarithmic derivative numerator

For f_(u,v)=X^u ((X+c)/(X-1))^v and
D=X(X+c)(X-1), the numerator of f'/f is
Q=u(X+c)(X-1)-v(c+1)X.

This module proves the algebraic cancellation of the two v terms.
It does not yet identify RatFunc.derivative, build a Wronskian,
or prove SpreadBound. Lean proof candidates until CI succeeds.
-/

noncomputable section

namespace ICFallDegrees.GrowthRationalLogNumerator

open Polynomial

variable {K : Type*} [Field K]

def denominator (c : K) : Polynomial K :=
  X * (X + C c) * (X - 1)

def logNumerator (c : K) (u v : ℕ) : Polynomial K :=
  C (u : K) * (X + C c) * (X - 1) +
    C (v : K) * X * (X - 1) -
    C (v : K) * X * (X + C c)

/-- Cancellation of the two terms arising from the logarithmic derivative
of the Möbius coordinate. -/
theorem logNumerator_eq (c : K) (u v : ℕ) :
    logNumerator c u v =
      C (u : K) * (X + C c) * (X - 1) -
        C (v : K) * C (c + 1) * X := by
  unfold logNumerator
  simp only [map_add, map_one]
  ring

/-- The degree bound required for the already developed recurrence:
the logarithmic numerator is a polynomial of degree at most two. -/
theorem logNumerator_natDegree_le (c : K) (u v : ℕ) :
    (logNumerator c u v).natDegree ≤ 2 := by
  rw [logNumerator_eq]
  have hlinear : (X + C c : Polynomial K).natDegree ≤ 1 := by
    simpa using natDegree_add_le_of_le (by simp : (X : Polynomial K).natDegree ≤ 1)
      (by simp : (C c : Polynomial K).natDegree ≤ 1)
  have hsub : (X - 1 : Polynomial K).natDegree ≤ 1 := by
    simpa using natDegree_sub_le_of_le (by simp : (X : Polynomial K).natDegree ≤ 1)
      (by simp : (1 : Polynomial K).natDegree ≤ 1)
  have hprod : (C (u : K) * (X + C c) * (X - 1) : Polynomial K).natDegree ≤ 2 := by
    have hconst : (C (u : K) : Polynomial K).natDegree ≤ 0 := by simp
    have hfirst := natDegree_mul_le_of_le hconst hlinear
    have hlast := natDegree_mul_le_of_le hfirst hsub
    omega
  have hsecond : (C (v : K) * C (c + 1) * X : Polynomial K).natDegree ≤ 2 := by
    have hc1 : (C (v : K) : Polynomial K).natDegree ≤ 0 := by simp
    have hc2 : (C (c + 1) : Polynomial K).natDegree ≤ 0 := by simp
    have hh := natDegree_mul_le_of_le hc1 hc2
    have hh2 := natDegree_mul_le_of_le hh (by simp : (X : Polynomial K).natDegree ≤ 1)
    omega
  exact natDegree_sub_le_of_le hprod hsecond


/-- The common logarithmic denominator of both rational coordinates
is a product of three linear factors. -/
theorem denominator_natDegree_le (c : K) :
    (denominator c).natDegree ≤ 3 := by
  unfold denominator
  have hx : (X : Polynomial K).natDegree ≤ 1 := by simp
  have hc : (C c : Polynomial K).natDegree ≤ 1 := by simp
  have hone : (1 : Polynomial K).natDegree ≤ 1 := by simp
  have ha : (X + C c : Polynomial K).natDegree ≤ 1 :=
    natDegree_add_le_of_le hx hc
  have hb : (X - 1 : Polynomial K).natDegree ≤ 1 :=
    natDegree_sub_le_of_le hx hone
  have h := natDegree_mul_le_of_le (natDegree_mul_le_of_le hx ha) hb
  omega

/-- The concrete rational-coordinate logarithmic data satisfy the
recurrence's degree budget, uniformly over all exponents. -/
theorem concrete_numeratorRec_natDegree_le (c : K) (u v s : ℕ) :
    (ICFallDegrees.GrowthRationalWronskianRec.numeratorRec
      (denominator c) (logNumerator c u v) s).natDegree ≤ 3 * s :=
  ICFallDegrees.GrowthRationalWronskianRec.numeratorRec_natDegree_le
    (denominator c) (logNumerator c u v)
    (denominator_natDegree_le c) (logNumerator_natDegree_le c u v) s

end ICFallDegrees.GrowthRationalLogNumerator
