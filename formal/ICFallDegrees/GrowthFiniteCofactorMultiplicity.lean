import ICFallDegrees.GrowthInfinityCofactorDegree
import Mathlib.Algebra.Polynomial.RingDivision
import Mathlib.Algebra.Polynomial.Degree.Domain
import Mathlib.Tactic

/-!
# N8: exact polynomial factors and root orders at X = 1

This is the finite-pole *polynomial* bridge needed for the true
rational cofactor comparison. It uses Mathlib's actual rootMultiplicity,
not an invented integer-valued order. No spread assumption is used.

The missing remaining theorem must equate cofactor-monomial quotients
from actual Wronskian maximal minors and transport that identity
through polynomial fractions to these root orders.
-/

noncomputable section
namespace ICFallDegrees.GrowthFiniteCofactorMultiplicity

open Polynomial
variable {K : Type*} [Field K]

/-- Every nonzero polynomial's root order at X = 1 is bounded
by its actual natural degree. -/
theorem rootMultiplicity_one_le_natDegree (A : Polynomial K) (hA : A ≠ 0) :
    rootMultiplicity (1 : K) A ≤ A.natDegree := by
  have hdvd :
      (X - C (1 : K)) ^ rootMultiplicity (1 : K) A ∣ A :=
    pow_rootMultiplicity_dvd A (1 : K)
  have hdeg := natDegree_le_of_dvd hdvd hA
  have hone : (X - (1 : Polynomial K)).natDegree = 1 := by
    simpa using (natDegree_X_sub_C (1 : K))
  simpa only [natDegree_pow, map_one, hone, mul_one] using hdeg

/-- The genuine N8 monomial numerator X^u (X+c)^v has no zero at
X=1 provided c+1 is nonzero. For the fixed curve c=8 this holds. -/
theorem numerator_rootMultiplicity_one_eq_zero
    (c : K) (hc : c + 1 ≠ 0) (u v : ℕ) :
    rootMultiplicity (1 : K)
      (X ^ u * (X + C c) ^ v : Polynomial K) = 0 := by
  have heval :
      (X ^ u * (X + C c) ^ v : Polynomial K).eval 1 =
        (1 + c) ^ v := by simp
  have hnotroot : ¬ IsRoot (X ^ u * (X + C c) ^ v : Polynomial K) (1 : K) := by
    intro h
    change (X ^ u * (X + C c) ^ v : Polynomial K).eval 1 = 0 at h
    rw [heval] at h
    exact pow_ne_zero v (by simpa [add_comm] using hc) h
  exact rootMultiplicity_eq_zero hnotroot

/-- The denominator of f_(u,v) has exactly v copies of X-1. -/
theorem denominator_rootMultiplicity_one (v : ℕ) :
    rootMultiplicity (1 : K) ((X - 1 : Polynomial K) ^ v) = v := by
  simpa using (rootMultiplicity_X_sub_C_pow (1 : K) v)

/-- Two genuine polynomial cofactors of degree at most B have a
root-multiplicity difference at X=1 of absolute value at most B. -/
theorem rootMultiplicity_gap_abs_le_budget
    (A B : Polynomial K) (budget : ℕ)
    (hA : A ≠ 0) (hB : B ≠ 0)
    (hAdeg : A.natDegree ≤ budget) (hBdeg : B.natDegree ≤ budget) :
    |(rootMultiplicity (1 : K) A : ℤ) -
      (rootMultiplicity (1 : K) B : ℤ)| ≤ budget := by
  have hAord := (rootMultiplicity_one_le_natDegree A hA).trans hAdeg
  have hBord := (rootMultiplicity_one_le_natDegree B hB).trans hBdeg
  have hAcast : (rootMultiplicity (1 : K) A : ℤ) ≤ budget := by
    exact_mod_cast hAord
  have hBcast : (rootMultiplicity (1 : K) B : ℤ) ≤ budget := by
    exact_mod_cast hBord
  apply abs_le.mpr
  constructor <;> omega

end ICFallDegrees.GrowthFiniteCofactorMultiplicity
