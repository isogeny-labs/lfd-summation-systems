import ICFallDegrees.GrowthRatFuncInfinityDegree
import Mathlib.Tactic

/-!
# N8: extract the first exponent gap from an *actual* cofactor ratio

If a constant times a ratio of nonzero genuine polynomial cofactors
equals a quotient of two rational N8 monomials, the pole order at
infinity identifies their first exponent difference with the difference
of cofactor degrees. This does not assume circuit spread or Wronskian
nonvanishing. The *cofactor identity* itself is an explicit hypothesis
that must be discharged using the genuine maximal-minor determinant
replacement theorem; it is not manufactured here.
-/

noncomputable section
namespace ICFallDegrees.GrowthInfinityCofactorDegree

open Polynomial
open ICFallDegrees.GrowthRatFuncInfinityDegree

variable {K : Type*} [Field K]

/-- An equality of actual rational functions, up to a nonzero
coefficient-field scalar, forces equality of their infinity degrees. -/
theorem exponent_gap_eq_cofactor_degree_gap
    (c gamma : K) (u v u' v' : ℕ) (A B : Polynomial K)
    (hA : A ≠ 0) (hB : B ≠ 0) (hgamma : gamma ≠ 0)
    (hratio :
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
        ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
        RatFunc.C gamma * RatFunc.mk A B) :
    (u : ℤ) - (u' : ℤ) = (A.natDegree : ℤ) - (B.natDegree : ℤ) := by
  have hCpoly : (Polynomial.C gamma : Polynomial K) ≠ 0 := by
    simpa using hgamma
  have hC : RatFunc.C gamma ≠ (0 : RatFunc K) := by
    rw [← RatFunc.algebraMap_C]
    exact RatFunc.algebraMap_ne_zero hCpoly
  have hAB : RatFunc.mk A B ≠ 0 := by
    rw [RatFunc.mk_eq_div]
    exact div_ne_zero (RatFunc.algebraMap_ne_zero hA)
      (RatFunc.algebraMap_ne_zero hB)
  calc
    (u : ℤ) - (u' : ℤ) =
        RatFunc.intDegree
          (ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
            ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v')) := by
              exact (monomial_ratio_intDegree c u v u' v').symm
    _ = RatFunc.intDegree (RatFunc.C gamma * RatFunc.mk A B) := by rw [hratio]
    _ = (A.natDegree : ℤ) - (B.natDegree : ℤ) := by
      rw [RatFunc.intDegree_mul hC hAB, RatFunc.intDegree_C,
        intDegree_mk A B hA hB]
      ring

/-- Once the actual cofactor ratio is established, the first-exponent
difference is bounded by the sharp polynomial cofactor budget. -/
theorem exponent_gap_abs_le_cofactor_budget
    (c gamma : K) (u v u' v' : ℕ) (A B : Polynomial K) (budget : ℕ)
    (hA : A ≠ 0) (hB : B ≠ 0) (hgamma : gamma ≠ 0)
    (hAdeg : A.natDegree ≤ budget) (hBdeg : B.natDegree ≤ budget)
    (hratio :
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
        ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
        RatFunc.C gamma * RatFunc.mk A B) :
    |(u : ℤ) - (u' : ℤ)| ≤ budget := by
  have hgap := exponent_gap_eq_cofactor_degree_gap
    c gamma u v u' v' A B hA hB hgamma hratio
  rw [hgap]
  have hAcast : (A.natDegree : ℤ) ≤ budget := by exact_mod_cast hAdeg
  have hBcast : (B.natDegree : ℤ) ≤ budget := by exact_mod_cast hBdeg
  apply abs_le.mpr
  constructor <;> omega

end ICFallDegrees.GrowthInfinityCofactorDegree
