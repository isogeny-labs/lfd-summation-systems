import ICFallDegrees.GrowthRatFuncRootMultiplicityCross
import ICFallDegrees.GrowthInfinityCofactorDegree
import Mathlib.Tactic

/-!
# N8: the genuine finite-pole exponent gap from a polynomial cofactor ratio

The earlier infinity bridge identifies u-u' with cofactor degree gaps.
This module extracts v-v' from *actual* X=1 polynomial root
multiplicities, using genuine equality in RatFunc K and the true N8
monomial fractions. No assumption about SpreadBound or a Wronskian
determinant occurs. The cofactor ratio is an explicit hypothesis:
proving it from omitted-column maximal minors is still required.
-/

noncomputable section
namespace ICFallDegrees.GrowthFiniteCofactorGap

open Polynomial
open ICFallDegrees.GrowthRatFuncRootMultiplicityCross
open ICFallDegrees.GrowthRatFuncInfinityDegree
open ICFallDegrees.GrowthFiniteCofactorMultiplicity

variable {K : Type*} [Field K]

/-- A quotient of actual rational-function fractions is the cross-product
fraction. This is valid with Lean's totalized division, even before
introducing nonzero hypotheses. -/
theorem mk_div_mk (P Q R S : Polynomial K) :
    RatFunc.mk P Q / RatFunc.mk R S =
      RatFunc.mk (P * S) (Q * R) := by
  simp only [RatFunc.mk_eq_div, map_mul, div_eq_mul_inv, mul_inv_rev, inv_inv]
  ring

/-- Actual rational monomial quotient in a common polynomial fraction. -/
theorem monomial_ratio_as_mk (c : K) (u v u' v' : ℕ) :
    ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
    RatFunc.mk
      ((X ^ u * (X + C c) ^ v) * ((X - 1 : Polynomial K) ^ v'))
      (((X - 1 : Polynomial K) ^ v) *
        (X ^ u' * (X + C c) ^ v')) := by
  unfold ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial
  exact mk_div_mk _ _ _ _

/-- In a *genuine* RatFunc cofactor-ratio equality, taking the actual
X=1 polynomial root multiplicity identifies the second-exponent gap.
This has not postulated an order at X=1. -/
theorem second_exponent_gap_eq_cofactor_root_gap
    (c gamma : K) (hc : c + 1 ≠ 0)
    (u v u' v' : ℕ) (A B : Polynomial K)
    (hA : A ≠ 0) (hB : B ≠ 0) (hgamma : gamma ≠ 0)
    (hratio :
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
        ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
      RatFunc.C gamma * RatFunc.mk A B) :
    (v' : ℤ) - (v : ℤ) =
      (rootMultiplicity (1 : K) A : ℤ) -
        (rootMultiplicity (1 : K) B : ℤ) := by
  obtain ⟨hN, hH⟩ := monomial_polynomials_ne_zero c u v
  obtain ⟨hN', hH'⟩ := monomial_polynomials_ne_zero c u' v'
  rw [monomial_ratio_as_mk] at hratio
  have hbal := rootMultiplicity_gap_of_scaled_mk_eq
    ((X ^ u * (X + C c) ^ v) * ((X - 1 : Polynomial K) ^ v'))
    (((X - 1 : Polynomial K) ^ v) * (X ^ u' * (X + C c) ^ v'))
    A B gamma
    (mul_ne_zero hN hH') (mul_ne_zero hH hN')
    hA hB hgamma hratio
  have hleft :
      rootMultiplicity (1 : K)
        ((X ^ u * (X + C c) ^ v) * ((X - 1 : Polynomial K) ^ v')) = v' := by
    rw [rootMultiplicity_mul (mul_ne_zero hN hH'),
      numerator_rootMultiplicity_one_eq_zero c hc u v,
      denominator_rootMultiplicity_one v', zero_add]
  have hright :
      rootMultiplicity (1 : K)
        (((X - 1 : Polynomial K) ^ v) * (X ^ u' * (X + C c) ^ v')) = v := by
    rw [rootMultiplicity_mul (mul_ne_zero hH hN'),
      denominator_rootMultiplicity_one v,
      numerator_rootMultiplicity_one_eq_zero c hc u' v', add_zero]
  rw [hleft, hright] at hbal
  exact hbal

/-- The second-exponent difference is controlled by the true
nonzero polynomial cofactor degree budget, using X=1 orders. -/
theorem second_exponent_gap_abs_le_cofactor_budget
    (c gamma : K) (hc : c + 1 ≠ 0)
    (u v u' v' : ℕ) (A B : Polynomial K) (budget : ℕ)
    (hA : A ≠ 0) (hB : B ≠ 0) (hgamma : gamma ≠ 0)
    (hAdeg : A.natDegree ≤ budget) (hBdeg : B.natDegree ≤ budget)
    (hratio :
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
        ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
      RatFunc.C gamma * RatFunc.mk A B) :
    |(v : ℤ) - (v' : ℤ)| ≤ budget := by
  have hgap := second_exponent_gap_eq_cofactor_root_gap
    c gamma hc u v u' v' A B hA hB hgamma hratio
  calc
    |(v : ℤ) - (v' : ℤ)| = |(v' : ℤ) - (v : ℤ)| := abs_sub_comm _ _
    _ = |(rootMultiplicity (1 : K) A : ℤ) -
      (rootMultiplicity (1 : K) B : ℤ)| := congrArg abs hgap
    _ ≤ budget := rootMultiplicity_gap_abs_le_budget A B budget hA hB hAdeg hBdeg

/-- Both *actual* pole computations combine into the 4B rather than 6B
L1 budget, conditional only on the genuine cofactor-ratio equality.
No formal SpreadBound or Wronskian-nonzero hypothesis is inserted. -/
theorem exponent_pair_gap_l1_le_four_budget
    (c gamma : K) (hc : c + 1 ≠ 0)
    (u v u' v' : ℕ) (A B : Polynomial K) (budget : ℕ)
    (hA : A ≠ 0) (hB : B ≠ 0) (hgamma : gamma ≠ 0)
    (hAdeg : A.natDegree ≤ budget) (hBdeg : B.natDegree ≤ budget)
    (hratio :
      ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u,v) /
        ICFallDegrees.GrowthRatFuncWronskianNormalization.monomial c (u',v') =
      RatFunc.C gamma * RatFunc.mk A B) :
    |2 * ((u : ℤ) - u')| + |2 * ((v : ℤ) - v')| ≤
      4 * (budget : ℤ) := by
  have hu := ICFallDegrees.GrowthInfinityCofactorDegree.exponent_gap_abs_le_cofactor_budget
    c gamma u v u' v' A B budget hA hB hgamma hAdeg hBdeg hratio
  have hv := second_exponent_gap_abs_le_cofactor_budget
    c gamma hc u v u' v' A B budget hA hB hgamma hAdeg hBdeg hratio
  rw [abs_mul, abs_mul]
  norm_num at *
  omega

end ICFallDegrees.GrowthFiniteCofactorGap
