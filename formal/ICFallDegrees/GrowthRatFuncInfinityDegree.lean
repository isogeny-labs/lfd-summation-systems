import ICFallDegrees.GrowthRatFuncWronskianNormalization
import Mathlib.FieldTheory.RatFunc.Degree
import Mathlib.Tactic

/-!
# N8: genuine rational-function pole order at infinity

The actual function-field degree of X^u ((X+c)/(X-1))^v
is precisely u, for every coefficient field and c. This is the
first of two pole computations required in the N8 circuit spread
theorem, not an assumed integer-valued functional.

The order at infinity is minus RatFunc.intDegree. The order at X=1
and the omitted-index cofactor ratio are separate open bridges.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncInfinityDegree

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K]

/-- Invert a nonzero rational function: its integer degree changes sign. -/
theorem intDegree_inv (f : RatFunc K) (hf : f ≠ 0) :
    RatFunc.intDegree f⁻¹ = - RatFunc.intDegree f := by
  have h := RatFunc.intDegree_mul hf (inv_ne_zero hf)
  rw [mul_inv_cancel₀ hf, RatFunc.intDegree_one] at h
  omega

/-- The integer degree of a genuine polynomial fraction with nonzero
numerator and denominator equals their ordinary degree difference. -/
theorem intDegree_mk (p q : Polynomial K)
    (hp : p ≠ 0) (hq : q ≠ 0) :
    RatFunc.intDegree (RatFunc.mk p q) =
      (p.natDegree : ℤ) - (q.natDegree : ℤ) := by
  rw [RatFunc.mk_eq_div, div_eq_mul_inv,
    RatFunc.intDegree_mul (RatFunc.algebraMap_ne_zero hp)
      (inv_ne_zero (RatFunc.algebraMap_ne_zero hq)),
    intDegree_inv _ (RatFunc.algebraMap_ne_zero hq)]
  simp only [RatFunc.intDegree_polynomial]
  ring

/-- The concrete numerator and denominator polynomials of every rational
two-torsion monomial are nonzero. No genericity on c is required. -/
theorem monomial_polynomials_ne_zero (c : K) (u v : ℕ) :
    (X ^ u * (X + C c) ^ v : Polynomial K) ≠ 0 ∧
      ((X - 1 : Polynomial K) ^ v) ≠ 0 := by
  have hx : (X : Polynomial K) ≠ 0 := X_ne_zero
  have hc : (X + C c : Polynomial K) ≠ 0 := by
    simpa using (X_sub_C_ne_zero (-c))
  have h1 : (X - 1 : Polynomial K) ≠ 0 := by
    simpa using (X_sub_C_ne_zero (1 : K))
  exact ⟨mul_ne_zero (pow_ne_zero _ hx) (pow_ne_zero _ hc),
    pow_ne_zero _ h1⟩

/-- Every concrete rational monomial is nonzero, uniformly in c,u,v. -/
theorem monomial_ne_zero (c : K) (u v : ℕ) :
    GrowthRatFuncWronskianNormalization.monomial c (u,v) ≠ 0 := by
  obtain ⟨hp,hq⟩ := monomial_polynomials_ne_zero c u v
  unfold GrowthRatFuncWronskianNormalization.monomial
  rw [RatFunc.mk_eq_div]
  exact div_ne_zero (RatFunc.algebraMap_ne_zero hp)
    (RatFunc.algebraMap_ne_zero hq)

/-- Exact pole-at-infinity formula for the *actual* N8 rational monomial:
the two equal-degree v factors cancel in intDegree. -/
theorem monomial_intDegree (c : K) (u v : ℕ) :
    RatFunc.intDegree (GrowthRatFuncWronskianNormalization.monomial c (u,v)) = (u : ℤ) := by
  obtain ⟨hp,hq⟩ := monomial_polynomials_ne_zero c u v
  have hdnum :
      (X ^ u * (X + C c) ^ v : Polynomial K).natDegree = u + v := by
    have hx : (X : Polynomial K) ≠ 0 := X_ne_zero
    have hc : (X + C c : Polynomial K) ≠ 0 := by
      simpa using (X_sub_C_ne_zero (-c))
    rw [Polynomial.natDegree_mul (pow_ne_zero _ hx) (pow_ne_zero _ hc)]
    simp [Polynomial.natDegree_pow, Polynomial.natDegree_X,
      Polynomial.natDegree_X_add_C]
  have hdden : ((X - 1 : Polynomial K) ^ v).natDegree = v := by
    have hsub : (X - 1 : Polynomial K).natDegree = 1 := by
      simpa using (Polynomial.natDegree_X_sub_C (1 : K))
    rw [Polynomial.natDegree_pow, hsub]
    omega
  unfold GrowthRatFuncWronskianNormalization.monomial
  rw [intDegree_mk _ _ hp hq, hdnum, hdden]
  push_cast
  omega

/-- The degree of a quotient of two nonzero N8 monomials is exactly
the difference of their first exponents. -/
theorem monomial_ratio_intDegree (c : K) (u v u' v' : ℕ) :
    RatFunc.intDegree
      (GrowthRatFuncWronskianNormalization.monomial c (u,v) / GrowthRatFuncWronskianNormalization.monomial c (u',v')) =
      (u : ℤ) - (u' : ℤ) := by
  have ha := monomial_ne_zero c u v
  have hb := monomial_ne_zero c u' v'
  rw [div_eq_mul_inv, RatFunc.intDegree_mul ha (inv_ne_zero hb),
    intDegree_inv _ hb, monomial_intDegree c u v,
    monomial_intDegree c u' v']
  ring

end ICFallDegrees.GrowthRatFuncInfinityDegree
