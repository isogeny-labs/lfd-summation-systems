import ICFallDegrees.GrowthFiniteCofactorMultiplicity
import Mathlib.Algebra.Polynomial.RingDivision
import Mathlib.Tactic

/-!
# N8: true rational-fraction equality preserves the finite-pole order

This module uses equality of actual `RatFunc K` elements, rather than
a formal valuation symbol. Cross-multiplication gives an equality of
polynomials, and actual Mathlib `Polynomial.rootMultiplicity` is
additive on the resulting nonzero products.

Nothing here assumes circuit spread, a cofactor comparison, or a
last-fall-degree bound. A separate bridge must identify the specific
Wronskian minor ratios with the N8 rational monomial ratios.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncRootMultiplicityCross

open Polynomial

variable {K : Type*} [Field K]

/-- Equality of nonzero fractions in the actual rational function
field preserves the sum of numerator and opposite-denominator root
multiplicities. This statement is valid without reduced fractions. -/
theorem rootMultiplicity_balance_of_mk_eq
    (A B C D : Polynomial K)
    (hA : A ≠ 0) (hB : B ≠ 0)
    (hC : C ≠ 0) (hD : D ≠ 0)
    (heq : RatFunc.mk A B = RatFunc.mk C D) :
    rootMultiplicity (1 : K) A + rootMultiplicity (1 : K) D =
      rootMultiplicity (1 : K) C + rootMultiplicity (1 : K) B := by
  let φ : Polynomial K →+* RatFunc K :=
    algebraMap (Polynomial K) (RatFunc K)
  have hφB : φ B ≠ 0 := RatFunc.algebraMap_ne_zero hB
  have hφD : φ D ≠ 0 := RatFunc.algebraMap_ne_zero hD
  have hdiv : φ A / φ B = φ C / φ D := by
    simpa only [RatFunc.mk_eq_div] using heq
  have hcross : φ (A * D) = φ (C * B) := by
    have h := (div_eq_div_iff hφB hφD).mp hdiv
    simpa only [map_mul] using h
  have hpoly : A * D = C * B := by
    by_contra hne
    have hdiff : A * D - C * B ≠ 0 := sub_ne_zero.mpr hne
    have hnonzero : φ (A * D - C * B) ≠ 0 :=
      RatFunc.algebraMap_ne_zero hdiff
    apply hnonzero
    rw [map_sub, hcross, sub_self]
  calc
    rootMultiplicity (1 : K) A + rootMultiplicity (1 : K) D =
      rootMultiplicity (1 : K) (A * D) := by
        exact (rootMultiplicity_mul (mul_ne_zero hA hD)).symm
    _ = rootMultiplicity (1 : K) (C * B) := by
      rw [hpoly]
    _ = rootMultiplicity (1 : K) C + rootMultiplicity (1 : K) B := by
      exact rootMultiplicity_mul (mul_ne_zero hC hB)

/-- Multiplication by a genuine nonzero coefficient-field constant
does not change the integer X=1 order of a nonzero rational fraction.
The equality is derived by cross-multiplying genuine RatFunc elements. -/
theorem rootMultiplicity_gap_of_scaled_mk_eq
    (A B C D : Polynomial K) (gamma : K)
    (hA : A ≠ 0) (hB : B ≠ 0)
    (hC : C ≠ 0) (hD : D ≠ 0)
    (hgamma : gamma ≠ 0)
    (heq : RatFunc.mk A B = RatFunc.C gamma * RatFunc.mk C D) :
    (rootMultiplicity (1 : K) A : ℤ) -
        (rootMultiplicity (1 : K) B : ℤ) =
      (rootMultiplicity (1 : K) C : ℤ) -
        (rootMultiplicity (1 : K) D : ℤ) := by
  have hpolyGamma : (Polynomial.C gamma : Polynomial K) ≠ 0 := by
    simpa using hgamma
  have hscaled :
      RatFunc.C gamma * RatFunc.mk C D =
        RatFunc.mk (Polynomial.C gamma * C) D := by
    simp only [RatFunc.mk_eq_div, map_mul, RatFunc.algebraMap_C]
    ring
  have heq' :
      RatFunc.mk A B = RatFunc.mk (Polynomial.C gamma * C) D :=
    heq.trans hscaled
  have hbal := rootMultiplicity_balance_of_mk_eq
    A B (Polynomial.C gamma * C) D hA hB
    (mul_ne_zero hpolyGamma hC) hD heq'
  have hconst : rootMultiplicity (1 : K) (Polynomial.C gamma) = 0 := by
    apply rootMultiplicity_eq_zero
    intro hroot
    apply hgamma
    change (Polynomial.C gamma : Polynomial K).eval (1 : K) = 0 at hroot
    simpa only [eval_C] using hroot
  have hmul :
      rootMultiplicity (1 : K) (Polynomial.C gamma * C) =
        rootMultiplicity (1 : K) (Polynomial.C gamma) +
          rootMultiplicity (1 : K) C :=
    rootMultiplicity_mul (mul_ne_zero hpolyGamma hC)
  rw [hmul, hconst, zero_add] at hbal
  have hbalZ : (rootMultiplicity (1 : K) A : ℤ) +
      (rootMultiplicity (1 : K) D : ℤ) =
      (rootMultiplicity (1 : K) C : ℤ) +
      (rootMultiplicity (1 : K) B : ℤ) := by
    exact_mod_cast hbal
  omega

end ICFallDegrees.GrowthRatFuncRootMultiplicityCross
