import ICFallDegrees.GrowthRationalCoordinates
import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum

/-!
# N8: actual rational-function-field coordinate model

Construct the two coordinate functions in the actual field K(X), not in
an arbitrary field with assumed x-genericity. The nonzero and on-curve
facts are algebraic; the intermediate field F^p and spread bound are NOT
part of this module. In particular this is not a CurveModel instance.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncModel

open ICFallDegrees.GrowthRationalCoordinates
open ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- The rational-function indeterminate differs from every scalar. -/
theorem X_ne_constant (c : K) :
    (RatFunc.X : RatFunc K) ≠ RatFunc.C c := by
  intro h
  have hpoly : (Polynomial.X : Polynomial K) = Polynomial.C c := by
    apply RatFunc.algebraMap_injective K
    simpa only [RatFunc.algebraMap_X, RatFunc.algebraMap_C] using h
  have hdegree := congrArg Polynomial.natDegree hpoly
  norm_num at hdegree

/-- The rational-function indeterminate is not the target t=1. -/
theorem X_ne_one : (RatFunc.X : RatFunc K) ≠ 1 := by
  simpa using X_ne_constant (1 : K)

/-- The second coordinate numerator is a nonzero rational function. -/
theorem numerator_ne_zero (a : K) :
    (RatFunc.X : RatFunc K) + 2 + RatFunc.C a ≠ 0 := by
  intro hz
  have hconst : (RatFunc.X : RatFunc K) = RatFunc.C (-2 - a) := by
    simp only [map_sub, map_neg, map_ofNat]
    linear_combination hz
  exact X_ne_constant (-2 - a) hconst

/-- The two actual rational functions for the two-torsion section. -/
def xi (a : K) : Fin 2 → RatFunc K :=
  fun i => if i = 0 then RatFunc.X else
    rationalZ (RatFunc.C a) (RatFunc.X : RatFunc K)

/-- The concrete algebraic portion of CurveModel for K(X): coordinate
nonvanishing and the exact Semaev relation over the original base field K.
The remaining missing hypothesis is SpreadBound over F^p. -/
theorem ratFunc_algebraic_data (a : K) :
    (∀ i : Fin 2, xi a i ≠ 0) ∧
      GrowthSpread.S3at a (-1 - a) (1 : K) (xi a) = 0 := by
  have data := rational_algebraic_data (RatFunc.C a) (RatFunc.X : RatFunc K)
    RatFunc.X_ne_zero X_ne_one (numerator_ne_zero a)
  constructor
  · simpa only [xi] using data.1
  · convert data.2 using 1 <;>
      simp [xi, GrowthSpread.S3at, RatFunc.algebraMap_eq_C, map_ofNat]

end ICFallDegrees.GrowthRatFuncModel
