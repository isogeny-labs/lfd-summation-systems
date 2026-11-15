import ICFallDegrees.GrowthRatFuncModel
import ICFallDegrees.GrowthActualReplacementSpread
import Mathlib.Tactic

/-!
# N8: actual two-coordinate monomial identification

This file transports the genuine rational-monomial results into the
`GrowthSpread` index type. It does not assert a circuit spread bound.
-/

noncomputable section

namespace ICFallDegrees.GrowthCircuitInterface

open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K]

/-- The two coordinates of a Fin 2 finitely supported exponent vector
determine that exponent vector. -/
theorem exponentPair_injective :
    Function.Injective (fun e : Fin 2 →₀ ℕ => (e 0, e 1)) := by
  intro e f h
  change (e 0, e 1) = (f 0, f 1) at h
  ext i
  fin_cases i
  · exact congrArg Prod.fst h
  · exact congrArg Prod.snd h

/-- In the actual rational section, the two-coordinate monomial used by
the growth argument is exactly the rational monomial used by the audited
Wronskian and cofactor proofs. -/
theorem mono_xi_eq_monomial (a : K) (e : Fin 2 →₀ ℕ) :
    ICFallDegrees.GrowthSpread.mono (xi a) e =
      monomial (a + 2) (e 0, e 1) := by
  change (RatFunc.X : RatFunc K) ^ e 0 *
      ((RatFunc.X + 2 + RatFunc.C a) / (RatFunc.X - 1)) ^ e 1 =
      RatFunc.mk
        (Polynomial.X ^ e 0 * (Polynomial.X + Polynomial.C (a + 2)) ^ e 1)
        ((Polynomial.X - 1 : Polynomial K) ^ e 1)
  rw [RatFunc.mk_eq_div]
  simp only [map_mul, map_pow, map_add, map_sub, map_one,
    RatFunc.algebraMap_X, RatFunc.algebraMap_C, map_ofNat, div_pow]
  ring

end ICFallDegrees.GrowthCircuitInterface
