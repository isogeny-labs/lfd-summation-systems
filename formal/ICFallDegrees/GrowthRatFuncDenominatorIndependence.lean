import ICFallDegrees.GrowthRatFuncCommonDenominator
import Mathlib.LinearAlgebra.LinearIndependent.Basic
import Mathlib.Tactic

/-!
# N8: common denominator clearing preserves linear independence

The exact common denominator H of the previous theorem gives polynomial
numerators P_j with (algebraMap H) f_j = algebraMap P_j. Over a field,
multiplication by the nonzero H and the injective polynomial embedding
both preserve K-linear independence. We prove the implication explicitly
using finite linear combinations, without postulating any nonzero Wronskian.

This is the last algebraic input before the characteristic-zero
polynomial Wronskian nonvanishing theorem can be applied.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncDenominatorIndependence

open ICFallDegrees.GrowthRatFuncCommonDenominator

variable {K : Type*} [Field K]

/-- An independent rational family remains independent after the
specified common polynomial denominator is cleared. -/
theorem polynomial_linearIndependent_of_clear
    (r : ℕ) (f : Fin r → RatFunc K)
    (H : Polynomial K) (p : Fin r → Polynomial K)
    (hH : H ≠ 0)
    (hp : ∀ j : Fin r,
      algebraMap (Polynomial K) (RatFunc K) H * f j =
        algebraMap (Polynomial K) (RatFunc K) (p j))
    (hf : LinearIndependent K f) :
    LinearIndependent K p := by
  classical
  have hHmap : algebraMap (Polynomial K) (RatFunc K) H ≠ 0 :=
    RatFunc.algebraMap_ne_zero hH
  apply Fintype.linearIndependent_iff.mpr
  intro c hc j
  apply (Fintype.linearIndependent_iff.mp hf) c ?_ j
  have hzero :
      algebraMap (Polynomial K) (RatFunc K) H *
        (∑ i : Fin r, c i • f i) = 0 := by
    calc
      algebraMap (Polynomial K) (RatFunc K) H *
          (∑ i : Fin r, c i • f i)
          = ∑ i : Fin r, c i •
              (algebraMap (Polynomial K) (RatFunc K) H * f i) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro i hi
              simp only [Algebra.smul_def]
              ring
      _ = ∑ i : Fin r, c i •
            (algebraMap (Polynomial K) (RatFunc K) (p i)) := by
          apply Finset.sum_congr rfl
          intro i hi
          rw [hp i]
      _ = algebraMap (Polynomial K) (RatFunc K)
            (∑ i : Fin r, c i • p i) := by
          rw [map_sum]
          apply Finset.sum_congr rfl
          intro i hi
          simp only [Algebra.smul_def, map_mul]
          rw [IsScalarTower.algebraMap_apply K (Polynomial K) (RatFunc K)]
      _ = 0 := by rw [hc, map_zero]
  exact (mul_eq_zero.mp hzero).resolve_left hHmap

/-- Every finite rational family independent over K admits one
nonzero common polynomial denominator and an independent polynomial
numerator family, simultaneously. -/
theorem exists_independent_polynomial_representatives
    (r : ℕ) (f : Fin r → RatFunc K)
    (hf : LinearIndependent K f) :
    ∃ (H : Polynomial K) (p : Fin r → Polynomial K),
      H ≠ 0 ∧
      (∀ j, algebraMap (Polynomial K) (RatFunc K) H * f j =
        algebraMap (Polynomial K) (RatFunc K) (p j)) ∧
      LinearIndependent K p := by
  obtain ⟨p, hp⟩ := exists_polynomial_numerators r f
  refine ⟨commonDenominator r f, p, commonDenominator_ne_zero r f, hp, ?_⟩
  exact polynomial_linearIndependent_of_clear r f
    (commonDenominator r f) p (commonDenominator_ne_zero r f) hp hf

end ICFallDegrees.GrowthRatFuncDenominatorIndependence
