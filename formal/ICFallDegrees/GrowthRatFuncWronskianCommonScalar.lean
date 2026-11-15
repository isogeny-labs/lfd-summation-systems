import ICFallDegrees.GrowthRatFuncScalarJetTriangular
import ICFallDegrees.GrowthRatFuncWronskianNormalization
import Mathlib.Tactic

/-!
# N8: common rational scalar factor in every Wronskian column

Prove the *actual* derivative matrix identity
 W(hf) = T(h) W(f)
using the already verified arbitrary-order rational Leibniz rule and
the independently verified lower triangular scalar jet T(h).
Its determinant gives det W(hf) = h^r det W(f).

This is a genuine polynomial/rational common-denominator bridge,
not a Wronskian nonvanishing assumption or an assumed SpreadBound.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncWronskianCommonScalar

open Finset
open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncHigherLeibniz
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncScalarJetTriangular

variable {K : Type*} [Field K]

/-- A finite sum indexed by the first i+1 natural numbers is exactly
the sum over Fin r with a cutoff at i, whenever i : Fin r. -/
private theorem sum_fin_cutoff (r : ℕ) (i : Fin r) (g : ℕ → RatFunc K) :
    (∑ j : Fin r, if j.val ≤ i.val then g j.val else 0) =
      ∑ k ∈ Finset.range (i.val + 1), g k := by
  classical
  rw [← Finset.sum_filter]
  symm
  apply Finset.sum_bij (fun k hk =>
    (⟨k, by
      have hk' : k ≤ i.val := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
      exact lt_of_le_of_lt hk' i.isLt⟩ : Fin r))
  · intro k hk
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
  · intro k hk l hl heq
    exact Fin.ext_iff.mp heq
  · intro j hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
    refine ⟨j.val, Finset.mem_range.mpr (Nat.lt_succ_of_le hj), ?_⟩
    apply Fin.ext
    rfl
  · intro k hk
    rfl

/-- The genuine rational Wronskian derivative matrix transforms via
the lower triangular scalar-jet operator, for arbitrary r. -/
theorem derivativeMatrix_mul_scalar
    (r : ℕ) (h : RatFunc K) (f : Fin r → RatFunc K) :
    derivativeMatrix r (fun j => h * f j) =
      scalarJet r h * derivativeMatrix r f := by
  classical
  ext i j
  simp only [derivativeMatrix, Matrix.of_apply, Matrix.mul_apply]
  rw [iterDerivative_mul]
  simp only [nsmul_eq_mul]
  rw [← sum_fin_cutoff r i
    (fun k => (i.val.choose k : RatFunc K) *
      (iterDerivative (i.val - k) h * iterDerivative k (f j)))]
  apply Finset.sum_congr rfl
  intro k hk
  by_cases hki : k.val ≤ i.val
  · simp only [scalarJet, Matrix.of_apply, if_pos hki]
    ring
  · simp only [scalarJet, Matrix.of_apply, if_neg hki, zero_mul, ite_false]

/-- A common nonzero (or zero) rational scalar multiplies
an r-column genuine Wronskian determinant by h^r. -/
theorem determinant_mul_scalar
    (r : ℕ) (h : RatFunc K) (f : Fin r → RatFunc K) :
    (derivativeMatrix r (fun j => h * f j)).det =
      h ^ r * (derivativeMatrix r f).det := by
  rw [derivativeMatrix_mul_scalar, Matrix.det_mul, scalarJet_det]

/-- Multiplication by a nonzero rational h does not change whether
the genuine Wronskian vanishes, at any size. -/
theorem determinant_mul_scalar_ne_zero_iff
    (r : ℕ) (h : RatFunc K) (f : Fin r → RatFunc K) (hh : h ≠ 0) :
    (derivativeMatrix r (fun j => h * f j)).det ≠ 0 ↔
      (derivativeMatrix r f).det ≠ 0 := by
  rw [determinant_mul_scalar, mul_ne_zero_iff]
  simp [hh]

end ICFallDegrees.GrowthRatFuncWronskianCommonScalar
