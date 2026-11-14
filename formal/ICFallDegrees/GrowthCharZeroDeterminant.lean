import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# N8: finite-free determinant reduction criterion

The group-algebra multiplier U becomes identity modulo two.  Once a
specific multiplication matrix on the residue monomial basis has been
constructed over Z, the following lemma turns this congruence into
a nonzero determinant over Z and hence into rational invertibility.

Proof candidate until the pinned Lean CI compile and axiom audit pass.
-/

namespace ICFallDegrees.GrowthCharZeroDeterminant

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Reducing a square integer matrix to the identity over F₂ forces
its determinant to reduce to one. -/
theorem det_cast_two_eq_one (M : Matrix ι ι ℤ)
    (hM : M.map (Int.castRingHom (ZMod 2)) = 1) :
    (M.det : ZMod 2) = 1 := by
  calc
    (M.det : ZMod 2) = (M.map (Int.castRingHom (ZMod 2))).det := by
      exact (Int.castRingHom (ZMod 2)).map_det M
    _ = 1 := by rw [hM]; simp

/-- Therefore the integer determinant is nonzero. -/
theorem det_ne_zero_of_mod_two_one (M : Matrix ι ι ℤ)
    (hM : M.map (Int.castRingHom (ZMod 2)) = 1) :
    M.det ≠ 0 := by
  intro h
  have hdet := det_cast_two_eq_one M hM
  rw [h] at hdet
  norm_num at hdet


/-- No explicit determinant computation is required: any matrix of the
form 2*T-I is invertible after extension of scalars to Q.  This is the
matrix shape of the integral U-multiplication operator. -/
theorem det_ne_zero_two_mul_sub_one (T : Matrix ι ι ℤ) :
    (2 • T - 1).det ≠ 0 := by
  apply det_ne_zero_of_mod_two_one
  ext i j
  simp only [Matrix.map_apply, Matrix.sub_apply, Matrix.smul_apply]
  have h2 : (2 : ZMod 2) = 0 := by decide
  have hneg : (-1 : ZMod 2) = 1 := by decide
  by_cases hij : i = j
  · simp [Matrix.one_apply, smul_eq_mul, hij, h2, hneg]
  · simp [Matrix.one_apply, smul_eq_mul, hij, h2, hneg]


/-- The integer matrix of the N8 characteristic-zero multiplier is a
linear combination of three translation permutation matrices and the
identity. No assumption on these matrices is needed for nonvanishing. -/
def multiplierMatrix (A B C : Matrix ι ι ℤ) : Matrix ι ι ℤ :=
  2 • A - 2 • B - 2 • C - 5 • (1 : Matrix ι ι ℤ)

/-- An exact additive decomposition that exposes the odd determinant. -/
theorem multiplierMatrix_eq (A B C : Matrix ι ι ℤ) :
    multiplierMatrix A B C =
      2 • (A - B - C - 2 • (1 : Matrix ι ι ℤ)) - 1 := by
  unfold multiplierMatrix
  abel

/-- Any integer matrix with the coefficient pattern of multiplication by
U=2xz-2z-2x-5 has nonzero determinant. This is uniform in the index set. -/
theorem multiplierMatrix_det_ne_zero (A B C : Matrix ι ι ℤ) :
    (multiplierMatrix A B C).det ≠ 0 := by
  rw [multiplierMatrix_eq]
  exact det_ne_zero_two_mul_sub_one _

end ICFallDegrees.GrowthCharZeroDeterminant
