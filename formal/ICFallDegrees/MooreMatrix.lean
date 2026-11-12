import Mathlib.FieldTheory.Finite.Basic
import Mathlib.LinearAlgebra.LinearIndependent.Lemmas
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

/-! Invertibility of the binary Moore matrix from an actual finite-field
basis. Dedekind independence of field homomorphisms and the proved exact
order of finite-field Frobenius supply the rank; no determinant certificate
or assumed matrix inverse is used. -/
noncomputable section
namespace ICFallDegrees.MooreMatrix
open scoped BigOperators
variable {K : Type*} [Field K] [Finite K] [Algebra (ZMod 2) K] {n : ℕ}

def matrix (basis : Basis (Fin n) (ZMod 2) K) : Matrix (Fin n) (Fin n) K :=
  fun i j => (basis j)^(2^(i : ℕ))

theorem rows_independent (basis : Basis (Fin n) (ZMod 2) K) :
    LinearIndependent K (matrix basis).row := by
  classical
  have dimension : Module.finrank (ZMod 2) K=n := by
    simpa using Module.finrank_eq_card_basis basis
  let f : Fin n → (K →ₐ[ZMod 2] K) := fun i => FiniteField.frobeniusAlgHom (ZMod 2) K^(i : ℕ)
  have distinct : Function.Injective f := by
    intro i j same
    let a : Fin (Module.finrank (ZMod 2) K) := ⟨i.val,by simpa only [dimension] using i.isLt⟩
    let b : Fin (Module.finrank (ZMod 2) K) := ⟨j.val,by simpa only [dimension] using j.isLt⟩
    have equal : a=b := (FiniteField.bijective_frobeniusAlgHom_pow (ZMod 2) K).injective same
    exact Fin.ext (congrArg (fun t : Fin (Module.finrank (ZMod 2) K) => t.val) equal)
  have homs : LinearIndependent K (fun i : Fin n => (f i).toLinearMap) :=
    (linearIndependent_algHom_toLinearMap (ZMod 2) K K).comp f distinct
  let coordinates : (K →ₗ[ZMod 2] K) ≃ₗ[K] (Fin n → K) := (basis.constr K).symm
  have mapped := homs.map' coordinates.toLinearMap (LinearMap.ker_eq_bot.mpr coordinates.injective)
  convert mapped using 1
  ext i j
  change (basis j)^(2^(i : ℕ))=(FiniteField.frobeniusAlgHom (ZMod 2) K^(i : ℕ)) (basis j)
  simp only [AlgHom.coe_pow,FiniteField.coe_frobeniusAlgHom,ZMod.card,pow_iterate]

theorem isUnit_matrix (basis : Basis (Fin n) (ZMod 2) K) : IsUnit (matrix basis) :=
  Matrix.linearIndependent_rows_iff_isUnit.mp (rows_independent basis)

theorem determinant_nonzero (basis : Basis (Fin n) (ZMod 2) K) : (matrix basis).det≠0 := by
  exact ((Matrix.isUnit_iff_isUnit_det _).mp (isUnit_matrix basis)).ne_zero

theorem inverse_identities (basis : Basis (Fin n) (ZMod 2) K) :
    matrix basis*(matrix basis)⁻¹=1 ∧ (matrix basis)⁻¹*matrix basis=1 := by
  have h := (Matrix.isUnit_iff_isUnit_det _).mp (isUnit_matrix basis)
  exact ⟨Matrix.mul_nonsing_inv _ h,Matrix.nonsing_inv_mul _ h⟩

end ICFallDegrees.MooreMatrix
