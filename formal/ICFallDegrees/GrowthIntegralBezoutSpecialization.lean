import ICFallDegrees.GrowthPrimeRationalSeed
import Mathlib.Tactic

/-!
# N8 integral Bézout specialization

A concrete identity with integer polynomial coefficients survives in any
field provided its constant coefficient remains nonzero after reduction.
This is a genuine ideal-membership theorem, but it does *not* construct
an integer Bézout certificate from the complex E6 unit theorem and does
not prove prime-field bounded-closure nonmembership.
-/

noncomputable section
namespace ICFallDegrees.GrowthIntegralBezoutSpecialization

open MvPolynomial

variable {σ : Type*} {n : ℕ}
variable {K : Type*} [Field K]

/-- An integer polynomial Bézout identity with nonzero constant c gives
unit-ideal membership after any field-valued coefficient specialization.

Unlike an assumed prime-field unit ideal, the hypothesis is an explicit
integer identity; producing it from the E6 complex certificate is a
separate remaining task. -/
theorem unit_of_integral_bezout
    (F G : Fin n → MvPolynomial σ ℤ) (c : ℤ)
    (hcert : (∑ i : Fin n, G i * F i) = C c)
    (hc : (c : K) ≠ 0) :
    (1 : MvPolynomial σ K) ∈
      Ideal.span (Set.range (fun i : Fin n =>
        MvPolynomial.map (Int.castRingHom K) (F i))) := by
  classical
  let S : Set (MvPolynomial σ K) :=
    Set.range (fun i : Fin n => MvPolynomial.map (Int.castRingHom K) (F i))
  have hmap :
      (∑ i : Fin n,
        MvPolynomial.map (Int.castRingHom K) (G i) *
        MvPolynomial.map (Int.castRingHom K) (F i)) =
      C (c : K) := by
    simpa only [map_sum, map_mul, map_C] using
      congrArg (MvPolynomial.map (Int.castRingHom K)) hcert
  have hmem : C (c : K) ∈ Ideal.span S := by
    rw [← hmap]
    apply Submodule.sum_mem
    intro i hi
    exact (Ideal.span S).mul_mem_left _
      (Ideal.subset_span (show
        MvPolynomial.map (Int.castRingHom K) (F i) ∈ S from ⟨i, rfl⟩))
  have hscaled := (Ideal.span S).mul_mem_left (C ((c : K)⁻¹)) hmem
  have hunit : C ((c : K)⁻¹) * C (c : K) = (1 : MvPolynomial σ K) := by
    rw [← map_mul, inv_mul_cancel₀ hc, map_one]
  simpa only [hunit] using hscaled

end ICFallDegrees.GrowthIntegralBezoutSpecialization
