import ICFallDegrees.GrowthRationalBezoutIntegerLift
import Mathlib.Tactic

/-!
# N8: a genuine rational Bézout identity yields an integral scalar identity

If a finite collection of integer polynomials generates the unit ideal
over the rationals, it contains a nonzero integer in the integral
polynomial ideal. The rational unit certificate itself is not assumed
for E6; it must still be constructed from the complex result.
-/

noncomputable section
namespace ICFallDegrees.GrowthRationalBezoutCertificate

open MvPolynomial
open ICFallDegrees.GrowthRationalBezoutIntegerLift

variable {σ : Type*} {n : ℕ}

/-- A finite explicit rational Bézout identity can be cleared to a
nonzero integer Bézout identity for the actual integer generators. -/
theorem integer_multiple_of_rational_bezout
    (F : Fin n → MvPolynomial σ ℤ)
    (H : Fin n → MvPolynomial σ ℚ)
    (hcert : (∑ i : Fin n,
      H i * (F i).map (Int.castRingHom ℚ)) = 1) :
    ∃ d : ℕ, 0 < d ∧ ∃ G : Fin n → MvPolynomial σ ℤ,
      (∑ i : Fin n, G i * F i) = C (d : ℤ) := by
  classical
  obtain ⟨d, hd, G, hG⟩ := clear_denominators_finite H
  refine ⟨d, hd, G, ?_⟩
  have hm : (∑ i : Fin n, G i * F i).map (Int.castRingHom ℚ) =
      (C (d : ℤ)).map (Int.castRingHom ℚ) := by
    calc
      (∑ i : Fin n, G i * F i).map (Int.castRingHom ℚ) =
          ∑ i : Fin n, (G i).map (Int.castRingHom ℚ) *
            (F i).map (Int.castRingHom ℚ) := by
              simp only [map_sum, map_mul]
      _ = ∑ i : Fin n, (C (d : ℚ) * H i) *
            (F i).map (Int.castRingHom ℚ) := by
              apply Finset.sum_congr rfl
              intro i _
              rw [hG i]
      _ = C (d : ℚ) * ∑ i : Fin n,
            H i * (F i).map (Int.castRingHom ℚ) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro i _
              ring
      _ = (C (d : ℤ)).map (Int.castRingHom ℚ) := by
              rw [hcert, mul_one]
              simp
  exact MvPolynomial.map_injective (Int.castRingHom ℚ)
    (by exact Int.cast_injective) hm

#print axioms integer_multiple_of_rational_bezout

end ICFallDegrees.GrowthRationalBezoutCertificate
