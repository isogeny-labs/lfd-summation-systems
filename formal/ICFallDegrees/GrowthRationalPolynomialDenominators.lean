import ICFallDegrees.GrowthIntegralBezoutSpecialization
import Mathlib.Tactic

/-!
# N8 rational polynomial denominator clearing

Every finite rational-coefficient multivariate polynomial is an integral
polynomial after multiplication by a nonzero positive integer.
This is the missing algebraic ingredient for lifting rational Bézout
identities to integer scalar identities. It does not imply a rational
Bézout certificate exists for the actual E6 family.
-/

noncomputable section
namespace ICFallDegrees.GrowthRationalPolynomialDenominators

open MvPolynomial

variable {σ : Type*}

/-- Every rational multivariate polynomial has a common positive integral
denominator. No fixed degree, number of variables or coefficient support
is assumed. -/
theorem clear_denominator (f : MvPolynomial σ ℚ) :
    ∃ d : ℕ, 0 < d ∧ ∃ g : MvPolynomial σ ℤ,
      C (d : ℚ) * f = g.map (Int.castRingHom ℚ) := by
  induction f using MvPolynomial.induction_on with
  | C a =>
      refine ⟨a.den, Rat.pos a, C a.num, ?_⟩
      have hd : (a.den : ℚ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt (Rat.pos a))
      have hnum : (a.den : ℚ) * a = (a.num : ℚ) := by
        calc
          (a.den : ℚ) * a =
              (a.den : ℚ) * ((a.num : ℚ) / (a.den : ℚ)) := by
                rw [Rat.num_div_den a]
          _ = (a.num : ℚ) := by field_simp [hd]
      simpa only [map_C, map_intCast, map_natCast, C_mul] using
        congrArg (C : ℚ →+* MvPolynomial σ ℚ) hnum
  | add p q ihp ihq =>
      obtain ⟨d, hd, g, hg⟩ := ihp
      obtain ⟨e, he, h, hh⟩ := ihq
      refine ⟨d * e, Nat.mul_pos hd he,
        C (e : ℤ) * g + C (d : ℤ) * h, ?_⟩
      calc
        C ((d * e : ℕ) : ℚ) * (p + q) =
            C (e : ℚ) * (C (d : ℚ) * p) +
            C (d : ℚ) * (C (e : ℚ) * q) := by
              simp only [Nat.cast_mul, map_mul]
              ring
        _ = C (e : ℚ) * g.map (Int.castRingHom ℚ) +
            C (d : ℚ) * h.map (Int.castRingHom ℚ) := by rw [hg, hh] <;> rfl
        _ = (C (e : ℤ) * g + C (d : ℤ) * h).map (Int.castRingHom ℚ) := by
              simp only [map_add, map_mul, map_C, map_intCast, map_natCast] <;> rfl
      all_goals rfl
  | mul_X p i ih =>
      obtain ⟨d, hd, g, hg⟩ := ih
      refine ⟨d, hd, g * X i, ?_⟩
      calc
        C (d : ℚ) * (p * X i) =
            (C (d : ℚ) * p) * X i := by ring
        _ = g.map (Int.castRingHom ℚ) * X i := by rw [hg] <;> rfl
        _ = (g * X i).map (Int.castRingHom ℚ) := by
              simp only [map_mul, map_X] <;> rfl
      all_goals rfl

#print axioms clear_denominator

end ICFallDegrees.GrowthRationalPolynomialDenominators
