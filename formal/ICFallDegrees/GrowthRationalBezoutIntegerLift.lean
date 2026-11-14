import ICFallDegrees.GrowthRationalPolynomialDenominators
import Mathlib.Tactic

/-!
# Finite rational Bézout coefficients admit one common integer denominator

This is a purely algebraic, finite-support interface. It does not deduce
a rational Bézout identity for the fixed E6 curve from the complex theorem;
that faithful scalar-extension descent is an independent open Lean gap.
-/

noncomputable section
namespace ICFallDegrees.GrowthRationalBezoutIntegerLift

open MvPolynomial
open ICFallDegrees.GrowthRationalPolynomialDenominators

variable {σ ι : Type*}

/-- A single integer denominator clears all rational polynomial coefficients
in an arbitrary finite index set. -/
theorem clear_denominators_on_finset (s : Finset ι)
    (f : ι → MvPolynomial σ ℚ) :
    ∃ d : ℕ, 0 < d ∧ ∃ g : ι → MvPolynomial σ ℤ,
      ∀ i ∈ s, C (d : ℚ) * f i = (g i).map (Int.castRingHom ℚ) := by
  classical
  induction s using Finset.induction_on with
  | empty =>
      refine ⟨1, by omega, fun _ => 0, ?_⟩
      intro i hi
      simp at hi
  | @insert a s ha ih =>
      obtain ⟨d, hd, g, hg⟩ := ih
      obtain ⟨e, he, ga, hga⟩ := clear_denominator (f a)
      refine ⟨d * e, Nat.mul_pos hd he,
        (fun i => if i = a then C (d : ℤ) * ga else C (e : ℤ) * g i), ?_⟩
      intro i hi
      by_cases hia : i = a
      · subst i
        simp only [ite_true]
        calc
          C (((d * e : ℕ) : ℚ)) * f a =
              C (d : ℚ) * (C (e : ℚ) * f a) := by
                simp only [Nat.cast_mul, map_mul]
                ring
          _ = C (d : ℚ) * ga.map (Int.castRingHom ℚ) := by
                rw [hga]
          _ = (C (d : ℤ) * ga).map (Int.castRingHom ℚ) := by
                simp only [map_mul, map_C, map_intCast, map_natCast] <;> rfl
        all_goals rfl
      · have his : i ∈ s :=
          (Finset.mem_insert.mp hi).resolve_left hia
        simp only [if_neg hia]
        calc
          C (((d * e : ℕ) : ℚ)) * f i =
              C (e : ℚ) * (C (d : ℚ) * f i) := by
                simp only [Nat.cast_mul, map_mul]
                ring
          _ = C (e : ℚ) * (g i).map (Int.castRingHom ℚ) := by
                rw [hg i his]
          _ = (C (e : ℤ) * g i).map (Int.castRingHom ℚ) := by
                simp only [map_mul, map_C, map_intCast, map_natCast] <;> rfl
        all_goals rfl

/-- A finite list of rational polynomial coefficients may be cleared by
one positive integer. -/
theorem clear_denominators_finite [Fintype ι] (f : ι → MvPolynomial σ ℚ) :
    ∃ d : ℕ, 0 < d ∧ ∃ g : ι → MvPolynomial σ ℤ,
      ∀ i, C (d : ℚ) * f i = (g i).map (Int.castRingHom ℚ) := by
  classical
  obtain ⟨d, hd, g, hg⟩ := clear_denominators_on_finset Finset.univ f
  exact ⟨d, hd, g, fun i => hg i (Finset.mem_univ i)⟩

#print axioms clear_denominators_on_finset
#print axioms clear_denominators_finite

end ICFallDegrees.GrowthRationalBezoutIntegerLift
