import ICFallDegrees.GrowthFixedE6IntegralInputs
import ICFallDegrees.GrowthRationalE6Unit
import ICFallDegrees.GrowthRationalBezoutCertificate
import ICFallDegrees.GrowthIntegralBezoutSpecialization
import Mathlib.Tactic

/-!
# Actual fixed E6 integer Bézout scalar from the complex unit theorem

The source unit theorem is genuinely about the normalized E6 chain in
characteristic zero, not an assumed rational certificate. Its Q-cofactor
vector is obtained by the field-linear retraction, then the common
rational denominator theorem produces a nonzero integer scalar identity
for the *exact* integer E6 input family.

This establishes the unit-ideal half for any characteristic not dividing
one fixed positive integer depending on k and depth. Finite-field
bounded-closure nonmembership is an independent, still-open Lean bridge.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6IntegralBezout

open MvPolynomial
open ICFallDegrees.GrowthFixedE6IntegralInputs
open ICFallDegrees.GrowthRationalBezoutCertificate
open ICFallDegrees.GrowthRationalE6Unit
open ICFallDegrees.GrowthIntegralBezoutSpecialization

/-- For every fixed depth, the actual normalized E6 chain contains an
explicit (existential) positive integer in its integral ideal. No
rational Bézout coefficient is assumed: it is extracted from the already
proven complex E6 unit ideal and coefficient retraction. -/
theorem exists_integer_e6_bezout (k : ℕ) (depth : 2 ≤ k) :
    ∃ d : ℕ, 0 < d ∧
      ∃ G : Fin (2 * k + 1) → MvPolynomial (Fin (2 * k)) ℤ,
        (∑ i : Fin (2 * k + 1), G i * e6InputsInt k depth i) =
          C (d : ℤ) := by
  classical
  let F : Fin (2 * k + 1) → MvPolynomial (Fin (2 * k)) ℚ :=
    ICFallDegrees.PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1
  have hunit : (1 : MvPolynomial (Fin (2 * k)) ℚ) ∈
      Ideal.span (Set.range F) := by
    exact rational_e6_unit k depth
  obtain ⟨H, hH⟩ := finite_span_explicit F hunit
  apply integer_multiple_of_rational_bezout (e6InputsInt k depth) H
  simpa only [F, map_e6_inputs] using hH.symm

/-- An actual integer E6 Bézout identity specializes to the unit ideal
over every field in which its scalar remains nonzero. -/
theorem exists_e6_unit_specialization_scalar
    (k : ℕ) (depth : 2 ≤ k) :
    ∃ d : ℕ, 0 < d ∧
      ∀ (K : Type*) [Field K], (d : K) ≠ 0 →
        (1 : MvPolynomial (Fin (2 * k)) K) ∈
          Ideal.span (Set.range (fun i : Fin (2 * k + 1) =>
            (e6InputsInt k depth i).map (Int.castRingHom K))) := by
  obtain ⟨d, hd, G, hG⟩ := exists_integer_e6_bezout k depth
  refine ⟨d, hd, ?_⟩
  intro K _ hnz
  apply unit_of_integral_bezout (e6InputsInt k depth) G (d : ℤ) hG
  simpa using hnz

#print axioms exists_integer_e6_bezout
#print axioms exists_e6_unit_specialization_scalar

end ICFallDegrees.GrowthFixedE6IntegralBezout
