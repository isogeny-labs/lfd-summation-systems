import ICFallDegrees.GrowthFieldCoefficientRetraction
import ICFallDegrees.GrowthFixedE6Existence
import ICFallDegrees.ScalarExtensionTransport
import ICFallDegrees.SelectorExtension
import Mathlib.Tactic

/-!
# N8: descend an actual complex E6 unit ideal to Q

A field extension admits a linear coefficient retraction. Because this
projection commutes with multiplication by lifted base-field polynomials,
an actual complex Bezout identity with rational generator coefficients
descends to a rational Bezout identity.

This proves the *rational unit-ideal half*, not any prime-field statement.
An independent integral scalar clearing lemma and finite-prime closure
specialization remain necessary.
-/

noncomputable section
namespace ICFallDegrees.GrowthRationalE6Unit

open MvPolynomial
open ICFallDegrees.GrowthFieldCoefficientRetraction
open ICFallDegrees.PrimeM2General

variable {R : Type*} [CommRing R] {n : ℕ}

/-- Actual ideal membership in a finite family has a finite indexed
polynomial Bézout coefficient vector. -/
theorem finite_span_explicit (F : Fin n → R) {z : R}
    (hz : z ∈ Ideal.span (Set.range F)) :
    ∃ G : Fin n → R, z = ∑ i : Fin n, G i * F i := by
  classical
  induction hz using Submodule.span_induction with
  | mem x hx =>
      obtain ⟨j, rfl⟩ := hx
      refine ⟨fun i => if i = j then 1 else 0, ?_⟩
      simp
  | zero =>
      exact ⟨fun _ => 0, by simp⟩
  | add x y _ _ hx hy =>
      obtain ⟨G, hG⟩ := hx
      obtain ⟨H, hH⟩ := hy
      refine ⟨fun i => G i + H i, ?_⟩
      simp only [add_mul, Finset.sum_add_distrib]
      rw [hG, hH]
  | smul a x _ hx =>
      obtain ⟨G, hG⟩ := hx
      refine ⟨fun i => a * G i, ?_⟩
      calc
        a • x = a * (∑ i : Fin n, G i * F i) := by
          simp only [smul_eq_mul, hG]
        _ = ∑ i : Fin n, (a * G i) * F i := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro i _
          ring

section Field

variable {K L σ : Type*} [Field K] [Field L] [Algebra K L]

/-- Project an explicit complex-field Bézout relation onto the base field.
All coefficients are derived from a *constructed* linear retraction, not
an assumed algebra map or generic multiplicativity of projection. -/
theorem descend_explicit_bezout
    (F : Fin n → MvPolynomial σ K)
    (H : Fin n → MvPolynomial σ L)
    (h : (∑ i : Fin n,
      H i * (F i).map (algebraMap K L)) = 1) :
    ∃ G : Fin n → MvPolynomial σ K,
      (∑ i : Fin n, G i * F i) = 1 := by
  let φ : L →ₗ[K] K := scalarRetraction (K := K) (L := L)
  refine ⟨fun i => coefficientProjection φ (H i), ?_⟩
  have hproject := congrArg (coefficientProjection φ) h
  have hone : coefficientProjection φ (1 : MvPolynomial σ L) = 1 := by
    have hunit := rationalRetraction_lift (K := K) (L := L)
      (σ := σ) (1 : MvPolynomial σ K)
    simpa only [map_one] using hunit
  rw [map_sum] at hproject
  simp only [coefficientProjection_mul_lift, hone] at hproject
  exact hproject

/-- The contraction of a finitely generated base-field polynomial ideal
from an arbitrary field extension contains the original ideal unit.

The output has no unit-ideal hypothesis: the complex one is consumed. -/
theorem unit_descends_from_field_extension
    (F : Fin n → MvPolynomial σ K)
    (h : (1 : MvPolynomial σ L) ∈
      Ideal.span (Set.range
        (fun i : Fin n => (F i).map (algebraMap K L)))) :
    (1 : MvPolynomial σ K) ∈ Ideal.span (Set.range F) := by
  obtain ⟨H, hH⟩ := finite_span_explicit
    (fun i : Fin n => (F i).map (algebraMap K L)) h
  obtain ⟨G, hG⟩ := descend_explicit_bezout F H hH.symm
  rw [← hG]
  apply Submodule.sum_mem
  intro i _
  exact (Ideal.span (Set.range F)).mul_mem_left (G i)
    (Ideal.subset_span (show F i ∈ Set.range F from ⟨i, rfl⟩))

end Field

/-- Actual rational E6 unit ideal, with no supplied rational certificate:
it follows from the existing *unconditional* complex E6 theorem and
the general coefficient-retraction descent. -/
theorem rational_e6_unit (k : ℕ) (depth : 2 ≤ k) :
    (1 : MvPolynomial (Fin (2 * k)) ℚ) ∈
      Ideal.span (Set.range
        (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1)) := by
  have hc := ICFallDegrees.GrowthFixedE6.unit k depth
  have hset :
      ScalarExtensionTransport.lift (K := ℂ) ''
          Set.range (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1) =
        Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1) := by
    simpa using (SelectorExtension.image_inputs (K := ℚ) (L := ℂ)
      k depth (6 : ℚ) (-7 : ℚ) (1 : ℚ))
  rw [← hset] at hc
  apply unit_descends_from_field_extension
    (K := ℚ) (L := ℂ)
    (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1)
  have hrange :
      ScalarExtensionTransport.lift (K := ℂ) ''
          Set.range (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1) =
        Set.range (fun i =>
          (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1 i).map
            (algebraMap ℚ ℂ)) := by
    ext q
    constructor
    · rintro ⟨_, ⟨i, rfl⟩, rfl⟩
      exact ⟨i, rfl⟩
    · rintro ⟨i, rfl⟩
      exact ⟨_, ⟨i, rfl⟩, rfl⟩
  rw [hrange] at hc
  exact hc

#print axioms finite_span_explicit
#print axioms unit_descends_from_field_extension
#print axioms rational_e6_unit

end ICFallDegrees.GrowthRationalE6Unit
