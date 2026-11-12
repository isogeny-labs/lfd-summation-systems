import ICFallDegrees.RawPresentation
import Mathlib.Algebra.Field.ZMod
import Mathlib.Tactic.NormNum

noncomputable section

namespace ICFallDegrees.PrimeCounterexample

open MvPolynomial

abbrev K := ZMod 5
instance : Fact (Nat.Prime 5) := ⟨by decide⟩

/-- Raw S3 at target (0,0), on E/F_5: y^2=x^3+2x. -/
def curveS3 : QuadraticFamily.P K :=
  ((QuadraticFamily.x : QuadraticFamily.P K)*QuadraticFamily.y-C (2 : K))^2

theorem nonsingular_discriminant : 4*(2 : K)^3+27*(0 : K)^2 ≠ 0 := by
  decide

theorem target_on_curve : (0 : K)^2=(0 : K)^3+2*(0 : K) := by simp

theorem raw_is_curveS3 : QuadraticFamily.raw (0 : K)=curveS3 := by
  have identity : curveS3-QuadraticFamily.raw (0 : K)=C (5 : K)*(1-QuadraticFamily.x*QuadraticFamily.y) := by
    rw [QuadraticFamily.raw_form]
    simp only [curveS3,QuadraticFamily.h,map_zero,sub_zero,map_ofNat]
    ring
  have five : C (5 : K)=(0 : QuadraticFamily.P K) := by
    have hfive : (5 : K)=0 := by decide
    rw [hfive,map_zero]
  have hz : curveS3-QuadraticFamily.raw (0 : K)=0 := by rw [identity,five,zero_mul]
  exact (sub_eq_zero.mp hz).symm

theorem normalized_lfd_three : IsLastFallDegree (QuadraticFamily.inputs (0 : K)) 3 :=
  QuadraticFamily.exact_lfd_three 0 (by norm_num)

theorem raw_lfd_four : IsLastFallDegree (QuadraticFamily.rawInputs (0 : K)) 4 :=
  QuadraticFamily.raw_exact_lfd_four 0 (by norm_num)

theorem presentations_same_ideal :
    Ideal.span (QuadraticFamily.inputs (0 : K))=Ideal.span (QuadraticFamily.rawInputs (0 : K)) := QuadraticFamily.same_ideal 0

/-- Explicit prime-field counterexample to presentation invariance of LFD.
This is not the larger FFD=4/LFD=5 counterexample, whose Lean port is pending. -/
theorem lfd_is_not_ideal_invariant :
    ∃ F G : Set (QuadraticFamily.P K), Ideal.span F=Ideal.span G ∧
      IsLastFallDegree F 3 ∧ IsLastFallDegree G 4 :=
  ⟨QuadraticFamily.inputs (0 : K),QuadraticFamily.rawInputs (0 : K),
    presentations_same_ideal,normalized_lfd_three,raw_lfd_four⟩

end ICFallDegrees.PrimeCounterexample
