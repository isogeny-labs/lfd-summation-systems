import ICFallDegrees.GrowthCharZeroComplexUnit
import Mathlib.Tactic

/-!
# A norm-separated open family of genuine N8 unit ideals over C

The parameter a is not fixed: for every a in C satisfying ||a|| > 5,
the two-torsion curve y^2=x^3+a*x-1-a has no two-summand
S3 root on any 2-power root grid, hence the full chain ideal is
the unit ideal for every k >= 2.

This result does NOT prove circuit spread, CurveModel or LFD growth.
The theorem uses the inherited pinned and audited grid Nullstellensatz.
-/

noncomputable section
namespace ICFallDegrees.GrowthCharZeroSeparatedFamily

open MvPolynomial PrimeM2General
open ICFallDegrees.GrowthSecant
open ICFallDegrees.GrowthRationalUnit
open ICFallDegrees.GrowthRationalChainUnit
open ICFallDegrees.GrowthCharZeroComplexUnit

/-- The reduction to a scalar a works uniformly on the unit-circle grid. -/
theorem s3_ne_zero_on_complex_root_grid_of_norm_gt
    (k : ℕ) (a x z : ℂ) (ha : 5 < ‖a‖)
    (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1) :
    s3 a (-1 - a) x z 1 ≠ 0 := by
  have hxn : ‖x‖ = 1 := norm_eq_one_of_root k x hx
  have hzn : ‖z‖ = 1 := norm_eq_one_of_root k z hz
  have hx1 : ‖x - 1‖ ≤ 2 := by
    calc
      ‖x - 1‖ ≤ ‖x‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ = 2 := by rw [hxn]; norm_num
  have hq : ‖(x - 1) * z - x - (2 : ℂ)‖ ≤ 5 := by
    calc
      ‖(x - 1) * z - x - (2 : ℂ)‖ ≤
          ‖(x - 1) * z - x‖ + ‖(2 : ℂ)‖ := norm_sub_le _ _
      _ ≤ (‖(x - 1) * z‖ + ‖x‖) + ‖(2 : ℂ)‖ := by
            gcongr
            exact norm_sub_le _ _
      _ = ‖x - 1‖ + 3 := by
        rw [norm_mul, hzn, hxn]
        norm_num <;> ring
      _ ≤ 5 := by linarith
  intro h
  have hp := parameter_of_s3_zero a x z h
  have hbound : ‖a‖ ≤ 5 := by
    rw [hp]
    exact hq
  linarith


/-- The chain unit ideal for every complex parameter outside the norm-five disk. -/
theorem unit_of_primitiveRoot_complex_of_norm_gt
    (k : ℕ) (depth : 2 ≤ k) (a : ℂ) (ha : 5 < ‖a‖)
    {ζ : ℂ} (root : IsPrimitiveRoot ζ (2 ^ k)) :
    (1 : P k ℂ) ∈
      Ideal.span (Set.range (inputs k depth a (-1 - a) 1)) := by
  classical
  have pos : 0 < 2 ^ k := by positivity
  refine UnitInstance.one_mem_of_no_grid_zero
    (fun _ => Polynomial.nthRootsFinset (2 ^ k) (1 : ℂ)) _ ?_ ?_
  · intro i
    have identity := congrArg (Polynomial.aeval (X i : P k ℂ))
      (Polynomial.X_pow_sub_one_eq_prod pos root)
    simp only [map_sub, map_pow, Polynomial.aeval_X, map_one, map_prod,
      Polynomial.aeval_C, MvPolynomial.algebraMap_eq] at identity
    rw [← identity]
    exact UnitInstance.var_pow_mem k depth a (-1 - a) 1 i
  · intro x hx
    have roots : ∀ i, x i ^ 2 ^ k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth a (-1 - a) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth a (-1 - a) 1
      (active k depth 0) (by simp [active]; omega)
    have second := UnitInstance.chain_step_mem k depth a (-1 - a) 1
      (active k depth 2) (by simp [active]; omega)
    have first_eq :
        (⟨(active k depth 0).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
          active k depth 1 := by apply Fin.ext; simp [active]
    have second_eq :
        (⟨(active k depth 2).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
          active k depth 3 := by apply Fin.ext; simp [active]
    rw [first_eq] at first
    rw [second_eq] at second
    by_cases c₁ : eval x
        (X (active k depth 0)^2 - X (active k depth 1) : P k ℂ) = 0
    swap
    · exact ⟨_, first, c₁⟩
    by_cases c₂ : eval x
        (X (active k depth 2)^2 - X (active k depth 3) : P k ℂ) = 0
    swap
    · exact ⟨_, second, c₂⟩
    exfalso
    have identity := congrArg (eval x)
      (raw_normalization_identity k depth a (-1 - a) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, eval_rawS3_rational] at identity
    exact s3_ne_zero_on_complex_root_grid_of_norm_gt k a
      (x (active k depth 0)) (x (active k depth 2))
      ha (roots (active k depth 0)) (roots (active k depth 2)) identity


/-- Uniform, explicit primitive-root removal over C. -/
theorem unit_complex_of_norm_gt (k : ℕ) (depth : 2 ≤ k)
    (a : ℂ) (ha : 5 < ‖a‖) :
    (1 : P k ℂ) ∈
      Ideal.span (Set.range (inputs k depth a (-1 - a) 1)) := by
  have hn : 2 ^ k ≠ 0 := by positivity
  exact unit_of_primitiveRoot_complex_of_norm_gt k depth a ha
    (Complex.isPrimitiveRoot_exp (2 ^ k) hn)

/-- This whole open family is nonsingular: the discriminant factors exactly. -/
theorem nonsingular_of_norm_gt (a : ℂ) (ha : 5 < ‖a‖) :
    4 * a ^ 3 + 27 * (-1 - a) ^ 2 ≠ 0 := by
  have h3 : a + 3 ≠ 0 := by
    intro h
    have heq : a = -(3 : ℂ) := by linear_combination h
    rw [heq] at ha
    norm_num at ha
  have h4 : 4 * a + 3 ≠ 0 := by
    intro h
    have heq : a = -(3/4 : ℂ) := by linear_combination (1/4 : ℂ) * h
    rw [heq] at ha
    norm_num at ha
  intro h
  have hf : 4 * a ^ 3 + 27 * (-1 - a) ^ 2 =
      (a + 3) ^ 2 * (4 * a + 3) := by ring
  rw [hf] at h
  rcases mul_eq_zero.mp h with hz | hz
  · exact h3 (pow_eq_zero hz)
  · exact h4 hz

end ICFallDegrees.GrowthCharZeroSeparatedFamily
