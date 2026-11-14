import ICFallDegrees.GrowthRationalChainUnit
import ICFallDegrees.GrowthCharZeroDeterminant
import Mathlib.Analysis.Complex.Basic
import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Tactic

/-!
# N8: characteristic-zero unit ideal from a separated complex parameter

For a=6, b=-7, t=1, the Semaev factor is
  S3(x,z,1) = ((x-1)*z - x - 8)^2.
If x^N=z^N=1 over C, both x,z have norm one, whereas
|(x-1)*z-x-2| <= 5.  Thus S3 never vanishes on the root grid.
The inherited finite-grid Nullstellensatz then gives an actual unit
ideal for the full two-summand chain at every depth, without a
group-algebra determinant and without assuming unit ideal as a premise.

This does NOT prove GrowthSpread.SpreadBound; growth is still conditional.
-/

noncomputable section
namespace ICFallDegrees.GrowthCharZeroComplexUnit

open MvPolynomial PrimeM2General
open ICFallDegrees.GrowthSecant
open ICFallDegrees.GrowthRationalUnit
open ICFallDegrees.GrowthRationalChainUnit

/-- A nonnegative real whose 2^k-th power is one is itself one. -/
theorem nonneg_eq_one_of_pow_two (k : ℕ) (r : ℝ)
    (hr : 0 ≤ r) (h : r ^ (2 ^ k) = 1) : r = 1 := by
  induction k generalizing r with
  | zero =>
      simpa using h
  | succ k ih =>
      have hs : (r ^ (2 ^ k)) ^ 2 = 1 := by
        rw [pow_succ, pow_mul] at h
        exact h
      have hn : 0 ≤ r ^ (2 ^ k) := pow_nonneg hr _
      have h' : r ^ (2 ^ k) = 1 := by nlinarith
      exact ih r hr h'

/-- The norm of a 2-power root of unity in C is one. -/
theorem norm_eq_one_of_root (k : ℕ) (x : ℂ)
    (hx : x ^ (2 ^ k) = 1) : ‖x‖ = 1 := by
  have hnorm : ‖x‖ ^ (2 ^ k) = (1 : ℝ) := by
    rw [← norm_pow, hx, norm_one]
  exact nonneg_eq_one_of_pow_two k ‖x‖ (norm_nonneg _) hnorm

/-- The reduced rational Semaev polynomial cannot vanish on a
2^k-root grid for the concrete rational parameter a=6. -/
theorem s3_ne_zero_on_complex_root_grid
    (k : ℕ) (x z : ℂ)
    (hx : x ^ (2 ^ k) = 1) (hz : z ^ (2 ^ k) = 1) :
    s3 (6 : ℂ) (-1 - 6) x z 1 ≠ 0 := by
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
  have hp := parameter_of_s3_zero (6 : ℂ) x z h
  have hbound : ‖(6 : ℂ)‖ ≤ 5 := by
    rw [hp]
    exact hq
  norm_num at hbound

/-- The actual chain ideal is the unit ideal over C, for every k >= 2,
for a nonsingular two-torsion family (a,b,t)=(6,-7,1). -/
theorem unit_of_primitiveRoot_complex
    (k : ℕ) (depth : 2 ≤ k)
    {ζ : ℂ} (root : IsPrimitiveRoot ζ (2 ^ k)) :
    (1 : P k ℂ) ∈
      Ideal.span (Set.range (inputs k depth (6 : ℂ) (-1 - 6) 1)) := by
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
    exact UnitInstance.var_pow_mem k depth (6 : ℂ) (-1 - 6) 1 i
  · intro x hx
    have roots : ∀ i, x i ^ 2 ^ k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth (6 : ℂ) (-1 - 6) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth (6 : ℂ) (-1 - 6) 1
      (active k depth 0) (by simp [active]; omega)
    have second := UnitInstance.chain_step_mem k depth (6 : ℂ) (-1 - 6) 1
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
      (raw_normalization_identity k depth (6 : ℂ) (-1 - 6) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, eval_rawS3_rational] at identity
    exact s3_ne_zero_on_complex_root_grid k
      (x (active k depth 0)) (x (active k depth 2))
      (roots (active k depth 0)) (roots (active k depth 2)) identity

/-- The explicit root of unity in C removes the remaining root-existence
hypothesis from the chain unit-ideal certificate. -/
theorem unit_complex (k : ℕ) (depth : 2 ≤ k) :
    (1 : P k ℂ) ∈
      Ideal.span (Set.range (inputs k depth (6 : ℂ) (-1 - 6) 1)) := by
  have hn : 2 ^ k ≠ 0 := by positivity
  exact unit_of_primitiveRoot_complex k depth
    (Complex.isPrimitiveRoot_exp (2 ^ k) hn)

/-- The exhibited rational curve is nonsingular. -/
theorem concrete_nonsingular :
    4 * (6 : ℂ) ^ 3 + 27 * (-1 - (6 : ℂ)) ^ 2 ≠ 0 := by
  norm_num

end ICFallDegrees.GrowthCharZeroComplexUnit
