import ICFallDegrees.GrowthRationalGridExclusion
import ICFallDegrees.UnitInstance

/-!
# N8: a genuine unit ideal for the two-torsion rational prime chain

Unlike the earlier Frobenius inverse in the abstract two-variable
quotient, this lemma uses the existing chain-ideal grid Nullstellensatz
and explicitly connects the normalized Semaev generator to the raw
S3 identity. It proves that the **chain system itself** generates
the unit ideal under the stated primitive-root and parameter hypotheses.

This does NOT establish the low-degree non-closure or the Wronskian
spread bound. It remains separate from the analytic CurveModel.
-/

noncomputable section

namespace ICFallDegrees.GrowthRationalChainUnit

open MvPolynomial PrimeM2General
open ICFallDegrees.GrowthRationalGridExclusion
open ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- The raw chain S3 is the scalar Semaev polynomial, including b≠0. -/
theorem eval_rawS3_rational (k : ℕ) (depth : 2 ≤ k) (a : K)
    (x : Fin (2*k) → K) :
    eval x (rawS3 k depth a (-1-a) 1) =
      s3 a (-1-a) (x (active k depth 0)) (x (active k depth 2)) 1 := by
  simp only [rawS3, s3, map_sub, map_add, map_mul, map_pow,
    eval_X, eval_C, map_ofNat]

/-- **Unconditional unit-ideal lemma** for the rational two-torsion
family at any depth k≥2. The target t=1 has y-coordinate 0, so this
does not follow from the old UnitInstance.unit_zmod lemma. -/
theorem unit_of_primitiveRoot_frobenius
    (p k : ℕ) [Fact p.Prime] [CharP K p]
    (depth : 2 ≤ k) (hdvd : 2^k ∣ p - 1)
    (a : K) (ha : a ^ p ≠ a)
    {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) :
    (1 : P k K) ∈
      Ideal.span (Set.range (inputs k depth a (-1-a) 1)) := by
  classical
  have pos : 0 < 2^k := by positivity
  refine UnitInstance.one_mem_of_no_grid_zero
    (fun _ => Polynomial.nthRootsFinset (2^k) (1 : K)) _ ?_ ?_
  · intro i
    have identity := congrArg (Polynomial.aeval (X i : P k K))
      (Polynomial.X_pow_sub_one_eq_prod pos root)
    simp only [map_sub, map_pow, Polynomial.aeval_X, map_one, map_prod,
      Polynomial.aeval_C, MvPolynomial.algebraMap_eq] at identity
    rw [← identity]
    exact UnitInstance.var_pow_mem k depth a (-1-a) 1 i
  · intro x hx
    have roots : ∀ i, x i^2^k = 1 :=
      fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth a (-1-a) 1) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := UnitInstance.chain_step_mem k depth a (-1-a) 1
      (active k depth 0) (by simp [active]; omega)
    have second := UnitInstance.chain_step_mem k depth a (-1-a) 1
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
        (X (active k depth 0)^2 - X (active k depth 1) : P k K) = 0
    swap
    · exact ⟨_, first, c₁⟩
    by_cases c₂ : eval x
        (X (active k depth 2)^2 - X (active k depth 3) : P k K) = 0
    swap
    · exact ⟨_, second, c₂⟩
    exfalso
    have identity := congrArg (eval x)
      (raw_normalization_identity k depth a (-1-a) 1)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero,
      mul_zero, add_zero, sub_zero, eval_rawS3_rational] at identity
    exact s3_ne_zero_on_root_grid p (2^k) hdvd a
      (x (active k depth 0)) (x (active k depth 2))
      ha (roots (active k depth 0)) (roots (active k depth 2)) identity

end ICFallDegrees.GrowthRationalChainUnit
