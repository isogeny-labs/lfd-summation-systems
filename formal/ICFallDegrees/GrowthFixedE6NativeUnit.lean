import ICFallDegrees.GrowthFixedE6IntegralBezout
import Mathlib.Tactic

/-!
# Actual E6 prime-field unit ideal outside finitely many bad primes

The integral normalized input family is mapped into arbitrary fields and
identified with the *actual* PrimeM2General inputs. Combined with the
integer Bézout scalar obtained from the complex E6 theorem, this yields a
native unit ideal over any field avoiding one nonzero integer, and hence
for all sufficiently large prime fields at each fixed (k, depth).

This theorem does NOT address low-degree InClosure nonmembership. It is
not a prime-field LFD lower bound or a fixed-prime sequence assertion.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6NativeUnit

open MvPolynomial
open ICFallDegrees.GrowthFixedE6IntegralInputs
open ICFallDegrees.GrowthFixedE6IntegralBezout
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GeneralQuadraticTop

variable {K : Type*} [Field K]

/-- All chain endpoint coefficients commute with integer reduction to
any field, not only the rationals. -/
theorem map_chain_tail_field (k : ℕ) (i : Fin (2 * k)) :
    (chainTailInt k i).map (Int.castRingHom K) =
      PrimeM2General.chainTail (K := K) k i := by
  unfold chainTailInt PrimeM2General.chainTail
  split_ifs <;> simp

/-- The integer E6 normalized quadratic maps to the exact normalized
quadratic of E6 in every characteristic. -/
theorem map_normalized_field (k : ℕ) (depth : 2 ≤ k) :
    (normalizedInt k depth).map (Int.castRingHom K) =
      PrimeM2General.normalized k depth (6 : K) (-7 : K) 1 := by
  unfold normalizedInt PrimeM2General.normalized PrimeM2General.tail
    GeneralQuadraticTop.top
  simp only [map_add, map_mul, map_C, map_X, map_sub, map_neg]
  simp only [map_ofNat, map_natCast, map_intCast]
  norm_num <;> simp only [map_ofNat] <;> ring

/-- Exact actual E6 native input identity over any field. -/
theorem map_e6_inputs_field (k : ℕ) (depth : 2 ≤ k)
    (i : Fin (2 * k + 1)) :
    (e6InputsInt k depth i).map (Int.castRingHom K) =
      PrimeM2General.inputs k depth (6 : K) (-7 : K) 1 i := by
  refine Fin.lastCases ?_ (fun i => ?_) i
  · simpa only [e6InputsInt, PrimeM2General.inputs, Fin.lastCases_last]
      using (map_normalized_field (K := K) k depth)
  · simpa only [e6InputsInt, PrimeM2General.inputs,
      Fin.lastCases_castSucc, map_add, map_pow, map_X] using
      congrArg (fun f : MvPolynomial (Fin (2 * k)) K =>
        (X i : MvPolynomial (Fin (2 * k)) K) ^ 2 + f)
        (map_chain_tail_field (K := K) k i)

/-- A genuine native E6 unit ideal for every field in which one
positive integral scalar d survives. No parameter- or generator-level
native unit assumption is made. -/
theorem exists_native_e6_unit_scalar (k : ℕ) (depth : 2 ≤ k) :
    ∃ d : ℕ, 0 < d ∧
      ∀ (L : Type*) [Field L], (d : L) ≠ 0 →
        (1 : MvPolynomial (Fin (2 * k)) L) ∈
          Ideal.span (Set.range
            (PrimeM2General.inputs k depth (6 : L) (-7 : L) 1)) := by
  obtain ⟨d, hd, hunit⟩ := exists_e6_unit_specialization_scalar k depth
  refine ⟨d, hd, ?_⟩
  intro L _ hnonzero
  have h := hunit L hnonzero
  simpa only [map_e6_inputs_field] using h

/-- For fixed (k,depth), all sufficiently large prime fields have the
actual normalized E6 unit ideal. The bound is existential, not evaluated. -/
theorem exists_actual_e6_prime_unit_bound (k : ℕ) (depth : 2 ≤ k) :
    ∃ B : ℕ, ∀ (p : ℕ) [Fact p.Prime], B < p →
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∈
        Ideal.span (Set.range
          (PrimeM2General.inputs k depth
            (6 : ZMod p) (-7 : ZMod p) 1)) := by
  obtain ⟨d, hd, hunit⟩ := exists_native_e6_unit_scalar k depth
  refine ⟨d, ?_⟩
  intro p _ hgt
  have hnonzero : (d : ZMod p) ≠ 0 := by
    intro hz
    have hdiv : p ∣ d :=
      (ZMod.natCast_zmod_eq_zero_iff_dvd d p).mp hz
    have hle : p ≤ d := Nat.le_of_dvd hd hdiv
    omega
  exact hunit (ZMod p) hnonzero

#print axioms map_chain_tail_field
#print axioms map_normalized_field
#print axioms map_e6_inputs_field
#print axioms exists_native_e6_unit_scalar
#print axioms exists_actual_e6_prime_unit_bound

end ICFallDegrees.GrowthFixedE6NativeUnit
