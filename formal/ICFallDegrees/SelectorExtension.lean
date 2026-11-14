import ICFallDegrees.SelectorBound
import ICFallDegrees.ScalarExtensionTransport
import ICFallDegrees.PrimeM2Lower
import Mathlib.NumberTheory.Cyclotomic.PrimitiveRoots

noncomputable section

/-!
# The selector bound without roots of unity in the base field

`SelectorBound.lfd_le_depth_add_two` needs a primitive `2^k`-th root of unity in the field.
The last fall degree is invariant under a finite extension of the coefficient field
(`ScalarExtensionTransport.bound_iff`), and the normalized presentation is defined over the
base field, so the bound descends from the `2^k`-th cyclotomic extension: it holds over every
field of characteristic other than two, for every nonsingular cubic and every target.

Over a field without the roots of unity the chain variables no longer parametrize a subgroup
of the base field; the statement is about the polynomial system only.
-/

namespace ICFallDegrees.SelectorExtension

open MvPolynomial
open PrimeM2General (inputs normalized tail chainTail)
open GeneralQuadraticTop (top)

variable {K L : Type*} [Field K] [Field L] [Algebra K L]

/-- The normalized presentation commutes with extension of the coefficient field. -/
theorem map_inputs (k : ℕ) (depth : 2 ≤ k) (a b t : K) (i : Fin (2*k+1)) :
    map (algebraMap K L) (inputs k depth a b t i) =
      inputs k depth (algebraMap K L a) (algebraMap K L b) (algebraMap K L t) i := by
  refine Fin.lastCases ?_ (fun j => ?_) i
  · simp only [inputs, Fin.lastCases_last, normalized, top, tail, map_add, map_mul, map_sub,
      map_neg, map_pow, map_ofNat, MvPolynomial.map_C, MvPolynomial.map_X]
  · simp only [inputs, Fin.lastCases_castSucc, chainTail]
    split_ifs <;> simp

theorem image_inputs (k : ℕ) (depth : 2 ≤ k) (a b t : K) :
    ScalarExtensionTransport.lift (K := L) '' Set.range (inputs k depth a b t) =
      Set.range (inputs k depth (algebraMap K L a) (algebraMap K L b) (algebraMap K L t)) := by
  ext p
  constructor
  · rintro ⟨_, ⟨i, rfl⟩, rfl⟩
    exact ⟨i, (map_inputs k depth a b t i).symm⟩
  · rintro ⟨i, rfl⟩
    exact ⟨_, ⟨i, rfl⟩, map_inputs k depth a b t i⟩

/-- Descent of the selector bound from a finite extension that contains the roots of unity. -/
theorem lfd_le_of_extension [FiniteDimensional K L] (k : ℕ) (depth : 2 ≤ k) {ζ : L}
    (root : IsPrimitiveRoot ζ (2^k)) (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    HasLFDBound (Set.range (inputs k depth a b t)) (k+2) := by
  have nonsingular' : 4*(algebraMap K L a)^3 + 27*(algebraMap K L b)^2 ≠ 0 := by
    intro h
    apply nonsingular
    apply (algebraMap K L).injective
    simp only [map_add, map_mul, map_pow, map_ofNat, map_zero]
    exact h
  have bound := SelectorBound.lfd_le_depth_add_two k depth root (algebraMap K L a)
    (algebraMap K L b) (algebraMap K L t) nonsingular'
  rw [← image_inputs] at bound
  exact (ScalarExtensionTransport.bound_iff (Module.finBasis K L) _ _).mp bound

/-- **Selector bound, every field of characteristic other than two.** For every depth
`k ≥ 2`, every nonsingular cubic `x^3 + a x + b` and every target `t`, the normalized
two-summand chain presentation has last fall degree at most `k + 2`. -/
theorem lfd_le_depth_add_two_of_two_ne (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k)
    (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    HasLFDBound (Set.range (inputs k depth a b t)) (k+2) := by
  haveI : NeZero (((2^k : ℕ+) : ℕ) : K) := ⟨by simpa using pow_ne_zero k two_ne⟩
  haveI : FiniteDimensional K (CyclotomicField (2^k) K) :=
    IsCyclotomicExtension.finite {(2^k : ℕ+)} K (CyclotomicField (2^k) K)
  have root : IsPrimitiveRoot
      (IsCyclotomicExtension.zeta (2^k : ℕ+) K (CyclotomicField (2^k) K)) (2^k) := by
    have h := IsCyclotomicExtension.zeta_spec (2^k : ℕ+) K (CyclotomicField (2^k) K)
    simpa using h
  exact lfd_le_of_extension k depth root a b t nonsingular

/-- Prime-field entry point: every prime `p ≥ 3`, every depth, nonsingular curve, every
target; no condition `2^k ∣ p - 1`. -/
theorem prime_field_lfd_le_depth_add_two (p : ℕ) [Fact p.Prime] (odd : p ≠ 2) (k : ℕ)
    (depth : 2 ≤ k) (a b t : ZMod p) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    HasLFDBound (Set.range (inputs k depth a b t)) (k+2) := by
  refine lfd_le_depth_add_two_of_two_ne ?_ k depth a b t nonsingular
  intro h
  have dvd : p ∣ 2 := (ZMod.natCast_zmod_eq_zero_iff_dvd 2 p).mp (by exact_mod_cast h)
  exact odd ((Nat.prime_dvd_prime_iff_eq Fact.out Nat.prime_two).mp dvd)

/-- Consistency with the lower bound: at depth two the selector bound is attained. For a
unit ideal with `3 t^2 + a ≠ 0` on a nonsingular cubic the last fall degree is exactly four. -/
theorem exact_lfd_four_at_depth_two (two_ne : (2 : K) ≠ 0) (a b t : K)
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (nondegenerate : 3*t^2 + a ≠ 0)
    (unit : (1 : MvPolynomial (Fin (2*2)) K) ∈ Ideal.span (Set.range (inputs 2 le_rfl a b t))) :
    IsLastFallDegree (Set.range (inputs 2 le_rfl a b t)) 4 :=
  ⟨lfd_le_depth_add_two_of_two_ne two_ne 2 le_rfl a b t nonsingular,
    PrimeM2Lower.lfd_ge_four_of_unit two_ne 2 le_rfl a b t nondegenerate unit⟩

end ICFallDegrees.SelectorExtension
