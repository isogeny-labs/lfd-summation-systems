import ICFallDegrees.GrowthConcreteClosureStages
import ICFallDegrees.GrowthIntegralBezoutTransfer
import ICFallDegrees.GrowthFixedE6Existence
import Mathlib.Tactic

/-!
# N8 native-field LFD assembly from actual closure and unit certificates

This proves the last logical step of the prime-field LFD program *without*
assuming any characteristic-p CurveModel, Wronskian or SpreadBound.

The input certificates are **native-field statements**: a unit ideal and
nonmembership of 1 in the actual inductive degree-D InClosure, where
D = floor((k-10)/14).  Neither certificate is assumed to follow merely
from the characteristic-zero theorem: their prime specialization remains
the substantive open problem.

No theorem in this file asserts that every p, or all sufficiently large
p uniformly in k, satisfies either certificate.
-/

noncomputable section

namespace ICFallDegrees.GrowthPrimeNativeWindow

open ICFallDegrees.GrowthMain
open ICFallDegrees.PrimeM2General

/-- Native-field transfer target: unit ideal plus actual degree-D
nonmembership gives the *least* last-fall degree and a sharp linear window,
with no auxiliary model or spread assumption.  This is not a transfer
of the two supplied certificates from characteristic zero. -/
theorem exists_lfd_window_of_native_certificates
    {K : Type*} [Field K]
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (a b t : K) (two_ne : (2 : K) ≠ 0)
    (nonsingular : 4 * a ^ 3 + 27 * b ^ 2 ≠ 0)
    (unit : (1 : MvPolynomial (Fin (2 * k)) K) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k depth a b t)))
    (outside : ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t))
      ((k - 10) / 14) (1 : MvPolynomial (Fin (2 * k)) K)) :
    ∃ d : ℕ,
      IsLastFallDegree (Set.range (PrimeM2General.inputs k depth a b t)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 := by
  classical
  let F : Set (MvPolynomial (Fin (2 * k)) K) :=
    Set.range (PrimeM2General.inputs k depth a b t)
  have hex : ∃ n : ℕ, HasLFDBound F n := by
    exact ⟨k + 2,
      SelectorExtension.lfd_le_depth_add_two_of_two_ne
        two_ne k depth a b t nonsingular⟩
  let d : ℕ := Nat.find hex
  have hd : IsLastFallDegree F d := by
    constructor
    · exact Nat.find_spec hex
    · intro e he hbound
      exact (Nat.find_min hex he) hbound
  have hlo : (k - 10) / 14 < d := by
    by_contra h
    have hle : d ≤ (k - 10) / 14 := by omega
    have hmem : InClosure F d (1 : MvPolynomial (Fin (2 * k)) K) :=
      one_of_hasLFDBound unit hd.1
    exact outside (hmem.mono hle)
  have hhi : d ≤ k + 2 := by
    by_contra h
    have hlt : k + 2 < d := by omega
    exact hd.2 (k + 2) hlt
      (SelectorExtension.lfd_le_depth_add_two_of_two_ne
        two_ne k depth a b t nonsingular)
  exact ⟨d, hd, hlo, hhi⟩

/-- The native certificates also imply explicit linear Theta constants.
The quantifiers are genuinely over a native field, but the certificates
still must be proved for the primes one wants to use. -/
theorem exists_lfd_theta_of_native_certificates
    {K : Type*} [Field K]
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (a b t : K) (two_ne : (2 : K) ≠ 0)
    (nonsingular : 4 * a ^ 3 + 27 * b ^ 2 ≠ 0)
    (unit : (1 : MvPolynomial (Fin (2 * k)) K) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k depth a b t)))
    (outside : ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t))
      ((k - 10) / 14) (1 : MvPolynomial (Fin (2 * k)) K)) :
    ∃ d : ℕ,
      IsLastFallDegree (Set.range (PrimeM2General.inputs k depth a b t)) d ∧
      k ≤ 28 * d ∧ d ≤ 2 * k := by
  obtain ⟨d, hd, hlo, hhi⟩ :=
    exists_lfd_window_of_native_certificates k depth hk a b t
      two_ne nonsingular unit outside
  exact ⟨d, hd, by omega, by omega⟩

/-- For p > 2187 the fixed E6 discriminant and 2 both remain nonzero.
This excludes the singular p=3 without asserting any unit ideal.
The exact E6 p=17,k=2 counterexample remains valid. -/
theorem fixed_e6_good_characteristic (p : ℕ) [Fact p.Prime]
    (hp : 2187 < p) :
    (2 : ZMod p) ≠ 0 ∧
      4 * (6 : ZMod p) ^ 3 + 27 * (-7 : ZMod p) ^ 2 ≠ 0 := by
  have htwo : (2 : ZMod p) ≠ 0 := by
    intro hz
    have hdiv : p ∣ 2 :=
      (ZMod.natCast_zmod_eq_zero_iff_dvd 2 p).mp hz
    have hle : p ≤ 2 := Nat.le_of_dvd (by norm_num) hdiv
    omega
  have hdisc : (2187 : ZMod p) ≠ 0 := by
    intro hz
    have hdiv : p ∣ 2187 :=
      (ZMod.natCast_zmod_eq_zero_iff_dvd 2187 p).mp hz
    have hle : p ≤ 2187 := Nat.le_of_dvd (by norm_num) hdiv
    omega
  refine ⟨htwo, ?_⟩
  have hid : 4 * (6 : ZMod p) ^ 3 + 27 * (-7 : ZMod p) ^ 2 =
      (2187 : ZMod p) := by ring
  exact hid ▸ hdisc

/-- Actual prime-field E6 target interface. Assumptions are precisely
the *two native certificates* still missing from the proof.
The bound p>2187 is a convenient sufficient condition for a nonsingular
E6 model, NOT a claimed exceptional-prime bound for N8 closure/unit. -/
theorem fixed_e6_prime_theta_of_native_certificates
    (p : ℕ) [Fact p.Prime] (hp : 2187 < p)
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k)
    (unit : (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∈
      Ideal.span (Set.range
        (PrimeM2General.inputs k depth (6 : ZMod p) (-7 : ZMod p) 1)))
    (outside : ¬ InClosure
      (Set.range
        (PrimeM2General.inputs k depth (6 : ZMod p) (-7 : ZMod p) 1))
      ((k - 10) / 14) (1 : MvPolynomial (Fin (2 * k)) (ZMod p))) :
    ∃ d : ℕ,
      IsLastFallDegree
        (Set.range (PrimeM2General.inputs k depth (6 : ZMod p) (-7 : ZMod p) 1)) d ∧
      (k - 10) / 14 < d ∧ d ≤ k + 2 ∧ k ≤ 28 * d ∧ d ≤ 2 * k := by
  obtain ⟨htwo, hdisc⟩ := fixed_e6_good_characteristic p hp
  obtain ⟨d, hd, hlo, hhi⟩ :=
    exists_lfd_window_of_native_certificates k depth hk
      (6 : ZMod p) (-7 : ZMod p) 1 htwo hdisc unit outside
  exact ⟨d, hd, hlo, hhi, by omega, by omega⟩

#print axioms exists_lfd_window_of_native_certificates
#print axioms exists_lfd_theta_of_native_certificates
#print axioms fixed_e6_good_characteristic
#print axioms fixed_e6_prime_theta_of_native_certificates

end ICFallDegrees.GrowthPrimeNativeWindow
