import ICFallDegrees.GrowthPrimeE6NativeNonmembership
import ICFallDegrees.GrowthFixedE6NativeUnit
import ICFallDegrees.GrowthPrimeNativeWindow
import ICFallDegrees.GrowthPrimeAdmissible
import Mathlib.Tactic

/-!
# Fixed-E6 N8 actual prime-field last-fall-degree growth

Combine genuine native E6 prime-characteristic box nonmembership, the
integral Bezout good-prime unit, the already verified native LFD window,
and the elementary infinite supply of primes congruent to 1 modulo
2^(k+1). The prime is allowed to depend on k and an arbitrary lower
bound; no claim is made for all admissible primes or one fixed p.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6PrimeTheta

open MvPolynomial
open ICFallDegrees.PrimeM2General

/-- Arbitrarily far out in the primes p ≡ 1 (mod 2^(k+1)), the
*actual* fixed E6 N8 system has a least LFD d with linear growth.
The lower prime threshold may depend on k. -/
theorem exists_actual_prime_theta_above
    (k B : ℕ) (hk : 24 ≤ k) :
    ∃ p : ℕ, Nat.Prime p ∧ B < p ∧
      p ≡ 1 [MOD 2 ^ (k + 1)] ∧
      ∃ d : ℕ,
        (∀ [Fact p.Prime], IsLastFallDegree
          (Set.range (PrimeM2General.inputs k (by omega : 2 ≤ k)
            (6 : ZMod p) (-7 : ZMod p) 1)) d) ∧
        (k - 10) / 14 < d ∧ d ≤ k + 2 ∧
        k ≤ 28 * d ∧ 28 * d ≤ 56 * k := by
  let D : ℕ := (k - 10) / 14
  have hD : 1 ≤ D := by dsimp [D]; omega
  have hkD : 14 * D + 10 ≤ k := by dsimp [D]; omega
  have hdepth : 2 ≤ k := by omega
  obtain ⟨Bunit, hunit⟩ :=
    ICFallDegrees.GrowthFixedE6NativeUnit.exists_actual_e6_prime_unit_bound
      k hdepth
  let C : ℕ := max B (max Bunit (max 2187 (2 * 2 ^ k)))
  have hB : B ≤ C := by
    dsimp [C]
    exact le_max_left _ _
  have hU : Bunit ≤ C := by
    dsimp [C]
    exact (le_max_left _ _).trans (le_max_right _ _)
  have hDisc : 2187 ≤ C := by
    dsimp [C]
    exact ((le_max_left _ _).trans (le_max_right _ _)).trans
      (le_max_right _ _)
  have hPow : 2 * 2 ^ k ≤ C := by
    dsimp [C]
    exact ((le_max_right _ _).trans (le_max_right _ _)).trans
      (le_max_right _ _)
  obtain ⟨p, hpprime, hpC, hcong⟩ :=
    ICFallDegrees.GrowthPrimeAdmissible.exists_admissible_prime_above k C
  haveI : Fact p.Prime := ⟨hpprime⟩
  have hpB : B < p := lt_of_le_of_lt hB hpC
  have hpUnit : Bunit < p := lt_of_le_of_lt hU hpC
  have hpDisc : 2187 < p := lt_of_le_of_lt hDisc hpC
  have hpLarge : 2 * 2 ^ k + 1 ≤ p := by omega
  have hnativeUnit :
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k hdepth
        (6 : ZMod p) (-7 : ZMod p) 1)) :=
    hunit p hpUnit
  have hnativeOutside :
      ¬ InClosure (Set.range (PrimeM2General.inputs k hdepth
        (6 : ZMod p) (-7 : ZMod p) 1)) D
        (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) :=
    ICFallDegrees.GrowthPrimeE6NativeNonmembership.e6_not_one_of_large_prime
      p k hdepth D hD hkD hpLarge
  obtain ⟨d, hlfd, hlo, hhi, hklo, hkhi⟩ :=
    ICFallDegrees.GrowthPrimeNativeWindow.fixed_e6_prime_theta_of_native_certificates
      p hpDisc k hdepth hk hnativeUnit
      (by simpa only [D] using hnativeOutside)
  refine ⟨p, hpprime, hpB, hcong, d, ?_, hlo, hhi, hklo, ?_⟩
  · intro _
    simpa only using hlfd
  · omega

#print axioms exists_actual_prime_theta_above

end ICFallDegrees.GrowthFixedE6PrimeTheta
