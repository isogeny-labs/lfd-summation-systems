import ICFallDegrees.GrowthFixedE6PrimeTheta
import Mathlib.Tactic

/-!
# Fixed-E6 N8: explicit increasing sequence of admissible prime fields

The cofinal-depth theorem provides a prime above any prescribed bound.
Use classical recursion to choose the next prime above the previous one,
at depth k=n+24. The resulting sequence is strictly increasing and
each field has a genuine least last fall degree with the Theta window.
There is no single fixed prime and no all-admissible-primes assertion.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6PrimeSequence

open MvPolynomial
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GrowthFixedE6PrimeTheta

/-- The prime at depth n+24, chosen larger than the previous prime. -/
private def growingPrime : ℕ → ℕ
  | 0 =>
      Classical.choose
        (exists_actual_prime_theta_above (0 + 24) 0 (by omega))
  | n + 1 =>
      Classical.choose
        (exists_actual_prime_theta_above ((n + 1) + 24)
          (growingPrime n) (by omega))

/-- Strict growth is built into the lower threshold at each step. -/
private theorem growingPrime_lt_succ (n : ℕ) :
    growingPrime n < growingPrime (n + 1) := by
  have h := Classical.choose_spec
    (exists_actual_prime_theta_above ((n + 1) + 24)
      (growingPrime n) (by omega))
  exact h.2.1

/-- Each selected prime carries the native, unconditional fixed-E6
last-fall-degree certificate and the required admissible congruence. -/
private theorem growingPrime_spec (n : ℕ) :
    Nat.Prime (growingPrime n) ∧
    growingPrime n ≡ 1 [MOD 2 ^ ((n + 24) + 1)] ∧
    ∃ d : ℕ,
      (∀ [Fact (growingPrime n).Prime],
        IsLastFallDegree
          (Set.range (PrimeM2General.inputs (n + 24)
            (by omega : 2 ≤ n + 24)
            (6 : ZMod (growingPrime n))
            (-7 : ZMod (growingPrime n)) 1)) d) ∧
      ((n + 24) - 10) / 14 < d ∧
      d ≤ (n + 24) + 2 ∧
      n + 24 ≤ 28 * d ∧
      28 * d ≤ 56 * (n + 24) := by
  cases n with
  | zero =>
      have h := Classical.choose_spec
        (exists_actual_prime_theta_above (0 + 24) 0 (by omega))
      exact ⟨h.1, h.2.2.1, h.2.2.2⟩
  | succ n =>
      have h := Classical.choose_spec
        (exists_actual_prime_theta_above ((n + 1) + 24)
          (growingPrime n) (by omega))
      exact ⟨h.1, h.2.2.1, h.2.2.2⟩

/-- Fully kernel-checkable growing-prime theorem for the actual
normalized E6 N8 system, indexed by n with k=n+24.
The LFD is the least actual last-fall degree, not a formal surrogate. -/
theorem exists_actual_growing_prime_theta_sequence :
    ∃ p : ℕ → ℕ,
      (∀ n, p n < p (n + 1)) ∧
      ∀ n, Nat.Prime (p n) ∧
        p n ≡ 1 [MOD 2 ^ ((n + 24) + 1)] ∧
        ∃ d : ℕ,
          (∀ [Fact (p n).Prime],
            IsLastFallDegree
              (Set.range (PrimeM2General.inputs (n + 24)
                (by omega : 2 ≤ n + 24)
                (6 : ZMod (p n)) (-7 : ZMod (p n)) 1)) d) ∧
          ((n + 24) - 10) / 14 < d ∧
          d ≤ (n + 24) + 2 ∧
          n + 24 ≤ 28 * d ∧
          28 * d ≤ 56 * (n + 24) := by
  exact ⟨growingPrime, growingPrime_lt_succ, growingPrime_spec⟩

#print axioms exists_actual_growing_prime_theta_sequence

end ICFallDegrees.GrowthFixedE6PrimeSequence
