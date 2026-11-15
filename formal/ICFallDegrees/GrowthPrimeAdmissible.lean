import Mathlib.NumberTheory.PrimesCongruentOne
import Mathlib.Tactic

/-!
# An unconditional supply of arbitrarily large admissible N8 primes

Use Mathlib's elementary cyclotomic theorem for primes congruent to 1
modulo a prescribed positive integer. The N8 modulus is 2^(k+1).

This supplies the growing-prime existence interface only. It DOES NOT
certify N8's finite exceptional sets or a finite-field LFD lower bound.
-/

namespace ICFallDegrees.GrowthPrimeAdmissible

/-- For every depth k and bound B, there is a prime exceeding B
congruent to 1 modulo 2^(k+1), without any use of Dirichlet's theorem. -/
theorem exists_admissible_prime_above (k B : ℕ) :
    ∃ p : ℕ, Nat.Prime p ∧ B < p ∧
      p ≡ 1 [MOD 2 ^ (k + 1)] := by
  have hmod : (2 : ℕ) ^ (k + 1) ≠ 0 := by positivity
  exact Nat.exists_prime_gt_modEq_one B hmod

#print axioms exists_admissible_prime_above

end ICFallDegrees.GrowthPrimeAdmissible
