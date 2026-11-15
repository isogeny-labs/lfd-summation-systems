import ICFallDegrees.GrowthUnitCertificate
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Tactic

/-!
# N8: Frobenius fixedness from the root-of-unity grid

These are small algebraic ingredients of the explicit unit-ideal certificate
for the rational two-torsion family. They must be proved in the quotient
from the chain equations, not assumed for the final LFD theorem.
-/

namespace ICFallDegrees.GrowthRationalFrobenius

variable {R : Type*} [CommRing R]

/-- Any N-th root of unity is fixed by exponentiation by p when
N divides p-1. This does not require a field or a characteristic. -/
theorem pow_fixed_of_root (N p : ℕ) (hp : 1 ≤ p)
    (hdvd : N ∣ p - 1) (x : R) (hx : x ^ N = 1) :
    x ^ p = x := by
  obtain ⟨m, hm⟩ := hdvd
  have hexp : p = N * m + 1 := by omega
  rw [hexp, pow_add, pow_mul, hx, one_pow, one_mul, pow_one]

/-- The affine rational-section numerator is Frobenius-fixed
whenever the two coordinates are N-th roots of unity. -/
theorem q_frobenius_fixed (p N : ℕ) [Fact p.Prime] [CharP R p]
    (hdvd : N ∣ p - 1) (x z : R)
    (hx : x ^ N = 1) (hz : z ^ N = 1) :
    ((x - 1) * z - x - 2) ^ p = (x - 1) * z - x - 2 := by
  have hp : 1 ≤ p := (Fact.out : p.Prime).one_le
  have hx' := pow_fixed_of_root N p hp hdvd x hx
  have hz' := pow_fixed_of_root N p hp hdvd z hz
  have htwo : (2 : R) ^ p = 2 := by
    have hh : (2 : R) = 1 + 1 := by ring
    rw [hh, add_pow_char]
    simp
  rw [sub_pow_char, sub_pow_char, mul_pow, sub_pow_char]
  simp only [hx', hz', one_pow, htwo]
 
end ICFallDegrees.GrowthRationalFrobenius
