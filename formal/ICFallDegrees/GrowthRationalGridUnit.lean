import ICFallDegrees.GrowthRationalFrobenius
import ICFallDegrees.GrowthUnitCertificate

/-!
# N8 rational unit family: explicit inverse on the N-th-root quotient

If both coordinates satisfy x^N=z^N=1, p is the characteristic,
N divides p-1, and a^p != a in the coefficient field, the specialized
Semaev polynomial is a unit. The proof produces the inverse explicitly.
This closes the *algebraic* unit-ideal part after transporting the
root equations from the actual chain quotient. It is NOT a low-degree
closure certificate and does not prove the Wronskian spread bound.
-/

namespace ICFallDegrees.GrowthRationalGridUnit

variable {K R : Type*} [Field K] [CommRing R] [Algebra K R]

/-- In the quotient by x^N-1,z^N-1, S3=(q-a)^2 has the displayed
explicit inverse; no Nullstellensatz required. -/
theorem square_unit_of_grid_roots (p N : ℕ) [Fact p.Prime]
    [CharP R p] (hdvd : N ∣ p - 1)
    (a : K) (ha : a - a ^ p ≠ 0)
    (x z : R) (hx : x ^ N = 1) (hz : z ^ N = 1) :
    let q : R := (x - 1) * z - x - 2
    (q - algebraMap K R a) ^ 2 *
      (((q - algebraMap K R a) ^ (p - 1) - 1) *
        algebraMap K R (a - a ^ p)⁻¹) ^ 2 = 1 := by
  dsimp
  let q : R := (x - 1) * z - x - 2
  have hp : p - 1 + 1 = p :=
    Nat.sub_add_cancel (Fact.out : p.Prime).one_le
  have hq : q ^ (p - 1 + 1) = q := by
    rw [hp]
    exact GrowthRationalFrobenius.q_frobenius_fixed p N hdvd x z hx hz
  have hFrob :
      (q - algebraMap K R a) ^ (p - 1 + 1) =
        q ^ (p - 1 + 1) - (algebraMap K R a) ^ (p - 1 + 1) := by
    rw [hp]
    exact sub_pow_char q (algebraMap K R a)
  have hcert := GrowthUnitCertificate.square_unit_of_frobenius_mismatch
    (p - 1) q a (by simpa [hp] using ha) hq hFrob
  simpa only [q, hp] using hcert

end ICFallDegrees.GrowthRationalGridUnit
