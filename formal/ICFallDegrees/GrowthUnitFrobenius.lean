import ICFallDegrees.GrowthUnitCertificate
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Tactic.Omega

/-!
# N8: explicit characteristic-p S3 unit certificate in the root quotient

The equations x^N=1, z^N=1 with N | p-1 imply x^p=x, z^p=z.
For q=(x-1)z-x-2, Frobenius yields q^p=q.
If a^p != a, then u=q-a satisfies u^p-u=a-a^p, a nonzero scalar.
The generic algebraic lemma in `GrowthUnitCertificate` exhibits an
explicit inverse of u^2=S3. No Nullstellensatz and no Wronskian is used.

All premises refer to original root equations, the characteristic and a
concrete parameter. They do not assume a unit ideal or CurveModel.
Candidate Lean 4.19.0 code; the pinned CI and #print axioms must pass.
-/

namespace ICFallDegrees.GrowthUnitFrobenius

open ICFallDegrees.GrowthUnitCertificate
open ICFallDegrees.GrowthRationalUnit

/-- In any multiplicative commutative ring, N-th roots are p-Frobenius
fixed if N divides p-1. This part does not require characteristic p. -/
theorem pow_fixed_of_root_of_unity
    {R : Type*} [CommRing R]
    (N p : ℕ) (x : R) (hp : 0 < p)
    (hdiv : N ∣ p - 1) (hx : x ^ N = 1) : x ^ p = x := by
  obtain ⟨m, hm⟩ := hdiv
  have heq : p = N * m + 1 := by omega
  calc
    x ^ p = x ^ (N * m + 1) := by rw [heq]
    _ = (x ^ N) ^ m * x := by rw [pow_add, pow_mul, pow_one]
    _ = x := by rw [hx]; simp

/-- The bivariate q from the reduced S3 curve is fixed by Frobenius
as soon as both coordinates are Frobenius-fixed. -/
theorem q_frobenius_fixed
    {R : Type*} [CommRing R] (p : ℕ)
    [Fact p.Prime] [CharP R p]
    (x z : R) (hx : x ^ p = x) (hz : z ^ p = z) :
    ((x - 1) * z - x - 2) ^ p = (x - 1) * z - x - 2 := by
  have htwo : (2 : R) ^ p = 2 := by
    simpa using (add_pow_char (1 : R) (1 : R) p)
  calc
    ((x - 1) * z - x - 2) ^ p =
        (((x - 1) * z - x) ^ p - (2 : R) ^ p) := by rw [sub_pow_char]
    _ = (((x - 1) * z) ^ p - x ^ p) - (2 : R) ^ p := by rw [sub_pow_char]
    _ = ((x ^ p - 1 ^ p) * z ^ p - x ^ p) - (2 : R) ^ p := by
      rw [mul_pow, sub_pow_char]
    _ = (x - 1) * z - x - 2 := by
      rw [hx, hz, htwo]
      ring

/-- The complete Bézout certificate for the actual rational S3
generator in any commutative K-algebra satisfying the N-th-root
equations. This is *ideal-unit* information only, not an LFD bound. -/
theorem s3_unit_certificate_of_root_equations
    {K R : Type*} [Field K] [CommRing R] [Algebra K R]
    (p N : ℕ) [Fact p.Prime] [CharP R p]
    (x z : R) (a : K) (hp : 0 < p)
    (hdiv : N ∣ p - 1) (hx : x ^ N = 1) (hz : z ^ N = 1)
    (ha : a ^ p ≠ a) :
    GrowthSecant.s3 (algebraMap K R a) (-1 - algebraMap K R a) x z 1 *
      (((((x - 1) * z - x - 2) - algebraMap K R a) ^ (p - 1) - 1) *
        algebraMap K R (a - a ^ p)⁻¹) ^ 2 = 1 := by
  let q : R := (x - 1) * z - x - 2
  have hxp := pow_fixed_of_root_of_unity N p x hp hdiv hx
  have hzp := pow_fixed_of_root_of_unity N p z hp hdiv hz
  have hq : q ^ p = q := q_frobenius_fixed p x z hxp hzp
  have hc : a - a ^ p ≠ 0 := sub_ne_zero.mpr ha.symm
  have hFrob : (q - algebraMap K R a) ^ p =
      q ^ p - (algebraMap K R a) ^ p := sub_pow_char _ _
  have hp' : p - 1 + 1 = p := by omega
  have hc' : a - a ^ (p - 1 + 1) ≠ 0 := by simpa only [hp'] using hc
  have hq' : q ^ (p - 1 + 1) = q := by simpa only [hp'] using hq
  have hFrob' : (q - algebraMap K R a) ^ (p - 1 + 1) =
      q ^ (p - 1 + 1) - (algebraMap K R a) ^ (p - 1 + 1) := by
    simpa only [hp'] using hFrob
  have cert := square_unit_of_frobenius_mismatch (p - 1) q a hc' hq' hFrob'
  have hs3 : GrowthSecant.s3 (algebraMap K R a)
      (-1 - algebraMap K R a) x z 1 =
      (q - algebraMap K R a) ^ 2 := by
    rw [s3_square_family]
    dsimp [q]
    congr 1
    ring
  rw [hs3]
  simpa only [hp', q] using cert

end ICFallDegrees.GrowthUnitFrobenius
