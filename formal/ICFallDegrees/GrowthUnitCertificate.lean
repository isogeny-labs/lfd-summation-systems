import ICFallDegrees.GrowthRationalUnit
import Mathlib.Tactic.Ring

/-!
# N8: explicit Bézout-style unit certificate from a Frobenius mismatch

This module isolates the algebraic unit-certificate step in a generic
commutative K-algebra R. Once q^p=q in the two-chain quotient and the
characteristic-p Frobenius identity have been proved, set u=q-a and
c=a-a^p. For c != 0, an explicit inverse of u (and u^2=S3) follows.

No `GrowthMain.CurveModel` or Wronskian assumptions are used.
The premises `hq` and `hFrob` of the final theorem must be discharged
from the field equations and characteristic p; they are *not* concealed
axioms. Lean build and axiom audit required before claiming kernel.
-/

namespace ICFallDegrees.GrowthUnitCertificate

/-- A difference of two powers equal to a nonzero scalar gives an
explicit two-sided inverse in any commutative algebra over a field. -/
theorem power_difference_unit_certificate
    {K R : Type*} [Field K] [CommRing R] [Algebra K R]
    (n : ℕ) (u : R) (c : K) (hc : c ≠ 0)
    (h : u ^ (n + 1) - u = algebraMap K R c) :
    u * ((u ^ n - 1) * algebraMap K R c⁻¹) = 1 := by
  calc
    u * ((u ^ n - 1) * algebraMap K R c⁻¹) =
        (u ^ (n + 1) - u) * algebraMap K R c⁻¹ := by
      rw [pow_succ]
      ring
    _ = algebraMap K R c * algebraMap K R c⁻¹ := by rw [h]
    _ = 1 := by
      rw [← map_mul]
      simp [hc]

/-- In particular, the S3 square has an explicit inverse. -/
theorem square_power_difference_unit_certificate
    {K R : Type*} [Field K] [CommRing R] [Algebra K R]
    (n : ℕ) (u : R) (c : K) (hc : c ≠ 0)
    (h : u ^ (n + 1) - u = algebraMap K R c) :
    u ^ 2 * ((u ^ n - 1) * algebraMap K R c⁻¹) ^ 2 = 1 := by
  calc
    u ^ 2 * ((u ^ n - 1) * algebraMap K R c⁻¹) ^ 2 =
        (u * ((u ^ n - 1) * algebraMap K R c⁻¹)) ^ 2 := by ring
    _ = 1 := by
      rw [power_difference_unit_certificate n u c hc h]
      ring

/-- Explicit unit certificate for u=q-a in a Frobenius-fixed quotient.
The assumptions q^(n+1)=q and Freshman's dream must be proved
separately from the N-th-root equations and CharP with n+1=p. -/
theorem square_unit_of_frobenius_mismatch
    {K R : Type*} [Field K] [CommRing R] [Algebra K R]
    (n : ℕ) (q : R) (a : K)
    (hc : a - a ^ (n + 1) ≠ 0)
    (hq : q ^ (n + 1) = q)
    (hFrob : (q - algebraMap K R a) ^ (n + 1) =
      q ^ (n + 1) - (algebraMap K R a) ^ (n + 1)) :
    (q - algebraMap K R a) ^ 2 *
      (((q - algebraMap K R a) ^ n - 1) *
        algebraMap K R (a - a ^ (n + 1))⁻¹) ^ 2 = 1 := by
  have hdiff :
      (q - algebraMap K R a) ^ (n + 1) -
          (q - algebraMap K R a) =
        algebraMap K R (a - a ^ (n + 1)) := by
    calc
      (q - algebraMap K R a) ^ (n + 1) -
          (q - algebraMap K R a) =
          (q ^ (n + 1) - (algebraMap K R a) ^ (n + 1)) -
            (q - algebraMap K R a) := by rw [hFrob]
      _ = algebraMap K R (a - a ^ (n + 1)) := by
        rw [hq, map_sub, map_pow]
        ring
  exact square_power_difference_unit_certificate n
    (q - algebraMap K R a) (a - a ^ (n + 1)) hc hdiff

end ICFallDegrees.GrowthUnitCertificate
