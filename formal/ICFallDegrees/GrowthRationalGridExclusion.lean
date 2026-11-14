import ICFallDegrees.GrowthRationalUnit
import ICFallDegrees.GrowthRationalFrobenius

/-!
# N8: no two-torsion Semaev zero on the Frobenius-fixed grid

This is a direct algebraic no-grid-zero criterion for the rational
unit-ideal family, valid over any field of characteristic p with
a^p != a and with N | p-1.

It is an intermediate step toward an unconditional unit-ideal proof:
it does NOT itself show the chain ideal is the unit ideal, nor does it
provide the rational-function spread / Wronskian estimate.
-/

namespace ICFallDegrees.GrowthRationalGridExclusion

open ICFallDegrees.GrowthRationalUnit
open ICFallDegrees.GrowthRationalFrobenius
open ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- If a is not Frobenius fixed, the rational two-torsion S3 cannot
vanish at any two N-th roots of unity, provided N divides p-1. -/
theorem s3_ne_zero_on_root_grid
    (p N : ℕ) [Fact p.Prime] [CharP K p]
    (hdvd : N ∣ p - 1) (a x z : K)
    (ha : a ^ p ≠ a)
    (hx : x ^ N = 1) (hz : z ^ N = 1) :
    s3 a (-1 - a) x z 1 ≠ 0 := by
  intro hzero
  have hq :
      ((x - 1) * z - x - 2) ^ p = (x - 1) * z - x - 2 :=
    q_frobenius_fixed p N hdvd x z hx hz
  have hp : a = (x - 1) * z - x - 2 :=
    parameter_of_s3_zero a x z hzero
  apply ha
  calc
    a ^ p = ((x - 1) * z - x - 2) ^ p := by rw [hp]
    _ = (x - 1) * z - x - 2 := hq
    _ = a := hp.symm

end ICFallDegrees.GrowthRationalGridExclusion
