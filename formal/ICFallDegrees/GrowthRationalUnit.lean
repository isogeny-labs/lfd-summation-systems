import ICFallDegrees.GrowthSecant

/-!
# N8: concrete unit-ideal family with a rational two-torsion curve model

The special family E_a : y² = x³ + a*x - 1 - a has the point (1,0)
of order two. S3 is a square of a bidegree-(1,1) polynomial.

These lemmas are unconditional algebra, not a construction of
`GrowthMain.CurveModel` and not a proof of the Wronskian bound.
They permit the unit-ideal and geometry problems to be separated.

Candidate for pinned Lean 4.19.0 build; DO NOT label kernel checked
before CI has compiled it and printed its axioms.
-/

namespace ICFallDegrees.GrowthRationalUnit

open ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- The target t=1 is always a two-torsion abscissa in this family. -/
theorem target_is_root (a : K) : cubic a (-1 - a) 1 = 0 := by
  unfold cubic
  ring

/-- Exact S3 factorization with target t=1 and b=-1-a. -/
theorem s3_square_family (a x z : K) :
    s3 a (-1 - a) x z 1 = ((x - 1) * z - (x + 2 + a)) ^ 2 := by
  unfold s3
  ring

/-- A root of the special S3 necessarily determines a from x and z. -/
theorem parameter_of_s3_zero (a x z : K)
    (h : s3 a (-1 - a) x z 1 = 0) :
    a = (x - 1) * z - x - 2 := by
  rw [s3_square_family] at h
  have hm : ((x - 1) * z - (x + 2 + a)) *
      ((x - 1) * z - (x + 2 + a)) = 0 := by
    simpa only [pow_two] using h
  have hlin : (x - 1) * z - (x + 2 + a) = 0 :=
    (mul_eq_zero.mp hm).elim id id
  have heq : (x - 1) * z = x + 2 + a := sub_eq_zero.mp hlin
  calc
    a = (x - 1) * z - x - 2 := by rw [heq]; ring

/-- If both candidate abscissae come from a subfield but the parameter
does not, S3 cannot vanish. In applications the subfield contains all
N-th roots of unity, and this becomes a concrete unit-ideal criterion. -/
theorem s3_ne_zero_of_parameter_not_in_subfield
    {B : Type*} [Field B] [Algebra B K]
    (a : K) (u v : B)
    (hout : ∀ c : B, algebraMap B K c ≠ a) :
    s3 a (-1 - a) (algebraMap B K u) (algebraMap B K v) 1 ≠ 0 := by
  intro h
  have hp := parameter_of_s3_zero a (algebraMap B K u)
    (algebraMap B K v) h
  have hmap : algebraMap B K ((u - 1) * v - u - 2) = a := by
    calc
      algebraMap B K ((u - 1) * v - u - 2) =
          (algebraMap B K u - 1) * algebraMap B K v -
            algebraMap B K u - 2 := by simp only [map_sub, map_mul, map_one, map_ofNat]
      _ = a := hp.symm
  exact hout _ hmap

/-- The nonsingularity discriminant factors in the rational family.
It vanishes only at the two exceptional parameters a=-3 and 4a=-3. -/
theorem discriminant_factor (a : K) :
    4 * a ^ 3 + 27 * (-1 - a) ^ 2 =
      (a + 3) ^ 2 * (4 * a + 3) := by
  ring

theorem nonsingular_of_parameter_conditions (a : K)
    (h3 : a + 3 ≠ 0) (h4 : 4 * a + 3 ≠ 0) :
    4 * a ^ 3 + 27 * (-1 - a) ^ 2 ≠ 0 := by
  rw [discriminant_factor]
  exact mul_ne_zero (pow_ne_zero _ h3) h4

end ICFallDegrees.GrowthRationalUnit
