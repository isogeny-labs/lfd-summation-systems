import ICFallDegrees.GrowthSecant

noncomputable section

/-!
# N8: algebraic double-cover identity for the Semaev curve

Pure algebra. Does not construct a function field or discharge CurveModel.
Candidate Lean 4.19.0 source: compilation and axiom audit required.
-/

namespace ICFallDegrees.GrowthSecantFactorization

open ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- S3 is a quadratic double cover after clearing the diagonal denominator.
The identity holds even when the target ordinate is zero. -/
theorem s3_square_factorization (a b x z t : K) :
    (x - t) ^ 2 * s3 a b x z t =
      ((x - t) ^ 2 * z - companion a b x t) ^ 2 -
        4 * cubic a b x * cubic a b t := by
  rw [s3_quadratic]
  calc
    (x - t) ^ 2 *
        ((x - t) ^ 2 * z ^ 2 -
          2 * companion a b x t * z + companionConstant a b x t) =
        ((x - t) ^ 2 * z - companion a b x t) ^ 2 -
          (companion a b x t ^ 2 -
            (x - t) ^ 2 * companionConstant a b x t) := by ring
    _ = ((x - t) ^ 2 * z - companion a b x t) ^ 2 -
        4 * cubic a b x * cubic a b t := by
      rw [companion_discriminant]

/-- Every S3 point satisfies a square equation. -/
theorem square_of_s3_zero (a b x z t : K)
    (h : s3 a b x z t = 0) :
    ((x - t) ^ 2 * z - companion a b x t) ^ 2 =
      4 * cubic a b x * cubic a b t := by
  have hfac := s3_square_factorization a b x z t
  rw [h, mul_zero] at hfac
  exact sub_eq_zero.mp hfac.symm

/-- The converse away from the diagonal. -/
theorem s3_zero_of_square (a b x z t : K) (hne : x ≠ t)
    (h : ((x - t) ^ 2 * z - companion a b x t) ^ 2 =
      4 * cubic a b x * cubic a b t) :
    s3 a b x z t = 0 := by
  have hfac := s3_square_factorization a b x z t
  rw [h, sub_self] at hfac
  have hpow : (x - t) ^ 2 ≠ 0 := pow_ne_zero _ (sub_ne_zero.mpr hne)
  exact (mul_eq_zero.mp hfac).resolve_left hpow

/-- The target 2-torsion case is a doubled section. It must not be
silently excluded when constructing the concrete function-field model. -/
theorem s3_square_of_target_root (a b x z t : K)
    (ht : cubic a b t = 0) :
    (x - t) ^ 2 * s3 a b x z t =
      ((x - t) ^ 2 * z - companion a b x t) ^ 2 := by
  rw [s3_square_factorization, ht, mul_zero, sub_zero]


/-- If the target is a 2-torsion point, S3 is *globally a perfect square*.
Thus the plane S3 scheme is nonreduced in this case: the geometric
CurveModel must take the reduced linear component rather than treating
S3 itself as an irreducible genus-one plane curve. -/
theorem s3_perfect_square_of_target_root (a b x z t : K)
    (ht : cubic a b t = 0) :
    s3 a b x z t =
      ((x - t) * z - (t * x + 2 * t ^ 2 + a)) ^ 2 := by
  have hb : b = -(t ^ 3 + a * t) := by
    dsimp [cubic] at ht
    calc
      b = (t ^ 3 + a * t + b) - (t ^ 3 + a * t) := by ring
      _ = 0 - (t ^ 3 + a * t) := by rw [ht]
      _ = -(t ^ 3 + a * t) := by ring
  rw [hb]
  unfold s3
  ring


/-- At a 2-torsion target on a nonsingular short Weierstrass cubic,
the linear factor of S3 has a nonvanishing numerator at x=t. This is
the algebraic input for the order-two pole of xi2 at R. -/
theorem two_torsion_derivative_ne_of_nonsingular (a b t : K)
    (ht : cubic a b t = 0)
    (hdisc : 4 * a ^ 3 + 27 * b ^ 2 ≠ 0) :
    3 * t ^ 2 + a ≠ 0 := by
  intro hder
  have ha : a = -(3 * t ^ 2) := by
    calc
      a = (3 * t ^ 2 + a) - 3 * t ^ 2 := by ring
      _ = 0 - 3 * t ^ 2 := by rw [hder]
      _ = -(3 * t ^ 2) := by ring
  have ht' : t ^ 3 + (-(3 * t ^ 2)) * t + b = 0 := by
    dsimp [cubic] at ht
    simpa only [ha] using ht
  have hb : b = 2 * t ^ 3 := by
    calc
      b = (t ^ 3 + (-(3 * t ^ 2)) * t + b) -
          (t ^ 3 + (-(3 * t ^ 2)) * t) := by ring
      _ = 0 - (t ^ 3 + (-(3 * t ^ 2)) * t) := by rw [ht']
      _ = 2 * t ^ 3 := by ring
  apply hdisc
  rw [ha, hb]
  ring

end ICFallDegrees.GrowthSecantFactorization
