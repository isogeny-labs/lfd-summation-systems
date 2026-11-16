import ICFallDegrees.GrowthSecantFactorization

noncomputable section

/-!
# Recover the elliptic coordinate from a non-two-torsion S3 target

This supplies the inverse of the algebraic secant coordinate away from x=t.
It isolates why the rational-function section used by the E6 proof cannot
simply be reused when the target ordinate is nonzero. No circuit-spread,
valuation, CurveModel, or general growth theorem is claimed here.
-/

namespace ICFallDegrees.GrowthNonTorsionBridge

open GrowthSecant GrowthSecantFactorization

variable {K : Type*} [Field K]

def recoverY (a b x z t Y : K) : K :=
  ((x - t) ^ 2 * z - companion a b x t) / (2 * Y)

/-- Every S3 point for a target with nonzero ordinate recovers a point on E. -/
theorem recoverY_square (a b x z t Y : K)
    (h2 : (2 : K) ≠ 0) (hY : Y ≠ 0)
    (hR : Y ^ 2 = cubic a b t) (hS : s3 a b x z t = 0) :
    recoverY a b x z t Y ^ 2 = cubic a b x := by
  have hs := square_of_s3_zero a b x z t hS
  rw [← hR] at hs
  unfold recoverY
  rw [div_pow]
  apply (div_eq_iff (pow_ne_zero 2 (mul_ne_zero h2 hY))).2
  calc
    _ = 4 * cubic a b x * Y ^ 2 := hs
    _ = cubic a b x * (2 * Y) ^ 2 := by ring

/-- Recovering the ordinate in the actual secant formula is exact. -/
theorem recoverY_secantX (a b x t y Y : K)
    (h2 : (2 : K) ≠ 0) (hY : Y ≠ 0) (hxt : x ≠ t) :
    recoverY a b x (secantX a b x t y Y) t Y = y := by
  unfold recoverY secantX secantNumerator
  field_simp [sub_ne_zero.mpr hxt, h2, hY] <;> ring

/-- In the other direction the secant formula recovers the S3 abscissa. -/
theorem secantX_recoverY (a b x z t Y : K)
    (h2 : (2 : K) ≠ 0) (hY : Y ≠ 0) (hxt : x ≠ t) :
    secantX a b x t (recoverY a b x z t Y) Y = z := by
  unfold secantX secantNumerator recoverY
  field_simp [sub_ne_zero.mpr hxt, h2, hY] <;> ring

/-- A rational section over a field where the cubic is nonsquare is impossible.
For the generic x in K(X), the remaining polynomial-degree argument must
establish this nonsquare hypothesis; it is not a spread assumption. -/
theorem no_section_of_nonsquare (a b x z t Y : K)
    (h2 : (2 : K) ≠ 0) (hY : Y ≠ 0)
    (hR : Y ^ 2 = cubic a b t)
    (hnonsquare : ∀ y : K, y ^ 2 ≠ cubic a b x) :
    s3 a b x z t ≠ 0 := by
  intro hs
  exact hnonsquare _ (recoverY_square a b x z t Y h2 hY hR hs)

/-- Nonzero constant term in the S3 quadratic forces its companion nonzero. -/
theorem secantX_ne_zero (a b x t y Y : K)
    (hP : y ^ 2 = cubic a b x) (hR : Y ^ 2 = cubic a b t)
    (hxt : x ≠ t) (hc : companionConstant a b x t ≠ 0) :
    secantX a b x t y Y ≠ 0 := by
  intro hz
  have hs := secant_s3 a b x t y Y hP hR hxt
  rw [s3_quadratic, hz] at hs
  simp only [zero_pow (by decide : 2 ≠ 0), mul_zero, sub_zero, zero_add] at hs
  exact hc hs

#print axioms recoverY_square
#print axioms recoverY_secantX
#print axioms secantX_recoverY
#print axioms no_section_of_nonsquare
#print axioms secantX_ne_zero

end ICFallDegrees.GrowthNonTorsionBridge
