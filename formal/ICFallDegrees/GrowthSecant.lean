import ICFallDegrees.GrowthSpread
import Mathlib.Tactic.FieldSimp

noncomputable section

/-!
# N8 growth: explicit algebraic Semaev companion-coordinate identity

This proves the algebraic part of the concrete function-field curve model.
For points P=(x,y), R=(t,Y) on the same short Weierstrass curve
with x ≠ t, the chord formula gives a companion abscissa x(R-P),
and the third summation polynomial vanishes at (x, x(R-P), t).

**Scope:** These are algebraic identities over an arbitrary field,
not a construction of the function field and not a proof of
`GrowthSpread.SpreadBound`. The full N8 growth theorem remains conditional
on `GrowthMain.CurveModel` until the geometric/valuation part is proved.

Verification: the identities were independently expanded symbolically and
tested on nonsingular curves over F_101, F_103, F_127, F_131, F_251.
The Lean file itself requires a pinned lake build/axiom audit before its
results can be called kernel-checked.
-/

namespace ICFallDegrees.GrowthSecant

variable {K : Type*} [Field K]

/-- The cubic on the right of the short Weierstrass equation. -/
def cubic (a b x : K) : K := x ^ 3 + a * x + b

/-- Mixed coefficient of Semaev's S3 as a quadratic in its middle variable. -/
def companion (a b x t : K) : K :=
  (x + t) * (x * t + a) + 2 * b

/-- Constant coefficient of that quadratic. -/
def companionConstant (a b x t : K) : K :=
  (x * t - a) ^ 2 - 4 * b * (x + t)

/-- The short Weierstrass third summation polynomial. -/
def s3 (a b x z t : K) : K :=
  (x - z) ^ 2 * t ^ 2 -
      2 * ((x + z) * (x * z + a) + 2 * b) * t +
    (x * z - a) ^ 2 - 4 * b * (x + z)

/-- S3 is symmetric, in particular quadratic in the middle abscissa. -/
theorem s3_quadratic (a b x z t : K) :
    s3 a b x z t =
      (x - t) ^ 2 * z ^ 2 -
        2 * companion a b x t * z + companionConstant a b x t := by
  unfold s3 companion companionConstant
  ring

/-- Its discriminant identity. -/
theorem companion_discriminant (a b x t : K) :
    companion a b x t ^ 2 -
      (x - t) ^ 2 * companionConstant a b x t =
    4 * cubic a b x * cubic a b t := by
  unfold companion companionConstant cubic
  ring

/-- The numerator in the formula for x(R-P). -/
def secantNumerator (a b x t y Y : K) : K :=
  companion a b x t + 2 * y * Y

/-- The companion abscissa away from the diagonal x=t. -/
def secantX (a b x t y Y : K) : K :=
  secantNumerator a b x t y Y / (x - t) ^ 2

/-- A polynomial identity under the two Weierstrass equations. -/
theorem secant_numerator (a b x t y Y : K)
    (hP : y ^ 2 = cubic a b x)
    (hR : Y ^ 2 = cubic a b t) :
    secantNumerator a b x t y Y ^ 2 -
      2 * companion a b x t * secantNumerator a b x t y Y +
      (x - t) ^ 2 * companionConstant a b x t = 0 := by
  calc
    _ = 4 * (y ^ 2 * Y ^ 2 - cubic a b x * cubic a b t) := by
      unfold secantNumerator companion companionConstant cubic
      ring
    _ = 0 := by rw [hP, hR]; ring

/-- **Concrete Semaev identity for the chord formula.**
This is the algebraic on-curve component needed by the N8 CurveModel
construction, but does not yet provide a generic function-field element. -/
theorem secant_s3 (a b x t y Y : K)
    (hP : y ^ 2 = cubic a b x)
    (hR : Y ^ 2 = cubic a b t)
    (hne : x ≠ t) :
    s3 a b x (secantX a b x t y Y) t = 0 := by
  let z : K := secantX a b x t y Y
  let L : K := secantNumerator a b x t y Y
  have hsub : x - t ≠ 0 := sub_ne_zero.mpr hne
  have hz : (x - t) ^ 2 * z = L := by
    dsimp [z, L, secantX]
    field_simp [hsub] <;> ring
  have hnum : L ^ 2 - 2 * companion a b x t * L +
      (x - t) ^ 2 * companionConstant a b x t = 0 := by
    simpa only [L] using secant_numerator a b x t y Y hP hR
  have hq :
      (x - t) ^ 2 *
        ((x - t) ^ 2 * z ^ 2 -
          2 * companion a b x t * z + companionConstant a b x t) = 0 := by
    calc
      _ = ((x - t) ^ 2 * z) ^ 2 -
        2 * companion a b x t * ((x - t) ^ 2 * z) +
        (x - t) ^ 2 * companionConstant a b x t := by ring
      _ = L ^ 2 - 2 * companion a b x t * L +
        (x - t) ^ 2 * companionConstant a b x t := by rw [hz]
      _ = 0 := hnum
  have hnz : (x - t) ^ 2 ≠ 0 := pow_ne_zero _ hsub
  have hzq :
      (x - t) ^ 2 * z ^ 2 -
        2 * companion a b x t * z + companionConstant a b x t = 0 :=
    (mul_eq_zero.mp hq).resolve_left hnz
  change s3 a b x z t = 0
  rw [s3_quadratic]
  exact hzq

/-! ## Bridge to the exact polynomial used by the N8 Lean development -/

/-- The scalar S3 expression is exactly `GrowthSpread.S3at` evaluated with
coordinates x and z. This is an algebraic identity, not a CurveModel instance. -/
theorem s3_eq_S3at (a b x z t : K) :
    GrowthSpread.S3at a b t (fun i : Fin 2 => if i = 0 then x else z) =
      s3 a b x z t := by
  simp [GrowthSpread.S3at, s3]

/-- The explicit companion coordinate satisfies the **same** S3at equation
required in `GrowthMain.CurveModel.on_curve`. The function-field construction,
genericity and spread estimate remain open. -/
theorem secant_on_curve (a b x t y Y : K)
    (hP : y ^ 2 = cubic a b x)
    (hR : Y ^ 2 = cubic a b t)
    (hne : x ≠ t) :
    GrowthSpread.S3at a b t
      (fun i : Fin 2 => if i = 0 then x else secantX a b x t y Y) = 0 := by
  rw [s3_eq_S3at]
  exact secant_s3 a b x t y Y hP hR hne

end ICFallDegrees.GrowthSecant
