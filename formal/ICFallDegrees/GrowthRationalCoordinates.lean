import ICFallDegrees.GrowthRationalUnit
import Mathlib.Tactic.FieldSimp

/-!
# N8 rational section: concrete on-curve coordinates

The t=1 two-torsion unit-ideal family admits the rational coordinate
z=(x+2+a)/(x-1). This proves the algebraic on-curve and nonzero fields
of a CurveModel at every generic x avoiding 0,1,-2-a.

It does NOT establish that x is a transcendental coordinate in K(X),
construct the rational-function field, or prove SpreadBound. Those are
separate formalization obligations.
-/

noncomputable section
namespace ICFallDegrees.GrowthRationalCoordinates

open ICFallDegrees.GrowthSecant
open ICFallDegrees.GrowthRationalUnit

variable {K : Type*} [Field K]

def rationalZ (a x : K) : K := (x + 2 + a) / (x - 1)

/-- Explicit rational section lies on the reduced Semaev curve. -/
theorem rational_s3_on_curve (a x : K) (hx : x ≠ 1) :
    s3 a (-1 - a) x (rationalZ a x) 1 = 0 := by
  rw [s3_square_family]
  have hsub : x - 1 ≠ 0 := sub_ne_zero.mpr hx
  have hz : (x - 1) * rationalZ a x = x + 2 + a := by
    dsimp [rationalZ]
    field_simp [hsub]
  rw [hz]
  ring

/-- The second rational coordinate is nonzero away from its zero/pole. -/
theorem rationalZ_ne_zero (a x : K)
    (hx : x ≠ 1) (hnum : x + 2 + a ≠ 0) :
    rationalZ a x ≠ 0 := by
  unfold rationalZ
  exact div_ne_zero hnum (sub_ne_zero.mpr hx)

/-- Both algebraic fields of CurveModel, without the analytic spread field.
This explicitly does not build a CurveModel. -/
theorem rational_algebraic_data (a x : K)
    (hzero : x ≠ 0) (hone : x ≠ 1) (hnum : x + 2 + a ≠ 0) :
    (∀ i : Fin 2, (if i = 0 then x else rationalZ a x) ≠ 0) ∧
      GrowthSpread.S3at a (-1 - a) 1
        (fun i : Fin 2 => if i = 0 then x else rationalZ a x) = 0 := by
  constructor
  · intro i
    fin_cases i
    · simpa using hzero
    · simpa using rationalZ_ne_zero a x hone hnum
  · rw [GrowthSecant.s3_eq_S3at]
    exact rational_s3_on_curve a x hone

end ICFallDegrees.GrowthRationalCoordinates
