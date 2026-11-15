import ICFallDegrees.GrowthRatFuncModel
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# N8: Möbius involution and distinguished rational poles

For the two-torsion target family E_a : y² = x³ + a*x - 1 - a,
the reduced Semaev component has z = (x + 2 + a)/(x - 1).

Here we show directly that the fractional-linear coordinate change is
an involution away from the singular parameter a = -3, and make its
simple-pole decomposition explicit.

These are genuine rational function identities needed to transport a
valuation at X=∞ to one at X=1; they do not construct valuations, prove
Wronskian estimates, or discharge GrowthMain.CurveModel. Do not label
the full growth theorem unconditional.
-/

noncomputable section

namespace ICFallDegrees.GrowthRationalMobius

open GrowthRationalCoordinates
open GrowthRatFuncModel

variable {K : Type*} [Field K]

/-- Clearing the denominator of the rational companion coordinate. -/
theorem mul_rationalZ (a x : K) (hx : x ≠ 1) :
    (x - 1) * rationalZ a x = x + 2 + a := by
  unfold rationalZ
  have h : x - 1 ≠ 0 := sub_ne_zero.mpr hx
  field_simp [h]

/-- A normalized local parameter at X=1: the second coordinate has
the exact principal part (a+3)/(X-1). -/
theorem rationalZ_sub_one (a x : K) (hx : x ≠ 1) :
    rationalZ a x - 1 = (a + 3) / (x - 1) := by
  unfold rationalZ
  have h : x - 1 ≠ 0 := sub_ne_zero.mpr hx
  field_simp [h]
  ring

/-- The distinguished target pole X=1 cannot also map to z=1
when a+3 is nonzero (the nonsingular branch). -/
theorem rationalZ_ne_one (a x : K) (hx : x ≠ 1)
    (ha : a + 3 ≠ 0) : rationalZ a x ≠ 1 := by
  intro heq
  have hmul := mul_rationalZ a x hx
  rw [heq] at hmul
  have hlin : x - 1 = x + 2 + a := by simpa using hmul
  apply ha
  calc
    a + 3 = (x + 2 + a) - (x - 1) := by ring
    _ = 0 := sub_eq_zero.mpr hlin.symm

/-- The Möbius companion transformation is its own inverse.
This explicitly rules out collapse of the two generic coordinates. -/
theorem rationalZ_involutive (a x : K) (hx : x ≠ 1)
    (ha : a + 3 ≠ 0) :
    rationalZ a (rationalZ a x) = x := by
  let z := rationalZ a x
  have hz : z ≠ 1 := rationalZ_ne_one a x hx ha
  have hfirst : (x - 1) * z = x + 2 + a :=
    mul_rationalZ a x hx
  have hsecond : (z - 1) * rationalZ a z = z + 2 + a :=
    mul_rationalZ a z hz
  have hswap : (z - 1) * x = z + 2 + a := by
    calc
      (z - 1) * x = (x - 1) * z + (z - x) := by ring
      _ = z + 2 + a := by rw [hfirst]; ring
  have hzero : (z - 1) * (rationalZ a z - x) = 0 := by
    calc
      (z - 1) * (rationalZ a z - x) =
          (z - 1) * rationalZ a z - (z - 1) * x := by ring
      _ = 0 := sub_eq_zero.mpr (hsecond.trans hswap.symm)
  have h := (mul_eq_zero.mp hzero).resolve_left (sub_ne_zero.mpr hz)
  exact sub_eq_zero.mp h

/-- The involution on the actual rational-function coordinate, not merely
on an abstract field satisfying an on-curve assumption. -/
theorem ratFunc_xi2_involutive (a : K)
    (ha : (RatFunc.C a : RatFunc K) + 3 ≠ 0) :
    rationalZ (RatFunc.C a) (xi a (1 : Fin 2)) = RatFunc.X := by
  simpa [xi] using
    (rationalZ_involutive (RatFunc.C a) (RatFunc.X : RatFunc K)
      X_ne_one ha)

end ICFallDegrees.GrowthRationalMobius
