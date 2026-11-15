import ICFallDegrees.GrowthSecantFactorization
import Mathlib.Tactic.FieldSimp

noncomputable section

/-!
# N8: the target-2-torsion Möbius model

For a nonsingular short Weierstrass curve and a target abscissa `t`
with `t^3 + a*t + b = 0`, the reduced Semaev curve is a rational
graph.  Its second abscissa can be expressed without adjoining a
generic ordinate:

  z = (t*x + 2*t^2 + a)/(x-t).

The map is an involution away from its pole when `3*t^2+a ≠ 0`.
This avoids treating the doubled Semaev plane scheme as a smooth
genus-one curve. These are algebraic lemmas only: they do NOT construct
the function field, prove valuations/Wronskians, or discharge
`GrowthMain.CurveModel`.

Lean 4.19 candidate; requires the pinned CI build and axiom audit.
-/

namespace ICFallDegrees.GrowthTwoTorsion

open ICFallDegrees.GrowthSecant
open ICFallDegrees.GrowthSecantFactorization

variable {K : Type*} [Field K]

/-- The companion abscissa for a target with ordinate zero. -/
def mobius (a t x : K) : K :=
  (t * x + 2 * t ^ 2 + a) / (x - t)

/-- The defining linear equation of the reduced Semaev component. -/
theorem mul_mobius (a t x : K) (hx : x ≠ t) :
    (x - t) * mobius a t x = t * x + 2 * t ^ 2 + a := by
  unfold mobius
  field_simp [sub_ne_zero.mpr hx]

/-- The Möbius coordinate satisfies the *actual* S3 equation,
not merely its squared resultant, even when the target is 2-torsion. -/
theorem mobius_on_s3 (a b t x : K)
    (ht : cubic a b t = 0) (hx : x ≠ t) :
    s3 a b x (mobius a t x) t = 0 := by
  have hlin :
      (x - t) * mobius a t x - (t * x + 2 * t ^ 2 + a) = 0 := by
    rw [mul_mobius a t x hx]
    ring
  rw [s3_perfect_square_of_target_root a b x (mobius a t x) t ht, hlin]
  ring

/-- Bridge to the exact N8 predicate.  The circuit spread bound
is not asserted here. -/
theorem mobius_on_S3at (a b t x : K)
    (ht : cubic a b t = 0) (hx : x ≠ t) :
    GrowthSpread.S3at a b t
      (fun i : Fin 2 => if i = 0 then x else mobius a t x) = 0 := by
  rw [s3_eq_S3at]
  exact mobius_on_s3 a b t x ht hx

/-- The pole and target are separated by the tangent coefficient.
The nonzero derivative is forced by nonsingularity when t is a root. -/
theorem mobius_target_difference (a t x : K) (hx : x ≠ t) :
    (x - t) * (mobius a t x - t) = 3 * t ^ 2 + a := by
  rw [mul_sub, mul_mobius a t x hx]
  ring

theorem mobius_ne_target (a t x : K) (hx : x ≠ t)
    (hder : 3 * t ^ 2 + a ≠ 0) : mobius a t x ≠ t := by
  intro heq
  have h := mobius_target_difference a t x hx
  rw [heq, sub_self, mul_zero] at h
  exact hder h.symm

/-- Translation by the 2-torsion point acts as an involution on
the abscissa, away from the pole. -/
theorem mobius_involutive (a t x : K) (hx : x ≠ t)
    (hder : 3 * t ^ 2 + a ≠ 0) :
    mobius a t (mobius a t x) = x := by
  let z : K := mobius a t x
  have hz : z ≠ t := mobius_ne_target a t x hx hder
  have h1 : (x - t) * z = t * x + 2 * t ^ 2 + a :=
    mul_mobius a t x hx
  have h2 : (z - t) * mobius a t z = t * z + 2 * t ^ 2 + a :=
    mul_mobius a t z hz
  have hzero : (z - t) * (mobius a t z - x) = 0 := by
    calc
      _ = (z - t) * mobius a t z - (z - t) * x := by ring
      _ = (t * z + 2 * t ^ 2 + a) - (z - t) * x := by rw [h2]
      _ = (t * x + 2 * t ^ 2 + a) - (x - t) * z := by ring
      _ = 0 := by rw [h1]; ring
  have hfinal : mobius a t z = x :=
    sub_eq_zero.mp ((mul_eq_zero.mp hzero).resolve_left (sub_ne_zero.mpr hz))
  simpa only [z] using hfinal

/-- Explicit fixed characteristic-zero nonsingular curve and target.
This does *not* claim that the chain ideal is a unit ideal. -/
theorem fixed_curve_target_root :
    cubic (10 : ℚ) (-11) 1 = 0 := by
  norm_num [cubic]

theorem fixed_curve_nonsingular :
    4 * (10 : ℚ) ^ 3 + 27 * (-11 : ℚ) ^ 2 ≠ 0 := by
  norm_num

theorem fixed_curve_s3 (x z : ℚ) :
    s3 (10 : ℚ) (-11) x z 1 = (x * z - z - x - 12) ^ 2 := by
  rw [s3_perfect_square_of_target_root (10 : ℚ) (-11) x z 1
    fixed_curve_target_root]
  ring

end ICFallDegrees.GrowthTwoTorsion
