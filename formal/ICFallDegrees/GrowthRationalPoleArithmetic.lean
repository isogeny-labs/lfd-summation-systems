import ICFallDegrees.GrowthPoleArithmetic
import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.Ring.Abs

/-!
# N8 rational two-torsion model: pole-order-to-spread bridge

The target-2-torsion construction in `GrowthTwoTorsion` is a rational
map on P1. At its two distinguished poles, the coordinate functions
have orders -1 and 0, rather than -2 and 0 on the elliptic curve.
The genus-zero four-place Wronskian estimate would give a stronger
height bound: `|o| + |r| <= 2*(m-1)*(m-2)`.

This file *only* formalizes the exact integer consequence of those
explicit hypotheses. It does not assert that the place orders or
Wronskian height estimate have been proved for the rational function
field. No CurveModel instance or analytic axiom is introduced.
-/

namespace ICFallDegrees.GrowthRationalPoleArithmetic

/-- A rational four-place Wronskian height bound and two simple pole
congruences imply the N8 `SpreadBound` inequality. The factor 6
required by N8 has slack: the hypotheses yield factor 4. -/
theorem spread_from_simple_poles
    (p m : ℕ) (hm : 2 ≤ m)
    (a a' b b' : ℕ) (o r co cr : ℤ)
    (ho : o = -((a : ℤ) - a') + (p : ℤ) * co)
    (hr : r = -((b : ℤ) - b') + (p : ℤ) * cr)
    (hheight : |o| + |r| ≤
      2 * ((m : ℤ) - 1) * ((m : ℤ) - 2)) :
    ∃ z₁ z₂ : ℤ,
      |2 * ((a : ℤ) - a') - (p : ℤ) * z₁| +
      |2 * ((b : ℤ) - b') - (p : ℤ) * z₂| ≤
        6 * ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
  have hmZ : (2 : ℤ) ≤ m := by exact_mod_cast hm
  have hprod : 0 ≤ ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
    apply mul_nonneg <;> omega
  have ho2 : 2 * o = -2 * ((a : ℤ) - a') + (p : ℤ) * (2 * co) := by
    rw [ho]
    ring
  have hr2 : 2 * r = -2 * ((b : ℤ) - b') + (p : ℤ) * (2 * cr) := by
    rw [hr]
    ring
  have htwice : |2 * o| + |2 * r| = 2 * (|o| + |r|) := by
    simp [abs_mul]
    ring
  have hheight2 : |2 * o| + |2 * r| ≤
      6 * ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
    rw [htwice]
    nlinarith [hheight, hprod]
  exact GrowthPoleArithmetic.spread_from_two_pole_orders
    p a a' b b' (6 * ((m : ℤ) - 1) * ((m : ℤ) - 2))
    (2 * o) (2 * r) (2 * co) (2 * cr) ho2 hr2 hheight2

end ICFallDegrees.GrowthRationalPoleArithmetic
