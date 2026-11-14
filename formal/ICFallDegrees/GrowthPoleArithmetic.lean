import Mathlib.Data.Int.Basic
import Mathlib.Tactic.Ring

/-!
# N8: exact integer bridge from pole orders to circuit spread

This isolates the final arithmetic step of the genus-one function-field
Wronskian argument (N8 Theorem 3.3). It is independent of the curve model:
the required pairwise *height bound* and the two pole valuations are
explicit inputs, rather than being hidden in a new structure.

The remaining analytic task is to construct these orders from the
function field and to prove the height bound; the present lemma does not
establish `GrowthSpread.SpreadBound` on its own.

The proof works with the modulus p=0 and with positive characteristics.
-/

namespace ICFallDegrees.GrowthPoleArithmetic

/-- If the order differences at the two distinguished places are congruent
to minus twice the exponent differences modulo p, and their combined
absolute value is at most H, then the exact distance-to-pZ bound follows.
The witnesses are the integers describing the p-divisible errors. -/
theorem spread_from_two_pole_orders
    (p : ℕ) (a a' b b' : ℕ) (H o r co cr : ℤ)
    (ho : o = -2 * ((a : ℤ) - a') + (p : ℤ) * co)
    (hr : r = -2 * ((b : ℤ) - b') + (p : ℤ) * cr)
    (hheight : |o| + |r| ≤ H) :
    ∃ z₁ z₂ : ℤ,
      |2 * ((a : ℤ) - a') - (p : ℤ) * z₁| +
      |2 * ((b : ℤ) - b') - (p : ℤ) * z₂| ≤ H := by
  refine ⟨co, cr, ?_⟩
  have hO : 2 * ((a : ℤ) - a') - (p : ℤ) * co = -o := by
    rw [ho]
    ring
  have hR : 2 * ((b : ℤ) - b') - (p : ℤ) * cr = -r := by
    rw [hr]
    ring
  rw [hO, hR, abs_neg, abs_neg]
  exact hheight

/-- The last step of Theorem 3.3, with the genus-one Wronskian bound.
The preceding geometric/Wronskian estimate is exactly `hheight`.
This does not claim that it is proved for the N8 curve. -/
theorem spread_from_genus_one_height
    (p m : ℕ) (a a' b b' : ℕ) (o r co cr : ℤ)
    (ho : o = -2 * ((a : ℤ) - a') + (p : ℤ) * co)
    (hr : r = -2 * ((b : ℤ) - b') + (p : ℤ) * cr)
    (hheight : |o| + |r| ≤ 6 * ((m : ℤ) - 1) * ((m : ℤ) - 2)) :
    ∃ z₁ z₂ : ℤ,
      |2 * ((a : ℤ) - a') - (p : ℤ) * z₁| +
      |2 * ((b : ℤ) - b') - (p : ℤ) * z₂| ≤
        6 * ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
  exact spread_from_two_pole_orders p a a' b b'
    (6 * ((m : ℤ) - 1) * ((m : ℤ) - 2)) o r co cr ho hr hheight

end ICFallDegrees.GrowthPoleArithmetic
