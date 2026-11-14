import ICFallDegrees.GrowthPoleArithmetic
import Mathlib.Tactic

/-!
# N8 characteristic-zero cofactor spread arithmetic

This isolates the purely integer part of the rational Wronskian argument.
No rational derivative, Wronskian nonvanishing, divisor, or valuation
identity is asserted here. Those mathematical bridges remain necessary
for a concrete `GrowthSpread.SpreadBound`.
-/

namespace ICFallDegrees.GrowthComplexCofactorArithmetic

/-- If two nonzero polynomial cofactors have degrees and multiplicities
at X=1 in the interval [0,B], and the exponent differences are
identified with the respective degree/multiplicity differences,
then their doubled L1 distance is at most 4B.

The degree/multiplicity identifications are *inputs*: constructing them
from genuine rational Wronskian cofactors remains an open Lean task. -/
theorem doubled_l1_le_four_budget
    (u u' v v' d d' o o' B : ℤ)
    (hd : 0 ≤ d ∧ d ≤ B)
    (hd' : 0 ≤ d' ∧ d' ≤ B)
    (ho : 0 ≤ o ∧ o ≤ B)
    (ho' : 0 ≤ o' ∧ o' ≤ B)
    (hu : u - u' = d - d')
    (hv : v - v' = o' - o) :
    2 * |u - u'| + 2 * |v - v'| ≤ 4 * B := by
  have hdu : -B ≤ u - u' ∧ u - u' ≤ B := by omega
  have hdv : -B ≤ v - v' ∧ v - v' ≤ B := by omega
  have hau : |u - u'| ≤ B := abs_le.mpr hdu
  have hav : |v - v'| ≤ B := abs_le.mpr hdv
  omega

/-- The specific cofactor budget B ≤ (m-1)(m-2) gives the stronger
4(m-1)(m-2) rational circuit-spread bound. This is only an
arithmetic implication, not a proof that such cofactors exist. -/
theorem doubled_l1_le_circuit_budget
    (m : ℕ) (hm : 2 ≤ m)
    (u u' v v' d d' o o' B : ℤ)
    (hd : 0 ≤ d ∧ d ≤ B)
    (hd' : 0 ≤ d' ∧ d' ≤ B)
    (ho : 0 ≤ o ∧ o ≤ B)
    (ho' : 0 ≤ o' ∧ o' ≤ B)
    (hu : u - u' = d - d')
    (hv : v - v' = o' - o)
    (hB : B ≤ ((m : ℤ) - 1) * ((m : ℤ) - 2)) :
    2 * |u - u'| + 2 * |v - v'| ≤
      4 * ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
  have h := doubled_l1_le_four_budget u u' v v' d d' o o' B
    hd hd' ho ho' hu hv
  nlinarith

/-- The same hypotheses imply the 6(m-1)(m-2) inequality
appearing in GrowthSpread.SpreadBound at characteristic zero. -/
theorem doubled_l1_le_required_budget
    (m : ℕ) (hm : 2 ≤ m)
    (u u' v v' d d' o o' B : ℤ)
    (hd : 0 ≤ d ∧ d ≤ B)
    (hd' : 0 ≤ d' ∧ d' ≤ B)
    (ho : 0 ≤ o ∧ o ≤ B)
    (ho' : 0 ≤ o' ∧ o' ≤ B)
    (hu : u - u' = d - d')
    (hv : v - v' = o' - o)
    (hB : B ≤ ((m : ℤ) - 1) * ((m : ℤ) - 2)) :
    2 * |u - u'| + 2 * |v - v'| ≤
      6 * ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
  have h := doubled_l1_le_circuit_budget m hm
    u u' v v' d d' o o' B hd hd' ho ho' hu hv hB
  have hm' : 0 ≤ ((m : ℤ) - 1) * ((m : ℤ) - 2) := by
    have hm'' : (2 : ℤ) ≤ m := by exact_mod_cast hm
    nlinarith
  nlinarith

end ICFallDegrees.GrowthComplexCofactorArithmetic
