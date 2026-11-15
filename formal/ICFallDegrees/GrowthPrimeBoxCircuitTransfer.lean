import ICFallDegrees.GrowthBoxSpread
import ICFallDegrees.GrowthActualCircuitSpread
import ICFallDegrees.GrowthRatFuncModel
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# N8 finite-box circuit specialization bridge

For the fixed E6 rational two-torsion model, a *one-way* transfer of
actual bounded minimal dependencies from characteristic p to Q suffices
to establish the finite-box spread bound over ZMod p.

The premise is deliberately the finite circuit-transfer condition, not
a disguised SpreadBound. It must ultimately be proved from the explicit
integer coefficient matrix for the N8 bounded monomials and the finite
nonzero-minor theorem. No characteristic-p Wronskian theorem is assumed.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeBoxCircuitTransfer

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthBoxSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthActualCircuitSpread

/-- If all minimal characteristic-p dependencies in the exponent box
survive as rational minimal dependencies, then the true rational E6
Wronskian inequality provides BoxSpreadBound over the prime field.

The only missing premise is a finite-matroid specialization statement
for the actual bounded rational monomials. -/
theorem fixed_e6_box_spread_of_circuit_transfer
    (p k : ℕ) [Fact p.Prime]
    (htransfer :
      ∀ J : Finset (Fin 2 →₀ ℕ),
        (∀ d ∈ J, d 0 < 2 ^ k ∧ d 1 < 2 ^ k) →
        Circuit (ZMod p)
          (mono (GrowthRatFuncModel.xi (6 : ZMod p))) J →
        Circuit ℚ
          (mono (GrowthRatFuncModel.xi (6 : ℚ))) J) :
    BoxSpreadBound (ZMod p)
      (GrowthRatFuncModel.xi (6 : ZMod p)) p k := by
  intro J hbox hc x hx y hy
  have hq : Circuit ℚ
      (mono (GrowthRatFuncModel.xi (6 : ℚ))) J :=
    htransfer J hbox hc
  have hspread :
      SpreadBound ℚ (GrowthRatFuncModel.xi (6 : ℚ)) 0 :=
    GrowthActualCircuitSpread.actual_spread_bound (6 : ℚ) (by norm_num)
  obtain ⟨z₁,z₂,h⟩ := hspread J hq x hx y hy
  refine ⟨0,0,?_⟩
  simpa only [Nat.cast_zero, zero_mul, mul_zero, sub_zero] using h

#print axioms fixed_e6_box_spread_of_circuit_transfer

end ICFallDegrees.GrowthPrimeBoxCircuitTransfer
