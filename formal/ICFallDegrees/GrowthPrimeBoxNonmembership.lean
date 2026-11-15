import ICFallDegrees.GrowthPrimeBoxCircuitTransfer
import ICFallDegrees.GrowthMain
import Mathlib.Tactic

/-!
# N8: native prime-field E6 bounded-closure nonmembership from finite-box spread

This is an actual prime-field InClosure conclusion for the normalized E6
chain. The only circuit premise is BoxSpreadBound in the exponent box.
The accompanying corollary uses the precise one-way finite-matroid
transfer condition instead. This is not an unconditional all-prime theorem:
the finite-matroid preservation condition remains to be established
outside a finite set of bad primes.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeBoxNonmembership

open MvPolynomial
open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthBoxSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthPrimeBoxCircuitTransfer
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GrowthReduction

/-- The true finite-box spread inequality on the fixed E6 rational model
implies native mod-p closure nonmembership at linear depth. This proof
uses no characteristic-p Wronskian or global SpreadBound hypothesis. -/
theorem e6_not_one_of_box_spread
    (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (D : ℕ) (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p)
    (hspread : BoxSpreadBound (ZMod p)
      (GrowthRatFuncModel.xi (6 : ZMod p)) p k) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth
        (6 : ZMod p) (-7 : ZMod p) 1)) D
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) := by
  let ξ : Fin 2 → RatFunc (ZMod p) := GrowthRatFuncModel.xi (6 : ZMod p)
  have hdata := GrowthRatFuncModel.ratFunc_algebraic_data (6 : ZMod p)
  have hcurve : GrowthSpread.S3at
      (6 : ZMod p) (-7 : ZMod p) 1 ξ = 0 := by
    simpa only [ξ, show (-1 - (6 : ZMod p)) = (-7 : ZMod p) by ring]
      using hdata.2
  have hlarge : 5 ≤ p := by
    have hpow : 2 ^ 2 ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) depth
    omega
  have htwo : (2 : ZMod p) ≠ 0 :=
    PrimeFieldConclusions.two_ne_zero p hlarge
  have hfei : FEI
      (RingHom.ker (aeval ξ :
        MvPolynomial (Fin 2) (ZMod p) →ₐ[ZMod p] RatFunc (ZMod p)))
      k (D - 1) := by
    apply GrowthBoxSpread.fei_of_box_spread
      (K := ZMod p) (F₀ := ZMod p) (F := RatFunc (ZMod p))
      ξ p k (D - 1) (13 * (D - 1) + 23) (Or.inr hp)
      hspread (GrowthCounting.lambdaBound_linear (D - 1))
      (by omega) (by omega)
  exact GrowthReduction.not_one_of_fei htwo k depth
    (6 : ZMod p) (-7 : ZMod p) 1
    (RingHom.ker (aeval ξ :
       MvPolynomial (Fin 2) (ZMod p) →ₐ[ZMod p] RatFunc (ZMod p)))
    (GrowthMain.X_not_mem_ker hdata.1)
    (GrowthMain.S3poly_mem_ker hcurve)
    D hfei

/-- The remaining specialized mathematical input is the one-way
preservation of all bounded minimal circuit supports from Fp to Q.
A finite integer coefficient matrix and good-prime minor theorem should
eventually eliminate this hypothesis outside finitely many primes. -/
theorem e6_not_one_of_bounded_circuit_transfer
    (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k)
    (D : ℕ) (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k)
    (hp : 2 * 2 ^ k + 1 ≤ p)
    (htransfer :
      ∀ J : Finset (Fin 2 →₀ ℕ),
        (∀ d ∈ J, d 0 < 2 ^ k ∧ d 1 < 2 ^ k) →
        Circuit (ZMod p)
          (mono (GrowthRatFuncModel.xi (6 : ZMod p))) J →
        Circuit ℚ (mono (GrowthRatFuncModel.xi (6 : ℚ))) J) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth
        (6 : ZMod p) (-7 : ZMod p) 1)) D
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) := by
  exact e6_not_one_of_box_spread p k depth D hD hk hp
    (fixed_e6_box_spread_of_circuit_transfer p k htransfer)

#print axioms e6_not_one_of_box_spread
#print axioms e6_not_one_of_bounded_circuit_transfer

end ICFallDegrees.GrowthPrimeBoxNonmembership
