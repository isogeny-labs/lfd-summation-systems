import ICFallDegrees.GrowthFixedE6Existence
import ICFallDegrees.ScalarExtensionTransport
import ICFallDegrees.SelectorExtension
import Mathlib.Tactic

/-!
# N8: unconditional rational E6 bounded-closure nonmembership

The fixed complex E6 proof is unconditional and kernel-checked. Since
the coefficients of the E6 input family are rational, and a rational
InClosure derivation would remain a valid derivation after coefficient
extension to C, its complex nonmembership theorem descends to Q.

This proves a missing exact *characteristic-zero rational* starting point
for finite-stage integer-matrix specialization. No prime-field
nonmembership, unit ideal or LFD claim is made. In particular it does
not exclude the exact admissible p=17,k=2 obstruction.

The descent is by the existing kernel theorem
ScalarExtensionTransport.map_closure and the exact input-image identity
SelectorExtension.image_inputs; no rational-complex equivalence hypothesis.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeRationalSeed

open ICFallDegrees.PrimeM2General
open ICFallDegrees.GrowthMain
open ICFallDegrees.GrowthRatFuncModel

/-- The actual rational E6 degree-D closure does not contain 1 when
k >= 14*D+10. This is a theorem about Q-polynomials, not a
statement about the reduction modulo *any* prime. -/
theorem rational_e6_not_one_linear
    (k : ℕ) (depth : 2 ≤ k) (D : ℕ)
    (hD : 1 ≤ D) (hk : 14 * D + 10 ≤ k) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1))
      D (1 : MvPolynomial (Fin (2 * k)) ℚ) := by
  have hcomplex :
      ¬ InClosure
        (Set.range (PrimeM2General.inputs k depth (6 : ℂ) (-7 : ℂ) 1))
        D (1 : MvPolynomial (Fin (2 * k)) ℂ) :=
    GrowthMain.growth_theorem_linear
      (K := ℂ) (F₀ := ℂ) (F := RatFunc ℂ)
      (by norm_num) k depth (6 : ℂ) (-7 : ℂ) (1 : ℂ)
      (GrowthRatFuncModel.xi (6 : ℂ)) 0 GrowthFixedE6.model
      (Or.inl rfl) D hD hk
  intro hrational
  have hlift := ScalarExtensionTransport.map_closure (K := ℂ)
    (Set.range (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1))
    D hrational
  rw [SelectorExtension.image_inputs (K := ℚ) (L := ℂ)
    k depth (6 : ℚ) (-7 : ℚ) (1 : ℚ)] at hlift
  apply hcomplex
  simpa [ScalarExtensionTransport.lift] using hlift

/-- A native rational nonmembership theorem with the final N8
linear threshold, ready for integer-matrix finite-stage specialization. -/
theorem rational_e6_not_one_at_linear_window
    (k : ℕ) (depth : 2 ≤ k) (hk : 24 ≤ k) :
    ¬ InClosure
      (Set.range (PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1))
      ((k - 10) / 14) (1 : MvPolynomial (Fin (2 * k)) ℚ) := by
  have hD : 1 ≤ (k - 10) / 14 := by omega
  have hthreshold : 14 * ((k - 10) / 14) + 10 ≤ k := by omega
  exact rational_e6_not_one_linear k depth _ hD hthreshold

#print axioms rational_e6_not_one_linear
#print axioms rational_e6_not_one_at_linear_window

end ICFallDegrees.GrowthPrimeRationalSeed
