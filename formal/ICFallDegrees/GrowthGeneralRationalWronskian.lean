import ICFallDegrees.GrowthPolynomialWronskianBasis
import ICFallDegrees.GrowthRatFuncDenominatorIndependence
import ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
import ICFallDegrees.GrowthRatFuncWronskianCommonScalar
import Mathlib.Tactic

/-!
# N8: general characteristic-zero rational Wronskian independence criterion

This composes *proved* interfaces, without assuming Wronskian nonvanishing:
* a K-linearly independent finite polynomial family is the tautological
  Basis.span of its span;
* the pinned degree-echelon basis theorem proves its Wronskian nonzero;
* common-denominator clearing preserves K-linear independence;
* the polynomial Wronskian embeds in the actual RatFunc derivation; and
* multiplying every column by one nonzero rational factor changes the
  Wronskian by its r-th power.

No new assumption about N8 circuit spread or last-fall degree is introduced.
-/

noncomputable section

namespace ICFallDegrees.GrowthGeneralRationalWronskian

open Polynomial
open ICFallDegrees.GrowthRatFuncDenominatorIndependence
open ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
open ICFallDegrees.GrowthRatFuncWronskianCommonScalar
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K] [CharZero K]

/-- The ordinary polynomial Wronskian of *any* finite linearly
independent polynomial family is nonzero in characteristic zero.
The family need not initially have pairwise distinct degrees. -/
theorem polynomial_wronskian_ne_zero_of_linearIndependent
    (r : ℕ) (p : Fin r → Polynomial K) (hp : LinearIndependent K p) :
    (polynomialDerivativeMatrix r p).det ≠ 0 := by
  classical
  let B : Submodule K (Polynomial K) := Submodule.span K (Set.range p)
  let bas : Basis (Fin r) K B := Basis.span hp
  have hbas (j : Fin r) : ((bas j : B) : Polynomial K) = p j := by
    exact Basis.span_apply hp j
  let k : ℕ := (Finset.univ.sup fun j : Fin r => (p j).natDegree) + 1
  have hdeg : ∀ j : Fin r, ((bas j : B) : Polynomial K).natDegree < k := by
    intro j
    rw [hbas j]
    have hle : (p j).natDegree ≤
        Finset.univ.sup (fun i : Fin r => (p i).natDegree) := by
      exact Finset.le_sup (s := (Finset.univ : Finset (Fin r)))
        (f := fun i : Fin r => (p i).natDegree) (Finset.mem_univ j)
    exact Nat.lt_succ_of_le hle
  have hwr := ICFallDegrees.GrowthPolynomialWronskianBasis.polynomial_wronskian_ne_zero_of_basis bas hdeg
  change (polynomialDerivativeMatrix r (fun j => ((bas j : B) : Polynomial K))).det ≠ 0 at hwr
  have heq : (fun j : Fin r => ((bas j : B) : Polynomial K)) = p := by
    funext j
    exact hbas j
  rw [heq] at hwr
  exact hwr

/-- A bona fide K-linearly independent family in the *actual* rational
function field K(X) has a nonzero all-column ordinary rational Wronskian.
No nonvanishing premise about any cofactor/circuit is assumed. -/
theorem rational_wronskian_ne_zero_of_linearIndependent
    (r : ℕ) (f : Fin r → RatFunc K) (hf : LinearIndependent K f) :
    (derivativeMatrix r f).det ≠ 0 := by
  obtain ⟨H, p, hH, hp, hind⟩ :=
    exists_independent_polynomial_representatives r f hf
  have hpoly : (polynomialDerivativeMatrix r p).det ≠ 0 :=
    polynomial_wronskian_ne_zero_of_linearIndependent r p hind
  have hrat :
      (derivativeMatrix r (fun j =>
        algebraMap (Polynomial K) (RatFunc K) (p j))).det ≠ 0 :=
    (determinant_polynomialMap_ne_zero_iff r p).2 hpoly
  let h : RatFunc K := algebraMap (Polynomial K) (RatFunc K) H
  have hh : h ≠ 0 := RatFunc.algebraMap_ne_zero hH
  have hscaled :
      (derivativeMatrix r (fun j => h * f j)).det ≠ 0 := by
    have heq : (fun j : Fin r => h * f j) =
        (fun j => algebraMap (Polynomial K) (RatFunc K) (p j)) := by
      funext j
      exact hp j
    rw [heq]
    exact hrat
  exact (determinant_mul_scalar_ne_zero_iff r h f hh).mp hscaled

end ICFallDegrees.GrowthGeneralRationalWronskian
