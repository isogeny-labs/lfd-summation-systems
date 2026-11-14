import ICFallDegrees.GrowthPrimeE6BoxClearing
import ICFallDegrees.GrowthRatFuncCofactorDegree
import Mathlib.Tactic

/-!
# N8: genuine E6 cofactor nonvanishing below the characteristic

For the actual E6 rational section over any field, finite box exponents
and independence give nonzero true rational Wronskians from
GrowthPrimeE6BoxClearing. The existing determinant normalization then
gives nonzero polynomial cofactors of the exact E6 monomials.

The polynomial cofactor degree estimate holds in all characteristics.
These are concrete inputs for a positive-characteristic replacement
ratio and spread theorem; no SpreadBound or circuit transfer is assumed.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6BoxCofactor

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthPrimeE6BoxClearing

variable {K : Type*} [Field K]

/-- Independence of any subfamily of actual E6 monomials with
exponents strictly below N implies a genuinely nonzero polynomial
Wronskian cofactor in every field of characteristic zero or >=2N. -/
theorem e6_cofactor_ne_zero_of_box
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (polynomialMatrix r (8 : K) e).det ≠ 0 := by
  have hW : (derivativeMatrix r (fun j => monomial (8 : K) (e j))).det ≠ 0 :=
    box_wronskian_ne_zero_of_box_char_bound r N e hbox hind hchar
  have hmap := embedded_cofactor_ne_zero_of_wronskian_ne_zero r (8 : K) e hW
  let φ : Polynomial K →+* RatFunc K :=
    algebraMap (Polynomial K) (RatFunc K)
  have hmatrix :
      (Matrix.of fun i j : Fin r => φ (polynomialMatrix r (8 : K) e i j)) =
        φ.mapMatrix (polynomialMatrix r (8 : K) e) := by
    ext i j
    rfl
  rw [hmatrix, ← RingHom.map_det] at hmap
  intro hzero
  exact hmap (by rw [hzero, map_zero])

/-- The same actual cofactor has a characteristic-independent tight
degree budget, alongside proved nonvanishing. -/
theorem e6_cofactor_nonzero_and_degree_of_box
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (polynomialMatrix r (8 : K) e).det ≠ 0 ∧
    (polynomialMatrix r (8 : K) e).det.natDegree ≤
      ∑ j : Fin r, 2 * j.val := by
  exact ⟨e6_cofactor_ne_zero_of_box r N e hbox hind hchar,
    concrete_cofactor_natDegree_le_two r (8 : K) e⟩

#print axioms e6_cofactor_ne_zero_of_box
#print axioms e6_cofactor_nonzero_and_degree_of_box

end ICFallDegrees.GrowthPrimeE6BoxCofactor
