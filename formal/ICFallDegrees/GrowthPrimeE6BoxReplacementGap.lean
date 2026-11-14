import ICFallDegrees.GrowthPrimeE6BoxCofactorRatio
import ICFallDegrees.GrowthFiniteCofactorGap
import Mathlib.Tactic

/-!
# N8: actual E6 positive-characteristic bounded replacement exponent gap

The characteristic-free genuine X=1 root-multiplicity and infinity-pole
gap theorems turn the newly proved bounded positive-characteristic E6
cofactor fraction into the sharp four-times-cofactor-degree exponent
spread. The concrete exponent box and ringChar>=2N are explicit.
Neither SpreadBound nor a circuit-transfer assumption occurs here.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6BoxReplacementGap

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthPrimeE6BoxCofactor
open ICFallDegrees.GrowthPrimeE6BoxCofactorRatio
open ICFallDegrees.GrowthFiniteCofactorGap

variable {K : Type*} [Field K]

/-- Actual independent bounded E6 monomials, an exact replacement
relation with a nonzero coefficient, and the large-characteristic
guard imply the true four-budget ℓ¹ exponent separation inequality. -/
theorem e6_box_replacement_l1_le_four_budget
    (r N : ℕ) (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial (8 : K) enew =
      ∑ j : Fin r, coeff j • monomial (8 : K) (e j))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K)
    (h9 : (8 : K) + 1 ≠ 0) :
    |2 * ((enew.1 : ℤ) - (e i).1)| +
      |2 * ((enew.2 : ℤ) - (e i).2)| ≤
        4 * ((∑ j : Fin r, 2 * j.val : ℕ) : ℤ) := by
  let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
  let A : Polynomial K := (polynomialMatrix r (8 : K) e).det
  let B : Polynomial K := (polynomialMatrix r (8 : K) er).det
  have hA : A ≠ 0 :=
    e6_cofactor_ne_zero_of_box r N e hbox hind hchar
  have hq := e6_cofactor_ratio_of_box_replacement r N e enew
    coeff i hbox hind hcoeff hrel hchar
  dsimp only at hq
  obtain ⟨hB, hratio⟩ := hq
  have hdA : A.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r (8 : K) e
  have hdB : B.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r (8 : K) er
  exact exponent_pair_gap_l1_le_four_budget
    (8 : K) (coeff i) h9 enew.1 enew.2 (e i).1 (e i).2
    A B (∑ j : Fin r, 2 * j.val)
    hA hB hcoeff hdA hdB hratio

#print axioms e6_box_replacement_l1_le_four_budget

end ICFallDegrees.GrowthPrimeE6BoxReplacementGap
