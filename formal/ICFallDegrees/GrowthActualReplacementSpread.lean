import ICFallDegrees.GrowthActualCofactorQuotient
import ICFallDegrees.GrowthFiniteCofactorGap
import Mathlib.Tactic

/-!
# N8: sharp unconditional exponent spread for a genuine replacement relation

No cofactor ratio, nonzero determinant, pole identities or global
SpreadBound appears among this theorem's hypotheses. All are supplied
by actual Wronskian and rational-function lemmas.

This theorem applies to any independent r-family plus a genuine
constant-linear relation representing a new monomial with a nonzero
coefficient at the chosen index. Reindexing an arbitrary GrowthSpread
circuit into such a family remains the next open interface.
-/

noncomputable section

namespace ICFallDegrees.GrowthActualReplacementSpread

open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorDegree
open ICFallDegrees.GrowthActualCofactorNonzero
open ICFallDegrees.GrowthActualCofactorQuotient
open ICFallDegrees.GrowthFiniteCofactorGap

variable {K : Type*} [Field K] [CharZero K]

/-- Genuine replacement relations have the sharp 4·Σ(2j) exponent
distance bound. The Wronskian cofactors, their nonvanishing, their
degree bounds, the fraction identity and both pole comparisons are
all obtained from previously proved results rather than hypothesized. -/
theorem replacement_l1_le_four_cofactor_budget
    (r : ℕ) (c : K) (hc : c + 1 ≠ 0)
    (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hind : LinearIndependent K (fun j => monomial c (e j)))
    (hcoeff : coeff i ≠ 0)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j)) :
    |2 * ((enew.1 : ℤ) - (e i).1)| +
      |2 * ((enew.2 : ℤ) - (e i).2)| ≤
        4 * ((∑ j : Fin r, 2 * j.val : ℕ) : ℤ) := by
  let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
  let A : Polynomial K := (polynomialMatrix r c e).det
  let B : Polynomial K := (polynomialMatrix r c er).det
  have hA : A ≠ 0 :=
    polynomial_cofactor_ne_zero_of_linearIndependent r c e hind
  have hq := cofactor_ratio_of_independent_replacement r c e enew
    coeff i hind hcoeff hrel
  dsimp only at hq
  obtain ⟨hB, hratio⟩ := hq
  have hdA : A.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r c e
  have hdB : B.natDegree ≤ ∑ j : Fin r, 2 * j.val :=
    concrete_cofactor_natDegree_le_two r c er
  exact exponent_pair_gap_l1_le_four_budget
    c (coeff i) hc enew.1 enew.2 (e i).1 (e i).2
    A B (∑ j : Fin r, 2 * j.val)
    hA hB hcoeff hdA hdB hratio

end ICFallDegrees.GrowthActualReplacementSpread
