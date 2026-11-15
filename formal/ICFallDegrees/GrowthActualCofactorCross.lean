import ICFallDegrees.GrowthActualMinorReplacement
import ICFallDegrees.GrowthRatFuncInfinityDegree
import Mathlib.Tactic

/-!
# N8: cancel genuine nonzero common Wronskian monomial factors

The prior exact determinant-replacement factorization contains the
product of all rational monomial columns. Here the products are compared
after changing one column, and the genuinely nonzero shared factors are
cancelled to obtain the cofactor identity for that column. A linear
relation among the actual rational functions is the only hypothesis.

A further quotient/circuit-enumeration bridge is still required before
GrowthSpread.SpreadBound, and none is taken as an assumption.
-/

noncomputable section

namespace ICFallDegrees.GrowthActualCofactorCross

open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthActualMinorReplacement
open ICFallDegrees.GrowthRatFuncInfinityDegree

variable {K : Type*} [Field K]

/-- A product of genuine rational monomials changes by exactly one
factor when one exponent column is replaced. This is true even if
some factors vanish, although actual monomials here never do. -/
theorem monomialProduct_replace (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ) (i : Fin r) :
    monomialProduct r c (fun j => if j = i then enew else e j) *
        monomial c (e i) =
      monomialProduct r c e * monomial c enew := by
  classical
  let rest : RatFunc K :=
    ∏ j ∈ (Finset.univ.erase i), monomial c (e j)
  have hrest :
      (∏ j ∈ (Finset.univ.erase i),
        monomial c (if j = i then enew else e j)) = rest := by
    dsimp [rest]
    apply Finset.prod_congr rfl
    intro j hj
    have hji : j ≠ i := (Finset.mem_erase.mp hj).1
    simp [hji]
  have hnew :
      monomialProduct r c (fun j => if j = i then enew else e j) =
        monomial c enew * rest := by
    unfold monomialProduct
    calc
      (∏ j : Fin r, monomial c (if j = i then enew else e j)) =
          monomial c enew *
            (∏ j ∈ (Finset.univ.erase i),
              monomial c (if j = i then enew else e j)) := by
                simpa using
                  (Finset.mul_prod_erase (Finset.univ : Finset (Fin r))
                    (fun j : Fin r =>
                      monomial c (if j = i then enew else e j))
                    (Finset.mem_univ i)).symm
      _ = monomial c enew * rest := by rw [hrest]
  have hold :
      monomialProduct r c e = monomial c (e i) * rest := by
    unfold monomialProduct
    exact (Finset.mul_prod_erase (Finset.univ : Finset (Fin r))
      (fun j : Fin r => monomial c (e j))
      (Finset.mem_univ i)).symm
  rw [hnew, hold]
  ring

/-- Every factor *except* the replaced one is nonzero, since the
original monomial is a nonzero rational function. -/
theorem other_monomials_ne_zero (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) (i : Fin r) :
    (∏ j ∈ (Finset.univ.erase i), monomial c (e j)) ≠ 0 := by
  apply Finset.prod_ne_zero_iff.mpr
  intro j hj
  exact monomial_ne_zero c (e j).1 (e j).2

/-- Exact comparison of the ACTUAL polynomial Wronskian cofactors
after cancelling the product of unchanged, nonzero rational columns.

No cofactor nonvanishing, no characteristic-zero Wronskian criterion,
no pole identity and no spread bound is assumed here. -/
theorem cofactor_cross_of_const_relation (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j)) :
    let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
    cofactorValue r c er * monomial c enew =
      algebraMap K (RatFunc K) (coeff i) *
        cofactorValue r c e * monomial c (e i) := by
  dsimp only
  let rest : RatFunc K := ∏ j ∈ (Finset.univ.erase i), monomial c (e j)
  have hrest : rest ≠ 0 := other_monomials_ne_zero r c e i
  have hnew :
      monomialProduct r c (fun j => if j = i then enew else e j) =
        monomial c enew * rest := by
    have hp := monomialProduct_replace r c e enew i
    have hi : monomial c (e i) ≠ 0 :=
      monomial_ne_zero c (e i).1 (e i).2
    have hold :
        monomialProduct r c e = monomial c (e i) * rest := by
      unfold monomialProduct
      exact (Finset.mul_prod_erase (Finset.univ : Finset (Fin r))
        (fun j : Fin r => monomial c (e j))
        (Finset.mem_univ i)).symm
    rw [hold] at hp
    have hcancel :
        (monomialProduct r c (fun j => if j = i then enew else e j)) =
          monomial c enew * rest := by
      apply (mul_right_cancel₀ hi)
      calc
        monomialProduct r c (fun j => if j = i then enew else e j) *
          monomial c (e i) =
            (monomial c (e i) * rest) * monomial c enew := hp
        _ = (monomial c enew * rest) * monomial c (e i) := by ring
    exact hcancel
  have hold :
      monomialProduct r c e = monomial c (e i) * rest := by
    unfold monomialProduct
    exact (Finset.mul_prod_erase (Finset.univ : Finset (Fin r))
      (fun j : Fin r => monomial c (e j))
      (Finset.mem_univ i)).symm
  have hprod := cofactor_replacement_product r c e enew coeff i hrel
  dsimp only at hprod
  rw [hnew, hold] at hprod
  apply (mul_right_cancel₀ hrest)
  calc
    (cofactorValue r c (fun j => if j = i then enew else e j) *
        monomial c enew) * rest =
      cofactorValue r c (fun j => if j = i then enew else e j) *
        (monomial c enew * rest) := by ring
    _ = algebraMap K (RatFunc K) (coeff i) *
      (cofactorValue r c e * (monomial c (e i) * rest)) := hprod
    _ = (algebraMap K (RatFunc K) (coeff i) *
      cofactorValue r c e * monomial c (e i)) * rest := by ring

end ICFallDegrees.GrowthActualCofactorCross
