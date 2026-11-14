import ICFallDegrees.GrowthRatFuncCofactorReplacement
import ICFallDegrees.GrowthRatFuncCofactorDegree
import Mathlib.Tactic

/-!
# N8 genuine Wronskian maximal-minor replacement

The columns here are actual rational functions and their *genuine*
iterated derivatives. Replacing one member of an indexed independent
family by a constant-linear combination has an exactly computable
Wronskian determinant. Normalization transfers the identity to the
polynomial cofactors of the two-torsion monomials.

This is a necessary cofactor comparison for arbitrary minimal circuits,
not yet a theorem enumerating all omitted circuit indices and not a
proof of GrowthSpread.SpreadBound.
-/

noncomputable section

namespace ICFallDegrees.GrowthActualMinorReplacement

open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncCofactorReplacement

variable {K : Type*} [Field K]

/-- Updating the rational-function family at an index is literally
updating the corresponding column of its genuine derivative matrix. -/
theorem derivativeMatrix_replace (r : ℕ) (f : Fin r → RatFunc K)
    (g : RatFunc K) (i : Fin r) :
    derivativeMatrix r (fun j => if j = i then g else f j) =
      (derivativeMatrix r f).updateCol i (fun k => iterDerivative k.val g) := by
  classical
  ext k j
  by_cases hj : j = i
  · subst j
    simp [derivativeMatrix]
  · simp [derivativeMatrix, hj]

/-- Exact all-rank omitted-column determinant comparison when the new
column is a genuine constant-linear combination of the original ones. -/
theorem determinant_replace_of_const_sum (r : ℕ)
    (f : Fin r → RatFunc K) (coeff : Fin r → K)
    (g : RatFunc K) (hg : g = ∑ j : Fin r, coeff j • f j)
    (i : Fin r) :
    (derivativeMatrix r (fun j => if j = i then g else f j)).det =
      algebraMap K (RatFunc K) (coeff i) * (derivativeMatrix r f).det := by
  rw [derivativeMatrix_replace]
  exact determinant_replaced_column_of_const_sum r f coeff g hg i

/-- Embed the actual polynomial Wronskian cofactor determinant. -/
def cofactorValue (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ) : RatFunc K :=
  (Matrix.of (fun i j =>
    algebraMap (Polynomial K) (RatFunc K) (polynomialMatrix r c e i j))).det

def monomialProduct (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ) : RatFunc K :=
  ∏ j : Fin r, monomial c (e j)

def rowDenominatorProduct (r : ℕ) (c : K) : RatFunc K :=
  ∏ i : Fin r,
    (algebraMap (Polynomial K) (RatFunc K)
      (ICFallDegrees.GrowthRationalLogNumerator.denominator c)) ^ i.val

/-- Actual cofactor factorization, with convenient names for the factors. -/
theorem cofactor_factorization (r : ℕ) (c : K) (e : Fin r → ℕ × ℕ) :
    rowDenominatorProduct r c *
        (derivativeMatrix r (fun j => monomial c (e j))).det =
      cofactorValue r c e * monomialProduct r c e := by
  exact concrete_determinant_factorization r c e

/-- The determinant-replacement law *after polynomial normalization*,
using actual rational monomials and no assumed Wronskian identity.
An arbitrary constant-linear relation is the only input. -/
theorem cofactor_replacement_product (r : ℕ) (c : K)
    (e : Fin r → ℕ × ℕ) (enew : ℕ × ℕ)
    (coeff : Fin r → K) (i : Fin r)
    (hrel : monomial c enew =
      ∑ j : Fin r, coeff j • monomial c (e j)) :
    let er : Fin r → ℕ × ℕ := fun j => if j = i then enew else e j
    cofactorValue r c er * monomialProduct r c er =
      algebraMap K (RatFunc K) (coeff i) *
        (cofactorValue r c e * monomialProduct r c e) := by
  dsimp only
  have hdet := determinant_replace_of_const_sum r
    (fun j => monomial c (e j)) coeff (monomial c enew) hrel i
  have hfactorNew := cofactor_factorization r c
    (fun j => if j = i then enew else e j)
  have hfactorOld := cofactor_factorization r c e
  have hfun :
      (fun j : Fin r => monomial c (if j = i then enew else e j)) =
        (fun j : Fin r => if j = i then monomial c enew else monomial c (e j)) := by
    funext j
    split_ifs <;> rfl
  have hsame :
      (derivativeMatrix r (fun j =>
          monomial c (if j = i then enew else e j))).det =
      (derivativeMatrix r (fun j =>
          if j = i then monomial c enew else monomial c (e j))).det := by
    rw [hfun]
  calc
    cofactorValue r c (fun j => if j = i then enew else e j) *
      monomialProduct r c (fun j => if j = i then enew else e j) =
        rowDenominatorProduct r c *
          (derivativeMatrix r (fun j =>
            monomial c (if j = i then enew else e j))).det := hfactorNew.symm
    _ = rowDenominatorProduct r c *
          (derivativeMatrix r (fun j =>
            if j = i then monomial c enew else monomial c (e j))).det := by rw [hsame]
    _ = rowDenominatorProduct r c *
          (algebraMap K (RatFunc K) (coeff i) *
            (derivativeMatrix r (fun j => monomial c (e j))).det) := by rw [hdet]
    _ = algebraMap K (RatFunc K) (coeff i) *
          (cofactorValue r c e * monomialProduct r c e) := by
      rw [← hfactorOld]
      ring

end ICFallDegrees.GrowthActualMinorReplacement
