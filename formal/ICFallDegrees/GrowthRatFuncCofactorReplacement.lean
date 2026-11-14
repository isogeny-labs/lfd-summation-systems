import ICFallDegrees.GrowthRatFuncWronskianNormalization
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

/-!
# N8: maximal Wronskian cofactor under a constant-linear column replacement

For a genuine rational Wronskian, replacing a column by a K-constant
linear combination of all original columns multiplies its determinant
by the corresponding coefficient. This is a determinant-level bridge
to the maximal-minor comparison for a minimal circuit.

No independence, valuation, pole or circuit spread assumption is used.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncCofactorReplacement

open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K]

/-- Exact determinant comparison when a Wronskian column is replaced by a
genuine rational function known to be a constant-linear combination of
the original columns. This holds for *any* finite number of columns and
for the actual iterated quotient-rule derivation, not a formal derivative
matrix with assumed identities. -/
theorem determinant_replaced_column_of_const_sum
    (r : ℕ) (f : Fin r → RatFunc K) (c : Fin r → K)
    (g : RatFunc K)
    (hg : g = ∑ j : Fin r, c j • f j)
    (i : Fin r) :
    ((derivativeMatrix r f).updateCol i
        (fun k => iterDerivative k.val g)).det =
      algebraMap K (RatFunc K) (c i) * (derivativeMatrix r f).det := by
  classical
  have hder (k : Fin r) :
      iterDerivative k.val g =
        ∑ j : Fin r, algebraMap K (RatFunc K) (c j) *
          (derivativeMatrix r f) k j := by
    rw [hg, iterDerivative_const_sum k.val Finset.univ c f]
    simp only [derivativeMatrix, Matrix.of_apply, Algebra.smul_def]
  have heq :
      (fun k : Fin r => iterDerivative k.val g) =
        (fun k : Fin r => ∑ j : Fin r,
          algebraMap K (RatFunc K) (c j) *
            (derivativeMatrix r f) k j) := by
    funext k
    exact hder k
  have hdet := Matrix.det_updateCol_sum (derivativeMatrix r f) i
    (fun j : Fin r => (algebraMap K (RatFunc K) (c j) : RatFunc K))
  rw [heq]
  simpa only [smul_eq_mul] using hdet

end ICFallDegrees.GrowthRatFuncCofactorReplacement
