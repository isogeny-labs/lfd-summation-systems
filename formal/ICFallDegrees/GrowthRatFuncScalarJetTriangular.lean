import ICFallDegrees.GrowthRatFuncHigherLeibniz
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Tactic

/-!
# N8: determinant of the scalar-multiplication jet matrix

For every h in RatFunc K, multiplying a family of rational functions by h
acts on their first r derivative rows by a *lower triangular* operator.
The higher Leibniz identity (proved in GrowthRatFuncHigherLeibniz)
supplies its entries. This module proves the independent matrix
structure: its determinant is h^r, in every characteristic.

The still-required connecting identity is
  derivativeMatrix r (fun j => h * f j) = scalarJet r h * derivativeMatrix r f
which would finish the common-scalar Wronskian determinant law.
No assertion about Wronskian nonvanishing or spread is made.
-/

noncomputable section

namespace ICFallDegrees.GrowthRatFuncScalarJetTriangular

open ICFallDegrees.GrowthRatFuncIteratedLinear

variable {K : Type*} [Field K]

/-- Lower triangular jet action of multiplication by a rational function h. -/
def scalarJet (r : ℕ) (h : RatFunc K) : Matrix (Fin r) (Fin r) (RatFunc K) :=
  Matrix.of (fun i j =>
    if j.val ≤ i.val then
      (i.val.choose j.val : RatFunc K) * iterDerivative (i.val - j.val) h
    else 0)

/-- The scalar-jet matrix has no entries strictly above its diagonal. -/
theorem scalarJet_blockTriangular (r : ℕ) (h : RatFunc K) :
    (scalarJet r h).BlockTriangular OrderDual.toDual := by
  intro i j hij
  change i < j at hij
  have hji : ¬ j.val ≤ i.val := by omega
  simp only [scalarJet, Matrix.of_apply, if_neg hji]

/-- Each diagonal entry is exactly h (since choose(i,i)=1). -/
theorem scalarJet_diag (r : ℕ) (h : RatFunc K) (i : Fin r) :
    scalarJet r h i i = h := by
  simp [scalarJet]

/-- The determinant is h^r without any nonzero assumption on h.
This is the determinant factor required to clear one common rational
denominator from all Wronskian columns. -/
theorem scalarJet_det (r : ℕ) (h : RatFunc K) :
    (scalarJet r h).det = h ^ r := by
  classical
  rw [Matrix.det_of_lowerTriangular (scalarJet r h)
    (scalarJet_blockTriangular r h)]
  simp [scalarJet_diag, Finset.prod_const]

end ICFallDegrees.GrowthRatFuncScalarJetTriangular
