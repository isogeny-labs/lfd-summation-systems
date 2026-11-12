import ICFallDegrees.Closure

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K]

/-- A span with no nonzero rows below the budget is already closed: only
constant multipliers can fit. This is a useful exact lower-bound criterion. -/
theorem closure_le_top_span (F : Set (MvPolynomial σ K)) (D : ℕ)
    (W : Submodule K (MvPolynomial σ K))
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → f ∈ W)
    (no_low : ∀ f ∈ W, f ≠ 0 → D ≤ f.totalDegree) :
    closureSpace F D ≤ W := by
  apply closure_minimal F D W inputs
  intro a f hf budget
  by_cases hz : f = 0
  · simp [hz]
  have hmin := no_low f hf hz
  have hdegree : a.totalDegree = 0 := by omega
  have ha := totalDegree_eq_zero_iff_eq_C.mp hdegree
  rw [ha, C_mul']
  exact W.smul_mem _ hf

/-- Evaluation is an independently checked exclusion invariant when every
input available at this budget vanishes at the specified point. -/
theorem eval_zero_of_closure (F : Set (MvPolynomial σ K)) (D : ℕ)
    (point : σ → K)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → eval point f = 0)
    {f} (hf : InClosure F D f) : eval point f = 0 := by
  induction hf with
  | input member hd => exact inputs _ member hd
  | zero => simp
  | add _ _ ih₁ ih₂ => simp [map_add, ih₁, ih₂]
  | smul c _ ih => simp [smul_eq_C_mul, ih]
  | mul a _ _ ih => simp [map_mul, ih]

end ICFallDegrees
