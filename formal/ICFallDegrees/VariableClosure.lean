import ICFallDegrees.Closure
import Mathlib.Data.Finsupp.Weight
import Mathlib.Tactic.Ring

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K] [Nontrivial K]

/-- Checking variable products of all rows below the budget suffices for
closure under every permitted polynomial multiplier. -/
theorem variable_closed_mul_mem
    (W : Submodule K (MvPolynomial σ K)) (D : ℕ)
    (variable_stable : ∀ i p, p∈W → p.totalDegree<D → X i*p∈W)
    (a p : MvPolynomial σ K) (hp : p∈W)
    (budget : a.totalDegree+p.totalDegree≤D) : a*p∈W := by
  classical
  have monomials : ∀ e : σ →₀ ℕ, e.degree+p.totalDegree≤D → monomial e 1*p∈W := by
    intro e
    induction e using Finsupp.induction with
    | zero => intro _; simpa using hp
    | single_add i n e _ _ ih =>
        intro bound
        have degrees : (Finsupp.single i n+e).degree=n+e.degree := by simp
        rw [degrees] at bound
        have base := ih (by omega)
        have powers : ∀ j : ℕ, j+e.degree+p.totalDegree≤D →
            X i^j*(monomial e 1*p)∈W := by
          intro j
          induction j with
          | zero => intro _; simpa using base
          | succ j hj =>
              intro hbound
              have before := hj (by omega)
              have degree : (X i^j*(monomial e 1*p)).totalDegree≤j+e.degree+p.totalDegree := by
                apply (totalDegree_mul _ _).trans
                have dp := totalDegree_mul (monomial e (1 : K)) p
                have de : (monomial e (1 : K)).totalDegree≤e.degree := totalDegree_monomial_le e (1 : K)
                have dx : (X i^j : MvPolynomial σ K).totalDegree=j := by simp
                simpa only [add_assoc] using
                  add_le_add dx.le (dp.trans (add_le_add de le_rfl))
              have derived := variable_stable i _ before (by omega)
              simpa [pow_succ,mul_assoc,mul_left_comm] using derived
        have derived := powers n bound
        simpa only [X_pow_eq_monomial,←mul_assoc,monomial_mul,mul_one] using derived
  rw [a.as_sum,Finset.sum_mul]
  apply W.sum_mem
  intro e he
  have edeg : e.degree≤a.totalDegree := le_totalDegree he
  have hm := monomials e (by omega)
  simpa only [smul_eq_C_mul,←mul_assoc,C_mul_monomial,mul_one] using W.smul_mem (coeff e a) hm

theorem closure_minimal_by_variables
    (F : Set (MvPolynomial σ K)) (D : ℕ)
    (W : Submodule K (MvPolynomial σ K))
    (inputs : ∀ f∈F, f.totalDegree≤D → f∈W)
    (variable_stable : ∀ i p, p∈W → p.totalDegree<D → X i*p∈W) :
    closureSpace F D≤W := by
  apply closure_minimal F D W inputs
  intro a p hp budget
  exact variable_closed_mul_mem W D variable_stable a p hp budget

end ICFallDegrees
