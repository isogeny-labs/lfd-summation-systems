import ICFallDegrees.Closure
import Mathlib.RingTheory.Ideal.Maps
import Mathlib.Tactic.Ring

/-!
# Degree-bounded reduction to a point

Once every coordinate difference is derived at budget D, every polynomial
minus its value at the point is derived at max(D, totalDegree p). This supplies
a full-ideal upper bound for proper point ideals, not just linear consequences.
-/

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K] [Nontrivial K]

omit [Nontrivial K] in
private theorem point_reduction_mul
    (F : Set (MvPolynomial σ K)) (point : σ → K) (D e f : ℕ)
    (p q : MvPolynomial σ K)
    (dp : p.totalDegree ≤ e) (dq : q.totalDegree ≤ f)
    (hp : InClosure F (max D e) (p-C (eval point p)))
    (hq : InClosure F (max D f) (q-C (eval point q))) :
    InClosure F (max D (e+f)) (p*q-C (eval point (p*q))) := by
  have hp' := hp.mono (max_le_max le_rfl (by omega : e ≤ e+f))
  have hq' := hq.mono (max_le_max le_rfl (by omega : f ≤ e+f))
  have dr : (p-C (eval point p)).totalDegree ≤ e :=
    (totalDegree_sub_C_le _ _).trans dp
  have budget : q.totalDegree+(p-C (eval point p)).totalDegree ≤ max D (e+f) := by
    have := Nat.le_max_right D (e+f)
    omega
  have product := InClosure.mul q hp' budget
  convert product.add (hq'.smul (eval point p)) using 1
  simp only [map_mul,smul_eq_C_mul]
  ring

theorem point_reduction
    (F : Set (MvPolynomial σ K)) (point : σ → K) (D : ℕ)
    (linear : ∀ i, InClosure F D (X i-C (point i)))
    (p : MvPolynomial σ K) :
    InClosure F (max D p.totalDegree) (p-C (eval point p)) := by
  classical
  have powers : ∀ i n, InClosure F (max D n)
      ((X i : MvPolynomial σ K)^n-C (eval point (X i^n))) := by
    intro i n
    induction n with
    | zero => simp only [pow_zero,map_one,sub_self]; exact .zero
    | succ n ih =>
        have hi : InClosure F (max D 1) (X i-C (eval point (X i))) := by
          simpa using (linear i).mono (Nat.le_max_left D 1)
        simpa only [pow_succ] using point_reduction_mul F point D n 1
          (X i^n) (X i) (by simp) (by simp) ih hi
  have monomials : ∀ e : σ →₀ ℕ, InClosure F (max D (e.sum fun _ n => n))
      (monomial e 1-C (eval point (monomial e 1))) := by
    intro e
    induction e using Finsupp.induction with
    | zero => simp only [monomial_zero',map_one,sub_self]; exact .zero
    | single_add i n e _ _ ih =>
        have hm := point_reduction_mul F point D n (e.sum fun _ n => n)
          (X i^n) (monomial e 1) (by simp) (totalDegree_monomial_le _ _) (powers i n) ih
        simpa [X_pow_eq_monomial,monomial_mul,Finsupp.sum_add_index] using hm
  have sum_mem : ∀ s : Finset (σ →₀ ℕ),
      (∀ e ∈ s, e.sum (fun _ n => n) ≤ p.totalDegree) →
      InClosure F (max D p.totalDegree)
        (∑ e ∈ s, (monomial e (coeff e p)-C (eval point (monomial e (coeff e p))))) := by
    intro s bounds
    apply (closureSpace F (max D p.totalDegree)).sum_mem
    intro e he
    have hm := ((monomials e).mono (max_le_max le_rfl (bounds e he))).smul (coeff e p)
    convert hm using 1
    simp [smul_eq_C_mul,mul_sub,eval_monomial,C_mul_monomial,map_mul]
  have hs := sum_mem p.support (fun _ he => le_totalDegree he)
  have identity :
      (∑ e ∈ p.support, (monomial e (coeff e p)-C (eval point (monomial e (coeff e p)))))
        = p-C (eval point p) := by
    rw [Finset.sum_sub_distrib,← map_sum,← map_sum,← p.as_sum]
  rwa [identity] at hs

/-- A full-ideal bound: all input polynomials vanish at the point, and all
coordinate differences have explicit degree-budgeted derivations. -/
theorem point_ideal_lfd_bound
    (F : Set (MvPolynomial σ K)) (point : σ → K) (D : ℕ)
    (vanish : ∀ f ∈ F, eval point f=0)
    (linear : ∀ i, InClosure F D (X i-C (point i))) : HasLFDBound F D := by
  intro p hp
  have subset : Ideal.span F ≤ RingHom.ker (eval point) := by
    apply Ideal.span_le.mpr
    intro f hf
    exact vanish f hf
  have hz : eval point p=0 := subset hp
  simpa [hz] using point_reduction F point D linear p

theorem point_ideal_eq_ker
    (F : Set (MvPolynomial σ K)) (point : σ → K) (D : ℕ)
    (vanish : ∀ f ∈ F, eval point f=0)
    (linear : ∀ i, InClosure F D (X i-C (point i))) :
    Ideal.span F=RingHom.ker (eval point) := by
  apply le_antisymm
  · apply Ideal.span_le.mpr
    intro f hf
    exact vanish f hf
  · intro p hp
    have hz : eval point p=0 := hp
    have hr := (point_reduction F point D linear p).sound
    simpa [hz] using hr

end ICFallDegrees
