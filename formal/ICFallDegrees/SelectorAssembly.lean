import ICFallDegrees.Closure
import Mathlib.Data.Finsupp.Weight
import Mathlib.Algebra.MvPolynomial.Degrees

/-!
# A projector criterion for full-ideal LFD bounds

`LocalReductionBound.lfd_bound_of_local` works with a quotient algebra and a section.
Here the same induction is run for a linear projector `N` on the polynomial ring itself:
if `N` has bounded degree, absorbs itself under variable products, kills the ideal, and the
corrections `1 - N 1` and `X i * N p - N (X i * N p)` are derivable at degree `D`, then the
last fall degree is at most `D`. This is the induction of Section 2.5 of
`docs/proofs/U1-upper-bounds.md`.
-/

noncomputable section
namespace ICFallDegrees.SelectorAssembly
open MvPolynomial
variable {σ K : Type*} [Field K]

private theorem variable_step
    (F : Set (MvPolynomial σ K)) (N : MvPolynomial σ K →ₗ[K] MvPolynomial σ K)
    (sigma D : ℕ) (room : sigma+1 ≤ D) (bounded : ∀ p, (N p).totalDegree ≤ sigma)
    (absorb : ∀ i p, N (X i * N p) = N (X i * p))
    (step : ∀ i p, InClosure F D (X i * N p - N (X i * N p)))
    (i : σ) (p : MvPolynomial σ K) (e : ℕ) (degree : p.totalDegree ≤ e)
    (previous : InClosure F (max D e) (p - N p)) :
    InClosure F (max D (e+1)) (X i * p - N (X i * p)) := by
  have correction := (step i p).mono (Nat.le_max_left D (e+1))
  have difference : (p - N p).totalDegree ≤ max e sigma :=
    (totalDegree_sub_le _ _).trans (max_le_max degree (bounded _))
  have lifted := previous.mono (max_le_max le_rfl (by omega : e ≤ e+1))
  have multiplied := lifted.mul (X i) (by
    have hx : (X i : MvPolynomial σ K).totalDegree ≤ 1 := by simp
    have hD := Nat.le_max_left D (e+1)
    have he := Nat.le_max_right D (e+1)
    omega)
  rw [absorb] at correction
  have identity : X i * p - N (X i * p) =
      X i * (p - N p) + (X i * N p - N (X i * p)) := by ring
  rw [identity]
  exact multiplied.add correction

theorem global_reduction
    (F : Set (MvPolynomial σ K)) (N : MvPolynomial σ K →ₗ[K] MvPolynomial σ K)
    (sigma D : ℕ) (room : sigma+1 ≤ D) (bounded : ∀ p, (N p).totalDegree ≤ sigma)
    (absorb : ∀ i p, N (X i * N p) = N (X i * p))
    (one : InClosure F D (1 - N 1))
    (step : ∀ i p, InClosure F D (X i * N p - N (X i * N p)))
    (p : MvPolynomial σ K) : InClosure F (max D p.totalDegree) (p - N p) := by
  classical
  have monomials : ∀ e : σ →₀ ℕ, InClosure F (max D e.degree)
      (monomial e 1 - N (monomial e 1)) := by
    intro e
    induction e using Finsupp.induction with
    | zero =>
        simpa using one
    | single_add i n e _ _ ih =>
        have powers : ∀ j : ℕ, InClosure F (max D (j+e.degree))
            (X i^j * monomial e 1 - N (X i^j * monomial e 1)) := by
          intro j
          induction j with
          | zero => simpa using ih
          | succ j hj =>
              have degree : (X i^j * (monomial e (1 : K))).totalDegree ≤ j+e.degree := by
                exact (totalDegree_mul _ _).trans
                  (Nat.add_le_add (by simp) (totalDegree_monomial_le _ _))
              have next := variable_step F N sigma D room bounded absorb step
                i (X i^j * monomial e 1) (j+e.degree) degree hj
              simpa only [pow_succ, mul_assoc, mul_left_comm, Nat.add_assoc, Nat.add_left_comm,
                Nat.add_comm] using next
        simpa only [X_pow_eq_monomial, ← mul_assoc, monomial_mul, mul_one,
          Finsupp.degree_add, Finsupp.degree_single] using powers n
  have summed : InClosure F (max D p.totalDegree)
      (∑ e ∈ p.support, (monomial e (coeff e p) - N (monomial e (coeff e p)))) := by
    apply (closureSpace F (max D p.totalDegree)).sum_mem
    intro e he
    have scaled := ((monomials e).mono (max_le_max le_rfl (le_totalDegree he))).smul (coeff e p)
    simpa only [smul_sub, ← map_smul, smul_monomial, smul_eq_mul, mul_one] using scaled
  have identity : (∑ e ∈ p.support, (monomial e (coeff e p) - N (monomial e (coeff e p)))) =
      p - N p := by
    rw [Finset.sum_sub_distrib, ← map_sum, ← p.as_sum]
  rw [identity] at summed
  exact summed

/-- A degree-bounded projector with derivable corrections gives a full-ideal bound. -/
theorem lfd_bound_of_projector
    (F : Set (MvPolynomial σ K)) (N : MvPolynomial σ K →ₗ[K] MvPolynomial σ K)
    (sigma D : ℕ) (room : sigma+1 ≤ D) (bounded : ∀ p, (N p).totalDegree ≤ sigma)
    (absorb : ∀ i p, N (X i * N p) = N (X i * p))
    (vanish : ∀ f ∈ Ideal.span F, N f = 0)
    (one : InClosure F D (1 - N 1))
    (step : ∀ i p, InClosure F D (X i * N p - N (X i * N p))) :
    HasLFDBound F D := by
  intro p member
  simpa only [vanish p member, sub_zero] using
    global_reduction F N sigma D room bounded absorb one step p

end ICFallDegrees.SelectorAssembly
