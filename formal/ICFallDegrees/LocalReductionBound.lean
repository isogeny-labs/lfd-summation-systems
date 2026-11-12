import ICFallDegrees.PointIdeal
import ICFallDegrees.NormalFormBound
import Mathlib.Data.Finsupp.Weight

/-! A finite low-degree reduction criterion for full-ideal LFD bounds.
This supports proper-ideal certificates: a bounded algebra section plus
checked local reductions suffices; enumeration of later falls is unnecessary. -/
noncomputable section
namespace ICFallDegrees.LocalReductionBound
open MvPolynomial
variable {σ K A : Type*} [Field K] [CommRing A] [Algebra K A]

private theorem variable_step
    (F : Set (MvPolynomial σ K)) (q : MvPolynomial σ K →ₐ[K] A)
    (s : A →ₗ[K] MvPolynomial σ K) (section_id : ∀ a, q (s a)=a)
    (sigma D : ℕ) (room : sigma+1≤D) (bounded : ∀ a, (s a).totalDegree≤sigma)
    (local_reductions : ∀ p : MvPolynomial σ K, p.totalDegree≤D → InClosure F D (p-s (q p)))
    (i : σ) (p : MvPolynomial σ K) (e : ℕ) (degree : p.totalDegree≤e)
    (previous : InClosure F (max D e) (p-s (q p))) :
    InClosure F (max D (e+1)) (X i*p-s (q (X i*p))) := by
  have low : (X i*s (q p)).totalDegree≤D := by
    have h := totalDegree_mul (X i) (s (q p))
    have hx : (X i : MvPolynomial σ K).totalDegree≤1 := by simp
    have hs := bounded (q p)
    omega
  have correction := (local_reductions _ low).mono (Nat.le_max_left D (e+1))
  have difference : (p-s (q p)).totalDegree≤max e sigma :=
    (totalDegree_sub_le _ _).trans (max_le_max degree (bounded _))
  have lifted := previous.mono (max_le_max le_rfl (by omega : e≤e+1))
  have multiplied := lifted.mul (X i) (by
    have hx : (X i : MvPolynomial σ K).totalDegree≤1 := by simp
    have hD := Nat.le_max_left D (e+1)
    have he := Nat.le_max_right D (e+1)
    omega)
  have same : q (X i*s (q p))=q (X i*p) := by rw [map_mul,section_id,map_mul]
  rw [same] at correction
  have identity : X i*p-s (q (X i*p))=
      X i*(p-s (q p))+(X i*s (q p)-s (q (X i*p))) := by ring
  rw [identity]
  exact multiplied.add correction

theorem global_reduction_of_local
    (F : Set (MvPolynomial σ K)) (q : MvPolynomial σ K →ₐ[K] A)
    (s : A →ₗ[K] MvPolynomial σ K) (section_id : ∀ a, q (s a)=a)
    (sigma D : ℕ) (room : sigma+1≤D) (bounded : ∀ a, (s a).totalDegree≤sigma)
    (local_reductions : ∀ p : MvPolynomial σ K, p.totalDegree≤D → InClosure F D (p-s (q p)))
    (p : MvPolynomial σ K) : InClosure F (max D p.totalDegree) (p-s (q p)) := by
  classical
  let N : MvPolynomial σ K →ₗ[K] MvPolynomial σ K := s.comp q.toLinearMap
  have monomials : ∀ e : σ →₀ ℕ, InClosure F (max D e.degree)
      (monomial e 1-N (monomial e 1)) := by
    intro e
    induction e using Finsupp.induction with
    | zero =>
        simpa [N] using local_reductions 1 (by simp)
    | single_add i n e _ _ ih =>
        have powers : ∀ j : ℕ, InClosure F (max D (j+e.degree))
            (X i^j*monomial e 1-N (X i^j*monomial e 1)) := by
          intro j
          induction j with
          | zero => simpa using ih
          | succ j hj =>
              have degree : (X i^j*(monomial e (1 : K))).totalDegree≤j+e.degree := by
                exact (totalDegree_mul _ _).trans
                  (Nat.add_le_add (by simp) (totalDegree_monomial_le _ _))
              have step := variable_step F q s section_id sigma D room bounded local_reductions
                i (X i^j*monomial e 1) (j+e.degree) degree hj
              simpa only [N,LinearMap.comp_apply,AlgHom.toLinearMap_apply,pow_succ,
                mul_assoc,mul_left_comm,Nat.add_assoc,Nat.add_left_comm,Nat.add_comm] using step
        simpa only [N,X_pow_eq_monomial,←mul_assoc,monomial_mul,mul_one,
          Finsupp.degree_add,Finsupp.degree_single] using powers n
  have summed : InClosure F (max D p.totalDegree)
      (∑ e∈p.support, (monomial e (coeff e p)-N (monomial e (coeff e p)))) := by
    apply (closureSpace F (max D p.totalDegree)).sum_mem
    intro e he
    have scaled := ((monomials e).mono (max_le_max le_rfl (le_totalDegree he))).smul (coeff e p)
    simpa only [smul_sub,←map_smul,smul_monomial,smul_eq_mul,mul_one] using scaled
  have identity : (∑ e∈p.support, (monomial e (coeff e p)-N (monomial e (coeff e p))))=p-N p := by
    rw [Finset.sum_sub_distrib,←map_sum,←p.as_sum]
  rw [identity] at summed
  exact summed

/-- Once the local slice is certified and representatives fit below D,
this is a full-ideal bound, including every proper-ideal consequence. -/
theorem lfd_bound_of_local
    (F : Set (MvPolynomial σ K)) (q : MvPolynomial σ K →ₐ[K] A)
    (s : A →ₗ[K] MvPolynomial σ K) (section_id : ∀ a, q (s a)=a)
    (sigma D : ℕ) (room : sigma+1≤D) (bounded : ∀ a, (s a).totalDegree≤sigma)
    (vanish : ∀ f∈F, q f=0)
    (local_reductions : ∀ p : MvPolynomial σ K, p.totalDegree≤D → InClosure F D (p-s (q p))) :
    HasLFDBound F D := by
  intro p member
  have inclusion : Ideal.span F≤RingHom.ker q.toRingHom := Ideal.span_le.mpr vanish
  have zero : q p=0 := inclusion member
  simpa only [zero,map_zero,sub_zero] using
    global_reduction_of_local F q s section_id sigma D room bounded local_reductions p

end ICFallDegrees.LocalReductionBound
