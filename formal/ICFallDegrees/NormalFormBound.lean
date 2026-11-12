import ICFallDegrees.Closure

/-!
# A bounded normal-form criterion for last fall degree

The hypotheses are explicit algebraic obligations for a quotient and a linear
repr. In particular, bounded-degree reduction into the input closure is a
hypothesis to be discharged for the concrete encoding; no Groebner basis or
low-degree preprocessing is supplied for free.
-/

namespace ICFallDegrees

open MvPolynomial

variable {σ K A : Type*} [CommRing K] [Nontrivial K] [CommRing A] [Algebra K A]

omit [Nontrivial K] in
private theorem image_mul_mem (q : MvPolynomial σ K →ₐ[K] A)
    (U : Submodule K A)
    (var_closed : ∀ i a, a ∈ U → q (X i) * a ∈ U)
    (p : MvPolynomial σ K) {a : A} (ha : a ∈ U) : q p * a ∈ U := by
  induction p using MvPolynomial.induction_on with
  | C c =>
      rw [← MvPolynomial.algebraMap_eq, q.commutes]
      simpa only [Algebra.smul_def] using U.smul_mem c ha
  | add p r hp hr =>
      simpa only [map_add, add_mul] using U.add_mem hp hr
  | mul_X p i hp =>
      have h := var_closed i (q p * a) hp
      simpa only [map_mul, mul_assoc, mul_left_comm] using h

/-- If every quotient element has a representative of degree at most s,
and ordinary degree-bounded reduction is derivable from the actual inputs,
then the LFD is at most D for D >= s+1 and D >= every input degree. -/
theorem bounded_normal_form_lfd
    (F : Set (MvPolynomial σ K))
    (q : MvPolynomial σ K →ₐ[K] A)
    (repr : A →ₗ[K] MvPolynomial σ K)
    (right_inverse : ∀ a, q (repr a) = a)
    (s D : ℕ)
    (representative_degree : ∀ a, (repr a).totalDegree ≤ s)
    (reduction : ∀ p, InClosure F p.totalDegree (p - repr (q p)))
    (step_budget : s + 1 ≤ D)
    (input_budget : ∀ f ∈ F, f.totalDegree ≤ D) : HasLFDBound F D := by
  let W := closureSpace F D
  let U : Submodule K A := W.comap repr
  have normal_mem {p : MvPolynomial σ K} (hp : InClosure F D p) :
      InClosure F D (repr (q p)) := by
    have hr := (reduction p).mono hp.degree_le
    simpa only [sub_sub_cancel] using hp.sub hr
  have variable_stable : ∀ i a, a ∈ U → q (X i) * a ∈ U := by
    intro i a ha
    change InClosure F D (repr (q (X i) * a))
    change InClosure F D (repr a) at ha
    have hbudget : (X i : MvPolynomial σ K).totalDegree +
        (repr a).totalDegree ≤ D := by
      calc
        _ ≤ 1+s := add_le_add (by simp) (representative_degree a)
        _ ≤ D := by omega
    have hp := InClosure.mul (X i) ha hbudget
    have hn := normal_mem hp
    simpa only [map_mul, right_inverse] using hn
  let I : Ideal A :=
    { carrier := U
      zero_mem' := U.zero_mem
      add_mem' := U.add_mem
      smul_mem' := by
        intro a b hb
        change a * b ∈ U
        rw [← right_inverse a]
        exact image_mul_mem q U variable_stable (repr a) hb }
  have generators : ∀ f ∈ F, q f ∈ I := by
    intro f hf
    exact normal_mem (.input hf (input_budget f hf))
  intro f hf
  have image_mem : q f ∈ I := by
    induction hf using Submodule.span_induction with
    | mem p hp => exact generators p hp
    | zero => simpa only [map_zero] using I.zero_mem
    | add p r _ _ hp hr => simpa only [map_add] using I.add_mem hp hr
    | smul a p _ hp =>
        simpa only [smul_eq_mul, map_mul] using I.mul_mem_left (q a) hp
  have hn : InClosure F D (repr (q f)) := image_mem
  have hred := (reduction f).mono (Nat.le_max_right D f.totalDegree)
  have hn' := hn.mono (Nat.le_max_left D f.totalDegree)
  simpa only [sub_add_cancel] using hred.add hn'

end ICFallDegrees
