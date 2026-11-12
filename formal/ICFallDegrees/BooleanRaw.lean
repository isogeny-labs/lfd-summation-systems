import ICFallDegrees.BooleanPair

noncomputable section

namespace ICFallDegrees.BooleanPair

open MvPolynomial

variable {K ι : Type*} [Field K]

def raw : P K := x^2*y^2-1
def rawInputs (c : ι → K) : Set (P K) :=
  fun f => f=cx ∨ f=cy ∨ ∃ i, f=C (c i)*raw

private def e4 : Fin 2 →₀ ℕ := Finsupp.single 0 2+Finsupp.single 1 2

theorem raw_monomial : (raw : P K)=monomial e4 1-1 := by
  simp only [raw,x,y,X_pow_eq_monomial,monomial_mul,mul_one]
  rfl

theorem degree_scaled_raw_le (a : K) : (C a*raw : P K).totalDegree≤4 := by
  have dr : (raw : P K).totalDegree≤4 := by
    rw [raw_monomial]
    exact (totalDegree_sub_le _ _).trans (max_le
      ((totalDegree_monomial_le _ _).trans (by simp [e4,Finsupp.sum_add_index])) (by simp))
  exact (totalDegree_mul _ _).trans (by simp only [totalDegree_C]; omega)

theorem degree_scaled_raw (a : K) (ha : a≠0) : (C a*raw : P K).totalDegree=4 := by
  apply le_antisymm (degree_scaled_raw_le a)
  have exponent_nonzero : (0 : Fin 2 →₀ ℕ)≠e4 := by
    intro eq
    have he := congrArg (fun f => f 0) eq
    simp [e4] at he
  have hc : coeff e4 (C a*raw : P K)=a := by
    simp [raw_monomial,coeff_C_mul,coeff_one,exponent_nonzero]
  have hs : e4∈(C a*raw : P K).support := by rw [mem_support_iff,hc]; exact ha
  have hd := le_totalDegree hs
  simpa [e4,Finsupp.sum_add_index] using hd

theorem raw_h_at_four (c : ι → K) (nonzero : ∃ i, c i≠0) :
    InClosure (rawInputs c) 4 (h : P K) := by
  obtain ⟨i,hi⟩ := nonzero
  have scaled : InClosure (rawInputs c) 4 (C (c i)*raw) :=
    .input (Or.inr (Or.inr ⟨i,rfl⟩)) (degree_scaled_raw_le _)
  have hr : InClosure (rawInputs c) 4 (raw : P K) := by
    convert scaled.smul (c i)⁻¹ using 1
    rw [smul_eq_C_mul,←mul_assoc,←map_mul,inv_mul_cancel₀ hi,map_one,one_mul]
  have hx : InClosure (rawInputs c) 4 (cx : P K) :=
    .input (Or.inl rfl) (degree_cx.trans (by omega))
  have hy : InClosure (rawInputs c) 4 (cy : P K) :=
    .input (Or.inr (Or.inl rfl)) (degree_cy.trans (by omega))
  have ysq : (y^2 : P K).totalDegree≤2 := by simp [y]
  have h1 := InClosure.mul (y^2) hx (by have := degree_cx (K:=K); omega)
  have h2 := InClosure.mul x hy (by have := degree_x (K:=K); have := degree_cy (K:=K); omega)
  convert (hr.sub h1).sub h2 using 1
  simp only [raw,h,cx,cy]
  ring

theorem raw_linear_at_four (c : ι → K) (nonzero : ∃ i, c i≠0) :
    ∀ i : Fin 2, InClosure (rawInputs c) 4 (X i-C (1 : K)) := by
  have available : ∀ p ∈ inputs c, p.totalDegree≤3 → InClosure (rawInputs c) 4 p := by
    intro p hp _
    rcases hp with rfl | rfl | ⟨i,rfl⟩
    · exact .input (Or.inl rfl) (degree_cx.trans (by omega))
    · exact .input (Or.inr (Or.inl rfl)) (degree_cy.trans (by omega))
    · simpa [smul_eq_C_mul] using (raw_h_at_four c nonzero).smul (c i)
  intro i
  have hi : i=0 ∨ i=1 := by omega
  rcases hi with rfl | rfl
  · simpa [x] using InClosure.transport (by omega : 3≤4) available (x_at_three c nonzero)
  · simpa [y] using InClosure.transport (by omega : 3≤4) available (y_at_three c nonzero)

theorem raw_inputs_vanish (c : ι → K) : ∀ f ∈ rawInputs c, eval (fun _ => 1) f=0 := by
  intro f hf
  rcases hf with rfl | rfl | ⟨i,rfl⟩ <;> simp [cx,cy,raw,x,y]

theorem raw_bound (c : ι → K) (nonzero : ∃ i, c i≠0) : HasLFDBound (rawInputs c) 4 :=
  point_ideal_lfd_bound (rawInputs c) (fun _ => 1) 4
    (raw_inputs_vanish c) (raw_linear_at_four c nonzero)

theorem raw_x_not_at_three (c : ι → K) : ¬ InClosure (rawInputs c) 3 (x-1 : P K) := by
  intro hx
  have vanish : ∀ f ∈ rawInputs c, f.totalDegree≤3 → eval (fun _ => 0) f=0 := by
    intro f hf hd
    rcases hf with rfl | rfl | ⟨i,rfl⟩
    · simp [cx,x]
    · simp [cy,y]
    · by_cases hi : c i=0
      · simp [hi]
      · rw [degree_scaled_raw _ hi] at hd
        omega
  have he := eval_zero_of_closure (rawInputs c) 3 (fun _ => 0) vanish hx
  simp [x] at he

theorem raw_exact_lfd_four (c : ι → K) (nonzero : ∃ i, c i≠0) :
    IsLastFallDegree (rawInputs c) 4 := by
  refine ⟨raw_bound c nonzero,?_⟩
  intro e he bound
  have linear : InClosure (rawInputs c) 4 (x-1 : P K) := by
    simpa [x] using raw_linear_at_four c nonzero 0
  have hx := bound (x-1) linear.sound
  have hd := degree_x_sub_one (K:=K)
  exact raw_x_not_at_three c (hx.mono (by omega))

theorem same_ideal (c : ι → K) (nonzero : ∃ i, c i≠0) :
    Ideal.span (inputs c)=Ideal.span (rawInputs c) := by
  have normalized : Ideal.span (inputs c)=RingHom.ker (eval (fun _ => (1 : K))) := by
    apply point_ideal_eq_ker (inputs c) (fun _ => 1) 3 (inputs_vanish c)
    intro i
    have hi : i=0 ∨ i=1 := by omega
    rcases hi with rfl | rfl
    · simpa [x] using x_at_three c nonzero
    · simpa [y] using y_at_three c nonzero
  have raw_ideal := point_ideal_eq_ker (rawInputs c) (fun _ => 1) 4
    (raw_inputs_vanish c) (raw_linear_at_four c nonzero)
  exact normalized.trans raw_ideal.symm

theorem proper_ideal (c : ι → K) : (1 : P K)∉Ideal.span (inputs c) := by
  intro hp
  have subset : Ideal.span (inputs c)≤RingHom.ker (eval (fun _ => (1 : K))) := by
    apply Ideal.span_le.mpr
    exact inputs_vanish c
  have he : eval (fun _ => (1 : K)) (1 : P K)=0 := subset hp
  simp at he

end ICFallDegrees.BooleanPair
