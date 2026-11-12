import ICFallDegrees.QuadraticFamily

noncomputable section

namespace ICFallDegrees.QuadraticFamily

open MvPolynomial

variable {K : Type*} [Field K]

private def e4 : Fin 2 →₀ ℕ := Finsupp.single 0 2+Finsupp.single 1 2

def raw (c : K) : P K := monomial e4 1+h c-1
def rawInputs (c : K) : Set (P K) := {cx,cy,raw c}

theorem raw_form (c : K) : raw c=(x^2*y^2 : P K)+h c-1 := by
  simp only [raw,x,y,X_pow_eq_monomial,monomial_mul,mul_one]
  rfl

private theorem e4_ne_zero : (0 : Fin 2 →₀ ℕ) ≠ e4 := by
  intro he
  have hc := congrArg (fun f => f 0) he
  simp [e4] at hc

private theorem xy_ne_e4 : Finsupp.single (0 : Fin 2) 1+Finsupp.single 1 1 ≠ e4 := by
  intro he
  have hc := congrArg (fun f => f 0) he
  simp [e4] at hc

private theorem coeff_raw (c : K) : coeff e4 (raw c)=1 := by
  classical
  simp [raw,h,x,y,X,monomial_mul,coeff_one,e4_ne_zero,xy_ne_e4]

theorem raw_degree (c : K) : (raw c).totalDegree=4 := by
  apply le_antisymm
  · apply (totalDegree_sub_le _ _).trans
    simp only [totalDegree_one,max_le_iff,Nat.zero_le,and_true]
    apply (totalDegree_add _ _).trans
    apply max_le
    · exact (totalDegree_monomial_le _ _).trans (by simp [e4,Finsupp.sum_add_index])
    · exact (degree_h c).trans (by omega)
  · have hs : e4 ∈ (raw c).support := by
      rw [mem_support_iff,coeff_raw]
      exact one_ne_zero
    have hd := le_totalDegree hs
    simpa [e4,Finsupp.sum_add_index] using hd

private theorem raw_cx_degree : (cx : P K).totalDegree ≤ 4 := degree_cx.trans (by omega)
private theorem raw_cy_degree : (cy : P K).totalDegree ≤ 4 := degree_cy.trans (by omega)

theorem normalized_at_four (c : K) : InClosure (rawInputs c) 4 (h c) := by
  have hr : InClosure (rawInputs c) 4 (raw c) := .input (by simp [rawInputs]) (by rw [raw_degree])
  have hx : InClosure (rawInputs c) 4 (cx : P K) := .input (by simp [rawInputs]) raw_cx_degree
  have hy : InClosure (rawInputs c) 4 (cy : P K) := .input (by simp [rawInputs]) raw_cy_degree
  have hysq : (y^2 : P K).totalDegree ≤ 2 := by
    exact (totalDegree_pow _ 2).trans (by simp [y])
  have product := InClosure.mul (y^2) hx (by have := degree_cx (K:=K); omega)
  convert (hr.sub product).sub hy using 1
  rw [raw_form]
  simp only [cx,cy]
  ring

theorem raw_one_at_four (c : K) (hc : 1-c*c ≠ 0) : InClosure (rawInputs c) 4 1 := by
  apply InClosure.transport (by omega : 3 ≤ 4) ?_ (one_at_three c hc)
  intro f hf _
  simp only [inputs,Set.mem_insert_iff,Set.mem_singleton_iff] at hf
  rcases hf with rfl | rfl | rfl
  · exact .input (by simp [rawInputs]) raw_cx_degree
  · exact .input (by simp [rawInputs]) raw_cy_degree
  · exact normalized_at_four c

theorem raw_no_one_at_three (c : K) : ¬ InClosure (rawInputs c) 3 1 := by
  intro hone
  have available : ∀ f ∈ rawInputs c, f.totalDegree ≤ 3 → eval (fun _ => (1 : K)) f=0 := by
    intro f hf hd
    simp only [rawInputs,Set.mem_insert_iff,Set.mem_singleton_iff] at hf
    rcases hf with rfl | rfl | rfl
    · simp [cx,x]
    · simp [cy,y]
    · rw [raw_degree] at hd
      omega
  have he := eval_zero_of_closure (rawInputs c) 3 (fun _ => (1 : K)) available hone
  simpa using he

theorem raw_exact_lfd_four (c : K) (hc : 1-c*c ≠ 0) :
    IsLastFallDegree (rawInputs c) 4 := by
  apply exactLFD_of_first_one (by omega) (raw_one_at_four c hc)
  simpa using raw_no_one_at_three c

theorem raw_in_normalized_at_four (c : K) : InClosure (inputs c) 4 (raw c) := by
  have hx : InClosure (inputs c) 4 (cx : P K) := .input (by simp [inputs]) raw_cx_degree
  have hy : InClosure (inputs c) 4 (cy : P K) := .input (by simp [inputs]) raw_cy_degree
  have hh : InClosure (inputs c) 4 (h c) := .input (by simp [inputs]) ((degree_h c).trans (by omega))
  have hysq : (y^2 : P K).totalDegree ≤ 2 := (totalDegree_pow _ 2).trans (by simp [y])
  have product := InClosure.mul (y^2) hx (by have := degree_cx (K:=K); omega)
  convert (product.add hy).add hh using 1
  rw [raw_form]
  simp only [cx,cy]
  ring

theorem same_ideal (c : K) : Ideal.span (inputs c)=Ideal.span (rawInputs c) := by
  apply le_antisymm
  · apply Ideal.span_le.mpr
    intro f hf
    simp only [inputs,Set.mem_insert_iff,Set.mem_singleton_iff] at hf
    rcases hf with rfl | rfl | rfl
    · exact Ideal.subset_span (by simp [rawInputs])
    · exact Ideal.subset_span (by simp [rawInputs])
    · exact (normalized_at_four c).sound
  · apply Ideal.span_le.mpr
    intro f hf
    simp only [rawInputs,Set.mem_insert_iff,Set.mem_singleton_iff] at hf
    rcases hf with rfl | rfl | rfl
    · exact Ideal.subset_span (by simp [inputs])
    · exact Ideal.subset_span (by simp [inputs])
    · exact (raw_in_normalized_at_four c).sound

end ICFallDegrees.QuadraticFamily
