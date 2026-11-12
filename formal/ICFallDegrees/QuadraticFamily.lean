import ICFallDegrees.LowerBounds
import ICFallDegrees.FirstFall
import Mathlib.Tactic.Ring

noncomputable section

/-!
# A sharp two-variable last-fall-degree family

For any field K and c^2 != 1, the input list
`{x^2-1, y^2-1, x*y-C c}` has exact last fall degree 3.
Its lower bound is a filtered-span argument, not an assumed rank from Python.
-/

namespace ICFallDegrees.QuadraticFamily

open MvPolynomial

variable {K : Type*} [Field K]

abbrev P (K : Type*) [CommRing K] := MvPolynomial (Fin 2) K

def x : P K := X 0
def y : P K := X 1
def cx : P K := x^2-1
def cy : P K := y^2-1
def h (c : K) : P K := x*y-C c
def inputs (c : K) : Set (P K) := {cx,cy,h c}

private def ex : Fin 2 →₀ ℕ := Finsupp.single 0 2
private def ey : Fin 2 →₀ ℕ := Finsupp.single 1 2
private def exy : Fin 2 →₀ ℕ := Finsupp.single 0 1 + Finsupp.single 1 1

private def combination (c a b d : K) : P K := C a*cx+C b*cy+C d*h c

private theorem exponent_ne : ex ≠ ey ∧ ex ≠ exy ∧ ey ≠ exy ∧
    (0 : Fin 2 →₀ ℕ) ≠ ex ∧ (0 : Fin 2 →₀ ℕ) ≠ ey ∧
    (0 : Fin 2 →₀ ℕ) ≠ exy := by
  refine ⟨?_,?_,?_,?_,?_,?_⟩
  all_goals intro eq
  · have he := congrArg (fun f => f 0) eq
    simp [ex,ey] at he
  · have he := congrArg (fun f => f 1) eq
    simp [ex,exy] at he
  · have he := congrArg (fun f => f 0) eq
    simp [ey,exy] at he
  · have he := congrArg (fun f => f 0) eq
    simp [ex] at he
  · have he := congrArg (fun f => f 1) eq
    simp [ey] at he
  · have he := congrArg (fun f => f 0) eq
    simp [exy] at he

private theorem cx_form : (cx : P K) = monomial ex 1-1 := by
  rw [cx,x,X_pow_eq_monomial]
  rfl

private theorem cy_form : (cy : P K) = monomial ey 1-1 := by
  rw [cy,y,X_pow_eq_monomial]
  rfl

private theorem h_form (c : K) : h c = monomial exy 1-C c := by
  simp [h,x,y,X,monomial_mul,exy]

private def quadSpan (c : K) : Submodule K (P K) where
  carrier := fun f => ∃ a b d, f = combination c a b d
  zero_mem' := ⟨0,0,0,by simp [combination]⟩
  add_mem' := by
    rintro f g ⟨a,b,d,rfl⟩ ⟨a',b',d',rfl⟩
    refine ⟨a+a',b+b',d+d',?_⟩
    simp only [combination, map_add]
    ring
  smul_mem' := by
    rintro a f ⟨b,c',d,rfl⟩
    refine ⟨a*b,a*c',a*d,?_⟩
    simp only [combination, smul_eq_C_mul, map_mul]
    ring

private theorem coeff_ex (c a b d : K) :
    coeff ex (combination c a b d) = a := by
  classical
  rcases exponent_ne with ⟨h₁,h₂,h₃,h₄,h₅,h₆⟩
  simp [combination,cx_form,cy_form,h_form,coeff_C_mul,coeff_one,h₁,h₂,h₃,
    h₁.symm,h₂.symm,h₃.symm,h₄,h₅,h₆]

private theorem coeff_ey (c a b d : K) :
    coeff ey (combination c a b d) = b := by
  classical
  rcases exponent_ne with ⟨h₁,h₂,h₃,h₄,h₅,h₆⟩
  simp [combination,cx_form,cy_form,h_form,coeff_C_mul,coeff_one,h₁,h₂,h₃,
    h₁.symm,h₂.symm,h₃.symm,h₄,h₅,h₆]

private theorem coeff_exy (c a b d : K) :
    coeff exy (combination c a b d) = d := by
  classical
  rcases exponent_ne with ⟨h₁,h₂,h₃,h₄,h₅,h₆⟩
  simp [combination,cx_form,cy_form,h_form,coeff_C_mul,coeff_one,h₁,h₂,h₃,
    h₁.symm,h₂.symm,h₃.symm,h₄,h₅,h₆]

private theorem no_low (c : K) :
    ∀ f ∈ quadSpan c, f ≠ 0 → 2 ≤ f.totalDegree := by
  rintro f ⟨a,b,d,rfl⟩ nonzero
  by_contra small
  have hd : (combination c a b d).totalDegree < 2 := by omega
  have ha : a = 0 := by
    rw [← coeff_ex c a b d]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < ex.sum (fun _ n => n)
    simpa [ex] using hd
  have hb : b = 0 := by
    rw [← coeff_ey c a b d]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < ey.sum (fun _ n => n)
    simpa [ey] using hd
  have hdd : d = 0 := by
    rw [← coeff_exy c a b d]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < exy.sum (fun _ n => n)
    simpa [exy,Finsupp.sum_add_index] using hd
  apply nonzero
  simp [combination,ha,hb,hdd]

private theorem closure_two_le_span (c : K) : closureSpace (inputs c) 2 ≤ quadSpan c := by
  have available : ∀ f ∈ inputs c, f.totalDegree ≤ 2 → f ∈ quadSpan c := by
    intro f hf _
    simp only [inputs,Set.mem_insert_iff,Set.mem_singleton_iff] at hf
    rcases hf with rfl | rfl | rfl
    · exact ⟨1,0,0,by simp [combination]⟩
    · exact ⟨0,1,0,by simp [combination]⟩
    · exact ⟨0,0,1,by simp [combination]⟩
  exact closure_le_top_span (inputs c) 2 (quadSpan c) available (no_low c)

theorem closure_two_no_low (c : K) {f : P K}
    (hf : InClosure (inputs c) 2 f) (degree : f.totalDegree < 2) : f=0 := by
  by_contra hz
  have hd := no_low c f (closure_two_le_span c hf) hz
  omega

theorem no_one_at_two (c : K) : ¬ InClosure (inputs c) 2 1 := by
  intro hone
  have hw := closure_two_le_span c hone
  have hd := no_low c 1 hw one_ne_zero
  simpa using hd

theorem degree_x : (x : P K).totalDegree ≤ 1 := by simp [x]
theorem degree_y : (y : P K).totalDegree ≤ 1 := by simp [y]
theorem degree_cx : (cx : P K).totalDegree ≤ 2 := by
  apply (totalDegree_sub_le _ _).trans
  simp only [totalDegree_one, max_le_iff, Nat.zero_le, and_true]
  exact (totalDegree_pow _ 2).trans (by simpa [x] using (le_refl 2))
theorem degree_cy : (cy : P K).totalDegree ≤ 2 := by
  apply (totalDegree_sub_le _ _).trans
  simp only [totalDegree_one, max_le_iff, Nat.zero_le, and_true]
  exact (totalDegree_pow _ 2).trans (by simpa [y] using (le_refl 2))
theorem degree_h (c : K) : (h c).totalDegree ≤ 2 := by
  apply (totalDegree_sub_le _ _).trans
  simp only [totalDegree_C, max_le_iff, Nat.zero_le, and_true]
  exact (totalDegree_mul _ _).trans (add_le_add degree_x degree_y)

theorem one_at_three (c : K) (hc : 1-c*c ≠ 0) :
    InClosure (inputs c) 3 1 := by
  have hx : InClosure (inputs c) 3 (cx : P K) :=
    .input (by simp [inputs]) (degree_cx.trans (by omega))
  have hy : InClosure (inputs c) 3 (cy : P K) :=
    .input (by simp [inputs]) (degree_cy.trans (by omega))
  have hh : InClosure (inputs c) 3 (h c) :=
    .input (by simp [inputs]) ((degree_h c).trans (by omega))
  have hA : InClosure (inputs c) 3 (x-C c*y) := by
    have h₁ := InClosure.mul y hh (by have := degree_y (K:=K); have := degree_h c; omega)
    have h₂ := InClosure.mul x hy (by have := degree_x (K:=K); have := degree_cy (K:=K); omega)
    convert h₁.sub h₂ using 1 <;> simp [h,cy] <;> ring
  have hB : InClosure (inputs c) 3 (y-C c*x) := by
    have h₁ := InClosure.mul x hh (by have := degree_x (K:=K); have := degree_h c; omega)
    have h₂ := InClosure.mul y hx (by have := degree_y (K:=K); have := degree_cx (K:=K); omega)
    convert h₁.sub h₂ using 1 <;> simp [h,cx] <;> ring
  have linear : InClosure (inputs c) 3 (C (1-c*c)*x) := by
    convert hA.add (hB.smul c) using 1 <;> simp [smul_eq_C_mul,map_sub,map_mul] <;> ring
  have xmem : InClosure (inputs c) 3 (x : P K) := by
    convert linear.smul (1-c*c)⁻¹ using 1
    rw [smul_eq_C_mul,← mul_assoc,← map_mul,inv_mul_cancel₀ hc,map_one,one_mul]
  have square := InClosure.mul x xmem (by have := degree_x (K:=K); omega)
  convert square.sub hx using 1 <;> simp [cx,pow_two]

theorem exact_lfd_three (c : K) (hc : 1-c*c ≠ 0) :
    IsLastFallDegree (inputs c) 3 := by
  apply exactLFD_of_first_one (by omega) (one_at_three c hc)
  simpa using no_one_at_two c

theorem first_strict_fall_three (c : K) (hc : 1-c*c ≠ 0) :
    IsFirstStrictFall (inputs c) 3 := by
  refine ⟨⟨1,one_at_three c hc,by simp,no_one_at_two c⟩,?_⟩
  intro e he
  rintro ⟨f,hf,degree,absent⟩
  have hbudget : e ≤ 2 := by omega
  have hz := closure_two_no_low c (hf.mono hbudget) (by omega)
  subst f
  exact absent InClosure.zero

end ICFallDegrees.QuadraticFamily
