import ICFallDegrees.LowerBounds
import ICFallDegrees.FirstFall
import ICFallDegrees.PointIdeal
import Mathlib.Tactic.Ring

noncomputable section

/-!
# Sharp Boolean two-coordinate systems

For any field and any coefficient family with a nonzero member, the Boolean
field equations together with coefficient multiples of xy-1 have exact LFD 3.
Replacing xy-1 by x²y²-1 gives exact LFD 4. Both ideals are proper point ideals.
The coefficient-family form allows actual descended coordinate equations.
-/

namespace ICFallDegrees.BooleanPair

open MvPolynomial

variable {K ι : Type*} [Field K]

abbrev P (K : Type*) [CommRing K] := MvPolynomial (Fin 2) K
def x : P K := X 0
def y : P K := X 1
def cx : P K := x^2-x
def cy : P K := y^2-y
def h : P K := x*y-1

def inputs (c : ι → K) : Set (P K) :=
  fun f => f=cx ∨ f=cy ∨ ∃ i, f=C (c i)*h

private def ex : Fin 2 →₀ ℕ := Finsupp.single 0 2
private def ey : Fin 2 →₀ ℕ := Finsupp.single 1 2
private def exy : Fin 2 →₀ ℕ := Finsupp.single 0 1+Finsupp.single 1 1

private def combination (a b d : K) : P K := C a*cx+C b*cy+C d*h

private def quadSpan : Submodule K (P K) where
  carrier := fun f => ∃ a b d, f=combination a b d
  zero_mem' := ⟨0,0,0,by simp [combination]⟩
  add_mem' := by
    rintro f g ⟨a,b,d,rfl⟩ ⟨a',b',d',rfl⟩
    refine ⟨a+a',b+b',d+d',?_⟩
    simp only [combination,map_add]
    ring
  smul_mem' := by
    rintro a f ⟨b,d,e,rfl⟩
    refine ⟨a*b,a*d,a*e,?_⟩
    simp only [combination,smul_eq_C_mul,map_mul]
    ring

private theorem top_coefficients (a b d : K) :
    coeff ex (combination a b d)=a ∧ coeff ey (combination a b d)=b ∧
      coeff exy (combination a b d)=d := by
  classical
  have hxy : ex≠ey := by
    intro eq
    have he := congrArg (fun f => f 0) eq
    simp [ex,ey] at he
  have hxz : ex≠exy := by
    intro eq
    have he := congrArg (fun f => f 1) eq
    simp [ex,exy] at he
  have hyz : ey≠exy := by
    intro eq
    have he := congrArg (fun f => f 0) eq
    simp [ey,exy] at he
  have degree_ex : ex.sum (fun _ n => n)=2 := by simp [ex]
  have degree_ey : ey.sum (fun _ n => n)=2 := by simp [ey]
  have degree_exy : exy.sum (fun _ n => n)=2 := by simp [exy,Finsupp.sum_add_index]
  have low_zero : ∀ f : P K, f.totalDegree≤1 →
      coeff ex f=0 ∧ coeff ey f=0 ∧ coeff exy f=0 := by
    intro f hf
    refine ⟨?_,?_,?_⟩
    all_goals apply coeff_eq_zero_of_totalDegree_lt
    · change f.totalDegree < ex.sum (fun _ n => n)
      omega
    · change f.totalDegree < ey.sum (fun _ n => n)
      omega
    · change f.totalDegree < exy.sum (fun _ n => n)
      omega
  have lx := low_zero (x : P K) (by simp [x])
  have ly := low_zero (y : P K) (by simp [y])
  have lone := low_zero (1 : P K) (by simp)
  have form : combination a b d =
      C a*(monomial ex 1-x)+C b*(monomial ey 1-y)+C d*(monomial exy 1-1) := by
    have fx : (cx : P K)=monomial ex 1-x := by rw [cx,x,X_pow_eq_monomial]; rfl
    have fy : (cy : P K)=monomial ey 1-y := by rw [cy,y,X_pow_eq_monomial]; rfl
    have fh : (h : P K)=monomial exy 1-1 := by simp [h,x,y,X,monomial_mul,exy]
    rw [combination,fx,fy,fh]
  rw [form]
  simp [coeff_C_mul,lx.1,lx.2.1,lx.2.2,ly.1,ly.2.1,ly.2.2,
    lone.1,lone.2.1,lone.2.2,hxy,hxz,hyz,hxy.symm,hxz.symm,hyz.symm]

private theorem no_low : ∀ f ∈ (quadSpan : Submodule K (P K)), f≠0 → 2≤f.totalDegree := by
  rintro f ⟨a,b,d,rfl⟩ nonzero
  by_contra small
  have hd : (combination a b d).totalDegree<2 := by omega
  rcases top_coefficients a b d with ⟨ha,hb,hc⟩
  have za : a=0 := by
    rw [←ha]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < ex.sum (fun _ n => n)
    simpa [ex] using hd
  have zb : b=0 := by
    rw [←hb]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < ey.sum (fun _ n => n)
    simpa [ey] using hd
  have zd : d=0 := by
    rw [←hc]
    apply coeff_eq_zero_of_totalDegree_lt
    change _ < exy.sum (fun _ n => n)
    simpa [exy,Finsupp.sum_add_index] using hd
  apply nonzero
  simp [combination,za,zb,zd]

theorem degree_x : (x : P K).totalDegree≤1 := by simp [x]
theorem degree_y : (y : P K).totalDegree≤1 := by simp [y]
theorem degree_cx : (cx : P K).totalDegree≤2 := by
  exact (totalDegree_sub_le _ _).trans (max_le
    ((totalDegree_pow _ _).trans (by simp [x])) (degree_x.trans (by omega)))
theorem degree_cy : (cy : P K).totalDegree≤2 := by
  exact (totalDegree_sub_le _ _).trans (max_le
    ((totalDegree_pow _ _).trans (by simp [y])) (degree_y.trans (by omega)))
theorem degree_h : (h : P K).totalDegree≤2 := by
  exact (totalDegree_sub_le _ _).trans (max_le
    ((totalDegree_mul _ _).trans (by have := degree_x (K:=K); have := degree_y (K:=K); omega))
    (by simp))
theorem degree_scaled_h (a : K) : (C a*h : P K).totalDegree≤2 := by
  exact (totalDegree_mul _ _).trans (by have := degree_h (K:=K); simp only [totalDegree_C]; omega)

theorem closure_two_no_low (c : ι → K) {f : P K}
    (hf : InClosure (inputs c) 2 f) (degree : f.totalDegree<2) : f=0 := by
  have containing : closureSpace (inputs c) 2 ≤ quadSpan := by
    apply closure_le_top_span (inputs c) 2 quadSpan ?_ no_low
    intro p hp _
    rcases hp with rfl | rfl | ⟨i,rfl⟩
    · exact ⟨1,0,0,by simp [combination]⟩
    · exact ⟨0,1,0,by simp [combination]⟩
    · exact ⟨0,0,c i,by simp [combination]⟩
  by_contra hz
  have := no_low f (containing hf) hz
  omega

theorem h_at_three (c : ι → K) (nonzero : ∃ i, c i≠0) : InClosure (inputs c) 3 (h : P K) := by
  obtain ⟨i,hi⟩ := nonzero
  have scaled : InClosure (inputs c) 3 (C (c i)*h) :=
    .input (Or.inr (Or.inr ⟨i,rfl⟩)) ((degree_scaled_h _).trans (by omega))
  convert scaled.smul (c i)⁻¹ using 1
  rw [smul_eq_C_mul,←mul_assoc,←map_mul,inv_mul_cancel₀ hi,map_one,one_mul]

theorem x_at_three (c : ι → K) (nonzero : ∃ i, c i≠0) :
    InClosure (inputs c) 3 (x-1 : P K) := by
  have hh := h_at_three c nonzero
  have hc : InClosure (inputs c) 3 (cx : P K) :=
    .input (Or.inl rfl) (degree_cx.trans (by omega))
  have h1 := InClosure.mul x hh (by have := degree_x (K:=K); have := degree_h (K:=K); omega)
  have h2 := InClosure.mul y hc (by have := degree_y (K:=K); have := degree_cx (K:=K); omega)
  convert (hh.sub h1).add h2 using 1
  simp only [h,cx]
  ring

theorem y_at_three (c : ι → K) (nonzero : ∃ i, c i≠0) :
    InClosure (inputs c) 3 (y-1 : P K) := by
  have hh := h_at_three c nonzero
  have hc : InClosure (inputs c) 3 (cy : P K) :=
    .input (Or.inr (Or.inl rfl)) (degree_cy.trans (by omega))
  have h1 := InClosure.mul y hh (by have := degree_y (K:=K); have := degree_h (K:=K); omega)
  have h2 := InClosure.mul x hc (by have := degree_x (K:=K); have := degree_cy (K:=K); omega)
  convert (hh.sub h1).add h2 using 1
  simp only [h,cy]
  ring

theorem inputs_vanish (c : ι → K) : ∀ f ∈ inputs c, eval (fun _ => 1) f=0 := by
  intro f hf
  rcases hf with rfl | rfl | ⟨i,rfl⟩ <;> simp [cx,cy,h,x,y]

theorem normalized_bound (c : ι → K) (nonzero : ∃ i, c i≠0) : HasLFDBound (inputs c) 3 := by
  apply point_ideal_lfd_bound (inputs c) (fun _ => 1) 3 (inputs_vanish c)
  intro i
  have hi : i=0 ∨ i=1 := by omega
  rcases hi with rfl | rfl
  · simpa [x] using x_at_three c nonzero
  · simpa [y] using y_at_three c nonzero

theorem degree_x_sub_one : (x-1 : P K).totalDegree≤1 := by
  exact (totalDegree_sub_le _ _).trans (by simp [x])

theorem x_sub_one_ne_zero : (x-1 : P K)≠0 := by
  intro hz
  have he := congrArg (eval (fun _ => (0 : K))) hz
  simp [x] at he

theorem exact_lfd_three (c : ι → K) (nonzero : ∃ i, c i≠0) : IsLastFallDegree (inputs c) 3 := by
  refine ⟨normalized_bound c nonzero,?_⟩
  intro e he bound
  have hx := bound (x-1) (x_at_three c nonzero).sound
  have hd := degree_x_sub_one (K:=K)
  have low := hx.mono (by omega : max e (x-1 : P K).totalDegree≤2)
  exact x_sub_one_ne_zero (closure_two_no_low c low (by omega))

theorem first_strict_fall_three (c : ι → K) (nonzero : ∃ i, c i≠0) :
    IsFirstStrictFall (inputs c) 3 := by
  have absent : ¬InClosure (inputs c) 2 (x-1 : P K) := by
    intro hx
    exact x_sub_one_ne_zero (closure_two_no_low c hx
      (lt_of_le_of_lt degree_x_sub_one (by omega)))
  refine ⟨⟨x-1,x_at_three c nonzero,
    lt_of_le_of_lt degree_x_sub_one (by omega),absent⟩,?_⟩
  intro e he
  rintro ⟨f,hf,degree,not_before⟩
  have zero := closure_two_no_low c (hf.mono (by omega : e≤2)) (by omega)
  subst f
  exact not_before InClosure.zero

end ICFallDegrees.BooleanPair
