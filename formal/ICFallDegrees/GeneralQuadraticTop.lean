import ICFallDegrees.LowDegreeCoefficients
import Mathlib.Algebra.BigOperators.Fin

noncomputable section
namespace ICFallDegrees.GeneralQuadraticTop
open MvPolynomial HomogeneousFirstFall PowerIdeals LowDegreeCoefficients
open scoped BigOperators
variable {r : ℕ} {K : Type*} [CommRing K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

def top (a : Fin 4 → Fin r) (α β : K) : P r K :=
  C α*(X (a 0)*X (a 2))+C β*(X (a 0)*X (a 3)+X (a 1)*X (a 2))+X (a 1)*X (a 3)
def q (a : Fin 4 → Fin r) (α β : K) : Fin (r+1) → P r K :=
  Fin.lastCases (top a α β) (fun i => X i^2)
def pivot (a : Fin 4 → Fin r) : Fin (r+1) → Fin r →₀ ℕ :=
  Fin.lastCases (pair (a 1) (a 3)) (fun i => Finsupp.single i 2)

theorem top_mul (a : Fin 4 → Fin r) (α β : K) (p : P r K) :
    p*top a α β=C α*(p*X (a 0)*X (a 2))+C β*(p*X (a 0)*X (a 3))+
      C β*(p*X (a 1)*X (a 2))+p*X (a 1)*X (a 3) := by unfold top; ring

theorem top_square_coefficient (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (i : Fin r) : coeff (Finsupp.single i 2) (top a α β)=0 := by
  simp only [top,coeff_add,coeff_C_mul]
  rw [coefficient_square_pair i (a 0) (a 2) (injective.ne (by decide)),
    coefficient_square_pair i (a 0) (a 3) (injective.ne (by decide)),
    coefficient_square_pair i (a 1) (a 2) (injective.ne (by decide)),
    coefficient_square_pair i (a 1) (a 3) (injective.ne (by decide))]
  ring

theorem top_cross_coefficient (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) : coeff (pair (a 1) (a 3)) (top a α β)=1 := by
  have h0 : pair (a 1) (a 3) (a 0)=0 := by simp [pair,Finsupp.single_apply,injective.eq_iff]
  have h2 : pair (a 1) (a 3) (a 2)=0 := by simp [pair,Finsupp.single_apply,injective.eq_iff]
  have first := coefficient_pair_absent_left (1 : P r K) _ (a 0) (a 2) h0
  have second := coefficient_pair_absent_left (1 : P r K) _ (a 0) (a 3) h0
  have third := coefficient_pair_absent_right (1 : P r K) _ (a 1) (a 2) h2
  have fourth := coefficient_pair_sum (1 : P r K) 0 (a 1) (a 3)
  simp only [one_mul,zero_add,coeff_zero_one] at first second third fourth
  simp only [top,coeff_add,coeff_C_mul,first,second,third,fourth]
  ring

theorem coefficient_pivots (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (i j : Fin (r+1)) : coeff (pivot a i) (q a α β j)=if i=j then 1 else 0 := by
  classical
  refine Fin.lastCases ?_ (fun i => ?_) i
  · refine Fin.lastCases ?_ (fun j => ?_) j
    · simpa [q,pivot] using top_cross_coefficient a injective α β
    · have member : (X j : P r K)^2∈powerIdeal 2 := Ideal.subset_span ⟨j,rfl⟩
      have zero := squarefree_coefficient_zero member (pair (a 1) (a 3))
        (pair_squarefree _ _ (injective.ne (by decide)))
      simpa [q,pivot,Ne.symm (Fin.castSucc_ne_last j)] using zero
  · refine Fin.lastCases ?_ (fun j => ?_) j
    · simpa [q,pivot] using top_square_coefficient a injective α β i
    · simp [q,pivot,coeff_single_X_pow,eq_comm]

theorem homogeneous_top (a : Fin 4 → Fin r) (α β : K) : (top a α β).IsHomogeneous 2 := by
  have hp (i j : Fin r) : (X i*X j : P r K).IsHomogeneous 2 :=
    (isHomogeneous_X K i).mul (isHomogeneous_X K j)
  exact (((hp (a 0) (a 2)).C_mul α).add
    (((hp (a 0) (a 3)).add (hp (a 1) (a 2))).C_mul β)).add (hp (a 1) (a 3))

theorem homogeneous_q (a : Fin 4 → Fin r) (α β : K) (i : Fin (r+1)) :
    (q a α β i).IsHomogeneous 2 := by
  refine Fin.lastCases ?_ (fun i => ?_) i
  · simpa only [q,Fin.lastCases_last] using homogeneous_top a α β
  · simpa only [q,Fin.lastCases_castSucc] using isHomogeneous_X_pow (R:=K) i 2

theorem other_in_squares (a : Fin 4 → Fin r) (α β : K) (i : Fin (r+1))
    (other : i≠Fin.last r) : q a α β i∈powerIdeal 2 := by
  revert other
  refine Fin.lastCases ?_ (fun i other => ?_) i
  · intro impossible
    exact False.elim (impossible rfl)
  · simpa only [q,Fin.lastCases_castSucc] using (Ideal.subset_span ⟨i,rfl⟩ : (X i : P r K)^2∈powerIdeal 2)

theorem coefficient_uvw (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (p : P r K) :
    coeff (triple (a 0) (a 1) (a 2)) (p*top a α β)=
      α*coeff (Finsupp.single (a 1) 1) p+β*coeff (Finsupp.single (a 0) 1) p := by
  have absent : triple (a 0) (a 1) (a 2) (a 3)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff]
  rw [top_mul]
  simp only [coeff_add,coeff_C_mul]
  rw [coefficient_triple_13,coefficient_triple_23,
    coefficient_pair_absent_right p _ (a 0) (a 3) absent,
    coefficient_pair_absent_right p _ (a 1) (a 3) absent]
  ring

theorem coefficient_uvz (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (p : P r K) :
    coeff (triple (a 0) (a 1) (a 3)) (p*top a α β)=
      β*coeff (Finsupp.single (a 1) 1) p+coeff (Finsupp.single (a 0) 1) p := by
  have absent : triple (a 0) (a 1) (a 3) (a 2)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff]
  rw [top_mul]
  simp only [coeff_add,coeff_C_mul]
  rw [coefficient_triple_13,coefficient_triple_23,
    coefficient_pair_absent_right p _ (a 0) (a 2) absent,
    coefficient_pair_absent_right p _ (a 1) (a 2) absent]
  ring

theorem coefficient_uwz (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (p : P r K) :
    coeff (triple (a 0) (a 2) (a 3)) (p*top a α β)=
      α*coeff (Finsupp.single (a 3) 1) p+β*coeff (Finsupp.single (a 2) 1) p := by
  have absent : triple (a 0) (a 2) (a 3) (a 1)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff]
  rw [top_mul]
  simp only [coeff_add,coeff_C_mul]
  rw [coefficient_triple_12,coefficient_triple_13,
    coefficient_pair_absent_left p _ (a 1) (a 2) absent,
    coefficient_pair_absent_left p _ (a 1) (a 3) absent]
  ring

theorem coefficient_vwz (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (p : P r K) :
    coeff (triple (a 1) (a 2) (a 3)) (p*top a α β)=
      β*coeff (Finsupp.single (a 3) 1) p+coeff (Finsupp.single (a 2) 1) p := by
  have absent : triple (a 1) (a 2) (a 3) (a 0)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff]
  rw [top_mul]
  simp only [coeff_add,coeff_C_mul]
  rw [coefficient_triple_12,coefficient_triple_13,
    coefficient_pair_absent_left p _ (a 0) (a 2) absent,
    coefficient_pair_absent_left p _ (a 0) (a 3) absent]
  ring

theorem coefficient_unused (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) (p : P r K) (i : Fin r) (unused : ∀ j, i≠a j) :
    coeff (triple (a 1) (a 3) i) (p*top a α β)=coeff (Finsupp.single i 1) p := by
  have absent0 : triple (a 1) (a 3) i (a 0)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff,unused 0]
  have absent2 : triple (a 1) (a 3) i (a 2)=0 := by
    simp [triple,pair,Finsupp.single_apply,injective.eq_iff,unused 2]
  rw [top_mul]
  simp only [coeff_add,coeff_C_mul]
  rw [coefficient_triple_12,
    coefficient_pair_absent_left p _ (a 0) (a 2) absent0,
    coefficient_pair_absent_left p _ (a 0) (a 3) absent0,
    coefficient_pair_absent_right p _ (a 1) (a 2) absent2]
  ring

end ICFallDegrees.GeneralQuadraticTop
