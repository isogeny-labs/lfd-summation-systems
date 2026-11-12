import ICFallDegrees.GeneralQuadraticTop
import Mathlib.Tactic.LinearCombination

/-! Exact ordinary Koszul first fall degree, with arbitrary ambient variable
count and four distinct active variables. Lower-degree tails do not matter. -/
noncomputable section
set_option maxHeartbeats 2000000
namespace ICFallDegrees.GeneralQuadraticClassification
open MvPolynomial HomogeneousFirstFall PowerIdeals LowDegreeCoefficients GeneralQuadraticTop
open scoped BigOperators
variable {r : ℕ} {K : Type*} [CommRing K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

def linearWitness (a : Fin 4 → Fin r) (β : K) : Fin (r+1) → P r K :=
  Fin.lastCases (X (a 1)-C β*X (a 0))
    (fun i => if i=a 0 then C (β^2)*(X (a 3)+C β*X (a 2))
      else if i=a 1 then -(X (a 3)+C β*X (a 2)) else 0)

def quadraticWitness (a : Fin 4 → Fin r) (α β : K) : Fin (r+1) → P r K :=
  Fin.lastCases (X (a 0)*X (a 1))
    (fun i => if i=a 0 then -(C α*X (a 1)*X (a 2)+C β*X (a 1)*X (a 3))
      else if i=a 1 then -(C β*X (a 0)*X (a 2)+X (a 0)*X (a 3)) else 0)

private theorem sum_two (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (A B : P r K) :
    (∑ i : Fin r, (if i=a 0 then A else if i=a 1 then B else 0)*X i^2)=
      A*X (a 0)^2+B*X (a 1)^2 := by
  classical
  have entry (i : Fin r) : (if i=a 0 then A else if i=a 1 then B else 0)=
      (if i=a 0 then A else 0)+(if i=a 1 then B else 0) := by
    by_cases h0 : i=a 0
    · subst i; simp [injective.eq_iff]
    · simp [h0]
  simp_rw [entry,add_mul,ite_mul,zero_mul]
  simp only [Finset.sum_add_distrib,Finset.sum_ite_eq',Finset.mem_univ,if_true]

theorem linear_witness_syzygy (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (β : K) : IsHomogeneousSyzygy (q a (β^2) β) (fun _ => 2) 3 (linearWitness a β) := by
  classical
  have linear (i j : Fin r) : (X i+C β*X j : P r K).IsHomogeneous 1 :=
    (isHomogeneous_X K i).add ((isHomogeneous_X K j).C_mul β)
  constructor
  · intro i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · simp only [linearWitness,Fin.lastCases_last]
      exact Or.inr ⟨by decide,(isHomogeneous_X K _).sub ((isHomogeneous_X K _).C_mul β)⟩
    · simp only [linearWitness,Fin.lastCases_castSucc]
      split_ifs
      · exact Or.inr ⟨by decide,(linear _ _).C_mul (β^2)⟩
      · exact Or.inr ⟨by decide,(linear _ _).neg⟩
      · exact Or.inl rfl
  · rw [Fin.sum_univ_castSucc]
    simp only [linearWitness,q,Fin.lastCases_castSucc,Fin.lastCases_last]
    rw [sum_two a injective]
    simp only [top,map_pow]
    ring

theorem quadratic_witness_syzygy (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) : IsHomogeneousSyzygy (q a α β) (fun _ => 2) 4 (quadraticWitness a α β) := by
  classical
  have pair (i j : Fin r) : (X i*X j : P r K).IsHomogeneous 2 :=
    (isHomogeneous_X K i).mul (isHomogeneous_X K j)
  constructor
  · intro i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · simp only [quadraticWitness,Fin.lastCases_last]
      exact Or.inr ⟨by decide,pair _ _⟩
    · simp only [quadraticWitness,Fin.lastCases_castSucc]
      split_ifs
      · right
        refine ⟨by decide,?_⟩
        simpa only [mul_assoc] using (((pair (a 1) (a 2)).C_mul α).add ((pair (a 1) (a 3)).C_mul β)).neg
      · right
        refine ⟨by decide,?_⟩
        simpa only [mul_assoc] using (((pair (a 0) (a 2)).C_mul β).add (pair (a 0) (a 3))).neg
      · exact Or.inl rfl
  · rw [Fin.sum_univ_castSucc]
    simp only [quadraticWitness,q,Fin.lastCases_castSucc,Fin.lastCases_last]
    rw [sum_two a injective]
    simp only [top]
    ring

variable [Nontrivial K]

theorem linear_witness_not_koszul (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (β : K) : linearWitness a β∉koszulSpan (q a (β^2) β) := by
  intro member
  have inside := koszul_coordinate_mem (q a (β^2) β) (Fin.last r) (powerIdeal 2)
    (other_in_squares a (β^2) β) member
  have zero := squarefree_coefficient_zero inside (Finsupp.single (a 1) 1)
    (fun i => by simp only [Finsupp.single_apply]; split_ifs <;> omega)
  have unit : coeff (Finsupp.single (a 1) 1) (linearWitness a β (Fin.last r))=1 := by
    simp [linearWitness,coeff_single_X,injective.eq_iff]
  exact one_ne_zero (unit.symm.trans zero)

theorem quadratic_witness_not_koszul (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : K) : quadraticWitness a α β∉koszulSpan (q a α β) := by
  intro member
  have inside := koszul_coordinate_mem (q a α β) (Fin.last r) (powerIdeal 2)
    (other_in_squares a α β) member
  have zero := squarefree_coefficient_zero inside (pair (a 0) (a 1))
    (pair_squarefree _ _ (injective.ne (by decide)))
  have unit : coeff (pair (a 0) (a 1)) (quadraticWitness a α β (Fin.last r))=1 := by
    simp only [quadraticWitness,Fin.lastCases_last]
    rw [pair_polynomial,coeff_monomial,if_pos rfl]
  exact one_ne_zero (unit.symm.trans zero)

section Field
variable {F : Type*} [Field F]

/-- Away from the determinant-zero branch there is no homogeneous linear
syzygy at all, a stronger statement than nonexistence modulo Koszul. -/
theorem degree_three_syzygy_zero (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : F) (nondegenerate : α-β^2≠0) {h : Fin (r+1) → P r F}
    (relation : IsHomogeneousSyzygy (q a α β) (fun _ => 2) 3 h) : h=0 := by
  classical
  have homogeneous (i : Fin (r+1)) : (h i).IsHomogeneous 1 := by
    rcases relation.1 i with zero | degree
    · rw [zero]; exact isHomogeneous_zero _ _ 1
    · exact degree.2
  have sum_zero : (∑ i : Fin r, h i.castSucc*X i^2)+h (Fin.last r)*top a α β=0 := by
    simpa only [Fin.sum_univ_castSucc,q,Fin.lastCases_castSucc,Fin.lastCases_last] using relation.2
  have sum_member : (∑ i : Fin r, h i.castSucc*X i^2)∈powerIdeal 2 := by
    apply Ideal.sum_mem
    intro i _
    exact (powerIdeal 2).mul_mem_left _ (Ideal.subset_span ⟨i,rfl⟩)
  have product_member : h (Fin.last r)*top a α β∈powerIdeal 2 := by
    have identity := eq_neg_of_add_eq_zero_right sum_zero
    rw [identity]
    exact (powerIdeal 2 : Ideal (P r F)).neg_mem sum_member
  have vanish (i j k : Fin r) (hij : i≠j) (hik : i≠k) (hjk : j≠k) :=
    squarefree_coefficient_zero product_member (triple i j k) (triple_squarefree i j k hij hik hjk)
  have e1 := vanish (a 0) (a 1) (a 2) (injective.ne (by decide))
    (injective.ne (by decide)) (injective.ne (by decide))
  have e2 := vanish (a 0) (a 1) (a 3) (injective.ne (by decide))
    (injective.ne (by decide)) (injective.ne (by decide))
  have e3 := vanish (a 0) (a 2) (a 3) (injective.ne (by decide))
    (injective.ne (by decide)) (injective.ne (by decide))
  have e4 := vanish (a 1) (a 2) (a 3) (injective.ne (by decide))
    (injective.ne (by decide)) (injective.ne (by decide))
  rw [coefficient_uvw a injective] at e1
  rw [coefficient_uvz a injective] at e2
  rw [coefficient_uwz a injective] at e3
  rw [coefficient_vwz a injective] at e4
  have c1 : coeff (Finsupp.single (a 1) 1) (h (Fin.last r))=0 := by
    apply (mul_eq_zero.mp (show (α-β^2)*coeff (Finsupp.single (a 1) 1) (h (Fin.last r))=0 by
      linear_combination e1-β*e2)).resolve_left nondegenerate
  have c0 : coeff (Finsupp.single (a 0) 1) (h (Fin.last r))=0 := by simpa [c1] using e2
  have c3 : coeff (Finsupp.single (a 3) 1) (h (Fin.last r))=0 := by
    apply (mul_eq_zero.mp (show (α-β^2)*coeff (Finsupp.single (a 3) 1) (h (Fin.last r))=0 by
      linear_combination e3-β*e4)).resolve_left nondegenerate
  have c2 : coeff (Finsupp.single (a 2) 1) (h (Fin.last r))=0 := by simpa [c3] using e4
  have last_zero : h (Fin.last r)=0 := by
    apply homogeneous_linear_eq_zero (homogeneous _)
    intro i
    by_cases h0 : i=a 0
    · simpa only [h0] using c0
    by_cases h1 : i=a 1
    · simpa only [h1] using c1
    by_cases h2 : i=a 2
    · simpa only [h2] using c2
    by_cases h3 : i=a 3
    · simpa only [h3] using c3
    have unused : ∀ j : Fin 4, i≠a j := by
      intro j
      fin_cases j <;> assumption
    have extra := vanish (a 1) (a 3) i (injective.ne (by decide)) (Ne.symm h1) (Ne.symm h3)
    rwa [coefficient_unused a injective α β _ i unused] at extra
  have squares_zero : (∑ i : Fin r, h i.castSucc*X i^2)=0 := by simpa [last_zero] using sum_zero
  funext i
  refine Fin.lastCases last_zero (fun i => ?_) i
  apply homogeneous_linear_eq_zero (homogeneous _)
  intro j
  have equation := congrArg (coeff (Finsupp.single i 2+Finsupp.single j 1)) squares_zero
  simpa only [coefficient_square_sum,coeff_zero] using equation

end Field

/-- Every number of unused variables is allowed, including all remaining
variables of a subgroup chain of any depth k≥2. -/
theorem first_fall_three (a : Fin 4 → Fin r) (injective : Function.Injective a) (β : K) :
    HasFall (q a (β^2) β) (fun _ => 2) 3 ∧
      ∀ D<3, ¬HasFall (q a (β^2) β) (fun _ => 2) D := by
  refine ⟨⟨linearWitness a β,linear_witness_syzygy a injective β,
    linear_witness_not_koszul a injective β⟩,?_⟩
  intro D small
  exact no_low_fall _ (pivot a) (coefficient_pivots a injective (β^2) β) 2 D (by omega)

theorem first_fall_four {F : Type*} [Field F] (a : Fin 4 → Fin r) (injective : Function.Injective a)
    (α β : F) (nondegenerate : α-β^2≠0) :
    HasFall (q a α β) (fun _ => 2) 4 ∧ ∀ D<4, ¬HasFall (q a α β) (fun _ => 2) D := by
  refine ⟨⟨quadraticWitness a α β,quadratic_witness_syzygy a injective α β,
    quadratic_witness_not_koszul a injective α β⟩,?_⟩
  intro D small ⟨h,relation,nontrivial⟩
  by_cases low : D≤2
  · exact no_low_fall _ (pivot a) (coefficient_pivots a injective α β) 2 D low ⟨h,relation,nontrivial⟩
  · have degree : D=3 := by omega
    subst D
    have zero := degree_three_syzygy_zero a injective α β nondegenerate relation
    exact nontrivial (zero.symm ▸ (koszulSpan _).zero_mem)

end ICFallDegrees.GeneralQuadraticClassification
