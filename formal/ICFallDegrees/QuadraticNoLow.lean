import ICFallDegrees.HomogeneousFirstFall
import ICFallDegrees.VariableClosure

/-! Excluding linear syzygies of quadratic tops excludes all degree-three
refutations in the full iterated closure, not merely a static Macaulay span. -/
noncomputable section
namespace ICFallDegrees.QuadraticNoLow
open MvPolynomial HomogeneousFirstFall
open scoped BigOperators
variable {r n : ℕ} {K : Type*} [Field K]
abbrev P (r : ℕ) (K : Type*) [Field K] := MvPolynomial (Fin r) K

def combine (f : Fin n → P r K) (c : Fin n → K) (l : Fin n → P r K) : P r K :=
  ∑ i, (C (c i)+l i)*f i

def affineSpan (f : Fin n → P r K) : Submodule K (P r K) where
  carrier := {p | ∃ c l, (∀ i, (l i).IsHomogeneous 1) ∧ p=combine f c l}
  zero_mem' := by
    refine ⟨0,0,fun _ => isHomogeneous_zero _ _ 1,?_⟩
    simp [combine]
  add_mem' := by
    rintro p q ⟨c,l,hl,rfl⟩ ⟨d,m,hm,rfl⟩
    refine ⟨c+d,l+m,fun i => (hl i).add (hm i),?_⟩
    simp only [combine,←Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _
    simp only [Pi.add_apply,map_add]
    ring
  smul_mem' := by
    rintro a p ⟨c,l,hl,rfl⟩
    refine ⟨fun i => a*c i,fun i => C a*l i,fun i => (hl i).C_mul a,?_⟩
    simp only [combine,smul_eq_C_mul,Finset.mul_sum,map_mul]
    apply Finset.sum_congr rfl
    intro i _
    ring

private theorem component_three (q t l : P r K) (c : K) (hq : q.IsHomogeneous 2)
    (ht : t.totalDegree≤1) (hl : l.IsHomogeneous 1) :
    homogeneousComponent 3 ((C c+l)*(q+t))=l*q := by
  have remainder : (C c*(q+t)+l*t).totalDegree<3 := by
    have hq2 := hq.totalDegree_le
    have hf2 : (q+t).totalDegree≤2 := (totalDegree_add _ _).trans (max_le hq2 (by omega))
    have h1 := totalDegree_mul (C c) (q+t)
    rw [totalDegree_C,zero_add] at h1
    have h2 := (totalDegree_mul l t).trans (Nat.add_le_add hl.totalDegree_le ht)
    have h := totalDegree_add (C c*(q+t)) (l*t)
    omega
  have identity : (C c+l)*(q+t)=l*q+(C c*(q+t)+l*t) := by ring
  rw [identity,map_add,homogeneousComponent_of_mem (hl.mul hq),if_pos rfl,
    homogeneousComponent_eq_zero _ _ remainder,add_zero]

private theorem linear_parts_zero (q t : Fin n → P r K)
    (hq : ∀ i, (q i).IsHomogeneous 2) (ht : ∀ i, (t i).totalDegree≤1)
    (no_linear : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 3 h → h=0)
    (c : Fin n → K) (l : Fin n → P r K) (hl : ∀ i, (l i).IsHomogeneous 1)
    (small : (combine (fun i => q i+t i) c l).totalDegree<3) : l=0 := by
  apply no_linear
  refine ⟨fun i => Or.inr ⟨by change 2≤3; decide,hl i⟩,?_⟩
  have zero := homogeneousComponent_eq_zero 3 _ small
  rw [combine,map_sum] at zero
  simpa only [component_three _ _ _ _ (hq _) (ht _) (hl _)] using zero

private theorem constants_zero (q t : Fin n → P r K)
    (hq : ∀ i, (q i).IsHomogeneous 2) (ht : ∀ i, (t i).totalDegree≤1)
    (no_constant : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 2 h → h=0)
    (c : Fin n → K) (small : (∑ i, C (c i)*(q i+t i)).totalDegree<2) : c=0 := by
  have relation : IsHomogeneousSyzygy q (fun _ => 2) 2 (fun i => C (c i)) := by
    refine ⟨fun i => Or.inr ⟨by change 2≤2; decide,isHomogeneous_C _ _⟩,?_⟩
    have zero := homogeneousComponent_eq_zero 2 _ small
    rw [map_sum] at zero
    have each (i : Fin n) : homogeneousComponent 2 (C (c i)*(q i+t i))=C (c i)*q i := by
      rw [homogeneousComponent_C_mul,map_add,homogeneousComponent_of_mem (hq i),if_pos rfl,
        homogeneousComponent_eq_zero _ _ (by have := ht i; omega),add_zero]
    simpa only [each] using zero
  have all_zero := no_constant _ relation
  funext i
  have zero := congrFun all_zero i
  have := congrArg constantCoeff zero
  simpa using this

/-- A uniform lower-bound tool: independent quadratic tops with no linear
syzygy admit no nonzero degree≤1 consequence at closure budget three. -/
theorem no_low_at_three (q t : Fin n → P r K)
    (hq : ∀ i, (q i).IsHomogeneous 2) (ht : ∀ i, (t i).totalDegree≤1)
    (no_constant : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 2 h → h=0)
    (no_linear : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 3 h → h=0)
    {p : P r K} (member : InClosure (Set.range (fun i => q i+t i)) 3 p)
    (small : p.totalDegree<2) : p=0 := by
  classical
  let f : Fin n → P r K := fun i => q i+t i
  have containment : closureSpace (Set.range f) 3≤affineSpan f := by
    apply closure_minimal_by_variables
    · rintro p ⟨i,rfl⟩ _
      refine ⟨fun j => if j=i then 1 else 0,0,fun _ => isHomogeneous_zero _ _ 1,?_⟩
      simp [combine]
    · intro v p hp degree
      obtain ⟨c,l,hl,rfl⟩ := hp
      have zero := linear_parts_zero q t hq ht no_linear c l hl degree
      have identity : combine f c l=∑ i, C (c i)*f i := by simp only [combine,zero,Pi.zero_apply,add_zero]
      refine ⟨0,fun i => C (c i)*X v,fun i => (isHomogeneous_X K v).C_mul (c i),?_⟩
      rw [identity,combine,Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      simp only [Pi.zero_apply,map_zero,zero_add]
      ring
  obtain ⟨c,l,hl,equal⟩ := containment member
  have zero_l : l=0 := linear_parts_zero q t hq ht no_linear c l hl (by rw [←equal]; omega)
  have identity : p=∑ i, C (c i)*(q i+t i) := by
    simpa only [combine,zero_l,Pi.zero_apply,add_zero] using equal
  have zero_c := constants_zero q t hq ht no_constant c (by rw [←identity]; exact small)
  simp only [identity,zero_c,Pi.zero_apply,map_zero,zero_mul,Finset.sum_const_zero]

theorem no_one_at_three (q t : Fin n → P r K)
    (hq : ∀ i, (q i).IsHomogeneous 2) (ht : ∀ i, (t i).totalDegree≤1)
    (no_constant : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 2 h → h=0)
    (no_linear : ∀ h, IsHomogeneousSyzygy q (fun _ => 2) 3 h → h=0) :
    ¬InClosure (Set.range (fun i => q i+t i)) 3 1 := by
  intro member
  exact one_ne_zero (no_low_at_three q t hq ht no_constant no_linear member (by simp))

end ICFallDegrees.QuadraticNoLow
