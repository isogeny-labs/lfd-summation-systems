import ICFallDegrees.BinaryTriple
import ICFallDegrees.AffineTransport
import ICFallDegrees.BinaryS3Bounds
import Mathlib.Algebra.CharP.Two
import Mathlib.RingTheory.MvPolynomial.Basic

/-! Explicit all-length specializations to the three-bit obstruction.
These prove the sharp lower bound4 for target1 from m3 onward and target0
from m4 onward. The separate written upper4 theorem is not assumed here. -/
noncomputable section
namespace ICFallDegrees.SplitChainLower
open MvPolynomial
abbrev Vars (r : ℕ) := Fin (r+3) ⊕ Fin (r+1)
abbrev P (r : ℕ) := MvPolynomial (Vars r) (ZMod 2)
abbrev Q := MvPolynomial (Fin 3) (ZMod 2)

def s3 {σ : Type*} (x y z : MvPolynomial σ (ZMod 2)) :=
  1+x*y+x*z+y*z+x*y*z

def factor (r : ℕ) (i : Fin (r+3)) : P r := X (Sum.inl i)
def auxiliary (r : ℕ) (i : Fin (r+1)) : P r := X (Sum.inr i)
def first (r : ℕ) := s3 (factor r 0) (factor r 1) (auxiliary r 0)
def middle (r : ℕ) (i : Fin r) :=
  s3 (auxiliary r ⟨i,by omega⟩) (factor r ⟨i+2,by omega⟩) (auxiliary r ⟨i+1,by omega⟩)
def last (r : ℕ) (t : ZMod 2) :=
  s3 (auxiliary r ⟨r,by omega⟩) (factor r ⟨r+2,by omega⟩) (C t)
def inputs (r : ℕ) (t : ZMod 2) : Set (P r) :=
  Set.range (fun v : Vars r => (X v)^2-X v) ∪
    {p | p=first r ∨ (∃ i, p=middle r i) ∨ p=last r t}
def witness (r : ℕ) : P r := factor r 0+factor r 1+auxiliary r 0

private theorem mid_one {σ : Type*} (p q : MvPolynomial σ (ZMod 2)) : s3 p 1 q=p+q+1 := by
  unfold s3
  ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero]
private theorem target_one {σ : Type*} (p q : MvPolynomial σ (ZMod 2)) : s3 p q 1=p+q+1 := by
  unfold s3
  ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero]
private theorem map_s3 {σ : Type*} (f : MvPolynomial σ (ZMod 2) →ₐ[ZMod 2] Q) (p q z) :
    f (s3 p q z)=s3 (f p) (f q) (f z) := by simp [s3]

private theorem bool_in_target (i : Fin 3) : InClosure BinaryTriple.inputs 3 ((X i : Q)^2-X i) := by
  have belongs : ((X i : Q)^2-X i)∈BinaryTriple.inputs := by
    left
    refine ⟨i,?_⟩
    simp [BinaryTriple.G,BinaryTriple.sets,GridQuotient.gridPolynomial]
    ring
  exact .input belongs ((totalDegree_sub_le _ _).trans (by simp))

private theorem affine_bool (i : Fin 3) (c : ZMod 2) :
    InClosure BinaryTriple.inputs 3 (((X i : Q)+C c)^2-(X i+C c)) := by
  have same : ((X i : Q)+C c)^2-(X i+C c)=X i^2-X i := by
    fin_cases c <;> simp only [map_zero,map_one,add_zero] <;>
      ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero] <;> ring
  rw [same]
  exact bool_in_target i

private theorem affine_degree (i : Fin 3) (c : ZMod 2) : ((X i : Q)+C c).totalDegree≤1 :=
  (totalDegree_add _ _).trans (by simp)

def valuesOne (r : ℕ) : Vars r → Q
  | Sum.inl i => if i.val=0 then X 0 else if i.val=1 then X 1 else
      if i.val=r+2 then 1+X 2+C (r : ZMod 2) else 1
  | Sum.inr i => X 2+C (i.val : ZMod 2)
def projectOne (r : ℕ) : P r →ₐ[ZMod 2] Q := aeval (valuesOne r)

private theorem projectOne_degree (r : ℕ) : ∀ v, ((projectOne r) (X v)).totalDegree≤1 := by
  intro v
  simp only [projectOne,aeval_X]
  cases v with
  | inl i =>
      dsimp [valuesOne]
      split_ifs
      · simp
      · simp
      · exact (totalDegree_add _ _).trans (max_le ((totalDegree_add _ _).trans (by simp))
          (by rw [totalDegree_C]; omega))
      · simp
  | inr i => exact affine_degree 2 _

private theorem projectOne_inputs (r : ℕ) :
    ∀ p∈inputs r 1, p.totalDegree≤3 → InClosure BinaryTriple.inputs 3 (projectOne r p) := by
  intro p hp _
  rcases hp with ⟨v,rfl⟩ | rfl | ⟨i,rfl⟩ | rfl
  · simp only [map_sub,map_pow,projectOne,aeval_X]
    cases v with
    | inl i =>
        dsimp [valuesOne]
        split_ifs
        · exact bool_in_target 0
        · exact bool_in_target 1
        · have same : 1+(X 2 : Q)+C (r : ZMod 2)=X 2+C (1+(r : ZMod 2)) := by simp only [map_add,map_one]; ring
          rw [same]
          exact affine_bool 2 _
        · simp only [one_pow,sub_self]; exact .zero
    | inr i => exact affine_bool 2 _
  · have same : projectOne r (first r)=BinaryTriple.h := by
      simp [first,map_s3,factor,auxiliary,projectOne,valuesOne,aeval_X,BinaryTriple.h_product,s3]
    rw [same]
    exact .input (Or.inr rfl) BinaryTriple.h_degree.le
  · have hi := i.isLt
    have n0 : i.val+2≠0 := by omega
    have n1 : i.val+2≠1 := by omega
    have nl : i.val+2≠r+2 := by omega
    have same : projectOne r (middle r i)=0 := by
      simp only [middle,map_s3,auxiliary,factor,projectOne,aeval_X,valuesOne,Fin.val_mk,
        n0,n1,nl,if_false]
      rw [mid_one]
      simp only [Nat.cast_add,Nat.cast_one,map_add,map_one]
      ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero]
    rw [same]; exact .zero
  · have nz : r+2≠0 := by omega
    have no : r+2≠1 := by omega
    have same : projectOne r (last r 1)=0 := by
      simp only [last,map_s3,auxiliary,factor,projectOne,aeval_X,valuesOne,Fin.val_mk,
        nz,no,if_false,if_pos rfl,map_one]
      rw [target_one]
      ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero]
    rw [same]; exact .zero

private theorem witness_degree (r : ℕ) : (witness r).totalDegree≤1 := by
  unfold witness factor auxiliary
  exact (totalDegree_add _ _).trans (max_le ((totalDegree_add _ _).trans (by simp)) (by simp))

private theorem witness_available (r : ℕ) (t : ZMod 2) : InClosure (inputs r t) 4 (witness r) := by
  let names : Fin 3 → P r := fun i => if i=0 then factor r 0 else if i=1 then factor r 1 else auxiliary r 0
  let lift : Q →ₐ[ZMod 2] P r := aeval names
  have linear : ∀ i, (lift (X i)).totalDegree≤1 := by
    intro i
    simp only [lift,aeval_X,names]
    split_ifs <;> simp [factor,auxiliary]
  have im0 : lift (X 0)=factor r 0 := by simp only [lift,aeval_X,names,if_pos rfl,ite_true]
  have im1 : lift (X 1)=factor r 1 := by
    simp only [lift,aeval_X,names,if_neg (by decide : (1 : Fin 3)≠0),if_pos rfl,ite_true]
  have im2 : lift (X 2)=auxiliary r 0 := by
    simp only [lift,aeval_X,names,if_neg (by decide : (2 : Fin 3)≠0),if_neg (by decide : (2 : Fin 3)≠1)]
  have base : InClosure BinaryTriple.inputs 4 BinaryTriple.parity := by
    have h := BinaryTriple.upper_bound BinaryTriple.parity BinaryTriple.parity_mem_ideal
    have degree : BinaryTriple.parity.totalDegree≤1 := by
      unfold BinaryTriple.parity
      exact (totalDegree_add _ _).trans (max_le ((totalDegree_add _ _).trans (by simp)) (by simp))
    exact h.mono (by omega)
  have mapped := AffineTransport.map_closure lift linear BinaryTriple.inputs (inputs r t) 4 ?_ base
  · simpa only [BinaryTriple.parity,map_add,im0,im1,im2,witness] using mapped
  · intro p hp _
    rcases hp with ⟨i,rfl⟩ | rfl
    · have same : lift (GridQuotient.gridPolynomial (BinaryTriple.sets i) i)=
          (lift (X i))^2-lift (X i) := by
        simp [GridQuotient.gridPolynomial,BinaryTriple.sets]
        ring
      rw [same]
      simp only [lift,aeval_X,names]
      split_ifs <;> exact .input (Or.inl ⟨_,rfl⟩) ((totalDegree_sub_le _ _).trans (by simp [factor,auxiliary]))
    · have same : lift BinaryTriple.h=first r := by
        rw [BinaryTriple.h_product]
        simp only [map_add,map_mul,map_one,im0,im1,im2,first,s3]
      rw [same]
      apply InClosure.input (Or.inr (Or.inl rfl))
      have h := AffineTransport.degree_le lift linear BinaryTriple.h
      rw [same] at h
      exact h.trans (BinaryTriple.h_degree.le.trans (by decide))

theorem target_one_lower (r : ℕ) : ¬HasLFDBound (inputs r 1) 3 := by
  intro bound
  have member := bound (witness r) (witness_available r 1).sound
  have degree : (witness r).totalDegree≤3 := (witness_degree r).trans (by decide)
  have small : InClosure (inputs r 1) 3 (witness r) := member.mono (max_le le_rfl degree)
  have projected := AffineTransport.map_closure (projectOne r) (projectOne_degree r)
    (inputs r 1) BinaryTriple.inputs 3 (projectOne_inputs r) small
  have same : projectOne r (witness r)=BinaryTriple.parity := by
    simp [witness,factor,auxiliary,projectOne,valuesOne,BinaryTriple.parity]
  rw [same] at projected
  exact BinaryTriple.parity_not_at_three projected



def valuesZero (s : ℕ) : Vars (s+1) → Q
  | Sum.inl i => if i.val=0 then X 0 else if i.val=1 then X 1 else
      if i.val=2 then 1+X 2 else if i.val=s+3 then 1 else 0
  | Sum.inr i => if i.val=0 then X 2 else 1
def projectZero (s : ℕ) : P (s+1) →ₐ[ZMod 2] Q := aeval (valuesZero s)

private theorem projectZero_degree (s : ℕ) : ∀ v, ((projectZero s) (X v)).totalDegree≤1 := by
  intro v
  simp only [projectZero,aeval_X]
  cases v with
  | inl i =>
      dsimp [valuesZero]
      split_ifs
      · simp
      · simp
      · exact (totalDegree_add _ _).trans (by simp)
      · simp
      · simp
  | inr i => dsimp [valuesZero]; split_ifs <;> simp

private theorem projectZero_inputs (s : ℕ) :
    ∀ p∈inputs (s+1) 0, p.totalDegree≤3 → InClosure BinaryTriple.inputs 3 (projectZero s p) := by
  intro p hp _
  rcases hp with ⟨v,rfl⟩ | rfl | ⟨i,rfl⟩ | rfl
  · simp only [map_sub,map_pow,projectZero,aeval_X]
    cases v with
    | inl i =>
        dsimp [valuesZero]
        split_ifs
        · exact bool_in_target 0
        · exact bool_in_target 1
        · simpa [map_one,add_comm] using affine_bool 2 1
        · simp only [one_pow,sub_self]; exact .zero
        · simp only [zero_pow (by decide : 2≠0),sub_self]; exact .zero
    | inr i =>
        dsimp [valuesZero]
        split_ifs
        · exact bool_in_target 2
        · simp only [one_pow,sub_self]; exact .zero
  · have same : projectZero s (first (s+1))=BinaryTriple.h := by
      simp [first,map_s3,factor,auxiliary,projectZero,valuesZero,aeval_X,BinaryTriple.h_product,s3]
    rw [same]
    exact .input (Or.inr rfl) BinaryTriple.h_degree.le
  · have hi := i.isLt
    by_cases iz : i.val=0
    · have same : projectZero s (middle (s+1) i)=0 := by
        simp only [middle,map_s3,auxiliary,factor,projectZero,aeval_X,valuesZero,Fin.val_mk,iz]
        norm_num
        rw [target_one]
        ring_nf <;> simp [CharTwo.two_eq_zero,CharTwo.neg_eq,CharTwo.add_self_eq_zero]
      rw [same]; exact .zero
    · have next : i.val+1≠0 := by omega
      have n0 : i.val+2≠0 := by omega
      have n1 : i.val+2≠1 := by omega
      have n2 : i.val+2≠2 := by omega
      have nl : i.val+2≠s+3 := by omega
      have last_ne : i.val≠s+1 := by omega
      have same : projectZero s (middle (s+1) i)=0 := by
        simp [middle,map_s3,auxiliary,factor,projectZero,aeval_X,valuesZero,iz,next,n0,n1,n2,nl,last_ne,s3,
          CharTwo.add_self_eq_zero]
      rw [same]; exact .zero
  · have nz : s+1≠0 := by omega
    have n0 : s+1+2≠0 := by omega
    have n1 : s+1+2≠1 := by omega
    have n2 : s+1+2≠2 := by omega
    have index : s+1+2=s+3 := by omega
    have same : projectZero s (last (s+1) 0)=0 := by
      simp [last,map_s3,auxiliary,factor,projectZero,aeval_X,valuesZero,nz,n0,n1,n2,index,s3,
        CharTwo.add_self_eq_zero]
    rw [same]; exact .zero

theorem target_zero_lower (s : ℕ) : ¬HasLFDBound (inputs (s+1) 0) 3 := by
  intro bound
  have member := bound (witness (s+1)) (witness_available (s+1) 0).sound
  have degree : (witness (s+1)).totalDegree≤3 := (witness_degree (s+1)).trans (by decide)
  have small : InClosure (inputs (s+1) 0) 3 (witness (s+1)) := member.mono (max_le le_rfl degree)
  have projected := AffineTransport.map_closure (projectZero s) (projectZero_degree s)
    (inputs (s+1) 0) BinaryTriple.inputs 3 (projectZero_inputs s) small
  have same : projectZero s (witness (s+1))=BinaryTriple.parity := by
    simp [witness,factor,auxiliary,projectZero,valuesZero,BinaryTriple.parity]
  rw [same] at projected
  exact BinaryTriple.parity_not_at_three projected

/-- Every chain with m=s+4 and either target bit needs at least budget4.
This is independent of the unformalized all-m upper bound. -/
theorem all_lengths_lower (s : ℕ) (t : ZMod 2) : ¬HasLFDBound (inputs (s+1) t) 3 := by
  fin_cases t
  · exact target_zero_lower s
  · exact target_one_lower (s+1)

theorem all_lengths_no_small_bound (s : ℕ) (t : ZMod 2) (D : ℕ) (small : D<4) :
    ¬HasLFDBound (inputs (s+1) t) D := by
  intro bound
  apply all_lengths_lower s t
  intro p hp
  exact (bound p hp).mono (max_le_max (by omega) le_rfl)

/-- The family uses the actual characteristic-two normalized S3 expression. -/
theorem normalized_s3_identity {q : ℕ} (x y z : MvPolynomial (Fin q) (ZMod 2)) :
    s3 x y z=BinaryS3Bounds.normalized x y z x y z 1 := by
  unfold s3 BinaryS3Bounds.normalized
  simp only [map_one]
  ring

end ICFallDegrees.SplitChainLower
