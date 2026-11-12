import ICFallDegrees.GeneralQuadraticClassification
import Mathlib.Data.Fin.VecNotation

/-! Any subgroup-chain depth k≥2, with two summands, one normalized S3
constraint, and no additional equations. The statements are algebraic over
any field of characteristic other than two; admissibility of an actual
prime-field curve/subgroup only restricts their parameters. -/
noncomputable section
open Classical
namespace ICFallDegrees.PrimeM2General
open MvPolynomial HomogeneousFirstFall GeneralQuadraticTop GeneralQuadraticClassification
variable {K : Type*} [Field K]
abbrev P (k : ℕ) (K : Type*) [Field K] := MvPolynomial (Fin (2*k)) K

def active (k : ℕ) (depth : 2≤k) : Fin 4 → Fin (2*k) :=
  ![⟨0,by omega⟩,⟨1,by omega⟩,⟨k,by omega⟩,⟨k+1,by omega⟩]

theorem active_injective (k : ℕ) (depth : 2≤k) : Function.Injective (active k depth) := by
  intro i j equal
  have values := congrArg Fin.val equal
  apply Fin.ext
  fin_cases i <;> fin_cases j <;> norm_num [active] at values ⊢ <;> omega

def chainTail (k : ℕ) (i : Fin (2*k)) : P k K :=
  if endBlock : i.val+1=k ∨ i.val+1=2*k then -1
  else -X ⟨i.val+1,by have := i.isLt; omega⟩

def tail (k : ℕ) (depth : 2≤k) (a b t : K) : P k K :=
  C (t^2)*(X (active k depth 1)+X (active k depth 3))+
    C (-2*a*t-4*b)*(X (active k depth 0)+X (active k depth 2))+C (a^2-4*b*t)

def normalized (k : ℕ) (depth : 2≤k) (a b t : K) : P k K :=
  top (active k depth) (-2*(t^2+a)) (-2*t)+tail k depth a b t

def inputs (k : ℕ) (depth : 2≤k) (a b t : K) : Fin (2*k+1) → P k K :=
  Fin.lastCases (normalized k depth a b t) (fun i => X i^2+chainTail k i)

def rawS3 (k : ℕ) (depth : 2≤k) (a b t : K) : P k K :=
  let u := X (active k depth 0)
  let w := X (active k depth 2)
  (u-w)^2*C (t^2)-2*((u+w)*(u*w+C a)+C (2*b))*C t+
    (u*w-C a)^2-C (4*b)*(u+w)

theorem raw_normalization_identity (k : ℕ) (depth : 2≤k) (a b t : K) :
    rawS3 k depth a b t-normalized k depth a b t=
      (X (active k depth 2)-C t)^2*
        (X (active k depth 0)^2-X (active k depth 1))+
      (X (active k depth 1)-C (2*t)*X (active k depth 0)+C (t^2))*
        (X (active k depth 2)^2-X (active k depth 3)) := by
  simp only [rawS3,normalized,top,tail,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring

theorem chain_tail_small (k : ℕ) (i : Fin (2*k)) : (chainTail k i : P k K).totalDegree<2 := by
  unfold chainTail
  split_ifs
  · simp
  · simp [ICFallDegrees.totalDegree_neg_eq]

theorem tail_small (k : ℕ) (depth : 2≤k) (a b t : K) : (tail k depth a b t).totalDegree<2 := by
  have first := ((isHomogeneous_X K (active k depth 1)).add
    (isHomogeneous_X K (active k depth 3))).C_mul (t^2)
  have second := ((isHomogeneous_X K (active k depth 0)).add
    (isHomogeneous_X K (active k depth 2))).C_mul (-2*a*t-4*b)
  have bound : (tail k depth a b t).totalDegree≤1 := by
    unfold tail
    apply (totalDegree_add _ _).trans
    apply max_le
    · exact (totalDegree_add _ _).trans (max_le first.totalDegree_le second.totalDegree_le)
    · rw [totalDegree_C]; decide
  omega

theorem degree_and_top (k : ℕ) (depth : 2≤k) (a b t : K) (i : Fin (2*k+1)) :
    (inputs k depth a b t i).totalDegree=2 ∧
      homogeneousComponent (inputs k depth a b t i).totalDegree (inputs k depth a b t i)=
        q (active k depth) (-2*(t^2+a)) (-2*t) i := by
  have nonzero (j : Fin (2*k+1)) : q (active k depth) (-2*(t^2+a)) (-2*t) j≠0 := by
    intro zero
    have equation := coefficient_pivots (active k depth) (active_injective k depth) (-2*(t^2+a)) (-2*t) j j
    rw [zero,coeff_zero,if_pos rfl] at equation
    exact zero_ne_one equation
  refine Fin.lastCases ?_ (fun j => ?_) i
  · simpa only [inputs,Fin.lastCases_last,normalized,q,Fin.lastCases_last] using
      degree_and_top_of_add (homogeneous_top (active k depth) (-2*(t^2+a)) (-2*t))
        (by simpa only [q,Fin.lastCases_last] using nonzero (Fin.last (2*k))) (tail_small k depth a b t)
  · simpa only [inputs,Fin.lastCases_castSucc,q,Fin.lastCases_castSucc] using
      degree_and_top_of_add (isHomogeneous_X_pow (R:=K) j 2)
        (by simpa only [q,Fin.lastCases_castSucc] using nonzero j.castSucc) (chain_tail_small (K:=K) k j)

/-- Exact ordinary homogeneous FFD for every k≥2. This is an FFD statement,
not an equality with LFD, solving degree, or relation-collection cost. -/
theorem exact_ordinary_ffd (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a b t : K) :
    IsKoszulFirstFall (inputs k depth a b t) (if 3*t^2+a=0 then 3 else 4) := by
  classical
  have tops : topForms (inputs k depth a b t)=q (active k depth) (-2*(t^2+a)) (-2*t) := by
    funext i; exact (degree_and_top k depth a b t i).2
  have degrees : (fun i => (inputs k depth a b t i).totalDegree)=(fun _ => 2) := by
    funext i; exact (degree_and_top k depth a b t i).1
  have determinant : -2*(t^2+a)-(-2*t)^2= -2*(3*t^2+a) := by ring
  change HasFall (topForms (inputs k depth a b t)) _ _ ∧ _
  simp only [HasKoszulFall,tops,degrees]
  by_cases zero : 3*t^2+a=0
  · rw [if_pos zero]
    have coefficient : -2*(t^2+a)=(-2*t)^2 := by
      apply sub_eq_zero.mp
      rw [determinant,zero,mul_zero]
    rw [coefficient]
    exact first_fall_three _ (active_injective k depth) (-2*t)
  · rw [if_neg zero]
    apply first_fall_four _ (active_injective k depth)
    rw [determinant]
    exact mul_ne_zero (neg_ne_zero.mpr two_ne) zero

end ICFallDegrees.PrimeM2General
