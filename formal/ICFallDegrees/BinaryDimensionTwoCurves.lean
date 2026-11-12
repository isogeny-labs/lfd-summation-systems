import ICFallDegrees.BinaryExtensionFields
import ICFallDegrees.BinaryDimensionTwoSystems
import ICFallDegrees.BinaryDescent
import Mathlib.Data.Fin.VecNotation

/-! Actual coefficient descents over the proved fields F8 and F16, with
explicit targets and non-Frobenius-invariant dimension-two coordinate bases. -/
noncomputable section
set_option maxHeartbeats 10000000
namespace ICFallDegrees.BinaryDimensionTwoCurves
open MvPolynomial BinaryExtensionFields BinaryDescent
open BinaryDimensionTwoSystems (fieldEquations ordinaryCoordinates koblitzCoordinates ordinaryInputs koblitzInputs)
open scoped BigOperators
abbrev P8 := MvPolynomial (Fin 4) K8
abbrev P16 := MvPolynomial (Fin 4) K16

def lift8 (p : BinaryDimensionTwoSystems.P) : P8 := map (algebraMap k K8) p
def lift16 (p : BinaryDimensionTwoSystems.P) : P16 := map (algebraMap k K16) p
def target8 : K8 := 1+alpha8^2
def target16 : K16 := 1+alpha16+alpha16^2

def lifted8 : P8 := lift8 (ordinaryCoordinates.get ⟨0,by decide⟩)+
  C alpha8*lift8 (ordinaryCoordinates.get ⟨1,by decide⟩)+
  (C alpha8)^2*lift8 (ordinaryCoordinates.get ⟨2,by decide⟩)
def lifted16 : P16 := lift16 (koblitzCoordinates.get ⟨0,by decide⟩)+
  C alpha16*lift16 (koblitzCoordinates.get ⟨1,by decide⟩)+
  (C alpha16)^2*lift16 (koblitzCoordinates.get ⟨2,by decide⟩)+
  (C alpha16)^3*lift16 (koblitzCoordinates.get ⟨3,by decide⟩)

/-- Explicit ring identity quotients; all identities below are checked by
ring normalization, including the integer multiple of two before CharP is used. -/
def correctionQ8 : P8 := let a : P8:=C alpha8;
  a^3*X 1+a^3*X 3+a*X 0+2*a*X 1*X 3+a*X 1+a*X 2+a*X 3+
    X 0*X 3+X 1*X 2-X 1-X 3
def correctionR8 : P8 := let a : P8:=C alpha8;
  -a^2*X 1*X 3-a*X 0-a*X 1*X 3-a*X 2+X 0*X 2-X 0*X 3-X 1*X 2

def correctionQ16 : P16 := let a : P16:=C alpha16;
  a^2*X 1+a^2*X 3+2*a*X 1+2*a*X 3+X 0+2*X 1*X 3+3*X 1+X 2+3*X 3
def correctionR16 : P16 := let a : P16:=C alpha16;
  a^3*X 0+a^3*X 2+a^2*X 0*X 3+a^2*X 0+a^2*X 1*X 2-a^2*X 1+
    a^2*X 2-a^2*X 3-a*X 1*X 3-3*a*X 1-3*a*X 3+
    X 0*X 2-X 1*X 3-2*X 1-2*X 3

theorem ordinary_normalized_coefficients : normalized alpha8 target8 alpha8=lifted8 := by
  have identity : normalized alpha8 target8 alpha8-lifted8=
      ((C alpha8)^3+C alpha8+1)*correctionQ8+2*correctionR8 := by
    simp [normalized,xbar,ybar,x,y,target8,lifted8,lift8,ordinaryCoordinates,
      correctionQ8,correctionR8,map_pow]
    ring
  have root : (C alpha8 : P8)^3+C alpha8+1=0 := by
    simpa only [map_add,map_pow,map_one,map_zero] using
      congrArg (C : K8 → P8) alpha8_relation
  rw [root,zero_mul,CharTwo.two_eq_zero (R:=P8),zero_mul,add_zero] at identity
  exact sub_eq_zero.mp identity

private theorem koblitz_last : koblitzCoordinates.get ⟨3,by decide⟩=
    (X 0*X 3+X 1*X 2+X 1*X 3+X 1+X 3 : BinaryDimensionTwoSystems.P) := rfl

theorem koblitz_normalized_coefficients : normalized alpha16 target16 1=lifted16 := by
  have identity : normalized alpha16 target16 1-lifted16=
      ((C alpha16)^4+C alpha16+1)*correctionQ16+2*correctionR16 := by
    unfold lifted16
    rw [koblitz_last]
    simp [normalized,xbar,ybar,x,y,target16,lift16,koblitzCoordinates,
      correctionQ16,correctionR16,map_pow]
    ring
  have root : (C alpha16 : P16)^4+C alpha16+1=0 := by
    simpa only [map_add,map_pow,map_one,map_zero] using
      congrArg (C : K16 → P16) alpha16_relation
  rw [root,zero_mul,CharTwo.two_eq_zero (R:=P16),zero_mul,add_zero] at identity
  exact sub_eq_zero.mp identity

theorem ordinary_expansion : normalized alpha8 target8 alpha8=
    ∑ i : Fin 3, C (basis8 i)*map (algebraMap k K8) (ordinaryCoordinates.get i) := by
  rw [ordinary_normalized_coefficients]
  simp [lifted8,lift8,Fin.sum_univ_succ,basis8_eq_pow,map_pow,add_assoc]

theorem koblitz_expansion : normalized alpha16 target16 1=
    ∑ i : Fin 4, C (basis16 i)*map (algebraMap k K16) (koblitzCoordinates.get i) := by
  rw [koblitz_normalized_coefficients]
  simp [lifted16,lift16,Fin.sum_univ_succ,basis16_eq_pow,map_pow,add_assoc]
  all_goals rfl

theorem ordinary_coefficient_descent (i : Fin 3) :
    component basis8 i (normalized alpha8 target8 alpha8)=ordinaryCoordinates.get i :=
  component_of_expansion basis8 _ _ ordinary_expansion i

theorem koblitz_coefficient_descent (i : Fin 4) :
    component basis16 i (normalized alpha16 target16 1)=koblitzCoordinates.get i :=
  component_of_expansion basis16 _ _ koblitz_expansion i

def ordinaryDescendedInputs : Set BinaryDimensionTwoSystems.P :=
  {p | p∈fieldEquations ∨ ∃ i : Fin 3, p=component basis8 i (normalized alpha8 target8 alpha8)}
def koblitzDescendedInputs : Set BinaryDimensionTwoSystems.P :=
  {p | p∈fieldEquations ∨ ∃ i : Fin 4, p=component basis16 i (normalized alpha16 target16 1)}

theorem ordinary_descended_input_semantics : ordinaryDescendedInputs=ordinaryInputs := by
  ext p
  simp only [ordinaryDescendedInputs,ordinaryInputs,Set.mem_setOf_eq,List.mem_append,
    ordinary_coefficient_descent]
  constructor
  · rintro (member | ⟨i,rfl⟩)
    · exact Or.inl member
    · exact Or.inr (List.get_mem ordinaryCoordinates i)
  · rintro (member | member)
    · exact Or.inl member
    · obtain ⟨i,equal⟩ := List.mem_iff_get.mp member
      exact Or.inr ⟨i,equal.symm⟩

theorem koblitz_descended_input_semantics : koblitzDescendedInputs=koblitzInputs := by
  ext p
  simp only [koblitzDescendedInputs,koblitzInputs,Set.mem_setOf_eq,List.mem_append,
    koblitz_coefficient_descent]
  constructor
  · rintro (member | ⟨i,rfl⟩)
    · exact Or.inl member
    · exact Or.inr (List.get_mem koblitzCoordinates i)
  · rintro (member | member)
    · exact Or.inl member
    · obtain ⟨i,equal⟩ := List.mem_iff_get.mp member
      exact Or.inr ⟨i,equal.symm⟩

theorem ordinary_exact_lfd_three : IsLastFallDegree ordinaryDescendedInputs 3 := by
  rw [ordinary_descended_input_semantics]
  exact BinaryDimensionTwoSystems.ordinary_exact_lfd_three

theorem koblitz_exact_lfd_three : IsLastFallDegree koblitzDescendedInputs 3 := by
  rw [koblitz_descended_input_semantics]
  exact BinaryDimensionTwoSystems.koblitz_exact_lfd_three

theorem ordinary_target_on_curve : (1 : K8)^2+target8*1=
    target8^3+0*target8^2+alpha8 := by
  have identity : (1 : K8)^2+target8*1-(target8^3+0*target8^2+alpha8)=
      (alpha8^3+alpha8+1)*(-alpha8^3-2*alpha8+1) := by unfold target8; ring
  apply sub_eq_zero.mp
  simpa only [alpha8_relation,zero_mul] using identity

theorem koblitz_target_on_curve : (0 : K16)^2+target16*0=
    target16^3+0*target16^2+1 := by
  have identity : target16^3+1=
      (alpha16^4+alpha16+1)*(alpha16^2+3*alpha16+6)+
        2*(3*alpha16^3+alpha16^2-3*alpha16-2) := by unfold target16; ring
  rw [alpha16_relation,zero_mul,CharTwo.two_eq_zero (R:=K16),zero_mul,add_zero] at identity
  simpa using identity.symm

theorem basis8_coord_power (i j : Fin 3) :
    basis8.coord i (alpha8^j.val)=if j=i then 1 else 0 := by
  rw [←basis8_eq_pow j]
  simp [Basis.coord_apply,Finsupp.single_apply,eq_comm]

theorem basis16_coord_power (i j : Fin 4) :
    basis16.coord i (alpha16^j.val)=if j=i then 1 else 0 := by
  rw [←basis16_eq_pow j]
  simp [Basis.coord_apply,Finsupp.single_apply,eq_comm]

theorem ordinary_curve_b_nonzero : alpha8≠0 := by
  intro zero
  have own : basis8.coord 1 alpha8=1 := by simpa using basis8_coord_power 1 1
  rw [zero,map_zero] at own
  exact zero_ne_one own

theorem koblitz_curve_b_nonzero : (1 : K16)≠0 := one_ne_zero

def factorBase8 : Submodule k K8 := Submodule.span k {1,alpha8}
def factorBase16 : Submodule k K16 := Submodule.span k {1,alpha16}

theorem ordinary_square_not_mem : alpha8^2∉factorBase8 := by
  have contained : factorBase8≤LinearMap.ker (basis8.coord 2) := by
    apply Submodule.span_le.mpr
    intro z hz
    rcases (by simpa using hz : z=1 ∨ z=alpha8) with rfl | rfl
    · simpa using basis8_coord_power 2 0
    · simpa using basis8_coord_power 2 1
  intro member
  have zero : basis8.coord 2 (alpha8^2)=0 := contained member
  have own : basis8.coord 2 (alpha8^2)=1 := by simpa using basis8_coord_power 2 2
  rw [own] at zero
  exact one_ne_zero zero

theorem koblitz_square_not_mem : alpha16^2∉factorBase16 := by
  have contained : factorBase16≤LinearMap.ker (basis16.coord 2) := by
    apply Submodule.span_le.mpr
    intro z hz
    rcases (by simpa using hz : z=1 ∨ z=alpha16) with rfl | rfl
    · simpa using basis16_coord_power 2 0
    · simpa using basis16_coord_power 2 1
  intro member
  have zero : basis16.coord 2 (alpha16^2)=0 := contained member
  have own : basis16.coord 2 (alpha16^2)=1 := by simpa using basis16_coord_power 2 2
  rw [own] at zero
  exact one_ne_zero zero

theorem ordinary_not_frobenius_invariant : ¬∀ a∈factorBase8, a^2∈factorBase8 := by
  intro invariant
  exact ordinary_square_not_mem (invariant alpha8 (Submodule.subset_span (by simp)))

theorem koblitz_not_frobenius_invariant : ¬∀ a∈factorBase16, a^2∈factorBase16 := by
  intro invariant
  exact koblitz_square_not_mem (invariant alpha16 (Submodule.subset_span (by simp)))

end ICFallDegrees.BinaryDimensionTwoCurves
