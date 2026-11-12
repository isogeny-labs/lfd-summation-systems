import ICFallDegrees.GeneratedBinaryOrdinary
import ICFallDegrees.GeneratedBinaryKoblitz
import Mathlib.RingTheory.MvPolynomial.Basic
import Mathlib.Algebra.CharP.Two
import Mathlib.Tactic.Ring

noncomputable section
set_option maxHeartbeats 10000000

namespace ICFallDegrees.BinaryDimensionTwoSystems
open MvPolynomial SparsePolynomial
abbrev k := ZMod 2
abbrev P := MvPolynomial (Fin 4) k

def fieldEquations : List P := [X 0^2+X 0,X 1^2+X 1,X 2^2+X 2,X 3^2+X 3]
def ordinaryCoordinates : List P := [
  X 0*X 3+X 0+X 1*X 2+X 1+X 2+X 3,
  X 0+X 2+1,
  X 0*X 2+X 0*X 3+X 0+X 1*X 2+X 1*X 3+X 2]
def koblitzCoordinates : List P := [
  X 1+X 3+1,
  X 0*X 2+X 0*X 3+X 0+X 1*X 2+X 1+X 2+X 3,
  X 0*X 2+X 0+X 1*X 3+X 2,
  X 0*X 3+X 1*X 2+X 1*X 3+X 1+X 3]
def ordinaryInputs : Set P := {p | p∈fieldEquations++ordinaryCoordinates}
def koblitzInputs : Set P := {p | p∈fieldEquations++koblitzCoordinates}

theorem field_equations_boolean : fieldEquations=
    [X 0^2-X 0,X 1^2-X 1,X 2^2-X 2,X 3^2-X 3] := by
  simp [fieldEquations,sub_eq_add_neg,CharTwo.neg_eq]

theorem toExponent_exponent (a b c d : ℕ) :
    toExponent (GeneratedBinaryOrdinary.exponent a b c d)=
      Finsupp.single 0 a+Finsupp.single 1 b+Finsupp.single 2 c+Finsupp.single 3 d := by
  ext i
  simp only [toExponent,Finsupp.onFinset_apply,Finsupp.add_apply,Finsupp.single_apply]
  fin_cases i <;> simp [GeneratedBinaryOrdinary.exponent]
  all_goals rfl

theorem ordinary_input_semantics : GeneratedBinaryOrdinary.inputs.map meaning=
    fieldEquations++ordinaryCoordinates := by
  simp only [GeneratedBinaryOrdinary.inputs,List.map_cons,List.map_nil,meaning,toExponent_exponent]
  simp [fieldEquations,ordinaryCoordinates,X_pow_eq_monomial,X,monomial_mul,monomial_pow,
    C_mul_monomial,add_assoc,add_comm,add_left_comm]

theorem koblitz_input_semantics : GeneratedBinaryKoblitz.inputs.map meaning=
    fieldEquations++koblitzCoordinates := by
  have exponent_equal : GeneratedBinaryKoblitz.exponent=GeneratedBinaryOrdinary.exponent := rfl
  simp only [GeneratedBinaryKoblitz.inputs,List.map_cons,List.map_nil,meaning,
    exponent_equal,toExponent_exponent]
  simp [fieldEquations,koblitzCoordinates,X_pow_eq_monomial,X,monomial_mul,monomial_pow,
    C_mul_monomial,add_assoc,add_comm,add_left_comm]

theorem ordinary_exact_lfd_three : IsLastFallDegree ordinaryInputs 3 := by
  have checked := GeneratedBinaryOrdinary.exact_lfd_three
  rw [ordinary_input_semantics] at checked
  exact checked

theorem koblitz_exact_lfd_three : IsLastFallDegree koblitzInputs 3 := by
  have checked := GeneratedBinaryKoblitz.exact_lfd_three
  rw [koblitz_input_semantics] at checked
  exact checked

end ICFallDegrees.BinaryDimensionTwoSystems
