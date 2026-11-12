import ICFallDegrees.GeneratedPrimeLower
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FinCases

noncomputable section

namespace ICFallDegrees.Prime257System

open MvPolynomial SparsePolynomial
abbrev K := GeneratedPrimeLower.K
abbrev P := MvPolynomial (Fin 6) K

def normalForm (a b t : K) (u v w z : P) : P :=
  v*z+C (-2*t)*(v*w+u*z)+C (-2*(t^2+a))*u*w+
    C (t^2)*(v+z)+C (-2*a*t-4*b)*(u+w)+C (a^2-4*b*t)

def curveS3 : P :=
  (X 0-X 3)^2*C ((42 : K)^2)-
    2*((X 0+X 3)*(X 0*X 3+C (105 : K))+C (2*(196 : K)))*C (42 : K)+
    (X 0*X 3-C (105 : K))^2-C (4*(196 : K))*(X 0+X 3)

def g : P := normalForm 105 196 42 (X 0) (X 1) (X 3) (X 4)
def chains : List P := [X 0^2-X 1,X 1^2-X 2,X 2^2-1,
  X 3^2-X 4,X 4^2-X 5,X 5^2-1]
def normalizedInputs : Set P := {p | p∈chains++[g]}

theorem target_on_curve : (56 : K)^2=(42 : K)^3+105*42+196 := by decide
theorem nonsingular_discriminant : 4*(105 : K)^3+27*(196 : K)^2≠0 := by decide
theorem derivative_zero : 3*(42 : K)^2+105=0 := by decide
theorem target_y_nonzero : (56 : K)≠0 := by decide
theorem target_outside_coordinate_base : (42 : K)^8≠1 := by decide

theorem raw_reduction_identity : curveS3-g=
    (X 3-C (42 : K))^2*(X 0^2-X 1)+
      (X 1-C (2*(42 : K))*X 0+C ((42 : K)^2))*(X 3^2-X 4) := by
  simp only [curveS3,g,normalForm,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring

theorem g_form : g=
    X 1*X 4+C (173 : K)*(X 1*X 3+X 0*X 4)+C (117 : K)*X 0*X 3+
      C (222 : K)*(X 1+X 4)+C (162 : K)*(X 0+X 3)+C (199 : K) := by
  have beta : -2*(42 : K)=173 := by decide
  have alpha : -2*((42 : K)^2+105)=117 := by decide
  have square : (42 : K)^2=222 := by decide
  have linear : -2*(105 : K)*42-4*196=162 := by decide
  have constant : (105 : K)^2-4*196*42=199 := by decide
  rw [g,normalForm,beta,alpha,square,linear,constant]

theorem toExponent_exponent (a b c d e f : ℕ) :
    toExponent (GeneratedPrimeLower.exponent a b c d e f)=
      Finsupp.single 0 a+Finsupp.single 1 b+Finsupp.single 2 c+
        Finsupp.single 3 d+Finsupp.single 4 e+Finsupp.single 5 f := by
  ext i
  simp only [toExponent,Finsupp.onFinset_apply,Finsupp.add_apply,Finsupp.single_apply]
  fin_cases i <;> simp [GeneratedPrimeLower.exponent]
  all_goals rfl

theorem input_list_semantics : GeneratedPrimeLower.inputs.map meaning=chains++[g] := by
  have negative : (256 : K)= -1 := by decide
  simp only [GeneratedPrimeLower.inputs,List.map_cons,List.map_nil,meaning,toExponent_exponent]
  simp [chains,g_form,negative,X_pow_eq_monomial,X,monomial_mul,monomial_pow,
    C_mul_monomial,mul_add,add_mul,sub_eq_add_neg,add_assoc]
  ring

theorem input_set_semantics : {p : P | p∈GeneratedPrimeLower.inputs.map meaning}=normalizedInputs := by
  rw [input_list_semantics]
  rfl

theorem one_not_at_three : ¬InClosure normalizedInputs 3 (1 : P) := by
  rw [←input_set_semantics]
  exact GeneratedPrimeLower.excludes_one

end ICFallDegrees.Prime257System
