import ICFallDegrees.HomogeneousFirstFall
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Ring
import Mathlib.Data.Fin.VecNotation

/-! A common-degree quadratic top system with an exact degree-three
non-Koszul syzygy. The lower-degree tails are arbitrary; actual S3
specialization is connected separately. -/
noncomputable section
namespace ICFallDegrees.QuadraticHomogeneousFamily
open MvPolynomial HomogeneousFirstFall
open scoped BigOperators
variable {K : Type*} [CommRing K] [Nontrivial K]
abbrev P (K : Type*) [CommRing K] := MvPolynomial (Fin 6) K

def top (β : K) : P K := (X 1+C β*X 0)*(X 4+C β*X 3)
def q (β : K) : Fin 7 → P K := ![X 0^2,X 1^2,X 2^2,X 3^2,X 4^2,X 5^2,top β]
def pivot : Fin 7 → Fin 6 →₀ ℕ := ![Finsupp.single 0 2,Finsupp.single 1 2,
  Finsupp.single 2 2,Finsupp.single 3 2,Finsupp.single 4 2,Finsupp.single 5 2,
  Finsupp.single 1 1+Finsupp.single 4 1]
def lower (tail : P K) : Fin 7 → P K := ![-X 1,-X 2,-1,-X 4,-X 5,-1,tail]
def system (β : K) (tail : P K) (i : Fin 7) : P K := q β i+lower tail i

def A (β : K) : P K := X 1-C β*X 0
def B (β : K) : P K := X 4+C β*X 3
def witness (β : K) : Fin 7 → P K := ![C (β^2)*B β,-B β,0,0,0,0,A β]

theorem coefficient_pivots (β : K) (i j : Fin 7) :
    coeff (pivot i) (q β j)=if i=j then 1 else 0 := by
  fin_cases i <;> fin_cases j
  all_goals simp [q,pivot,top,add_mul,mul_add]
  all_goals try simp only [X_pow_eq_monomial]
  all_goals simp [X,monomial_mul,C_mul_monomial,coeff_monomial,
    Finsupp.ext_iff,Fin.forall_fin_succ,Finsupp.single_apply]

theorem q_homogeneous (β : K) (i : Fin 7) : (q β i).IsHomogeneous 2 := by
  have linear (a b : Fin 6) : (X a+C β*X b : P K).IsHomogeneous 1 :=
    (isHomogeneous_X K a).add ((isHomogeneous_X K b).C_mul β)
  fin_cases i
  · exact isHomogeneous_X_pow (R:=K) 0 2
  · exact isHomogeneous_X_pow (R:=K) 1 2
  · exact isHomogeneous_X_pow (R:=K) 2 2
  · exact isHomogeneous_X_pow (R:=K) 3 2
  · exact isHomogeneous_X_pow (R:=K) 4 2
  · exact isHomogeneous_X_pow (R:=K) 5 2
  · exact (linear 1 0).mul (linear 4 3)

theorem q_nonzero (β : K) (i : Fin 7) : q β i≠0 := by
  intro zero
  have equal := coefficient_pivots β i i
  rw [zero,coeff_zero,if_pos rfl] at equal
  exact zero_ne_one equal

theorem degree_and_top (β : K) (tail : P K) (small : tail.totalDegree<2) (i : Fin 7) :
    (system β tail i).totalDegree=2 ∧
      homogeneousComponent (system β tail i).totalDegree (system β tail i)=q β i := by
  apply degree_and_top_of_add (q_homogeneous β i) (q_nonzero β i)
  fin_cases i <;> simp [lower,ICFallDegrees.totalDegree_neg_eq,small]

theorem witness_homogeneous (β : K) :
    IsHomogeneousSyzygy (q β) (fun _ => 2) 3 (witness β) := by
  have hB : (B β).IsHomogeneous 1 :=
    (isHomogeneous_X K 4).add ((isHomogeneous_X K 3).C_mul β)
  have hA : (A β).IsHomogeneous 1 :=
    (isHomogeneous_X K 1).sub ((isHomogeneous_X K 0).C_mul β)
  constructor
  · intro i
    fin_cases i
    all_goals first | left; rfl | right; refine ⟨by decide,?_⟩
    · exact hB.C_mul (β^2)
    · exact hB.neg
    · exact hA
  · simp [witness,q,top,A,B,Fin.sum_univ_succ,map_pow]
    ring

theorem origin_values (β : K) (i : Fin 7) : constantCoeff (q β i)=0 := by
  fin_cases i <;> simp [q,top]

theorem origin_derivatives (β : K) (i : Fin 7) :
    constantCoeff (pderiv 1 (q β i))=0 := by
  fin_cases i <;> simp [q,top,pderiv_mul,pderiv_pow]

theorem witness_not_koszul (β : K) : witness β∉koszulSpan (q β) := by
  intro member
  have zero := (koszulSpan_jet_zero (q β) 1 (origin_values β)
    (origin_derivatives β) member 6).2
  have unit : constantCoeff (pderiv 1 (witness β 6))=1 := by simp [witness,A]
  rw [unit] at zero
  exact one_ne_zero zero

theorem exact_first_fall_three (β : K) (tail : P K) (small : tail.totalDegree<2) :
    IsKoszulFirstFall (system β tail) 3 := by
  have tops : topForms (system β tail)=q β := by
    funext i
    exact (degree_and_top β tail small i).2
  have degrees : (fun i => (system β tail i).totalDegree)=(fun _ => 2) := by
    funext i
    exact (degree_and_top β tail small i).1
  constructor
  · change HasFall (topForms (system β tail)) (fun i => (system β tail i).totalDegree) 3
    rw [tops,degrees]
    exact ⟨witness β,witness_homogeneous β,witness_not_koszul β⟩
  · intro d hd
    change ¬HasFall (topForms (system β tail)) (fun i => (system β tail i).totalDegree) d
    rw [tops,degrees]
    exact no_low_fall (q β) pivot (coefficient_pivots β) 2 d (by omega)

end ICFallDegrees.QuadraticHomogeneousFamily
