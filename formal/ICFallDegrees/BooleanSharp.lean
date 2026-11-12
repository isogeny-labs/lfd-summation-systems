import ICFallDegrees.GridQuotient
import ICFallDegrees.QuotientLowerBound
import Mathlib.Tactic.Ring

noncomputable section

/-!
# A growing-dimension sharp Boolean family

For every field and r≥2, the Boolean field equations together with
1-∏(1-X_i) have exact LFD r+1. The ideal is proper. This is a general
Boolean-system family, not an asserted elliptic-curve summation family.
-/

namespace ICFallDegrees.BooleanSharp

open MvPolynomial MonomialOrder FiniteQuotient GridQuotient
open scoped BigOperators

variable {r : ℕ} {K : Type*} [Field K]

def sets (r : ℕ) (K : Type*) [Field K] : Fin r → Finset K := by
  classical
  exact fun _ => {0,1}
def cubeProduct (r : ℕ) (K : Type*) [Field K] : P r K := ∏ i : Fin r, (1-X i)
def h (r : ℕ) (K : Type*) [Field K] : P r K := 1-cubeProduct r K
def inputs (r : ℕ) (K : Type*) [Field K] : Set (P r K) :=
  gridInputs (sets r K)∪{h r K}

theorem field_polynomial (i : Fin r) :
    gridPolynomial (sets r K i) i=(X i : P r K)^2-X i := by
  classical
  simp [gridPolynomial,sets]
  ring

private theorem factor_degree (i : Fin r) : (1-X i : P r K).totalDegree=1 := by
  have identity : (1-X i : P r K)=-(X i-C 1) := by simp only [map_one]; ring
  rw [identity,totalDegree_neg,←degree_degLexDegree,degLex.degree_X_sub_C]
  simp

private theorem factor_nonzero (i : Fin r) : (1-X i : P r K)≠0 := by
  intro zero
  have he := congrArg (eval (fun _ => (0 : K))) zero
  simp at he

theorem cubeProduct_degree : (cubeProduct r K).totalDegree=r := by
  classical
  have product_nonzero : ∀ s : Finset (Fin r), (∏ i∈s, (1-X i : P r K))≠0 := by
    intro s zero
    have he := congrArg (eval (fun _ => (0 : K))) zero
    simp at he
  have all : ∀ s : Finset (Fin r), (∏ i∈s, (1-X i : P r K)).totalDegree=s.card := by
    intro s
    induction s using Finset.induction_on with
    | empty => simp
    | @insert i s hi ih =>
        rw [Finset.prod_insert hi,totalDegree_mul_eq (factor_nonzero i)
          (product_nonzero s),
          factor_degree,ih,Finset.card_insert_of_not_mem hi]
        omega
  simpa [cubeProduct] using all Finset.univ

theorem h_degree : (h r K).totalDegree=r := by
  have upper : (h r K).totalDegree≤r := by
    exact (totalDegree_sub_le _ _).trans (by simp [cubeProduct_degree])
  apply le_antisymm upper
  have identity : cubeProduct r K=1-h r K := by simp only [h]; ring
  have lower := totalDegree_sub_le (1 : P r K) (h r K)
  rw [←identity,cubeProduct_degree] at lower
  simpa using lower

theorem h_degreeOf (i : Fin r) : (h r K).degreeOf i≤1 := by
  classical
  have one_degree : degreeOf i (1 : P r K)=0 := degreeOf_C 1 i
  have product : (cubeProduct r K).degreeOf i≤1 := by
    calc
      _ ≤ ∑ j : Fin r, (1-X j : P r K).degreeOf i := degreeOf_prod_le _ _ _
      _ ≤ ∑ j : Fin r, if i=j then 1 else 0 := by
        apply Finset.sum_le_sum
        intro j _
        simpa only [one_degree,degreeOf_X,Nat.zero_max] using degreeOf_sub_le i (1 : P r K) (X j)
      _ = 1 := by simp
  exact (degreeOf_sub_le i (1 : P r K) (cubeProduct r K)).trans
    (max_le (by rw [one_degree]; omega) product)

theorem h_reduced : h r K∈reducedSpace (gridInputs (sets r K)) := by
  intro e he g hg
  rcases hg with ⟨i,rfl⟩
  rw [leading_gridPolynomial]
  intro divides
  have low := (degreeOf_le_iff.mp (h_degreeOf i)) e he
  have high := divides i
  simp [sets,zero_ne_one] at high
  omega

theorem upper_bound (dimension : 2≤r) : HasLFDBound (inputs r K) (r+1) := by
  apply normalized_boolean_lfd_bound (inputs r K)
  · intro i
    exact Or.inl ⟨i,field_polynomial i⟩
  · intro p hp
    rcases hp with ⟨i,rfl⟩ | rfl
    · have degree : (gridPolynomial (sets r K i) i).totalDegree=2 := by
        rw [←degree_degLexDegree,leading_gridPolynomial]
        simp [sets,zero_ne_one]
      exact degree.le.trans dimension
    · exact h_degree.le

theorem variable_mem_ideal (i : Fin r) : (X i : P r K)∈Ideal.span (inputs r K) := by
  classical
  let other : P r K := ∏ j∈Finset.univ.erase i, (1-X j)
  have split : cubeProduct r K=(1-X i)*other :=
    (Finset.mul_prod_erase Finset.univ (fun j : Fin r => (1-X j : P r K)) (Finset.mem_univ i)).symm
  have identity : (X i : P r K)=X i*h r K-other*(X i^2-X i) := by
    rw [h,split]
    ring
  rw [identity]
  apply (Ideal.span (inputs r K)).sub_mem
  · exact (Ideal.span (inputs r K)).mul_mem_left _ (Ideal.subset_span (Or.inr rfl))
  · exact (Ideal.span (inputs r K)).mul_mem_left _ (Ideal.subset_span (Or.inl ⟨i,field_polynomial i⟩))

theorem variable_not_at_r (dimension : 2≤r) (i : Fin r) :
    ¬InClosure (inputs r K) r (X i : P r K) := by
  let G := gridInputs (sets r K)
  have nonempty : ∀ i, (sets r K i).Nonempty := by intro i; simp [sets]
  have zero := reduced_grid_ideal_trivial (sets r K) nonempty
  have closed : closureSpace (inputs r K) r≤quotientSpan (Ideal.span G) (h r K) := by
    apply closure_le_quotient_span (inputs r K) G (gridInputs_nonzero _) zero
      (h r K) h_reduced r h_degree
    intro p hp _
    rcases hp with hp | rfl
    · refine ⟨0,?_⟩
      simpa using (Ideal.subset_span hp : p∈Ideal.span G)
    · exact ⟨1,by simp⟩
  intro hx
  have kernel := low_quotient_span_mem_ideal G (gridInputs_nonzero _) zero
    (h r K) h_reduced r h_degree (closed hx) (by simp; omega)
  have containing : Ideal.span G≤RingHom.ker (eval (fun _ => (1 : K))) := by
    apply Ideal.span_le.mpr
    rintro p ⟨j,rfl⟩
    change eval (fun _ => (1 : K)) (gridPolynomial (sets r K j) j)=0
    rw [field_polynomial]
    simp
  have he : eval (fun _ => (1 : K)) (X i : P r K)=0 := containing kernel
  simp at he

theorem exact_lfd (dimension : 2≤r) : IsLastFallDegree (inputs r K) (r+1) := by
  refine ⟨upper_bound dimension,?_⟩
  intro d hd bound
  let i : Fin r := ⟨0,by omega⟩
  have hx := bound (X i) (variable_mem_ideal i)
  exact variable_not_at_r dimension i (hx.mono (by simp; omega))

theorem proper_ideal : (1 : P r K)∉Ideal.span (inputs r K) := by
  intro hp
  have containing : Ideal.span (inputs r K)≤RingHom.ker (eval (fun _ => (0 : K))) := by
    apply Ideal.span_le.mpr
    intro p hp
    rcases hp with ⟨i,rfl⟩ | rfl
    · change eval (fun _ => (0 : K)) (gridPolynomial (sets r K i) i)=0
      rw [field_polynomial]
      simp
    · simp [h,cubeProduct]
  have he : eval (fun _ => (0 : K)) (1 : P r K)=0 := containing hp
  simp at he

end ICFallDegrees.BooleanSharp
