import ICFallDegrees.VariableClosure
import Mathlib.LinearAlgebra.Finsupp.LinearCombination
import Mathlib.Algebra.Module.Submodule.Range
import Mathlib.Data.Finsupp.Weight

noncomputable section

namespace ICFallDegrees.FilteredSpan

open MvPolynomial
open scoped BigOperators

variable {r n : ℕ} {K : Type*} [CommRing K] [Nontrivial K]

abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

def rowMap (b : Fin n → P r K) : (Fin n → K) →ₗ[K] P r K :=
  Fintype.linearCombination K b

def rowSpace (b : Fin n → P r K) : Submodule K (P r K) :=
  LinearMap.range (rowMap b)

theorem coefficient_rowMap (b : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ)
    (rref : ∀ i j, coeff (pivot i) (b j)=if i=j then 1 else 0)
    (c : Fin n → K) (i : Fin n) : coeff (pivot i) (rowMap b c)=c i := by
  classical
  simp [rowMap,Fintype.linearCombination_apply,coeff_sum,smul_eq_C_mul,coeff_C_mul,rref,mul_ite]

/-- The pivot identities certify the filtration. No rank is assumed: the
coefficient of each pivot in a represented row is its representation scalar. -/
theorem rowSpace_variable_stable
    (b : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ) (D : ℕ)
    (rref : ∀ i j, coeff (pivot i) (b j)=if i=j then 1 else 0)
    (products : ∀ i, (pivot i).degree<D → ∀ v, X v*b i∈rowSpace b) :
    ∀ v p, p∈rowSpace b → p.totalDegree<D → X v*p∈rowSpace b := by
  classical
  intro v p hp degree
  obtain ⟨c,hc⟩ := hp
  have identity : X v*p=∑ i : Fin n, c i • (X v*b i) := by
    rw [←hc]
    simp [rowMap,Fintype.linearCombination_apply,Finset.mul_sum,smul_eq_C_mul,mul_left_comm,mul_assoc]
  rw [identity]
  apply (rowSpace b).sum_mem
  intro i _
  by_cases zero : c i=0
  · simp [zero]
  · have nonzero : coeff (pivot i) p≠0 := by
      rw [←hc,coefficient_rowMap b pivot rref]
      exact zero
    have support : pivot i∈p.support := mem_support_iff.mpr nonzero
    have smaller : (pivot i).degree<D := lt_of_le_of_lt (le_totalDegree support) degree
    exact (rowSpace b).smul_mem _ (products i smaller v)

theorem closure_le_rowSpace
    (F : Set (P r K)) (b : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ) (D : ℕ)
    (rref : ∀ i j, coeff (pivot i) (b j)=if i=j then 1 else 0)
    (products : ∀ i, (pivot i).degree<D → ∀ v, X v*b i∈rowSpace b)
    (inputs : ∀ p∈F, p∈rowSpace b) : closureSpace F D≤rowSpace b := by
  apply closure_minimal_by_variables F D (rowSpace b) (fun p hp _ => inputs p hp)
  exact rowSpace_variable_stable b pivot D rref products

theorem one_not_mem_rowSpace
    (b : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ)
    (rref : ∀ i j, coeff (pivot i) (b j)=if i=j then 1 else 0)
    (positive : ∀ i, 0<(pivot i).degree) : (1 : P r K)∉rowSpace b := by
  intro member
  obtain ⟨c,hc⟩ := member
  have zero : ∀ i, c i=0 := by
    intro i
    rw [←coefficient_rowMap b pivot rref c i,hc]
    apply coeff_eq_zero_of_totalDegree_lt
    simpa using positive i
  simp [rowMap,Fintype.linearCombination_apply,zero] at hc

theorem one_not_in_closure
    (F : Set (P r K)) (b : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ) (D : ℕ)
    (rref : ∀ i j, coeff (pivot i) (b j)=if i=j then 1 else 0)
    (positive : ∀ i, 0<(pivot i).degree)
    (products : ∀ i, (pivot i).degree<D → ∀ v, X v*b i∈rowSpace b)
    (inputs : ∀ p∈F, p∈rowSpace b) : ¬InClosure F D (1 : P r K) := by
  intro derived
  exact one_not_mem_rowSpace b pivot rref positive
    (closure_le_rowSpace F b pivot D rref products inputs derived)

end ICFallDegrees.FilteredSpan
