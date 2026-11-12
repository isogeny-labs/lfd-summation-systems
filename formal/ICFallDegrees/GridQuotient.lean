import ICFallDegrees.FiniteQuotient
import Mathlib.Combinatorics.Nullstellensatz

noncomputable section

/-!
# Explicit grid and Boolean field equations

The reduced-ideal hypothesis is proved using mathlib's grid vanishing theorem.
Thus the Boolean bound has no assumed GB, normal form, closure rank, or
solution count. It applies to all additional equations with the stated input
degree, including actual normalized binary descents.
-/

namespace ICFallDegrees.GridQuotient

open MvPolynomial MonomialOrder
open scoped BigOperators MonomialOrder
open FiniteQuotient

variable {r : ℕ} {K : Type*} [Field K]

def gridPolynomial (S : Finset K) (i : Fin r) : P r K :=
  ∏ a∈S, (X i-C a)

def gridInputs (S : Fin r → Finset K) : Set (P r K) :=
  Set.range (fun i => gridPolynomial (S i) i)

theorem leading_gridPolynomial (S : Finset K) (i : Fin r) :
    degLex.degree (gridPolynomial S i)=Finsupp.single i S.card := by
  classical
  simp only [gridPolynomial]
  rw [degree_prod_of_regular]
  · simp [Finset.sum_congr rfl (fun a _ => degLex.degree_X_sub_C i a)]
  · intro a _
    rw [degLex.monic_X_sub_C]
    exact isRegular_one

theorem gridPolynomial_nonzero (S : Finset K) (i : Fin r) : gridPolynomial S i≠0 := by
  have monic : degLex.Monic (gridPolynomial S i) :=
    MonomialOrder.Monic.prod (fun a _ => degLex.monic_X_sub_C i a)
  exact monic.ne_zero

theorem gridInputs_nonzero (S : Fin r → Finset K) : (0 : P r K)∉gridInputs S := by
  rintro ⟨i,zero⟩
  exact gridPolynomial_nonzero (S i) i zero

theorem standard_coordinate_lt (S : Fin r → Finset K) (e : Fin r →₀ ℕ)
    (standard : Standard (gridInputs S) e) (i : Fin r) : e i<(S i).card := by
  have hd := standard (gridPolynomial (S i) i) ⟨i,rfl⟩
  rw [leading_gridPolynomial,Finsupp.single_le_iff] at hd
  omega

theorem standard_degree_bound (S : Fin r → Finset K) (e : Fin r →₀ ℕ)
    (standard : Standard (gridInputs S) e) :
    e.degree≤∑ i : Fin r, ((S i).card-1) := by
  have sum : e.degree=∑ i : Fin r, e i := by
    exact Finsupp.sum_fintype e (fun _ n => n) (fun _ => rfl)
  rw [sum]
  apply Finset.sum_le_sum
  intro i _
  have := standard_coordinate_lt S e standard i
  omega

theorem reduced_grid_ideal_trivial (S : Fin r → Finset K)
    (nonempty : ∀ i, (S i).Nonempty) : ReducedIdealTrivial (gridInputs S) := by
  classical
  intro p hp reduced
  apply eq_zero_of_eval_zero_at_prod_finset p S
  · intro i
    rw [degreeOf_eq_sup,Finset.sup_lt_iff (Finset.card_pos.mpr (nonempty i))]
    intro e he
    exact standard_coordinate_lt S e (reduced e he) i
  · intro point belongs
    have containing : Ideal.span (gridInputs S)≤RingHom.ker (eval point) := by
      apply Ideal.span_le.mpr
      rintro f ⟨i,rfl⟩
      change eval point (gridPolynomial (S i) i)=0
      simp only [gridPolynomial,map_prod]
      apply Finset.prod_eq_zero (belongs i)
      simp
    exact containing hp

/-- Exact explicit field-equation corollary, including all input degrees. -/
theorem grid_lfd_bound (S : Fin r → Finset K)
    (nonempty : ∀ i, (S i).Nonempty)
    (F : Set (P r K)) (included : gridInputs S⊆F)
    (delta : ℕ) (input_degree : ∀ f∈F, f.totalDegree≤delta) :
    HasLFDBound F (max ((∑ i : Fin r, ((S i).card-1))+1) delta) := by
  exact finite_quotient_lfd_of_reduced_zero F (gridInputs S) included (gridInputs_nonzero S)
    (reduced_grid_ideal_trivial S nonempty) _ delta (standard_degree_bound S) input_degree

/-- All Boolean field equations are explicit inputs; the other equations
may be arbitrary polynomials of degree at most delta. -/
theorem boolean_lfd_bound (F : Set (P r K))
    (field_equations : ∀ i, ((X i : P r K)^2-X i)∈F)
    (delta : ℕ) (input_degree : ∀ f∈F, f.totalDegree≤delta) :
    HasLFDBound F (max (r+1) delta) := by
  classical
  let S : Fin r → Finset K := fun _ => {0,1}
  have included : gridInputs S⊆F := by
    rintro f ⟨i,rfl⟩
    have identity : gridPolynomial (S i) i=(X i : P r K)^2-X i := by
      simp [gridPolynomial,S]
      ring
    change gridPolynomial (S i) i∈F
    rw [identity]
    exact field_equations i
  have bound := grid_lfd_bound S (by intro i; simp [S]) F included delta input_degree
  simpa [S,zero_ne_one] using bound

/-- The full-ideal r+1 bound for multilinear normalized Boolean systems. -/
theorem normalized_boolean_lfd_bound (F : Set (P r K))
    (field_equations : ∀ i, ((X i : P r K)^2-X i)∈F)
    (input_degree : ∀ f∈F, f.totalDegree≤r) : HasLFDBound F (r+1) := by
  simpa using boolean_lfd_bound F field_equations r input_degree

end ICFallDegrees.GridQuotient
