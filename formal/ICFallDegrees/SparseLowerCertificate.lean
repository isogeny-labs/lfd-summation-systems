import ICFallDegrees.SparseCoefficient
import ICFallDegrees.FilteredSpan
import Mathlib.Algebra.BigOperators.Fin

namespace ICFallDegrees.SparseLowerCertificate

open SparsePolynomial MvPolynomial
open scoped BigOperators

variable {r n : ℕ} {K : Type*} [CommRing K] [Nontrivial K] [DecidableEq K]

def weight (e : Exponent r) : ℕ := ∑ i : Fin r, e i

def project (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r) (p : Sparse r K) : Sparse r K :=
  SparsePolynomial.sum (List.ofFn fun i => scale (coefficient (pivot i) p) (basis i))

def contains (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r) (p : Sparse r K) : Bool :=
  (add p (scale (-1) (project basis pivot p))).isEmpty

def rowCheck (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r) (D : ℕ) (i : Fin n) : Bool :=
  decide (0<weight (pivot i)) &&
  decide (degree (basis i)≤weight (pivot i)) &&
  decide (weight (pivot i)≤D) &&
  (List.ofFn (fun j : Fin n => j)).all (fun j =>
    decide (coefficient (pivot i) (basis j)=if i=j then 1 else 0)) &&
  (if weight (pivot i)<D then
    (List.ofFn (fun v : Fin r => v)).all (fun v => contains basis pivot (mul (variableData v) (basis i)))
   else true)

/-- Every input is represented, not skipped using an unverified input-degree
estimate. Rows, pivots, all low-row variable products, and input inclusions
are checked on exact sparse data; no externally reported rank is accepted. -/
def check (inputs : List (Sparse r K)) (basis : Fin n → Sparse r K)
    (pivot : Fin n → Exponent r) (D : ℕ) : Bool :=
  (List.ofFn (fun i : Fin n => i)).all (rowCheck basis pivot D) &&
  inputs.all (contains basis pivot)

theorem weight_toExponent (e : Exponent r) : (toExponent e).degree=weight e := by
  exact Finsupp.sum_fintype (toExponent e) (fun _ n => n) (fun _ => rfl)

theorem meaning_project (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r)
    (p : Sparse r K) : meaning (project basis pivot p)=
      FilteredSpan.rowMap (fun i => meaning (basis i)) (fun i => coefficient (pivot i) p) := by
  simp [project,meaning_sum,List.map_ofFn,List.sum_ofFn,meaning_scale,
    FilteredSpan.rowMap,Fintype.linearCombination_apply]

theorem contains_sound (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r)
    (p : Sparse r K) (checked : contains basis pivot p=true) :
    meaning p∈FilteredSpan.rowSpace (fun i => meaning (basis i)) := by
  have empty : add p (scale (-1) (project basis pivot p))=[] := by
    simpa only [contains,List.isEmpty_iff] using checked
  have equal := congrArg meaning empty
  rw [meaning_add,meaning_scale] at equal
  have identity : meaning p=meaning (project basis pivot p) := by
    simpa only [meaning,neg_one_smul,←sub_eq_add_neg,sub_eq_zero] using equal
  rw [identity,meaning_project]
  exact ⟨_,rfl⟩

theorem rowCheck_spec (basis : Fin n → Sparse r K) (pivot : Fin n → Exponent r)
    (D : ℕ) (i : Fin n) (checked : rowCheck basis pivot D i=true) :
    (0<weight (pivot i)) ∧ (degree (basis i)≤weight (pivot i)) ∧ (weight (pivot i)≤D) ∧
    (∀ j, coefficient (pivot i) (basis j)=if i=j then 1 else 0) ∧
    (weight (pivot i)<D → ∀ v, contains basis pivot (mul (variableData v) (basis i))=true) := by
  by_cases degree : weight (pivot i)<D
  all_goals simpa [rowCheck,List.all_eq_true,List.forall_mem_ofFn_iff,and_assoc,degree] using checked

/-- A successful check proves exclusion from the full iterated ordinary
V-space, including all multipliers of newly fallen rows. -/
theorem check_excludes_one (inputs : List (Sparse r K)) (basis : Fin n → Sparse r K)
    (pivot : Fin n → Exponent r) (D : ℕ)
    (checked : check inputs basis pivot D=true) :
    ¬InClosure {p | p∈inputs.map meaning} D (1 : MvPolynomial (Fin r) K) := by
  have checks : (∀ i, rowCheck basis pivot D i=true) ∧
      (∀ p∈inputs, contains basis pivot p=true) := by
    simpa [check,List.all_eq_true,List.forall_mem_ofFn_iff] using checked
  have rows := fun i => rowCheck_spec basis pivot D i (checks.1 i)
  apply FilteredSpan.one_not_in_closure {p | p∈inputs.map meaning}
    (fun i => meaning (basis i)) (fun i => toExponent (pivot i)) D
  · intro i j
    rw [meaning_coefficient]
    exact (rows i).2.2.2.1 j
  · intro i
    rw [weight_toExponent]
    exact (rows i).1
  · intro i degree v
    rw [weight_toExponent] at degree
    have accepted := (rows i).2.2.2.2 degree v
    have member := contains_sound basis pivot (mul (variableData v) (basis i)) accepted
    simpa only [meaning_mul,meaning_variable] using member
  · intro p hp
    obtain ⟨data,member,rfl⟩ := List.mem_map.mp hp
    exact contains_sound basis pivot data (checks.2 data member)

end ICFallDegrees.SparseLowerCertificate
