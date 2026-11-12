import ICFallDegrees.Closure
import Mathlib.Algebra.BigOperators.Finsupp.Basic

/-!
# Executable sparse arithmetic and its semantic bridge

The data path uses finite exponent vectors and coefficient operations. Only
the interpretation into mathlib's `MvPolynomial` is noncomputable. Combining
equal terms allows genuine degree falls; no field-equation reduction occurs.
-/

namespace ICFallDegrees.SparsePolynomial

open MvPolynomial
open scoped BigOperators

abbrev Exponent (r : ℕ) := Fin r → ℕ
abbrev Sparse (r : ℕ) (K : Type*) := List (Exponent r × K)

variable {r : ℕ} {K : Type*} [CommRing K] [DecidableEq K]

def sameExponent (e f : Exponent r) : Bool := List.ofFn e==List.ofFn f

def insert (e : Exponent r) (c : K) : Sparse r K → Sparse r K
  | [] => if c=0 then [] else [(e,c)]
  | (f,d)::p =>
      if c=0 then (f,d)::p
      else if sameExponent e f then
        if c+d=0 then p else (f,c+d)::p
      else (f,d)::insert e c p

def normalize : Sparse r K → Sparse r K
  | [] => []
  | (e,c)::p => insert e c (normalize p)

def add (p q : Sparse r K) : Sparse r K := normalize (p++q)
def scale (c : K) (p : Sparse r K) : Sparse r K :=
  normalize (p.map fun (e,d) => (e,c*d))
def mul (p q : Sparse r K) : Sparse r K :=
  normalize (p.flatMap fun (e,c) => q.map fun (f,d) => ((fun i => e i+f i),c*d))

def degree : Sparse r K → ℕ
  | [] => 0
  | (e,_)::p => max (∑ i : Fin r, e i) (degree p)

noncomputable def toExponent (e : Exponent r) : Fin r →₀ ℕ :=
  Finsupp.onFinset Finset.univ e (by intro i _; exact Finset.mem_univ i)

noncomputable def meaning : Sparse r K → MvPolynomial (Fin r) K
  | [] => 0
  | (e,c)::p => monomial (toExponent e) c+meaning p

theorem toExponent_add (e f : Exponent r) :
    toExponent (fun i => e i+f i)=toExponent e+toExponent f := by ext i; rfl

@[simp] theorem sameExponent_eq_true (e f : Exponent r) : sameExponent e f=true ↔ e=f := by
  simp [sameExponent]

theorem meaning_insert (e : Exponent r) (c : K) (p : Sparse r K) :
    meaning (insert e c p)=monomial (toExponent e) c+meaning p := by
  classical
  induction p with
  | nil => by_cases hc : c=0 <;> simp [insert,hc,meaning]
  | cons term p ih =>
      obtain ⟨f,d⟩ := term
      by_cases hc : c=0
      · simp [insert,hc,meaning]
      · by_cases he : e=f
        · subst f
          by_cases zero : c+d=0
          · simp only [insert,if_neg hc,if_pos rfl,if_pos zero]
            rw [meaning,←add_assoc,←map_add,zero,monomial_zero,zero_add]
            simp
          · simp [insert,hc,zero,meaning,map_add,add_assoc]
        · simp [insert,hc,he,meaning,ih,add_assoc,add_left_comm]

theorem meaning_normalize (p : Sparse r K) : meaning (normalize p)=meaning p := by
  induction p with
  | nil => rfl
  | cons term p ih => obtain ⟨e,c⟩ := term; simp [normalize,meaning_insert,meaning,ih]

omit [DecidableEq K] in
theorem meaning_append (p q : Sparse r K) : meaning (p++q)=meaning p+meaning q := by
  induction p with
  | nil => simp [meaning]
  | cons term p ih => obtain ⟨e,c⟩ := term; simp [meaning,ih,add_assoc]

theorem meaning_add (p q : Sparse r K) : meaning (add p q)=meaning p+meaning q := by
  rw [add,meaning_normalize,meaning_append]

theorem meaning_scale (c : K) (p : Sparse r K) : meaning (scale c p)=c • meaning p := by
  rw [scale,meaning_normalize]
  induction p with
  | nil => simp [meaning]
  | cons term p ih =>
      obtain ⟨e,d⟩ := term
      simp only [List.map_cons,meaning,smul_add,smul_eq_C_mul,C_mul_monomial] at ih ⊢
      rw [ih]

omit [DecidableEq K] in
theorem meaning_term_mul (e : Exponent r) (c : K) (q : Sparse r K) :
    meaning (q.map fun (f,d) => ((fun i => e i+f i),c*d))=
      monomial (toExponent e) c*meaning q := by
  induction q with
  | nil => simp [meaning]
  | cons term q ih =>
      obtain ⟨f,d⟩ := term
      simp [meaning,toExponent_add,monomial_mul,mul_add,ih]

theorem meaning_mul (p q : Sparse r K) : meaning (mul p q)=meaning p*meaning q := by
  rw [mul,meaning_normalize]
  induction p with
  | nil => simp [meaning]
  | cons term p ih =>
      obtain ⟨e,c⟩ := term
      simp [List.flatMap_cons,meaning_append,meaning_term_mul,meaning,add_mul,ih]

omit [DecidableEq K] in
/-- A conservative degree computed entirely on data. Normalization can lower
it by exact coefficient cancellation before a multiplication is charged. -/
theorem meaning_degree_le (p : Sparse r K) : (meaning p).totalDegree≤degree p := by
  induction p with
  | nil => simp [meaning,degree]
  | cons term p ih =>
      obtain ⟨e,c⟩ := term
      apply (totalDegree_add _ _).trans
      apply max_le_max ?_ ih
      have sum : (toExponent e).sum (fun _ n => n)=∑ i : Fin r, e i :=
        Finsupp.sum_fintype _ _ (fun _ => rfl)
      exact (totalDegree_monomial_le _ _).trans_eq sum

end ICFallDegrees.SparsePolynomial
