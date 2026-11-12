import ICFallDegrees.PowerIdeals
import Mathlib.Tactic.Ring

noncomputable section
namespace ICFallDegrees.LowDegreeCoefficients
open MvPolynomial
open scoped BigOperators
variable {r : ℕ} {K : Type*} [CommRing K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

def pair (i j : Fin r) : Fin r →₀ ℕ := Finsupp.single i 1+Finsupp.single j 1
def triple (i j k : Fin r) : Fin r →₀ ℕ := pair i j+Finsupp.single k 1

theorem degree_one_single (e : Fin r →₀ ℕ) (degree : e.degree=1) :
    ∃ i, e=Finsupp.single i 1 := by
  classical
  induction e using Finsupp.induction with
  | zero => simp at degree
  | single_add i n e hi hn ih =>
      rw [Finsupp.degree_add,Finsupp.degree_single] at degree
      have n1 : n=1 := by omega
      have zero : e=0 := (Finsupp.degree_eq_zero_iff e).mp (by omega)
      exact ⟨i,by simp [n1,zero]⟩

theorem homogeneous_linear_eq_zero {p : P r K} (homogeneous : p.IsHomogeneous 1)
    (coefficients : ∀ i, coeff (Finsupp.single i 1) p=0) : p=0 := by
  apply MvPolynomial.eq_zero_iff.mpr
  intro e
  by_cases degree : e.degree=1
  · obtain ⟨i,rfl⟩ := degree_one_single e degree
    exact coefficients i
  · exact homogeneous.coeff_eq_zero degree

theorem pair_polynomial (i j : Fin r) :
    (X i : P r K)*X j=monomial (pair i j) 1 := by
  simp only [X,monomial_mul,one_mul,pair]

theorem coefficient_pair_sum (p : P r K) (e : Fin r →₀ ℕ) (i j : Fin r) :
    coeff (e+pair i j) (p*X i*X j)=coeff e p := by
  rw [mul_assoc,pair_polynomial,coeff_mul_monomial,mul_one]

theorem coefficient_triple_12 (p : P r K) (i j k : Fin r) :
    coeff (triple i j k) (p*X i*X j)=coeff (Finsupp.single k 1) p := by
  have order : triple i j k=Finsupp.single k 1+pair i j := by unfold triple; exact add_comm _ _
  rw [order,coefficient_pair_sum]

theorem coefficient_triple_13 (p : P r K) (i j k : Fin r) :
    coeff (triple i j k) (p*X i*X k)=coeff (Finsupp.single j 1) p := by
  have order : triple i j k=Finsupp.single j 1+pair i k := by unfold triple pair; abel
  rw [order,coefficient_pair_sum]

theorem coefficient_triple_23 (p : P r K) (i j k : Fin r) :
    coeff (triple i j k) (p*X j*X k)=coeff (Finsupp.single i 1) p := by
  have order : triple i j k=Finsupp.single i 1+pair j k := by unfold triple pair; abel
  rw [order,coefficient_pair_sum]

theorem coefficient_pair_absent_left (p : P r K) (e : Fin r →₀ ℕ) (i j : Fin r)
    (absent : e i=0) : coeff e (p*X i*X j)=0 := by
  classical
  rw [mul_assoc,pair_polynomial,coeff_mul_monomial']
  have no : ¬pair i j≤e := by
    intro bound
    have at_i := bound i
    simp only [pair,Finsupp.add_apply,Finsupp.single_eq_same,absent] at at_i
    omega
  rw [if_neg no]

theorem coefficient_pair_absent_right (p : P r K) (e : Fin r →₀ ℕ) (i j : Fin r)
    (absent : e j=0) : coeff e (p*X i*X j)=0 := by
  have identity : p*X i*X j=p*X j*X i := by ring
  rw [identity]
  exact coefficient_pair_absent_left p e j i absent

theorem pair_squarefree (i j : Fin r) (distinct : i≠j) : ∀ s, pair i j s≤1 := by
  intro s
  simp only [pair,Finsupp.add_apply,Finsupp.single_apply]
  split_ifs <;> simp_all

theorem triple_squarefree (i j k : Fin r) (hij : i≠j) (hik : i≠k) (hjk : j≠k) :
    ∀ s, triple i j k s≤1 := by
  intro s
  simp only [triple,pair,Finsupp.add_apply,Finsupp.single_apply]
  split_ifs <;> simp_all

theorem coefficient_mul_power_absent (p : P r K) (e : Fin r →₀ ℕ) (i : Fin r) (d : ℕ)
    (small : e i<d) : coeff e (p*X i^d)=0 := by
  classical
  rw [X_pow_eq_monomial,coeff_mul_monomial']
  have no : ¬Finsupp.single i d≤e := by
    rw [Finsupp.single_le_iff]
    omega
  rw [if_neg no]

theorem coefficient_square_pair (i j k : Fin r) (distinct : j≠k) :
    coeff (Finsupp.single i 2) ((X j : P r K)*X k)=0 := by
  classical
  by_cases equal : i=j
  · subst j
    have absent : (Finsupp.single i 2 : Fin r →₀ ℕ) k=0 := by exact Finsupp.single_eq_of_ne distinct
    simpa only [one_mul] using coefficient_pair_absent_right (1 : P r K) _ i k absent
  · have absent : (Finsupp.single i 2 : Fin r →₀ ℕ) j=0 := by exact Finsupp.single_eq_of_ne equal
    simpa only [one_mul] using coefficient_pair_absent_left (1 : P r K) _ j k absent

theorem coefficient_square_sum (h : Fin r → P r K) (i j : Fin r) :
    coeff (Finsupp.single i 2+Finsupp.single j 1) (∑ s, h s*X s^2)=
      coeff (Finsupp.single j 1) (h i) := by
  classical
  rw [coeff_sum,Finset.sum_eq_single i]
  · have reordered : Finsupp.single i 2+Finsupp.single j 1=
        Finsupp.single j 1+Finsupp.single i 2 := add_comm _ _
    rw [reordered,X_pow_eq_monomial,coeff_mul_monomial,mul_one]
  · intro s hs distinct
    apply coefficient_mul_power_absent
    simp only [Finsupp.add_apply,Finsupp.single_apply]
    split_ifs <;> simp_all
  · simp

end ICFallDegrees.LowDegreeCoefficients
