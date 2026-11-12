import Mathlib.Algebra.Field.ZMod
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.Algebra.Polynomial.SpecificDegree
import Mathlib.Tactic.FinCases
import Mathlib.Algebra.CharP.Algebra
import Mathlib.Tactic.ComputeDegree
import Mathlib.Tactic.Ring

/-! Concrete extension fields for the dimension-two descent examples.
Irreducibility is proved, not supplied as an assumption. -/
noncomputable section
namespace ICFallDegrees.BinaryExtensionFields
open Polynomial
abbrev k := ZMod 2

def modulus8 : Polynomial k := X^3+X+1
def modulus16 : Polynomial k := X^4+X+1

theorem modulus8_monic : modulus8.Monic := by unfold modulus8; monicity <;> norm_num
theorem modulus16_monic : modulus16.Monic := by unfold modulus16; monicity <;> norm_num
theorem modulus8_degree : modulus8.natDegree=3 := by unfold modulus8; compute_degree <;> norm_num
theorem modulus16_degree : modulus16.natDegree=4 := by unfold modulus16; compute_degree <;> norm_num

theorem modulus8_irreducible : Irreducible modulus8 := by
  apply irreducible_of_degree_le_three_of_not_isRoot
  · rw [modulus8_degree]
    decide
  · intro x
    fin_cases x <;> norm_num [modulus8,Polynomial.IsRoot] <;> decide

private theorem degree_two_form (q : Polynomial k) (degree : q.natDegree≤2) :
    q=C (q.coeff 0)+C (q.coeff 1)*X+C (q.coeff 2)*X^2 := by
  ext i
  by_cases small : i≤2
  · interval_cases i <;> simp
  · have zero : q.coeff i=0 := coeff_eq_zero_of_natDegree_lt (by omega)
    rw [zero]
    have h0 : i≠0 := by omega
    have h1 : i≠1 := by omega
    have h2 : i≠2 := by omega
    simp [coeff_add,coeff_C_mul,coeff_X_pow,coeff_X,coeff_C,h0,h1,h2,eq_comm]

private theorem two_polynomial : (2 : Polynomial k)=0 := by
  have two : (2 : k)=0 := by decide
  calc (2 : Polynomial k)=Polynomial.C (2 : k) := (map_ofNat Polynomial.C 2).symm
       _=0 := by rw [two,map_zero]

private theorem unit_of_bezout {q : Polynomial k} (divides : q∣modulus16)
    (a b : Polynomial k) (identity : a*q+b*modulus16=1) : IsUnit q := by
  apply isUnit_of_dvd_one
  rw [←identity]
  exact dvd_add (dvd_mul_left _ _) (dvd_mul_of_dvd_right divides _)

theorem modulus16_irreducible : Irreducible modulus16 := by
  apply (modulus16_monic.irreducible_iff_degree_lt ?_).2
  · intro q degree divides
    have bound : q.natDegree≤2 := by
      apply natDegree_le_of_degree_le
      simpa only [modulus16_degree] using degree
    have form := degree_two_form q bound
    have bits (c : k) : c=0 ∨ c=1 := by fin_cases c <;> simp
    rcases bits (q.coeff 0) with h0 | h0 <;>
      rcases bits (q.coeff 1) with h1 | h1 <;>
      rcases bits (q.coeff 2) with h2 | h2
    all_goals rw [form,h0,h1,h2] at divides ⊢
    all_goals simp only [map_zero,map_one,zero_mul,zero_add,add_zero,one_mul] at divides ⊢
    · exact False.elim (modulus16_monic.ne_zero (by simpa only [zero_dvd_iff] using divides))
    · apply unit_of_bezout divides (X^3+X^2+1) (X+1)
      unfold modulus16
      ring_nf
      simp [two_polynomial]
    · apply unit_of_bezout divides (X^3+1) 1
      unfold modulus16
      ring_nf
      simp [two_polynomial]
    · apply unit_of_bezout divides (X^2+X+1) 1
      unfold modulus16
      ring_nf
      simp [two_polynomial]
    · exact isUnit_one
    · apply unit_of_bezout divides (X^3+X+1) X
      unfold modulus16
      ring_nf
      simp [two_polynomial]
    · apply unit_of_bezout divides (X^3+X^2+X) 1
      unfold modulus16
      ring_nf
      simp [two_polynomial]
    · apply unit_of_bezout divides (X^2+X) 1
      unfold modulus16
      ring_nf
      simp [two_polynomial]
  · intro equal
    have degree := congrArg Polynomial.natDegree equal
    rw [modulus16_degree,natDegree_one] at degree
    contradiction

instance : Fact (Irreducible modulus8) := ⟨modulus8_irreducible⟩
instance : Fact (Irreducible modulus16) := ⟨modulus16_irreducible⟩
abbrev K8 := AdjoinRoot modulus8
abbrev K16 := AdjoinRoot modulus16
def alpha8 : K8 := AdjoinRoot.root modulus8
def alpha16 : K16 := AdjoinRoot.root modulus16

theorem alpha8_relation : alpha8^3+alpha8+1=0 := by
  have relation := AdjoinRoot.eval₂_root modulus8
  simpa [modulus8,alpha8] using relation

theorem alpha16_relation : alpha16^4+alpha16+1=0 := by
  have relation := AdjoinRoot.eval₂_root modulus16
  simpa [modulus16,alpha16] using relation

def basis8 : Basis (Fin 3) k K8 := (AdjoinRoot.powerBasis' modulus8_monic).basis.reindex
  (Fin.castOrderIso modulus8_degree).toEquiv
def basis16 : Basis (Fin 4) k K16 := (AdjoinRoot.powerBasis' modulus16_monic).basis.reindex
  (Fin.castOrderIso modulus16_degree).toEquiv

theorem basis8_eq_pow (i : Fin 3) : basis8 i=alpha8^i.val := by
  rw [basis8,Basis.reindex_apply]
  rw [(AdjoinRoot.powerBasis' modulus8_monic).basis_eq_pow]
  rfl

theorem basis16_eq_pow (i : Fin 4) : basis16 i=alpha16^i.val := by
  rw [basis16,Basis.reindex_apply]
  rw [(AdjoinRoot.powerBasis' modulus16_monic).basis_eq_pow]
  rfl

instance : Fintype K8 := Fintype.ofEquiv (Fin 3 → k) basis8.equivFun.symm.toEquiv
instance : Fintype K16 := Fintype.ofEquiv (Fin 4 → k) basis16.equivFun.symm.toEquiv
instance : CharP K8 2 := charP_of_injective_algebraMap' k K8 2
instance : CharP K16 2 := charP_of_injective_algebraMap' k K16 2

theorem card8 : Fintype.card K8=8 := by
  calc Fintype.card K8=Fintype.card (Fin 3 → k) := Fintype.card_congr basis8.equivFun.toEquiv
       _=8 := by decide

theorem card16 : Fintype.card K16=16 := by
  calc Fintype.card K16=Fintype.card (Fin 4 → k) := Fintype.card_congr basis16.equivFun.toEquiv
       _=16 := by decide

end ICFallDegrees.BinaryExtensionFields
