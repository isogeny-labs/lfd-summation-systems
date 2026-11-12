import ICFallDegrees.BinaryDescent
import ICFallDegrees.Closure
import Mathlib.Algebra.CharP.Frobenius
import Mathlib.Data.ZMod.Basic

/-! Arbitrary-dimensional ordinary polynomial identities for the inverse
coordinate bridge. Boolean corrections are explicitly charged at degree two
or four. These are actual polynomial identities, not identities of Boolean
functions with free reduction. -/
noncomputable section
namespace ICFallDegrees.BinaryInverseCoordinates
open MvPolynomial
open scoped BigOperators
variable {σ K : Type*} {n : ℕ} [CommRing K] [Nontrivial K] [CharP K 2]

def form (beta : Fin n → K) (v : Fin n → σ) : MvPolynomial σ K :=
  ∑ j, C (beta j)*X (v j)

def boolean (i : σ) : MvPolynomial σ K := X i^2-X i

def conjugate (beta : Fin n → K) (v : Fin n → σ) (i : ℕ) : MvPolynomial σ K :=
  form (fun j => beta j^(2^i)) v

theorem form_degree (beta : Fin n → K) (v : Fin n → σ) :
    (form beta v).totalDegree≤1 := by
  apply totalDegree_finsetSum_le
  intro j _
  exact (totalDegree_mul _ _).trans (by simp)

theorem boolean_degree (i : σ) : (boolean i : MvPolynomial σ K).totalDegree≤2 := by
  exact (ICFallDegrees.totalDegree_sub_le _ _).trans
    (max_le ((totalDegree_pow _ _).trans (by simp)) (by simp))

theorem square_form (beta : Fin n → K) (v : Fin n → σ) :
    (form beta v)^2-form (fun j => beta j^2) v=
      ∑ j, C (beta j^2)*boolean (v j) := by
  unfold form
  rw [CharTwo.sum_sq,←Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro j _
  simp only [mul_pow,←map_pow,boolean]
  ring

theorem square_form_at_two (F : Set (MvPolynomial σ K))
    (beta : Fin n → K) (v : Fin n → σ)
    (fields : ∀ j, boolean (v j)∈F) :
    InClosure F 2 ((form beta v)^2-form (fun j => beta j^2) v) := by
  rw [square_form]
  apply (closureSpace F 2).sum_mem
  intro j _
  simpa only [smul_eq_C_mul] using
    (InClosure.input (fields j) (boolean_degree (v j))).smul (beta j^2)

theorem conjugate_square_at_two (F : Set (MvPolynomial σ K))
    (beta : Fin n → K) (v : Fin n → σ)
    (fields : ∀ j, boolean (v j)∈F) (i : ℕ) :
    InClosure F 2 ((conjugate beta v i)^2-conjugate beta v (i+1)) := by
  have exponents : ∀ j, (beta j^(2^i))^2=beta j^(2^(i+1)) := by
    intro j
    rw [←pow_mul,pow_succ]
  have h := square_form_at_two F (fun j => beta j^(2^i)) v fields
  rw [funext exponents] at h
  exact h

theorem conjugate_cycle (beta : Fin n → K) (v : Fin n → σ)
    (fixed : ∀ j, beta j^(2^n)=beta j) : conjugate beta v n=conjugate beta v 0 := by
  simp only [conjugate,pow_zero,pow_one,fixed]

theorem target_zero_s3 (x y : MvPolynomial σ K) :
    (x*y)^2+1=(x*y+1)^2 := by
  rw [CharTwo.add_sq,one_pow]

theorem raw_inverse_correction (beta gamma : Fin n → K) (v w : Fin n → σ) :
    (form beta v*form gamma w)^2+1-
      (form (fun j => beta j^2) v*form (fun j => gamma j^2) w+1)=
    (form gamma w)^2*(∑ j, C (beta j^2)*boolean (v j))+
      form (fun j => beta j^2) v*(∑ j, C (gamma j^2)*boolean (w j)) := by
  rw [←square_form beta v,←square_form gamma w]
  ring

theorem raw_inverse_at_four (F : Set (MvPolynomial σ K))
    (beta gamma : Fin n → K) (v w : Fin n → σ)
    (fields_v : ∀ j, boolean (v j)∈F) (fields_w : ∀ j, boolean (w j)∈F)
    (raw : (form beta v*form gamma w)^2+1∈F) :
    InClosure F 4 (form (fun j => beta j^2) v*form (fun j => gamma j^2) w+1) := by
  have dv := form_degree beta v
  have dw := form_degree gamma w
  have dv' := form_degree (fun j => beta j^2) v
  have dw' := form_degree (fun j => gamma j^2) w
  have sqw : ((form gamma w)^2).totalDegree≤2 :=
    (totalDegree_pow _ _).trans (by omega)
  have first := (square_form_at_two F beta v fields_v).mono (by omega : 2≤4)
  have second := (square_form_at_two F gamma w fields_w).mono (by omega : 2≤4)
  have first_degree : ((form beta v)^2-form (fun j => beta j^2) v).totalDegree≤2 :=
    (ICFallDegrees.totalDegree_sub_le _ _).trans
      (max_le ((totalDegree_pow _ _).trans (by omega)) (dv'.trans (by omega)))
  have second_degree : ((form gamma w)^2-form (fun j => gamma j^2) w).totalDegree≤2 :=
    (ICFallDegrees.totalDegree_sub_le _ _).trans
      (max_le ((totalDegree_pow _ _).trans (by omega)) (dw'.trans (by omega)))
  have a := first.mul ((form gamma w)^2) (by omega)
  have b := second.mul (form (fun j => beta j^2) v) (by omega)
  have identity := raw_inverse_correction beta gamma v w
  rw [←square_form beta v,←square_form gamma w] at identity
  have budget : ((form beta v*form gamma w)^2+1).totalDegree≤4 := by
    apply (totalDegree_add _ _).trans
    apply max_le
    · have product_degree : (form beta v*form gamma w).totalDegree≤2 :=
        (totalDegree_mul _ _).trans (Nat.add_le_add dv dw)
      exact (totalDegree_pow _ _).trans (Nat.mul_le_mul_left 2 product_degree)
    · simp
  have input : InClosure F 4 ((form beta v*form gamma w)^2+1) := .input raw budget
  have correction := a.add b
  rw [←identity] at correction
  simpa only [sub_sub_cancel] using input.sub correction

section Coordinates
variable [Algebra (ZMod 2) K]

theorem component_reconstruction (basis : Basis (Fin n) (ZMod 2) K)
    (p : MvPolynomial σ K) :
    (∑ j, C (basis j)*map (algebraMap (ZMod 2) K) (BinaryDescent.component basis j p))=p := by
  classical
  ext e
  simp only [coeff_sum,coeff_C_mul,coeff_map]
  change (∑ j, basis j*algebraMap (ZMod 2) K (basis.repr (coeff e p) j))=coeff e p
  simpa only [Algebra.smul_def,mul_comm] using basis.sum_repr (coeff e p)

theorem prime_field_fixed (i : ℕ) (p : MvPolynomial σ (ZMod 2)) :
    map (iterateFrobenius K 2 i) (map (algebraMap (ZMod 2) K) p)=
      map (algebraMap (ZMod 2) K) p := by
  have identity : (iterateFrobenius K 2 i).comp (algebraMap (ZMod 2) K)=algebraMap (ZMod 2) K :=
    RingHom.ext_zmod _ _
  rw [map_map,identity]

theorem form_frobenius (beta : Fin n → K) (v : Fin n → σ) (i : ℕ) :
    map (iterateFrobenius K 2 i) (form beta v)=conjugate beta v i := by
  simp [form,conjugate,iterateFrobenius_def]

theorem conjugate_reconstruction (basis : Basis (Fin n) (ZMod 2) K)
    (p : MvPolynomial σ K) (i : ℕ) :
    (∑ j, C ((basis j)^(2^i))*map (algebraMap (ZMod 2) K) (BinaryDescent.component basis j p))=
      map (iterateFrobenius K 2 i) p := by
  have identity := congrArg (map (iterateFrobenius K 2 i)) (component_reconstruction basis p)
  simpa only [map_sum,map_mul,map_C,iterateFrobenius_def,prime_field_fixed] using identity

theorem inverse_conjugate_identity (basis : Basis (Fin n) (ZMod 2) K)
    (beta gamma : Fin n → K) (v w : Fin n → σ) (i : ℕ) :
    (∑ j, C ((basis j)^(2^i))*map (algebraMap (ZMod 2) K)
      (BinaryDescent.component basis j (form beta v*form gamma w+1)))=
      conjugate beta v i*conjugate gamma w i+1 := by
  rw [conjugate_reconstruction]
  simp only [map_add,map_mul,map_one,form_frobenius]

theorem mapped_component_degree (basis : Basis (Fin n) (ZMod 2) K)
    (p : MvPolynomial σ K) (j : Fin n) :
    (map (algebraMap (ZMod 2) K) (BinaryDescent.component basis j p)).totalDegree≤p.totalDegree := by
  classical
  apply Finset.sup_le
  intro e he
  apply le_totalDegree
  by_contra absent
  have zero : coeff e p=0 := not_mem_support_iff.mp absent
  have mapped_zero : coeff e (map (algebraMap (ZMod 2) K) (BinaryDescent.component basis j p))=0 := by
    rw [coeff_map]
    change algebraMap (ZMod 2) K (basis.coord j (coeff e p))=0
    rw [zero,map_zero,map_zero]
  exact (mem_support_iff.mp he) mapped_zero

theorem inverse_conjugate_at_two (F : Set (MvPolynomial σ K))
    (basis : Basis (Fin n) (ZMod 2) K)
    (beta gamma : Fin n → K) (v w : Fin n → σ)
    (coordinates : ∀ j, map (algebraMap (ZMod 2) K)
      (BinaryDescent.component basis j (form beta v*form gamma w+1))∈F) (i : ℕ) :
    InClosure F 2 (conjugate beta v i*conjugate gamma w i+1) := by
  rw [←inverse_conjugate_identity basis beta gamma v w i]
  apply (closureSpace F 2).sum_mem
  intro j _
  have degree : (form beta v*form gamma w+1).totalDegree≤2 :=
    (totalDegree_add _ _).trans
      (max_le ((totalDegree_mul _ _).trans (Nat.add_le_add (form_degree beta v) (form_degree gamma w))) (by simp))
  have input := InClosure.input (coordinates j) ((mapped_component_degree basis _ j).trans degree)
  simpa only [smul_eq_C_mul] using input.smul ((basis j)^(2^i))

end Coordinates

end ICFallDegrees.BinaryInverseCoordinates
