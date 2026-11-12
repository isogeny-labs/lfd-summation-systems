import ICFallDegrees.FilteredSpan
import Mathlib.RingTheory.MvPolynomial.Homogeneous
import Mathlib.Algebra.MvPolynomial.PDeriv

/-! Ordinary-polynomial homogeneous syzygies modulo the actual Koszul
submodule. This is separate from affine strict falls and from truncated
finite-field quotient conventions. -/
noncomputable section

namespace ICFallDegrees.HomogeneousFirstFall
open MvPolynomial
open scoped BigOperators

abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K
variable {r n : ℕ} {K : Type*} [CommRing K]

def koszulRow (q : Fin n → P r K) (i j : Fin n) : Fin n → P r K :=
  fun k => (if k=i then q j else 0)-(if k=j then q i else 0)

def koszulSpan (q : Fin n → P r K) : Submodule (P r K) (Fin n → P r K) :=
  Submodule.span (P r K) (Set.range fun ij : Fin n × Fin n => koszulRow q ij.1 ij.2)

/-- Zero components are allowed in every degree. Every nonzero component has
exact weighted degree D, without saturating negative degrees to zero. -/
def IsHomogeneousSyzygy (q : Fin n → P r K) (degrees : Fin n → ℕ)
    (D : ℕ) (h : Fin n → P r K) : Prop :=
  (∀ i, h i=0 ∨ degrees i≤D ∧ (h i).IsHomogeneous (D-degrees i)) ∧
  ∑ i, h i*q i=0

def HasFall (q : Fin n → P r K) (degrees : Fin n → ℕ) (D : ℕ) : Prop :=
  ∃ h, IsHomogeneousSyzygy q degrees D h ∧ h∉koszulSpan q

def topForms (f : Fin n → P r K) : Fin n → P r K :=
  fun i => homogeneousComponent (f i).totalDegree (f i)

def HasKoszulFall (f : Fin n → P r K) (D : ℕ) : Prop :=
  HasFall (topForms f) (fun i => (f i).totalDegree) D

def IsKoszulFirstFall (f : Fin n → P r K) (D : ℕ) : Prop :=
  HasKoszulFall f D ∧ ∀ d<D, ¬HasKoszulFall f d

/-- Every Koszul combination is a syzygy; the span allows arbitrary
polynomial coefficients, not just coefficients of the witness degree. -/
theorem koszulSpan_syzygy (q : Fin n → P r K) {h : Fin n → P r K}
    (member : h∈koszulSpan q) : ∑ i, h i*q i=0 := by
  change h∈Submodule.span _ _ at member
  induction member using Submodule.span_induction with
  | mem p hp =>
      obtain ⟨⟨i,j⟩,rfl⟩ := hp
      simp only [koszulRow,sub_mul,Finset.sum_sub_distrib,ite_mul,zero_mul,
        Finset.sum_ite_eq',Finset.mem_univ,if_true]
      ring
  | zero => simp
  | add p q _ _ hp hq =>
      simp only [Pi.add_apply,add_mul,Finset.sum_add_distrib,hp,hq,add_zero]
  | smul a p _ hp =>
      simp only [Pi.smul_apply,smul_eq_mul,mul_assoc,←Finset.mul_sum,hp,mul_zero]

/-- Vanishing of values and first derivatives at the origin is preserved by
arbitrary polynomial multiples and sums, hence by all Koszul combinations. -/
theorem koszulSpan_jet_zero (q : Fin n → P r K) (v : Fin r)
    (values : ∀ i, constantCoeff (q i)=0)
    (derivatives : ∀ i, constantCoeff (pderiv v (q i))=0)
    {h : Fin n → P r K} (member : h∈koszulSpan q) :
    ∀ i, constantCoeff (h i)=0 ∧
      constantCoeff (pderiv v (h i))=0 := by
  change h∈Submodule.span _ _ at member
  induction member using Submodule.span_induction with
  | mem p hp =>
      obtain ⟨⟨i,j⟩,rfl⟩ := hp
      intro k
      simp only [koszulRow]
      split_ifs <;> simp [values,derivatives]
  | zero => simp
  | add p q _ _ hp hq =>
      intro i
      simp [Pi.add_apply,map_add,(hp i).1,(hp i).2,(hq i).1,(hq i).2]
  | smul a p _ hp =>
      intro i
      simp [Pi.smul_apply,smul_eq_mul,pderiv_mul,(hp i).1,(hp i).2]

variable [Nontrivial K]

/-- Independent coefficient pivots rule out every homogeneous relation of
weighted degree at most the common input degree. -/
theorem low_syzygy_zero (q : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ)
    (rref : ∀ i j, coeff (pivot i) (q j)=if i=j then 1 else 0)
    (δ D : ℕ) (budget : D≤δ) {h : Fin n → P r K}
    (relation : IsHomogeneousSyzygy q (fun _ => δ) D h) : h=0 := by
  have constants : ∀ i, h i=C (coeff 0 (h i)) := by
    intro i
    rcases relation.1 i with zero | degree
    · simp [zero]
    · apply totalDegree_eq_zero_iff_eq_C.mp
      have bound := degree.2.totalDegree_le
      have subtract : D-δ=0 := Nat.sub_eq_zero_of_le budget
      rw [subtract] at bound
      exact Nat.eq_zero_of_le_zero bound
  have zero : FilteredSpan.rowMap q (fun i => coeff 0 (h i))=0 := by
    simpa only [FilteredSpan.rowMap,Fintype.linearCombination_apply,smul_eq_C_mul,
      ←constants] using relation.2
  funext i
  have scalar_zero : coeff 0 (h i)=0 := by
    have equal := congrArg (coeff (pivot i)) zero
    rw [FilteredSpan.coefficient_rowMap q pivot rref] at equal
    simpa using equal
  change h i=0
  rw [constants i,scalar_zero,map_zero]

theorem no_low_fall (q : Fin n → P r K) (pivot : Fin n → Fin r →₀ ℕ)
    (rref : ∀ i j, coeff (pivot i) (q j)=if i=j then 1 else 0)
    (δ D : ℕ) (budget : D≤δ) : ¬HasFall q (fun _ => δ) D := by
  rintro ⟨h,relation,nontrivial⟩
  have zero := low_syzygy_zero q pivot rref δ D budget relation
  exact nontrivial (zero.symm ▸ (koszulSpan q).zero_mem)

/-- Identify both the literal input degree and its highest homogeneous part
from a nonzero homogeneous summand and a strictly lower-degree tail. -/
theorem degree_and_top_of_add {q l : P r K} {d : ℕ}
    (homogeneous : q.IsHomogeneous d) (nonzero : q≠0) (lower : l.totalDegree<d) :
    (q+l).totalDegree=d ∧ homogeneousComponent (q+l).totalDegree (q+l)=q := by
  have degree : (q+l).totalDegree=d := by
    rw [totalDegree_add_eq_left_of_totalDegree_lt]
    · exact homogeneous.totalDegree nonzero
    · simpa only [homogeneous.totalDegree nonzero] using lower
  refine ⟨degree,?_⟩
  rw [degree,map_add,homogeneousComponent_of_mem homogeneous,
    if_pos rfl,homogeneousComponent_eq_zero _ _ lower,add_zero]

end ICFallDegrees.HomogeneousFirstFall
