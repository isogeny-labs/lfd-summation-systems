import ICFallDegrees.GeneralQuadraticClassification
import Mathlib.Algebra.CharP.Quotient
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Algebra.MvPolynomial.CommRing
import Mathlib.RingTheory.MvPolynomial.Basic

/-! Truncated finite-field FFD in B=K[z]/(z_i^p), including both Koszul
and Frobenius rows among trivial syzygies. Homogeneous polynomial
representatives specify the grading; no affine degree-fall proxy is used. -/
noncomputable section
set_option maxHeartbeats 1000000
namespace ICFallDegrees.TruncatedFirstFall
open MvPolynomial HomogeneousFirstFall PowerIdeals LowDegreeCoefficients GeneralQuadraticTop
open GeneralQuadraticClassification
open scoped BigOperators
variable {r n : ℕ} {K : Type*} [CommRing K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K
abbrev B (r : ℕ) (K : Type*) [CommRing K] (p : ℕ) := P r K ⧸ powerIdeal p

def image (p : ℕ) : P r K →+* B r K p :=
  Ideal.Quotient.mk (powerIdeal p : Ideal (P r K))

def koszulRow (p : ℕ) (q : Fin n → P r K) (i j : Fin n) : Fin n → B r K p :=
  fun k => (if k=i then image p (q j) else 0)-(if k=j then image p (q i) else 0)

def frobeniusRow (p : ℕ) (q : Fin n → P r K) (i : Fin n) : Fin n → B r K p :=
  fun j => if j=i then (image p (q i))^(p-1) else 0

def trivialSpan (p : ℕ) (q : Fin n → P r K) : Submodule (B r K p) (Fin n → B r K p) :=
  Submodule.span _ (Set.range (fun ij : Fin n × Fin n => koszulRow p q ij.1 ij.2) ∪
    Set.range (frobeniusRow p q))

def IsSyzygy (p : ℕ) (q : Fin n → P r K) (degrees : Fin n → ℕ)
    (D : ℕ) (h : Fin n → P r K) : Prop :=
  (∀ i, h i=0 ∨ degrees i≤D ∧ (h i).IsHomogeneous (D-degrees i)) ∧
  ∑ i, image p (h i)*image p (q i)=0

def HasFall (p : ℕ) (q : Fin n → P r K) (degrees : Fin n → ℕ) (D : ℕ) : Prop :=
  ∃ h, IsSyzygy p q degrees D h ∧ (fun i => image p (h i))∉trivialSpan p q

def IsFirstFall (p : ℕ) (f : Fin n → P r K) (D : ℕ) : Prop :=
  HasFall p (topForms f) (fun i => (f i).totalDegree) D ∧
    ∀ d<D, ¬HasFall p (topForms f) (fun i => (f i).totalDegree) d

/-- Truncation preserves the coefficient characteristic: no nonzero
constant belongs to a positive-degree pure-power ideal. -/
theorem quotient_charP (p : ℕ) [Fact p.Prime] [CharP K p] : CharP (B r K p) p := by
  apply CharP.quotient' p (powerIdeal p : Ideal (P r K))
  intro c member
  apply degree_lt_eq_zero p member
  have positive := (Fact.out : p.Prime).pos
  have cast : (c : P r K)=C (c : K) := (map_natCast (C : K →+* P r K) c).symm
  rw [cast,totalDegree_C]
  exact positive

/-- Freshman's dream in the actual quotient, proved for all polynomials. -/
theorem image_pow (p : ℕ) [Fact p.Prime] [CharP K p] (f : P r K) :
    (image p f)^p=image p (C ((constantCoeff f)^p)) := by
  letI := quotient_charP (r:=r) (K:=K) p
  have xp (i : Fin r) : (image p (X i : P r K))^p=0 := by
    rw [←map_pow]
    exact Ideal.Quotient.eq_zero_iff_mem.mpr (Ideal.subset_span ⟨i,rfl⟩)
  induction f using MvPolynomial.induction_on with
  | C c => simp only [constantCoeff_C,map_pow]
  | add f g hf hg =>
      simp only [map_add,add_pow_char,hf,hg]
  | mul_X f i hf =>
      simp only [map_mul,mul_pow,xp,mul_zero,constantCoeff_X,zero_pow
        (Fact.out : p.Prime).ne_zero,map_zero]

/-- The quotient-induced rows really are syzygies, for every homogeneous
quadratic input; their removal is part of the definition, not an assumption. -/
theorem frobenius_row_syzygy (p : ℕ) [Fact p.Prime] [CharP K p]
    (q : Fin n → P r K) (constant_zero : ∀ i, constantCoeff (q i)=0) (i : Fin n) :
    ∑ j, frobeniusRow p q i j*image p (q j)=0 := by
  classical
  simp only [frobeniusRow,ite_mul,zero_mul,Finset.sum_ite_eq',Finset.mem_univ,if_true]
  rw [←pow_succ]
  have positive := (Fact.out : p.Prime).pos
  rw [Nat.sub_add_cancel (by omega : 1≤p),image_pow,constant_zero,zero_pow
    (Fact.out : p.Prime).ne_zero,map_zero,map_zero]

/-- Low degree quotient syzygies lift exactly, because the truncation ideal
contains no nonzero polynomial below its generator degree. -/
theorem lift_syzygy (p D : ℕ) (below : D<p) (q : Fin n → P r K)
    (degrees : Fin n → ℕ) (homogeneous : ∀ i, (q i).IsHomogeneous (degrees i))
    {h : Fin n → P r K} (relation : IsSyzygy p q degrees D h) :
    IsHomogeneousSyzygy q degrees D h := by
  refine ⟨relation.1,degree_lt_eq_zero p ?_ ?_⟩
  · apply Ideal.Quotient.eq_zero_iff_mem.mp
    change image p (∑ i, h i*q i)=0
    simpa only [image,map_sum,map_mul] using relation.2
  · apply lt_of_le_of_lt ?_ below
    apply totalDegree_finsetSum_le
    intro i _
    rcases relation.1 i with zero | degree
    · simp [zero]
    · exact (totalDegree_mul _ _).trans ((Nat.add_le_add degree.2.totalDegree_le
        (homogeneous i).totalDegree_le).trans (by omega))

theorem ordinary_to_quotient (p : ℕ) (q : Fin n → P r K) (degrees : Fin n → ℕ)
    (D : ℕ) {h : Fin n → P r K} (relation : IsHomogeneousSyzygy q degrees D h) :
    IsSyzygy p q degrees D h := by
  refine ⟨relation.1,?_⟩
  have equation := congrArg (image p) relation.2
  simpa only [image,map_sum,map_mul,map_zero] using equation

/-- Coordinate obstruction valid against arbitrary B-valued combinations
of both kinds of trivial rows. -/
theorem trivial_coordinate_mem (p : ℕ) (q : Fin n → P r K) (j : Fin n)
    (I : Ideal (B r K p)) (other : ∀ i, i≠j → image p (q i)∈I)
    (powers : (image p (q j))^(p-1)∈I) {h : Fin n → B r K p}
    (member : h∈trivialSpan p q) : h j∈I := by
  change h∈Submodule.span _ _ at member
  induction member using Submodule.span_induction with
  | mem v hv =>
      rcases hv with koszul | frobenius
      · obtain ⟨⟨a,b⟩,rfl⟩ := koszul
        by_cases ha : a=j
        · subst a
          by_cases hb : b=j
          · subst b; simp [koszulRow]
          · simpa [koszulRow,hb,Ne.symm hb] using other b hb
        · by_cases hb : b=j
          · subst b; simpa [koszulRow,ha,Ne.symm ha] using I.neg_mem (other a ha)
          · simp [koszulRow,ha,hb,Ne.symm ha,Ne.symm hb]
      · obtain ⟨i,rfl⟩ := frobenius
        by_cases equal : j=i
        · subst i; simpa [frobeniusRow] using powers
        · simp [frobeniusRow,equal]
  | zero => exact I.zero_mem
  | add v w _ _ hv hw => exact I.add_mem hv hw
  | smul a v _ hv => exact I.mul_mem_left a hv

/-- The active bilinear form has cube zero modulo the square ideal,
independently of characteristic or the determinant branch. -/
theorem top_cube_mem_squares (a : Fin 4 → Fin r) (α β : K) :
    (top a α β)^3∈powerIdeal 2 := by
  let u : P r K := X (a 0)
  let v : P r K := X (a 1)
  let A : P r K := C α*X (a 2)+C β*X (a 3)
  let Q : P r K := C β*X (a 2)+X (a 3)
  have identity : (top a α β)^3=
      (u*A^3+3*v*A^2*Q)*u^2+(3*u*A*Q^2+v*Q^3)*v^2 := by
    dsimp [u,v,A,Q,top]
    ring
  rw [identity]
  exact (powerIdeal 2).add_mem
    ((powerIdeal 2).mul_mem_left _ (Ideal.subset_span ⟨a 0,rfl⟩))
    ((powerIdeal 2).mul_mem_left _ (Ideal.subset_span ⟨a 1,rfl⟩))

theorem top_power_mem_squares (a : Fin 4 → Fin r) (α β : K) (e : ℕ) (large : 3≤e) :
    (top a α β)^e∈powerIdeal 2 := by
  have identity : (top a α β)^e=(top a α β)^(e-3)*(top a α β)^3 := by
    rw [←pow_add]; congr 1; omega
  rw [identity]
  exact (powerIdeal 2).mul_mem_left _ (top_cube_mem_squares a α β)

theorem last_of_trivial_in_squares (p : ℕ) (large : 5≤p) (a : Fin 4 → Fin r) (α β : K)
    {h : Fin (r+1) → P r K} (member : (fun i => image p (h i))∈trivialSpan p (q a α β)) :
    h (Fin.last r)∈powerIdeal 2 := by
  let I : Ideal (B r K p) := (powerIdeal 2).map (Ideal.Quotient.mk (powerIdeal p))
  have other (i : Fin (r+1)) (different : i≠Fin.last r) : image p (q a α β i)∈I :=
    Ideal.mem_map_of_mem _ (other_in_squares a α β i different)
  have powers : (image p (q a α β (Fin.last r)))^(p-1)∈I := by
    rw [←map_pow]
    apply Ideal.mem_map_of_mem
    simpa only [q,Fin.lastCases_last] using top_power_mem_squares a α β (p-1) (by omega)
  have inside := trivial_coordinate_mem p (q a α β) (Fin.last r) I other powers member
  exact (Ideal.mem_quotient_iff_mem (powerIdeal_antitone (by omega : 2≤p))).mp inside

variable [Nontrivial K]

theorem linear_nontrivial (p : ℕ) (large : 5≤p) (a : Fin 4 → Fin r)
    (injective : Function.Injective a) (β : K) :
    (fun i => image p (linearWitness a β i))∉trivialSpan p (q a (β^2) β) := by
  intro member
  have inside := last_of_trivial_in_squares p large a (β^2) β member
  have zero := squarefree_coefficient_zero inside (Finsupp.single (a 1) 1)
    (fun i => by simp only [Finsupp.single_apply]; split_ifs <;> omega)
  have unit : coeff (Finsupp.single (a 1) 1) (linearWitness a β (Fin.last r))=1 := by
    simp [linearWitness,coeff_single_X,injective.eq_iff]
  exact one_ne_zero (unit.symm.trans zero)

theorem quadratic_nontrivial (p : ℕ) (large : 5≤p) (a : Fin 4 → Fin r)
    (injective : Function.Injective a) (α β : K) :
    (fun i => image p (quadraticWitness a α β i))∉trivialSpan p (q a α β) := by
  intro member
  have inside := last_of_trivial_in_squares p large a α β member
  have zero := squarefree_coefficient_zero inside (pair (a 0) (a 1))
    (pair_squarefree _ _ (injective.ne (by decide)))
  have unit : coeff (pair (a 0) (a 1)) (quadraticWitness a α β (Fin.last r))=1 := by
    simp only [quadraticWitness,Fin.lastCases_last]
    rw [pair_polynomial,coeff_monomial,if_pos rfl]
  exact one_ne_zero (unit.symm.trans zero)

/-- Exact truncated value in the determinant-zero branch, with all
Frobenius trivial rows included. -/
theorem first_fall_three (p : ℕ) (large : 5≤p) (a : Fin 4 → Fin r)
    (injective : Function.Injective a) (β : K) :
    HasFall p (q a (β^2) β) (fun _ => 2) 3 ∧
      ∀ D<3, ¬HasFall p (q a (β^2) β) (fun _ => 2) D := by
  refine ⟨⟨linearWitness a β,ordinary_to_quotient _ _ _ _ (linear_witness_syzygy a injective β),
    linear_nontrivial p large a injective β⟩,?_⟩
  intro D small ⟨h,relation,nontrivial⟩
  have ordinary := lift_syzygy p D (by omega) _ _ (homogeneous_q a (β^2) β) relation
  have zero := low_syzygy_zero _ (pivot a) (coefficient_pivots a injective (β^2) β) 2 D (by omega) ordinary
  apply nontrivial
  have images : (fun i => image p (h i))=0 := by funext i; simp only [zero,Pi.zero_apply,map_zero]
  rw [images]
  exact (trivialSpan _ _).zero_mem

section Field
variable {F : Type*} [Field F]

theorem first_fall_four (p : ℕ) (large : 5≤p) (a : Fin 4 → Fin r)
    (injective : Function.Injective a) (α β : F) (nondegenerate : α-β^2≠0) :
    HasFall p (q a α β) (fun _ => 2) 4 ∧ ∀ D<4, ¬HasFall p (q a α β) (fun _ => 2) D := by
  refine ⟨⟨quadraticWitness a α β,ordinary_to_quotient _ _ _ _ (quadratic_witness_syzygy a injective α β),
    quadratic_nontrivial p large a injective α β⟩,?_⟩
  intro D small ⟨h,relation,nontrivial⟩
  have ordinary := lift_syzygy p D (by omega) _ _ (homogeneous_q a α β) relation
  have zero : h=0 := by
    by_cases low : D≤2
    · exact low_syzygy_zero _ (pivot a) (coefficient_pivots a injective α β) 2 D low ordinary
    · have degree : D=3 := by omega
      subst D
      exact degree_three_syzygy_zero a injective α β nondegenerate ordinary
  apply nontrivial
  have images : (fun i => image p (h i))=0 := by funext i; simp only [zero,Pi.zero_apply,map_zero]
  rw [images]
  exact (trivialSpan _ _).zero_mem

end Field
end ICFallDegrees.TruncatedFirstFall
