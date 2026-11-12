import ICFallDegrees.BooleanSharp
import Mathlib.Algebra.Field.ZMod
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

/-! The three-bit S3 constraint has exact LFD4. This is the obstruction
used by affine specializations of arbitrarily long split chains. -/
noncomputable section
namespace ICFallDegrees.BinaryTriple
open MvPolynomial FiniteQuotient GridQuotient
open scoped BigOperators
abbrev P := MvPolynomial (Fin 3) (ZMod 2)
def xy : Fin 3 →₀ ℕ := Finsupp.single 0 1+Finsupp.single 1 1
def xz : Fin 3 →₀ ℕ := Finsupp.single 0 1+Finsupp.single 2 1
def yz : Fin 3 →₀ ℕ := Finsupp.single 1 1+Finsupp.single 2 1
def xyz : Fin 3 →₀ ℕ := xy+Finsupp.single 2 1
def h : P := 1+monomial xy 1+monomial xz 1+monomial yz 1+monomial xyz 1
def parity : P := X 0+X 1+X 2
def sets : Fin 3 → Finset (ZMod 2) := fun _ => {0,1}
def G : Set P := gridInputs sets
def inputs : Set P := G∪{h}

theorem h_product : h=1+X 0*X 1+X 0*X 2+X 1*X 2+X 0*X 1*X 2 := by
  simp [h,xy,xz,yz,xyz,X,monomial_mul]

private theorem monomial_reduced (e : Fin 3 →₀ ℕ) (squarefree : ∀ i, e i≤1) :
    (monomial e (1 : ZMod 2) : P)∈reducedSpace G := by
  intro d hd g hg
  have same : d=e := by simpa [support_monomial] using hd
  subst d
  obtain ⟨i,rfl⟩ := hg
  rw [leading_gridPolynomial]
  intro divides
  have high := divides i
  simp [sets,Finsupp.single_apply] at high
  have low := squarefree i
  omega

theorem h_reduced : h∈reducedSpace G := by
  have ex : ∀ i, xy i≤1 := by intro i; fin_cases i <;> simp [xy,Finsupp.single_apply]
  have ez : ∀ i, xz i≤1 := by intro i; fin_cases i <;> simp [xz,Finsupp.single_apply]
  have ey : ∀ i, yz i≤1 := by intro i; fin_cases i <;> simp [yz,Finsupp.single_apply]
  have et : ∀ i, xyz i≤1 := by intro i; fin_cases i <;> simp [xyz,xy,Finsupp.single_apply]
  have one : (1 : P)∈reducedSpace G := by simpa using monomial_reduced 0 (by simp)
  have h1 := (reducedSpace G).add_mem one (monomial_reduced xy ex)
  have h2 := (reducedSpace G).add_mem h1 (monomial_reduced xz ez)
  have h3 := (reducedSpace G).add_mem h2 (monomial_reduced yz ey)
  exact (reducedSpace G).add_mem h3 (monomial_reduced xyz et)

private theorem exponent_degrees : xy.degree=2 ∧ xz.degree=2 ∧ yz.degree=2 ∧ xyz.degree=3 := by
  simp [xy,xz,yz,xyz,Finsupp.degree_add,Finsupp.degree_single]

private theorem top_coefficient : coeff xyz h=1 := by
  have nonzero : xyz≠0 := by intro eq; have := congrArg (fun e => e (2 : Fin 3)) eq; simp [xyz,xy,Finsupp.single_apply] at this
  have nxy : xyz≠xy := by intro eq; have := congrArg (fun e => e (2 : Fin 3)) eq; simp [xyz,xy,Finsupp.single_apply] at this
  have nxz : xyz≠xz := by intro eq; have := congrArg (fun e => e (1 : Fin 3)) eq; simp [xyz,xy,xz,Finsupp.single_apply] at this
  have nyz : xyz≠yz := by intro eq; have := congrArg (fun e => e (0 : Fin 3)) eq; simp [xyz,xy,yz,Finsupp.single_apply] at this
  simp [h,coeff_add,coeff_monomial,coeff_one,nonzero,Ne.symm nonzero,nxy,Ne.symm nxy,nxz,Ne.symm nxz,nyz,Ne.symm nyz]

theorem h_degree : h.totalDegree=3 := by
  have ed := exponent_degrees
  have mono (e : Fin 3 →₀ ℕ) (degree : e.degree≤3) : (monomial e (1 : ZMod 2)).totalDegree≤3 :=
    (totalDegree_monomial_le e 1).trans degree
  have add {p q : P} (hp : p.totalDegree≤3) (hq : q.totalDegree≤3) : (p+q).totalDegree≤3 :=
    (totalDegree_add _ _).trans (max_le hp hq)
  have upper : h.totalDegree≤3 := add (add (add (add (by simp)
    (mono xy (by omega))) (mono xz (by omega))) (mono yz (by omega))) (mono xyz (by omega))
  have member : xyz∈h.support := mem_support_iff.mpr (by rw [top_coefficient]; exact one_ne_zero)
  have lower := le_totalDegree member
  change xyz.degree≤h.totalDegree at lower
  rw [ed.2.2.2] at lower
  omega

private theorem grid_mem_of_vanishing (p : P)
    (vanish : ∀ point : Fin 3 → ZMod 2, eval point p=0) : p∈Ideal.span G := by
  obtain ⟨c,_,identity⟩ := combinatorial_nullstellensatz_exists_linearCombination sets
    (by intro i; simp [sets]) p (fun point _ => vanish point)
  rw [identity,Finsupp.linearCombination_apply,Finsupp.sum]
  apply (Ideal.span G).sum_mem
  intro i _
  simp only [smul_eq_mul]
  exact (Ideal.span G).mul_mem_left _ (Ideal.subset_span ⟨i,rfl⟩)

theorem parity_mem_ideal : parity∈Ideal.span inputs := by
  have error : parity-parity*h∈Ideal.span G := by
    apply grid_mem_of_vanishing
    intro point
    have bit : ∀ z : ZMod 2, z=0 ∨ z=1 := by intro z; fin_cases z <;> simp
    rcases bit (point 0) with h0 | h0 <;>
      rcases bit (point 1) with h1 | h1 <;>
      rcases bit (point 2) with h2 | h2 <;>
      norm_num [parity,h_product,h0,h1,h2] <;> decide
  have included : Ideal.span G≤Ideal.span inputs := Ideal.span_mono (fun _ hp => Or.inl hp)
  have multiple : parity*h∈Ideal.span inputs := (Ideal.span inputs).mul_mem_left _ (Ideal.subset_span (Or.inr rfl))
  have result := (Ideal.span inputs).add_mem (included error) multiple
  simpa using result

theorem parity_not_at_three : ¬InClosure inputs 3 parity := by
  have nonempty : ∀ i, (sets i).Nonempty := by intro i; simp [sets]
  have zero := reduced_grid_ideal_trivial sets nonempty
  have closed : closureSpace inputs 3≤quotientSpan (Ideal.span G) h := by
    apply closure_le_quotient_span inputs G (gridInputs_nonzero sets) zero h h_reduced 3 h_degree
    intro p hp _
    rcases hp with hp | rfl
    · exact ⟨0,by simpa using (Ideal.subset_span hp : p∈Ideal.span G)⟩
    · exact ⟨1,by simp⟩
  intro member
  have small : parity.totalDegree<3 := by
    have : parity.totalDegree≤1 := by
      unfold parity
      exact (totalDegree_add _ _).trans (max_le ((totalDegree_add _ _).trans (by simp)) (by simp))
    omega
  have in_grid := low_quotient_span_mem_ideal G (gridInputs_nonzero sets) zero h h_reduced 3 h_degree
    (closed member) small
  let point : Fin 3 → ZMod 2 := fun i => if i=0 then 1 else 0
  have containing : Ideal.span G≤RingHom.ker (eval point) := by
    apply Ideal.span_le.mpr
    rintro p ⟨i,rfl⟩
    change eval point (gridPolynomial (sets i) i)=0
    simp [gridPolynomial,sets]
    fin_cases i <;> norm_num [point]
  have bad : eval point parity=0 := containing in_grid
  have value : eval point parity=1 := by norm_num [parity,point] <;> decide
  rw [value] at bad
  exact one_ne_zero bad

theorem upper_bound : HasLFDBound inputs 4 := by
  have bound := grid_lfd_bound sets (by intro i; simp [sets]) inputs
    (fun _ hp => Or.inl hp) 3 (by
      intro p hp
      rcases hp with ⟨i,rfl⟩ | rfl
      · have id : gridPolynomial (sets i) i=(X i : P)^2-X i := by
          simp [gridPolynomial,sets]
          ring
        change (gridPolynomial (sets i) i).totalDegree≤3
        rw [id]
        exact (totalDegree_sub_le _ _).trans (by simp)
      · exact h_degree.le)
  simpa [sets,Fin.sum_univ_three] using bound

theorem exact_lfd_four : IsLastFallDegree inputs 4 := by
  refine ⟨upper_bound,?_⟩
  intro d hd bound
  have derived := bound parity parity_mem_ideal
  have degree : parity.totalDegree≤1 := by
    unfold parity
    exact (totalDegree_add _ _).trans (max_le ((totalDegree_add _ _).trans (by simp)) (by simp))
  exact parity_not_at_three (derived.mono (by omega))

end ICFallDegrees.BinaryTriple
