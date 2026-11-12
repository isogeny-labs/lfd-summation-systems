import ICFallDegrees.NormalFormBound
import Mathlib.RingTheory.MvPolynomial.Groebner
import Mathlib.RingTheory.MvPolynomial.MonomialOrder.DegLex
import Mathlib.RingTheory.Ideal.Quotient.Operations

noncomputable section

/-!
# Finite-standard-monomial last-fall bound

The normal-form and degree-bounded reduction hypotheses of the earlier
criterion are constructed here from graded-lex division and the mathematical
Groebner-basis property. No normal form, rank, or reduction budget is assumed.
-/

namespace ICFallDegrees.FiniteQuotient

open MvPolynomial MonomialOrder
open scoped MonomialOrder

variable {r : ℕ} {K : Type*} [Field K]

abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

/-- The leading-monomial characterization of a Groebner basis, with zero
generators removed. The inequality is coordinatewise divisibility. -/
def IsGroebner (G : Set (P r K)) : Prop :=
  (0 : P r K)∉G ∧ ∀ p ∈ Ideal.span G, p≠0 →
    ∃ g ∈ G, degLex.degree g≤degLex.degree p

def Standard (G : Set (P r K)) (e : Fin r →₀ ℕ) : Prop :=
  ∀ g ∈ G, ¬degLex.degree g≤e

def reducedSpace (G : Set (P r K)) : Submodule K (P r K) where
  carrier := fun p => ∀ e∈p.support, Standard G e
  zero_mem' := by intro e he; simp at he
  add_mem' := by
    intro p q hp hq e he
    rcases Finset.mem_union.mp (MvPolynomial.support_add he) with he | he
    · exact hp e he
    · exact hq e he
  smul_mem' := by
    intro a p hp e he
    exact hp e (MvPolynomial.support_smul he)

/-- The zero-remainder characterization, useful for explicit grid ideals. -/
def ReducedIdealTrivial (G : Set (P r K)) : Prop :=
  ∀ p∈Ideal.span G, p∈reducedSpace G → p=0

theorem totalDegree_mul_eq {p q : P r K} (hp : p≠0) (hq : q≠0) :
    (p*q).totalDegree=p.totalDegree+q.totalDegree := by
  rw [←degree_degLexDegree,degLex.degree_mul hp hq,Finsupp.degree_add,
    degree_degLexDegree,degree_degLexDegree]

/-- Division has the literal ordinary-degree multiplier budget. -/
theorem exists_reduced_remainder (G : Set (P r K)) (nonzero : (0 : P r K)∉G)
    (p : P r K) : ∃ z : P r K, z∈reducedSpace G ∧
      InClosure G p.totalDegree (p-z) := by
  classical
  have units : ∀ g : G, IsUnit (degLex.leadingCoeff (g : P r K)) := by
    intro g
    apply (MonomialOrder.isUnit_leadingCoeff (m:=degLex)).mpr
    intro zero
    exact nonzero (zero ▸ g.property)
  obtain ⟨multipliers,z,identity,bounds,reduced⟩ := degLex.div units p
  refine ⟨z,?_,?_⟩
  · intro e he b hb
    exact reduced e he ⟨b,hb⟩
  have summed : InClosure G p.totalDegree
      (Finsupp.linearCombination (P r K) (fun g : G => (g : P r K)) multipliers) := by
    rw [Finsupp.linearCombination_apply,Finsupp.sum]
    apply (closureSpace G p.totalDegree).sum_mem
    intro g hg
    have hmul : multipliers g≠0 := Finsupp.mem_support_iff.mp hg
    have hgen : (g : P r K)≠0 := by
      intro zero
      exact nonzero (zero ▸ g.property)
    have hd : (g : P r K).totalDegree+(multipliers g).totalDegree≤p.totalDegree := by
      rw [←totalDegree_mul_eq hgen hmul]
      exact degLex_totalDegree_monotone (bounds g)
    have input : InClosure G p.totalDegree (g : P r K) :=
      .input g.property (by omega)
    simpa only [smul_eq_mul] using InClosure.mul (multipliers g) input (by omega)
  have rewrite : p-z=Finsupp.linearCombination (P r K) (fun g : G => (g : P r K)) multipliers := by
    rw [identity,add_sub_cancel_right]
  rwa [rewrite]

theorem reduced_unique (G : Set (P r K)) (groebner : IsGroebner G)
    {p q : P r K} (hp : p∈reducedSpace G) (hq : q∈reducedSpace G)
    (difference : p-q∈Ideal.span G) : p=q := by
  by_contra different
  have nonzero : p-q≠0 := sub_ne_zero.mpr different
  obtain ⟨g,hg,divides⟩ := groebner.2 (p-q) difference nonzero
  have hs : degLex.degree (p-q)∈(p-q).support := by
    rw [mem_support_iff]
    exact (MonomialOrder.coeff_degree_ne_zero_iff (m:=degLex)).mpr nonzero
  exact (reducedSpace G).sub_mem hp hq _ hs g hg divides

theorem reduced_unique_of_reduced_zero (G : Set (P r K)) (zero : ReducedIdealTrivial G)
    {p q : P r K} (hp : p∈reducedSpace G) (hq : q∈reducedSpace G)
    (difference : p-q∈Ideal.span G) : p=q :=
  sub_eq_zero.mp (zero (p-q) difference ((reducedSpace G).sub_mem hp hq))

/-- General finite-quotient theorem for an explicit degree-compatible GB
subsystem. A bound on all standard monomial degrees is the finite-quotient
hypothesis; all input degrees are retained in the final max. -/
theorem finite_quotient_lfd_of_reduced_zero
    (F G : Set (P r K)) (included : G⊆F) (nonzero : (0 : P r K)∉G)
    (zero : ReducedIdealTrivial G)
    (s delta : ℕ)
    (standard_degree : ∀ e, Standard G e → e.degree≤s)
    (input_degree : ∀ f∈F, f.totalDegree≤delta) :
    HasLFDBound F (max (s+1) delta) := by
  classical
  let A := (P r K) ⧸ Ideal.span G
  let q : P r K →ₐ[K] A := Ideal.Quotient.mkₐ K (Ideal.span G)
  have exists_repr : ∀ a : A, ∃ z : P r K, z∈reducedSpace G ∧ q z=a := by
    intro a
    obtain ⟨p,hp⟩ := Ideal.Quotient.mkₐ_surjective K (Ideal.span G) a
    obtain ⟨z,hz,hr⟩ := exists_reduced_remainder G nonzero p
    have same : q p=q z := Ideal.Quotient.eq.mpr hr.sound
    exact ⟨z,hz,same.symm.trans hp⟩
  choose repr reduced image using exists_repr
  have unique : ∀ (a : A) (p : P r K), p∈reducedSpace G → q p=a → p=repr a := by
    intro a p hp ha
    apply reduced_unique_of_reduced_zero G zero hp (reduced a)
    apply Ideal.Quotient.eq.mp
    exact ha.trans (image a).symm
  let liftRepr : A →ₗ[K] P r K :=
    { toFun := repr
      map_add' := by
        intro a b
        apply Eq.symm
        apply unique (a+b) (repr a+repr b) ((reducedSpace G).add_mem (reduced a) (reduced b))
        rw [map_add,image,image]
      map_smul' := by
        intro c a
        apply Eq.symm
        apply unique (c • a) (c • repr a) ((reducedSpace G).smul_mem c (reduced a))
        rw [map_smul,image] }
  have representative_degree : ∀ a, (liftRepr a).totalDegree≤s := by
    intro a
    apply Finset.sup_le
    intro e he
    exact standard_degree e (reduced a e he)
  have reduction : ∀ p : P r K, InClosure F p.totalDegree (p-liftRepr (q p)) := by
    intro p
    obtain ⟨z,hz,hr⟩ := exists_reduced_remainder G nonzero p
    have same : q z=q p := (Ideal.Quotient.eq.mpr hr.sound).symm
    have normal : z=liftRepr (q p) := unique (q p) z hz same
    rw [←normal]
    exact InClosure.transport le_rfl (fun f hf hd => .input (included hf) hd) hr
  exact bounded_normal_form_lfd F q liftRepr image s (max (s+1) delta)
    representative_degree reduction (Nat.le_max_left _ _)
    (fun f hf => (input_degree f hf).trans (Nat.le_max_right _ _))

/-- The finite-quotient theorem in its usual Groebner-basis formulation. -/
theorem finite_quotient_lfd
    (F G : Set (P r K)) (included : G⊆F) (groebner : IsGroebner G)
    (s delta : ℕ)
    (standard_degree : ∀ e, Standard G e → e.degree≤s)
    (input_degree : ∀ f∈F, f.totalDegree≤delta) :
    HasLFDBound F (max (s+1) delta) := by
  apply finite_quotient_lfd_of_reduced_zero F G included groebner.1 ?_ s delta
    standard_degree input_degree
  intro p hp reduced
  exact reduced_unique G groebner reduced (reducedSpace G).zero_mem (by simpa using hp)

end ICFallDegrees.FiniteQuotient
