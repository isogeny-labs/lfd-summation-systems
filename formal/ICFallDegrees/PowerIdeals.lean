import ICFallDegrees.HomogeneousFirstFall
import Mathlib.RingTheory.MvPolynomial.Ideal
import Mathlib.RingTheory.Ideal.Quotient.Operations

noncomputable section
namespace ICFallDegrees.PowerIdeals
open MvPolynomial HomogeneousFirstFall
variable {r n : ℕ} {K : Type*} [CommRing K]
abbrev P (r : ℕ) (K : Type*) [CommRing K] := MvPolynomial (Fin r) K

def powerIdeal (d : ℕ) : Ideal (P r K) :=
  Ideal.span (Set.range fun i : Fin r => (X i : P r K)^d)

theorem mem_powerIdeal (d : ℕ) (p : P r K) :
    p∈powerIdeal d ↔ ∀ e∈p.support, ∃ i, d≤e i := by
  classical
  have generators : Set.range (fun i : Fin r => (X i : P r K)^d)=
      (fun e => monomial e (1 : K)) '' Set.range (fun i : Fin r => Finsupp.single i d) := by
    ext p
    constructor
    · rintro ⟨i,rfl⟩
      exact ⟨Finsupp.single i d,⟨i,rfl⟩,(X_pow_eq_monomial (n := i) (e := d)).symm⟩
    · rintro ⟨e,⟨i,rfl⟩,rfl⟩
      exact ⟨i,X_pow_eq_monomial (n := i) (e := d)⟩
  rw [powerIdeal,generators,mem_ideal_span_monomial_image]
  simp only [Set.mem_range,exists_exists_eq_and,Finsupp.single_le_iff]

theorem coefficient_zero (d : ℕ) {p : P r K} (member : p∈powerIdeal d)
    (e : Fin r →₀ ℕ) (small : ∀ i, e i<d) : coeff e p=0 := by
  by_contra nonzero
  obtain ⟨i,large⟩ := (mem_powerIdeal d p).mp member e (mem_support_iff.mpr nonzero)
  exact Nat.not_le_of_lt (small i) large

theorem degree_lt_eq_zero (d : ℕ) {p : P r K}
    (member : p∈powerIdeal d) (small : p.totalDegree<d) : p=0 := by
  apply MvPolynomial.eq_zero_iff.mpr
  intro e
  by_contra nonzero
  have support : e∈p.support := mem_support_iff.mpr nonzero
  obtain ⟨i,large⟩ := (mem_powerIdeal d p).mp member e support
  have bound := (Finsupp.le_degree i e).trans (le_totalDegree support)
  omega

theorem powerIdeal_antitone {d e : ℕ} (smaller : d≤e) :
    (powerIdeal e : Ideal (P r K))≤powerIdeal d := by
  apply Ideal.span_le.mpr
  rintro _ ⟨i,rfl⟩
  have identity : (X i : P r K)^e=X i^(e-d)*X i^d := by
    rw [←pow_add]
    congr 1
    omega
  change (X i : P r K)^e∈powerIdeal d
  rw [identity]
  exact (powerIdeal d).mul_mem_left _ (Ideal.subset_span ⟨i,rfl⟩)

theorem squarefree_coefficient_zero {p : P r K} (member : p∈powerIdeal 2)
    (e : Fin r →₀ ℕ) (squarefree : ∀ i, e i≤1) : coeff e p=0 :=
  coefficient_zero 2 member e (fun i => by have := squarefree i; omega)

/-- Any coordinate of a Koszul combination lies in every ideal containing
all generators other than that coordinate's own generator. -/
theorem koszul_coordinate_mem (q : Fin n → P r K) (j : Fin n) (I : Ideal (P r K))
    (other : ∀ i, i≠j → q i∈I) {h : Fin n → P r K} (member : h∈koszulSpan q) : h j∈I := by
  change h∈Submodule.span _ _ at member
  induction member using Submodule.span_induction with
  | mem p hp =>
      obtain ⟨⟨a,b⟩,rfl⟩ := hp
      by_cases ha : a=j
      · subst a
        by_cases hb : b=j
        · subst b; simp [koszulRow]
        · simpa [koszulRow,hb,Ne.symm hb] using other b hb
      · by_cases hb : b=j
        · subst b
          simpa [koszulRow,ha,Ne.symm ha] using I.neg_mem (other a ha)
        · simp [koszulRow,ha,hb,Ne.symm ha,Ne.symm hb]
  | zero => exact I.zero_mem
  | add p q _ _ hp hq => exact I.add_mem hp hq
  | smul a p _ hp => exact I.mul_mem_left a hp

end ICFallDegrees.PowerIdeals
