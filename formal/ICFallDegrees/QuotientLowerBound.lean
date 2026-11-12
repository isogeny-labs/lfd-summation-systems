import ICFallDegrees.FiniteQuotient
import Mathlib.Tactic.Ring

noncomputable section

/-!
# Closed quotient spans as lower-bound witnesses

An ordinary input-closure invariant is proved by showing that every row of
degree below D in the proposed quotient span has zero quotient image. This
accounts for low-degree field equations rather than pretending they are absent.
-/

namespace ICFallDegrees.FiniteQuotient

open MvPolynomial

variable {r : ℕ} {K : Type*} [Field K]

def quotientSpan (J : Ideal (P r K)) (h : P r K) : Submodule K (P r K) where
  carrier := fun p => ∃ c : K, p-C c*h∈J
  zero_mem' := ⟨0,by simp⟩
  add_mem' := by
    rintro p q ⟨c,hc⟩ ⟨d,hd⟩
    refine ⟨c+d,?_⟩
    have identity : p+q-C (c+d)*h=(p-C c*h)+(q-C d*h) := by simp only [map_add]; ring
    rw [identity]
    exact J.add_mem hc hd
  smul_mem' := by
    rintro a p ⟨c,hc⟩
    refine ⟨a*c,?_⟩
    have identity : a • p-C (a*c)*h=C a*(p-C c*h) := by
      simp only [smul_eq_C_mul,map_mul]
      ring
    rw [identity]
    exact J.mul_mem_left _ hc

theorem remainder_degree_le {G : Set (P r K)} {p z : P r K}
    (reduction : InClosure G p.totalDegree (p-z)) : z.totalDegree≤p.totalDegree := by
  have identity : z=p-(p-z) := by ring
  rw [identity]
  exact (totalDegree_sub_le _ _).trans (max_le le_rfl reduction.degree_le)

/-- Low-degree rows in a scalar span of a reduced top-degree generator have
zero quotient image. This proves the multiplication part of the invariant. -/
theorem low_quotient_span_mem_ideal
    (G : Set (P r K)) (nonzero : (0 : P r K)∉G) (zero : ReducedIdealTrivial G)
    (h : P r K) (reduced : h∈reducedSpace G) (D : ℕ) (degree_h : h.totalDegree=D)
    {p : P r K} (hp : p∈quotientSpan (Ideal.span G) h) (degree_p : p.totalDegree<D) :
    p∈Ideal.span G := by
  obtain ⟨c,hc⟩ := hp
  obtain ⟨z,hz,hr⟩ := exists_reduced_remainder G nonzero p
  have difference : z-C c*h∈Ideal.span G := by
    have identity : z-C c*h=(p-C c*h)-(p-z) := by ring
    rw [identity]
    exact (Ideal.span G).sub_mem hc hr.sound
  have scaled_reduced : C c*h∈reducedSpace G := by
    simpa only [smul_eq_C_mul] using (reducedSpace G).smul_mem c reduced
  have same := reduced_unique_of_reduced_zero G zero hz scaled_reduced difference
  have dz := remainder_degree_le hr
  have scalar_zero : c=0 := by
    by_contra hc0
    have hh0 : h≠0 := by
      intro hh
      simp only [hh,totalDegree_zero] at degree_h
      omega
    have ds : (C c*h : P r K).totalDegree=D := by
      rw [totalDegree_mul_eq (C_ne_zero.mpr hc0) hh0,totalDegree_C,degree_h,zero_add]
    rw [same,ds] at dz
    omega
  simpa [scalar_zero] using hc

theorem closure_le_quotient_span
    (F G : Set (P r K)) (nonzero : (0 : P r K)∉G) (zero : ReducedIdealTrivial G)
    (h : P r K) (reduced : h∈reducedSpace G) (D : ℕ) (degree_h : h.totalDegree=D)
    (inputs : ∀ p∈F, p.totalDegree≤D → p∈quotientSpan (Ideal.span G) h) :
    closureSpace F D≤quotientSpan (Ideal.span G) h := by
  apply closure_minimal F D (quotientSpan (Ideal.span G) h) inputs
  intro a p hp budget
  by_cases constant : a.totalDegree=0
  · rw [totalDegree_eq_zero_iff_eq_C.mp constant,C_mul']
    exact (quotientSpan (Ideal.span G) h).smul_mem _ hp
  · have low : p.totalDegree<D := by omega
    have in_ideal := low_quotient_span_mem_ideal G nonzero zero h reduced D degree_h hp low
    refine ⟨0,?_⟩
    simpa using (Ideal.span G).mul_mem_left a in_ideal

end ICFallDegrees.FiniteQuotient
