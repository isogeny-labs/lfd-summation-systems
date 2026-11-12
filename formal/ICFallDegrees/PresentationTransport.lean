import ICFallDegrees.AffineTransport

/-! Exact transport of the full iterated closure. The input replacement
criterion charges each input at its own ordinary degree. Affine equivalences
are checked in both directions, so exact LFD values, not just upper bounds,
can be transported. No equality-of-ideals shortcut is used. -/
noncomputable section
namespace ICFallDegrees.PresentationTransport
open MvPolynomial
variable {σ τ K : Type*} [CommRing K]

theorem ideal_map_mem (f : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (F : Set (MvPolynomial σ K)) {p} (member : p∈Ideal.span F) :
    f p∈Ideal.span (f '' F) := by
  induction member using Submodule.span_induction with
  | mem p hp => exact Ideal.subset_span ⟨p,hp,rfl⟩
  | zero => simpa using (Ideal.span (f '' F)).zero_mem
  | add p q _ _ hp hq => simpa only [map_add] using (Ideal.span (f '' F)).add_mem hp hq
  | smul a p _ hp =>
      simpa only [smul_eq_mul,map_mul] using (Ideal.span (f '' F)).mul_mem_left (f a) hp

theorem affine_degree_eq (e : MvPolynomial σ K ≃ₐ[K] MvPolynomial τ K)
    (forward : ∀ i, (e (X i)).totalDegree≤1)
    (backward : ∀ i, (e.symm (X i)).totalDegree≤1) (p : MvPolynomial σ K) :
    (e p).totalDegree=p.totalDegree := by
  apply le_antisymm (AffineTransport.degree_le e.toAlgHom forward p)
  simpa using AffineTransport.degree_le e.symm.toAlgHom backward (e p)

theorem affine_closure_iff (e : MvPolynomial σ K ≃ₐ[K] MvPolynomial τ K)
    (forward : ∀ i, (e (X i)).totalDegree≤1)
    (backward : ∀ i, (e.symm (X i)).totalDegree≤1)
    (F : Set (MvPolynomial σ K)) (D : ℕ) (p : MvPolynomial σ K) :
    InClosure (e '' F) D (e p) ↔ InClosure F D p := by
  constructor
  · intro hp
    have mapped := AffineTransport.map_closure e.symm.toAlgHom backward (e '' F) F D
      (by
        rintro _ ⟨q,hq,rfl⟩ degree
        have dq : q.totalDegree≤D := by simpa using
          (AffineTransport.degree_le e.symm.toAlgHom backward (e q)).trans degree
        simpa using (InClosure.input hq dq : InClosure F D q)) hp
    simpa using mapped
  · intro hp
    exact AffineTransport.map_closure e.toAlgHom forward F (e '' F) D
      (fun q hq degree => .input ⟨q,hq,rfl⟩
        ((AffineTransport.degree_le e.toAlgHom forward q).trans degree)) hp

theorem affine_ideal_iff (e : MvPolynomial σ K ≃ₐ[K] MvPolynomial τ K)
    (F : Set (MvPolynomial σ K)) (p : MvPolynomial σ K) :
    e p∈Ideal.span (e '' F) ↔ p∈Ideal.span F := by
  constructor
  · intro hp
    have mapped := ideal_map_mem e.symm.toAlgHom (e '' F) hp
    have images : e.symm '' (e '' F)=F := by
      ext q
      constructor
      · rintro ⟨_,⟨r,hr,rfl⟩,rfl⟩
        simpa using hr
      · intro hq
        exact ⟨e q,⟨q,hq,rfl⟩,by simp⟩
    simpa only [AlgEquiv.coe_algHom,AlgEquiv.symm_apply_apply,images] using mapped
  · exact ideal_map_mem e.toAlgHom F

theorem affine_bound_iff (e : MvPolynomial σ K ≃ₐ[K] MvPolynomial τ K)
    (forward : ∀ i, (e (X i)).totalDegree≤1)
    (backward : ∀ i, (e.symm (X i)).totalDegree≤1)
    (F : Set (MvPolynomial σ K)) (d : ℕ) :
    HasLFDBound (e '' F) d ↔ HasLFDBound F d := by
  constructor
  · intro bound p hp
    have h := bound (e p) ((affine_ideal_iff e F p).mpr hp)
    rw [affine_degree_eq e forward backward p] at h
    exact (affine_closure_iff e forward backward F _ p).mp h
  · intro bound q hq
    let p := e.symm q
    have image : e p=q := e.apply_symm_apply q
    have hp : p∈Ideal.span F := (affine_ideal_iff e F p).mp (image.symm ▸ hq)
    have h := (affine_closure_iff e forward backward F _ p).mpr (bound p hp)
    have degree : q.totalDegree=p.totalDegree := image ▸ affine_degree_eq e forward backward p
    simpa only [image,←degree] using h

theorem affine_exact_iff (e : MvPolynomial σ K ≃ₐ[K] MvPolynomial τ K)
    (forward : ∀ i, (e (X i)).totalDegree≤1)
    (backward : ∀ i, (e.symm (X i)).totalDegree≤1)
    (F : Set (MvPolynomial σ K)) (d : ℕ) :
    IsLastFallDegree (e '' F) d ↔ IsLastFallDegree F d := by
  simp only [IsLastFallDegree,affine_bound_iff e forward backward F]

theorem replace_closure (F G : Set (MvPolynomial σ K))
    (inputs : ∀ p∈F, InClosure G p.totalDegree p) (D : ℕ)
    {p} (member : InClosure F D p) : InClosure G D p :=
  InClosure.transport le_rfl (fun q hq degree => (inputs q hq).mono degree) member

theorem replace_ideal (F G : Set (MvPolynomial σ K))
    (inputs : ∀ p∈F, InClosure G p.totalDegree p) : Ideal.span F≤Ideal.span G := by
  apply Ideal.span_le.mpr
  intro p hp
  exact (inputs p hp).sound

theorem same_budget_closure_iff (F G : Set (MvPolynomial σ K))
    (forward : ∀ p∈F, InClosure G p.totalDegree p)
    (backward : ∀ p∈G, InClosure F p.totalDegree p) (D : ℕ) (p : MvPolynomial σ K) :
    InClosure F D p ↔ InClosure G D p :=
  ⟨replace_closure F G forward D,replace_closure G F backward D⟩

theorem same_budget_bound_iff (F G : Set (MvPolynomial σ K))
    (forward : ∀ p∈F, InClosure G p.totalDegree p)
    (backward : ∀ p∈G, InClosure F p.totalDegree p) (d : ℕ) :
    HasLFDBound F d ↔ HasLFDBound G d := by
  constructor
  · intro bound p hp
    exact replace_closure F G forward _ (bound p (replace_ideal G F backward hp))
  · intro bound p hp
    exact replace_closure G F backward _ (bound p (replace_ideal F G forward hp))

theorem same_budget_exact_iff (F G : Set (MvPolynomial σ K))
    (forward : ∀ p∈F, InClosure G p.totalDegree p)
    (backward : ∀ p∈G, InClosure F p.totalDegree p) (d : ℕ) :
    IsLastFallDegree F d ↔ IsLastFallDegree G d := by
  simp only [IsLastFallDegree,same_budget_bound_iff F G forward backward]

end ICFallDegrees.PresentationTransport
