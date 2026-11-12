import ICFallDegrees.BinaryDescent
import ICFallDegrees.Closure

/-! Faithful transport of full ordinary-degree closures under a finite
coefficient-field extension. Coefficient projection of an arbitrary allowed
product is expanded over the field basis with all degrees retained. -/
noncomputable section
namespace ICFallDegrees.ScalarExtensionTransport
open MvPolynomial
open scoped BigOperators
variable {σ k K : Type*} {n : ℕ} [Field k] [Field K] [Algebra k K]

abbrev lift (p : MvPolynomial σ k) : MvPolynomial σ K := map (algebraMap k K) p

theorem map_degree (p : MvPolynomial σ k) : (lift (K:=K) p).totalDegree=p.totalDegree := by
  classical
  apply le_antisymm
  · apply Finset.sup_le
    intro e he
    apply le_totalDegree
    rw [mem_support_iff] at *
    intro zero
    exact he (by simp [lift,coeff_map,zero])
  · apply Finset.sup_le
    intro e he
    apply le_totalDegree
    rw [mem_support_iff] at *
    simp only [lift,coeff_map]
    change algebraMap k K (coeff e p)≠0
    intro zero
    exact he ((algebraMap k K).injective (zero.trans (map_zero _).symm))

theorem component_degree (basis : Basis (Fin n) k K) (i : Fin n) (p : MvPolynomial σ K) :
    (BinaryDescent.component basis i p).totalDegree≤p.totalDegree := by
  classical
  apply Finset.sup_le
  intro e he
  apply le_totalDegree
  by_contra absent
  have zero : coeff e p=0 := not_mem_support_iff.mp absent
  have projected : coeff e (BinaryDescent.component basis i p)=0 := by
    change basis.coord i (coeff e p)=0
    rw [zero,map_zero]
  exact (mem_support_iff.mp he) projected

theorem reconstruct (basis : Basis (Fin n) k K) (p : MvPolynomial σ K) :
    (∑ j, C (basis j)*lift (BinaryDescent.component basis j p))=p := by
  classical
  ext e
  simp only [coeff_sum,coeff_C_mul,lift,coeff_map]
  change (∑ j, basis j*algebraMap k K (basis.repr (coeff e p) j))=coeff e p
  simpa only [Algebra.smul_def,mul_comm] using basis.sum_repr (coeff e p)

theorem component_mul_X (basis : Basis (Fin n) k K) (i : Fin n)
    (p : MvPolynomial σ K) (v : σ) :
    BinaryDescent.component basis i (p*X v)=BinaryDescent.component basis i p*X v := by
  classical
  ext e
  change basis.coord i (coeff e (p*X v))=coeff e (BinaryDescent.component basis i p*X v)
  rw [coeff_mul_X',coeff_mul_X']
  split_ifs
  · rfl
  · exact map_zero _

theorem component_mul_lift (basis : Basis (Fin n) k K) (i : Fin n)
    (a : MvPolynomial σ K) (p : MvPolynomial σ k) :
    BinaryDescent.component basis i (a*lift p)=BinaryDescent.component basis i a*p := by
  induction p using MvPolynomial.induction_on with
  | C c =>
      have scalar : a*C (algebraMap k K c)=c • a := by
        ext e
        rw [mul_comm a (C _),coeff_C_mul,coeff_smul,Algebra.smul_def]
      simpa only [lift,map_C,scalar,map_smul,smul_eq_C_mul,mul_comm] using
        (BinaryDescent.component basis i).map_smul c a
  | add p q hp hq => simp only [lift,map_add,mul_add,map_add,hp,hq]
  | mul_X p v hp =>
      simp only [lift,map_mul,map_X,←mul_assoc,component_mul_X,hp]

theorem component_product (basis : Basis (Fin n) k K) (i : Fin n)
    (a p : MvPolynomial σ K) :
    BinaryDescent.component basis i (a*p)=
      ∑ j, BinaryDescent.component basis i (a*C (basis j))*BinaryDescent.component basis j p := by
  conv_lhs => rw [←reconstruct basis p]
  rw [Finset.mul_sum,map_sum]
  apply Finset.sum_congr rfl
  intro j _
  rw [←mul_assoc,component_mul_lift]

theorem project_product (basis : Basis (Fin n) k K)
    (F : Set (MvPolynomial σ k)) (D : ℕ) (a p : MvPolynomial σ K)
    (previous : ∀ j, InClosure F D (BinaryDescent.component basis j p))
    (budget : a.totalDegree+p.totalDegree≤D) (i : Fin n) :
    InClosure F D (BinaryDescent.component basis i (a*p)) := by
  rw [component_product]
  apply (closureSpace F D).sum_mem
  intro j _
  have first : (BinaryDescent.component basis i (a*C (basis j))).totalDegree≤a.totalDegree :=
    (component_degree basis i _).trans ((totalDegree_mul _ _).trans (by simp))
  exact (previous j).mul _ ((Nat.add_le_add first (component_degree basis j p)).trans budget)

theorem map_closure (F : Set (MvPolynomial σ k)) (D : ℕ) {p}
    (member : InClosure F D p) : InClosure (lift (K:=K) '' F) D (lift p) := by
  induction member with
  | input belongs degree => exact .input ⟨_,belongs,rfl⟩ ((map_degree _).le.trans degree)
  | zero => simpa only [lift,map_zero] using (InClosure.zero : InClosure (lift (K:=K) '' F) D 0)
  | add _ _ ih ih' => simpa only [lift,map_add] using ih.add ih'
  | smul c _ ih => simpa only [smul_eq_C_mul,lift,map_mul,map_C] using ih.smul (algebraMap k K c)
  | mul a hp budget ih =>
      have h := ih.mul (lift (K:=K) a) (by simpa only [map_degree] using budget)
      simpa only [lift,map_mul] using h

theorem project_closure (basis : Basis (Fin n) k K) (F : Set (MvPolynomial σ k)) (D : ℕ) {p}
    (member : InClosure (lift (K:=K) '' F) D p) :
    ∀ i, InClosure F D (BinaryDescent.component basis i p) := by
  induction member with
  | input belongs degree =>
      obtain ⟨q,hq,rfl⟩ := belongs
      intro i
      have available : InClosure F D q := .input hq (by simpa only [map_degree] using degree)
      have identity : BinaryDescent.component basis i (lift (K:=K) q)=C (basis.coord i 1)*q := by
        simpa only [one_mul,C_1] using BinaryDescent.component_scaled_map basis i 1 q
      rw [identity]
      simpa only [smul_eq_C_mul] using available.smul (basis.coord i 1)
  | zero => intro i; simpa using (InClosure.zero : InClosure F D 0)
  | add _ _ ih ih' => intro i; simpa only [map_add] using (ih i).add (ih' i)
  | smul c hp ih =>
      intro i
      rw [smul_eq_C_mul]
      exact project_product basis F D (C c) _ ih (by simpa using hp.degree_le) i
  | mul a hp budget ih => exact project_product basis F D a _ ih budget

theorem closure_iff (basis : Basis (Fin n) k K) (F : Set (MvPolynomial σ k))
    (D : ℕ) (p : MvPolynomial σ k) :
    InClosure (lift (K:=K) '' F) D (lift p) ↔ InClosure F D p := by
  constructor
  · intro hp
    have nonzero : ∃ i, basis.coord i (1 : K)≠0 := by
      by_contra absent
      have zero : ∀ i, basis.coord i (1 : K)=0 := by simpa using not_exists.mp absent
      have sum := basis.sum_repr (1 : K)
      simp only [Basis.coord_apply] at zero
      simp only [zero,zero_smul,Finset.sum_const_zero] at sum
      exact one_ne_zero sum.symm
    obtain ⟨i,hi⟩ := nonzero
    have h := project_closure basis F D hp i
    have identity : BinaryDescent.component basis i (lift (K:=K) p)=C (basis.coord i 1)*p := by
      simpa only [one_mul,C_1] using BinaryDescent.component_scaled_map basis i 1 p
    rw [identity] at h
    have scaled := h.smul ((basis.coord i (1 : K))⁻¹)
    simpa only [smul_eq_C_mul,←mul_assoc,←map_mul,inv_mul_cancel₀ hi,C_1,one_mul] using scaled
  · exact map_closure F D

theorem ideal_has_budget (F : Set (MvPolynomial σ k)) {p}
    (member : p∈Ideal.span F) : ∃ D, InClosure F D p := by
  induction member using Submodule.span_induction with
  | mem p hp => exact ⟨p.totalDegree,.input hp le_rfl⟩
  | zero => exact ⟨0,.zero⟩
  | add p q _ _ hp hq =>
      obtain ⟨D,hD⟩ := hp
      obtain ⟨E,hE⟩ := hq
      exact ⟨max D E,(hD.mono (Nat.le_max_left _ _)).add (hE.mono (Nat.le_max_right _ _))⟩
  | smul a p _ hp =>
      obtain ⟨D,hD⟩ := hp
      refine ⟨a.totalDegree+D,?_⟩
      exact (hD.mono (by omega)).mul a (Nat.add_le_add_left hD.degree_le _)

theorem bound_iff (basis : Basis (Fin n) k K) (F : Set (MvPolynomial σ k)) (d : ℕ) :
    HasLFDBound (lift (K:=K) '' F) d ↔ HasLFDBound F d := by
  constructor
  · intro bound p hp
    obtain ⟨D,hD⟩ := ideal_has_budget F hp
    have h := bound (lift p) (map_closure F D hD).sound
    rw [map_degree] at h
    exact (closure_iff basis F _ p).mp h
  · intro bound p hp
    obtain ⟨D,hD⟩ := ideal_has_budget (lift (K:=K) '' F) hp
    have projected := project_closure basis F D hD
    have recombined : InClosure (lift (K:=K) '' F) (max d p.totalDegree)
        (∑ j, C (basis j)*lift (BinaryDescent.component basis j p)) := by
      apply (closureSpace (lift (K:=K) '' F) (max d p.totalDegree)).sum_mem
      intro j _
      have derived := map_closure (K:=K) F _ (bound _ (projected j).sound)
      have lifted := derived.mono (max_le_max le_rfl (component_degree basis j p))
      simpa only [smul_eq_C_mul] using lifted.smul (basis j)
    simpa only [reconstruct] using recombined

theorem exact_iff (basis : Basis (Fin n) k K) (F : Set (MvPolynomial σ k)) (d : ℕ) :
    IsLastFallDegree (lift (K:=K) '' F) d ↔ IsLastFallDegree F d := by
  simp only [IsLastFallDegree,bound_iff basis F]

end ICFallDegrees.ScalarExtensionTransport
