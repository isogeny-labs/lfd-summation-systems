/-
Portions adapted from ArkLib/Data/Polynomial/ClassicalWronskian.lean
Copyright (c) 2026 ArkLib Contributors (Alexander Hicks).
Apache-2.0 license, https://github.com/Verified-zkEVM/ArkLib/blob/main/LICENSE
Adaptation: use pinned Mathlib 4.19 and N8 polynomial Wronskian determinant.
-/
import ICFallDegrees.GrowthPolynomialWronskianDistinctDegree
import Mathlib.Algebra.Polynomial.Degree.TrailingDegree
import Mathlib.LinearAlgebra.Basis.Fin
import Mathlib.LinearAlgebra.Dual.Lemmas
import Mathlib.LinearAlgebra.Matrix.Basis
import Mathlib.Tactic

/-!
# N8: Polynomial Wronskian of a basis of any finite-dimensional subspace

A finite-dimensional polynomial subspace has a basis with distinct degrees.
By the previous Vandermonde theorem, its polynomial Wronskian is nonzero.
Changing basis multiplies the Wronskian by an invertible constant determinant.
No spread bound or rational-function cofactor properties are assumed.
-/

noncomputable section
namespace ICFallDegrees.GrowthPolynomialWronskianBasis

open Polynomial Module Matrix
open scoped Function BigOperators

private def classicalWronskian {F : Type*} [CommRing F]
    (sigma : ℕ) (P : Fin sigma → F[X]) : F[X] :=
  (Matrix.of fun i j : Fin sigma => derivative^[i.val] (P j)).det

private theorem exists_basis_natDegree_injective_aux
    {F M : Type*} [Field F] [AddCommGroup M] [Module F M]
    [FiniteDimensional F M] (e : M →ₗ[F] F[X]) (he : Function.Injective e)
    (k : ℕ) (hdeg : ∀ x, x ≠ 0 → (e x).natDegree < k) :
    ∃ b : Basis (Fin (finrank F M)) F M,
      Function.Injective (fun i => (e (b i)).natDegree) := by
  induction k generalizing M with
  | zero =>
      have hM : ∀ x : M, x = 0 := by
        intro x
        by_contra hx
        exact (Nat.not_lt_zero _ (hdeg x hx))
      have hrank : finrank F M = 0 :=
        Module.finrank_zero_iff.mpr ⟨fun x y => (hM x).trans (hM y).symm⟩
      let b := Module.finBasis F M
      refine ⟨b, ?_⟩
      intro i j hij
      exact Fin.elim0 (Fin.cast hrank i)
  | succ k ih =>
      let f : M →ₗ[F] F := (lcoeff F k).comp e
      by_cases hf : f = 0
      · have hdeg' : ∀ x, x ≠ 0 → (e x).natDegree < k := by
          intro x hx
          have hle : (e x).natDegree ≤ k := Nat.lt_succ_iff.mp (hdeg x hx)
          apply lt_of_le_of_ne hle
          intro heq
          have hcoeff : (e x).coeff k = 0 := by
            have := LinearMap.congr_fun hf x
            simpa [f, lcoeff_apply] using this
          have hex : e x ≠ 0 := fun h => hx (he (h.trans (map_zero e).symm))
          rw [← heq, coeff_natDegree, leadingCoeff_eq_zero] at hcoeff
          exact hex hcoeff
        exact ih e he hdeg'
      · have hdegKer : ∀ x : LinearMap.ker f, x ≠ 0 →
            (e.comp (LinearMap.ker f).subtype x).natDegree < k := by
          intro x hx
          have hxM : (x : M) ≠ 0 := fun h => hx (Subtype.ext h)
          have hle : (e x).natDegree ≤ k := Nat.lt_succ_iff.mp (hdeg x hxM)
          apply lt_of_le_of_ne hle
          intro heq
          have hcoeff : (e x).coeff k = 0 := by
            have hxker : f (x : M) = 0 := LinearMap.mem_ker.mp x.2
            change (lcoeff F k) (e x) = 0 at hxker
            change (e x).coeff k = 0
            simpa only [lcoeff_apply] using hxker
          have hex : e x ≠ 0 := fun h => hxM (he (h.trans (map_zero e).symm))
          rw [← heq, coeff_natDegree, leadingCoeff_eq_zero] at hcoeff
          exact hex hcoeff
        have heKer : Function.Injective (e.comp (LinearMap.ker f).subtype) :=
          he.comp (LinearMap.ker f).subtype_injective
        obtain ⟨b, hb⟩ := ih (e.comp (LinearMap.ker f).subtype) heKer hdegKer
        obtain ⟨y, hy⟩ : ∃ y : M, f y ≠ 0 := by
          simpa [LinearMap.ext_iff] using hf
        have hli : ∀ (c : F), ∀ x ∈ LinearMap.ker f, c • y + x = 0 → c = 0 := by
          intro c x hxN hzero
          have hfx : f x = 0 := LinearMap.mem_ker.mp hxN
          have hmap : f (c • y + x) = f 0 := congr_arg f hzero
          rw [map_add, map_smul, map_zero, hfx, add_zero] at hmap
          have hmul : c * f y = 0 := by simpa [smul_eq_mul] using hmap
          exact (mul_eq_zero.mp hmul).resolve_right hy
        have hsp : ∀ z : M, ∃ c : F, z + c • y ∈ LinearMap.ker f := by
          intro z
          refine ⟨-(f z / f y), ?_⟩
          simp only [LinearMap.mem_ker, map_add, map_smul]
          simp [smul_eq_mul, hy]
        let b' : Basis (Fin (finrank F (LinearMap.ker f) + 1)) F M :=
          Basis.mkFinCons y b hli hsp
        have hrank : finrank F (LinearMap.ker f) + 1 = finrank F M :=
          Module.Dual.finrank_ker_add_one_of_ne_zero hf
        let b'' : Basis (Fin (finrank F M)) F M := b'.reindex (finCongr hrank)
        refine ⟨b'', ?_⟩
        have hydeg : (e y).natDegree = k := by
          have hy0 : y ≠ 0 := by
            intro h
            apply hy
            rw [h, map_zero]
          apply le_antisymm (Nat.lt_succ_iff.mp (hdeg y hy0))
          apply le_of_not_gt
          intro hlt
          have : (e y).coeff k = 0 := coeff_eq_zero_of_natDegree_lt hlt
          exact hy (by simpa [f, lcoeff_apply] using this)
        have hsmall : ∀ i, (e ((b i : LinearMap.ker f) : M)).natDegree < k := by
          intro i
          exact hdegKer (b i) (b.ne_zero i)
        have hcons : Function.Injective
            (fun i : Fin (finrank F (LinearMap.ker f) + 1) => (e (b' i)).natDegree) := by
          have hb' : (b' : Fin (_ + 1) → M) = Fin.cons y ((↑) ∘ b) := by
            simp [b']
          have heq : (fun i => (e (b' i)).natDegree) =
              Fin.cons k (fun i => (e ((b i : LinearMap.ker f) : M)).natDegree) := by
            funext i
            rw [hb']
            cases i using Fin.cases with
            | zero => simp [hydeg]
            | succ j => simp [Function.comp_apply]
          rw [heq]
          intro i j hij
          cases i using Fin.cases with
          | zero =>
              cases j using Fin.cases with
              | zero => rfl
              | succ j =>
                  have hj : k = (e ((b j : LinearMap.ker f) : M)).natDegree := by
                    exact hij
                  exact False.elim ((Nat.ne_of_lt (hsmall j)) hj.symm)
          | succ i =>
              cases j using Fin.cases with
              | zero =>
                  have hi : (e ((b i : LinearMap.ker f) : M)).natDegree = k := by
                    exact hij
                  exact False.elim ((Nat.ne_of_lt (hsmall i)) hi)
              | succ j =>
                  exact congrArg Fin.succ (hb hij)
        simp only [b'', Basis.reindex_apply]
        exact hcons.comp (finCongr hrank).symm.injective

/-- Replacing the polynomials by constant linear combinations right-multiplies the
Wronskian matrix by the coefficient matrix. -/
theorem classicalWronskian_of_linearComb {F : Type*} [Field F] {sigma : ℕ}
    (P c : Fin sigma → F[X]) (U : Matrix (Fin sigma) (Fin sigma) F)
    (hc : ∀ j, c j = ∑ i, U i j • P i) :
    classicalWronskian sigma c =
      classicalWronskian sigma P * C U.det := by
  classical
  have hM : (Matrix.of fun i j : Fin sigma => derivative^[i.val] (c j))
      = (Matrix.of fun i j : Fin sigma => derivative^[i.val] (P j)) *
        ((C : F →+* F[X]).mapMatrix U) := by
    refine Matrix.ext fun i j => ?_
    simp only [Matrix.of_apply, Matrix.mul_apply, RingHom.mapMatrix_apply, Matrix.map_apply]
    rw [hc j, iterate_derivative_sum]
    exact Finset.sum_congr rfl fun i' _ => by
      rw [iterate_derivative_smul]
      simp [smul_eq_C_mul, mul_comm]
  unfold classicalWronskian
  rw [hM, Matrix.det_mul, ← RingHom.map_det]

/-- If a finite-dimensional polynomial subspace has a basis of polynomials of degree `< k`,
its classical Wronskian is nonzero in characteristic zero, and in positive characteristic at
least `k`.

The proof changes to a basis in degree-echelon form, where the top-degree coefficient of the
Wronskian is a nonzero Vandermonde determinant; changing back multiplies the Wronskian by a
nonzero constant. -/
theorem polynomial_wronskian_ne_zero_of_basis_char_bound {F : Type*} [Field F] {sigma k : ℕ}
    {B : Submodule F F[X]} (bas : Basis (Fin sigma) F B)
    (hdeg : ∀ j, ((bas j : B) : F[X]).natDegree < k)
    (hchar : ringChar F = 0 ∨ k ≤ ringChar F) :
    classicalWronskian sigma (fun j => ((bas j : B) : F[X])) ≠ 0 := by
  classical
  letI : Module.Finite F B := Module.Finite.of_basis bas
  have hrank : finrank F B = sigma := by
    rw [Module.finrank_eq_card_basis bas, Fintype.card_fin]
  have hbound : ∀ x : B, x ≠ 0 → (B.subtype x).natDegree < k := by
    intro x hx
    have hmem : B.subtype x ∈ degreeLT F k := by
      rw [← bas.sum_repr x]
      rw [map_sum]
      exact Submodule.sum_mem _ fun j _ => Submodule.smul_mem _ _ (by
        rw [mem_degreeLT, degree_eq_natDegree, Nat.cast_lt]
        · exact hdeg j
        · exact fun h => bas.ne_zero j (Subtype.ext h))
    rw [mem_degreeLT, degree_eq_natDegree, Nat.cast_lt] at hmem
    · exact hmem
    · exact fun h => hx (B.subtype_injective (h.trans (map_zero B.subtype).symm))
  obtain ⟨cb₀, hcb₀⟩ :=
    exists_basis_natDegree_injective_aux B.subtype B.subtype_injective k hbound
  let cb : Basis (Fin sigma) F B := cb₀.reindex (finCongr hrank)
  have hcbinj : Function.Injective
      (fun j => (((cb j : B) : F[X])).natDegree) := by
    simp only [cb, Basis.reindex_apply]
    exact hcb₀.comp (finCongr hrank).symm.injective
  have hcbdeg : ∀ j, (((cb j : B) : F[X])).natDegree < k := fun j =>
    hbound (cb j) (cb.ne_zero j)
  have hcbW : classicalWronskian sigma (fun j => ((cb j : B) : F[X])) ≠ 0 :=
    ICFallDegrees.GrowthPolynomialWronskianDistinctDegree.polynomial_wronskian_ne_zero_of_degree_injective_char_bound _
      (fun j => fun h => cb.ne_zero j (Subtype.ext h)) hcbdeg hcbinj hchar
  let U : Matrix (Fin sigma) (Fin sigma) F := bas.toMatrix (⇑cb)
  have hcomb : ∀ j, ((cb j : B) : F[X]) =
      ∑ i, U i j • ((bas i : B) : F[X]) := by
    intro j
    have h1 : ∑ i, U i j • bas i = cb j :=
      Basis.sum_toMatrix_smul_self bas (⇑cb) j
    have h2 : B.subtype (∑ i, U i j • bas i) = B.subtype (cb j) := by rw [h1]
    rw [map_sum] at h2
    simp only [map_smul, Submodule.coe_subtype] at h2
    exact h2.symm
  have hdetU : U.det ≠ 0 := by
    have h := congrArg Matrix.det (Basis.toMatrix_mul_toMatrix_flip bas cb)
    rw [Matrix.det_mul, Matrix.det_one] at h
    intro h0
    rw [h0, zero_mul] at h
    exact zero_ne_one h
  have hW := classicalWronskian_of_linearComb
    (fun j => ((bas j : B) : F[X])) (fun j => ((cb j : B) : F[X])) U hcomb
  intro hzero
  rw [hzero, zero_mul] at hW
  exact hcbW hW


/-- Backwards-compatible characteristic-zero Wronskian basis theorem.
The genuine positive-characteristic extension has the explicit degree
guard `k ≤ ringChar F` in `polynomial_wronskian_ne_zero_of_basis_char_bound`. -/
theorem polynomial_wronskian_ne_zero_of_basis {F : Type*}
    [Field F] [CharZero F] {sigma k : ℕ}
    {B : Submodule F F[X]} (bas : Basis (Fin sigma) F B)
    (hdeg : ∀ j, ((bas j : B) : F[X]).natDegree < k) :
    classicalWronskian sigma (fun j => ((bas j : B) : F[X])) ≠ 0 := by
  have hzero : ringChar F = 0 :=
    (CharP.ringChar_zero_iff_CharZero F).2 inferInstance
  exact polynomial_wronskian_ne_zero_of_basis_char_bound bas hdeg (Or.inl hzero)

#print axioms polynomial_wronskian_ne_zero_of_basis_char_bound
#print axioms polynomial_wronskian_ne_zero_of_basis


end ICFallDegrees.GrowthPolynomialWronskianBasis
