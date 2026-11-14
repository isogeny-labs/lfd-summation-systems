/-
Portions adapted from ArkLib.Data.Polynomial.ClassicalWronskian,
copyright (c) 2026 ArkLib Contributors (Alexander Hicks),
licensed under the Apache License, Version 2.0.

This is a standalone Mathlib v4.19 port using no ArkLib dependency.
-/
import ICFallDegrees.GrowthPolynomialWronskianMonomials
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Data.Nat.Factorial.NatCast
import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.RingTheory.Polynomial.Pochhammer
import Mathlib.Tactic

/-!
# N8: nonzero polynomial Wronskian of distinct-degree polynomials

For nonzero polynomials of pairwise distinct degrees, the leading coefficient
of their Wronskian equals a nonzero leading-coefficient product times a
Vandermonde determinant in the degrees, in characteristic zero.
This is an actual general r-column nonvanishing theorem under a
degree-injectivity hypothesis. The degree-echelon reduction for an arbitrary
independent polynomial family is not yet proved here.
-/

noncomputable section
namespace ICFallDegrees.GrowthPolynomialWronskianDistinctDegree

open Polynomial Module Matrix
open scoped Function BigOperators

variable {K : Type*} [Field K]

private def classicalWronskian (sigma : ℕ) (P : Fin sigma → Polynomial K) :
    Polynomial K :=
  (Matrix.of fun i j : Fin sigma => (Polynomial.derivative^[i.val]) (P j)).det

private lemma isUnit_natCast_descFactorial {F : Type*} [Field F]
    {d i k : ℕ} (hi : i ≤ d) (hd : d < k)
    (hk : ringChar F = 0 ∨ k ≤ ringChar F) :
    IsUnit (d.descFactorial i : F) := by
  have hfac : IsUnit (d.factorial : F) := by
    rcases hk with hk0 | hkpos
    · letI : CharZero F := (CharP.ringChar_zero_iff_CharZero F).mp hk0
      exact IsUnit.natCast_factorial_of_algebra F d
    · letI : NeZero (ringChar F) :=
        ⟨Nat.ne_zero_of_lt ((Nat.zero_le d).trans_lt (hd.trans_le hkpos))⟩
      letI : Fact (ringChar F).Prime := CharP.char_is_prime_of_pos F _
      exact (IsUnit.natCast_factorial_iff_of_charP (ringChar F)).2 (hd.trans_le hkpos)
  have hmul : IsUnit (((d - i).factorial : F) * (d.descFactorial i : F)) := by
    rw [← Nat.cast_mul, Nat.factorial_mul_descFactorial hi]
    exact hfac
  exact (IsUnit.mul_iff.mp hmul).2

private lemma natDegree_iterate_derivative_eq {F : Type*} [Field F]
    {p : F[X]} (hp : p ≠ 0) {i k : ℕ} (hi : i ≤ p.natDegree)
    (hdeg : p.natDegree < k) (hk : ringChar F = 0 ∨ k ≤ ringChar F) :
    (derivative^[i] p).natDegree = p.natDegree - i := by
  apply le_antisymm (natDegree_iterate_derivative p i)
  apply le_of_not_gt
  intro hlt
  have hcoeff := coeff_eq_zero_of_natDegree_lt hlt
  rw [coeff_iterate_derivative] at hcoeff
  have hadd : p.natDegree - i + i = p.natDegree := Nat.sub_add_cancel hi
  rw [hadd, coeff_natDegree, nsmul_eq_mul] at hcoeff
  exact (mul_ne_zero
    (isUnit_natCast_descFactorial hi hdeg hk).ne_zero
    (leadingCoeff_ne_zero.mpr hp)) hcoeff

private lemma leadingCoeff_iterate_derivative {F : Type*} [Field F]
    {p : F[X]} (hp : p ≠ 0) {i k : ℕ} (hi : i ≤ p.natDegree)
    (hdeg : p.natDegree < k) (hk : ringChar F = 0 ∨ k ≤ ringChar F) :
    (derivative^[i] p).leadingCoeff =
      (p.natDegree.descFactorial i : F) * p.leadingCoeff := by
  rw [leadingCoeff, natDegree_iterate_derivative_eq hp hi hdeg hk,
    coeff_iterate_derivative, Nat.sub_add_cancel hi, nsmul_eq_mul, coeff_natDegree]

/-- The classical Wronskian of nonzero polynomials with pairwise distinct degrees is nonzero
in characteristic zero, and remains so in positive characteristic when their common strict
degree bound is no greater than the characteristic. -/
private theorem classicalWronskian_ne_zero_of_natDegree_injective
    {F : Type*} [Field F] {sigma k : ℕ}
    (P : Fin sigma → F[X]) (hP0 : ∀ j, P j ≠ 0)
    (hPdeg : ∀ j, (P j).natDegree < k)
    (hPinj : Function.Injective (fun j => (P j).natDegree))
    (hk : ringChar F = 0 ∨ k ≤ ringChar F) : classicalWronskian sigma P ≠ 0 := by
  classical
  let d : Fin sigma → ℕ := fun j => (P j).natDegree
  let D : ℕ := (∑ j, d j) - ∑ i : Fin sigma, i.val
  let A : Matrix (Fin sigma) (Fin sigma) F :=
    Matrix.of fun i j => (d j).descFactorial i.val
  have hcastinj : Function.Injective (fun j => (d j : F)) := by
    intro i j hij
    apply hPinj
    rcases hk with hk0 | hkpos
    · letI : CharZero F := (CharP.ringChar_zero_iff_CharZero F).mp hk0
      exact Nat.cast_injective hij
    · exact CharP.natCast_injOn_Iio F (ringChar F)
        ((hPdeg i).trans_le hkpos) ((hPdeg j).trans_le hkpos) hij
  have hAdet : A.det ≠ 0 := by
    have hv := (Matrix.det_vandermonde_ne_zero_iff (R := F)
      (v := fun j => (d j : F))).2 hcastinj
    have heval := Matrix.det_eval_matrixOfPolynomials_eq_det_vandermonde
      (v := fun j => (d j : F)) (p := fun i : Fin sigma => descPochhammer F i.val)
      (fun i => descPochhammer_natDegree F i.val)
      (fun i => monic_descPochhammer F i.val)
    have heq : (Matrix.vandermonde (fun j => (d j : F))).det = A.det := by
      rw [heval]
      rw [← Matrix.det_transpose]
      congr 1
      ext i j
      simp [A, d, descPochhammer_eval_eq_descFactorial]
    rwa [heq] at hv
  let L : Fin sigma → F := fun j => (P j).leadingCoeff
  have hL : ∀ j, L j ≠ 0 := fun j => leadingCoeff_ne_zero.mpr (hP0 j)
  let B : Matrix (Fin sigma) (Fin sigma) F :=
    Matrix.of fun i j => (d j).descFactorial i.val * L j
  have hBdet : B.det ≠ 0 := by
    have hBA : B = A * Matrix.diagonal L := by
      ext i j
      simp [B, A]
    rw [hBA, Matrix.det_mul, Matrix.det_diagonal]
    exact mul_ne_zero hAdet (Finset.prod_ne_zero_iff.mpr fun j _ => hL j)
  apply fun hW => hBdet ?_
  have hcoeff : (classicalWronskian sigma P).coeff D = B.det := by
    unfold classicalWronskian
    rw [Matrix.det_apply', Polynomial.finset_sum_coeff, Matrix.det_apply']
    apply Finset.sum_congr rfl
    intro g hg
    simp only [Matrix.of_apply]
    simp only [coeff_intCast_mul, mul_eq_mul_left_iff]
    apply Or.inl
    by_cases hall : ∀ j : Fin sigma, (g j).val ≤ d j
    · have hder0 : ∀ j : Fin sigma, derivative^[((g j).val)] (P j) ≠ 0 := by
        intro j hzero
        have hlc := leadingCoeff_iterate_derivative (hP0 j) (hall j) (hPdeg j) hk
        rw [hzero, leadingCoeff_zero] at hlc
        exact (mul_ne_zero
          (isUnit_natCast_descFactorial (hall j) (hPdeg j) hk).ne_zero
          (leadingCoeff_ne_zero.mpr (hP0 j))) (by simpa [d] using hlc.symm)
      have hD : (∑ j, (d j - (g j).val)) = D := by
        simp only [D]
        rw [Finset.sum_tsub_distrib Finset.univ (fun j _ => hall j)]
        congr 1
        exact Equiv.sum_comp g (fun i : Fin sigma => i.val)
      have hproddeg : (∏ j : Fin sigma, derivative^[((g j).val)] (P j)).natDegree = D := by
        rw [natDegree_prod Finset.univ _ (fun j _ => hder0 j)]
        simp_rw [natDegree_iterate_derivative_eq (hP0 _) (hall _) (hPdeg _) hk]
        exact hD
      rw [← hproddeg, coeff_natDegree, leadingCoeff_prod]
      simp only [B, Matrix.of_apply]
      apply Finset.prod_congr rfl
      intro j hj
      rw [leadingCoeff_iterate_derivative (hP0 j) (hall j) (hPdeg j) hk]
    · simp only [not_forall, not_le] at hall
      obtain ⟨j, hj⟩ := hall
      have hzero : derivative^[((g j).val)] (P j) = 0 :=
        iterate_derivative_eq_zero hj
      rw [Finset.prod_eq_zero (Finset.mem_univ j) hzero, coeff_zero]
      simp only [B, Matrix.of_apply]
      symm
      apply Finset.prod_eq_zero (Finset.mem_univ j)
      change ((d j).descFactorial (g j).val : F) * L j = 0
      rw [Nat.descFactorial_eq_zero_iff_lt.mpr hj, Nat.cast_zero, zero_mul]
  rw [hW, coeff_zero] at hcoeff
  exact hcoeff.symm



/-- General characteristic-bounded polynomial Wronskian nonvanishing.
Pairwise distinct degrees strictly below k are enough, provided either
the characteristic is zero or it is at least k. This exposes the
positive-characteristic case already established by the internal
descending-factorial/Vandermonde proof. -/
theorem polynomial_wronskian_ne_zero_of_degree_injective_char_bound
    {r k : ℕ} (P : Fin r → Polynomial K)
    (hP0 : ∀ j, P j ≠ 0)
    (hPdeg : ∀ j, (P j).natDegree < k)
    (hPinj : Function.Injective (fun j => (P j).natDegree))
    (hchar : ringChar K = 0 ∨ k ≤ ringChar K) :
    (ICFallDegrees.GrowthRatFuncPolynomialWronskianMap.polynomialDerivativeMatrix r P).det ≠ 0 :=
  classicalWronskian_ne_zero_of_natDegree_injective P hP0 hPdeg hPinj hchar

/-- A genuine r-column polynomial Wronskian is nonzero when the input
polynomials are nonzero and have pairwise distinct degrees in characteristic
zero. No Wronskian nonvanishing premise is used. -/
theorem polynomial_wronskian_ne_zero_of_degree_injective [CharZero K]
    {r k : ℕ} (P : Fin r → Polynomial K)
    (hP0 : ∀ j, P j ≠ 0)
    (hPdeg : ∀ j, (P j).natDegree < k)
    (hPinj : Function.Injective (fun j => (P j).natDegree)) :
    (ICFallDegrees.GrowthRatFuncPolynomialWronskianMap.polynomialDerivativeMatrix r P).det ≠ 0 := by
  have hchar : ringChar K = 0 := (CharP.ringChar_zero_iff_CharZero K).2 inferInstance
  exact classicalWronskian_ne_zero_of_natDegree_injective P hP0 hPdeg hPinj
    (Or.inl hchar)

end ICFallDegrees.GrowthPolynomialWronskianDistinctDegree
