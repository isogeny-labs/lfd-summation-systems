import ICFallDegrees.GrowthPrimeBoundedRationalWronskian
import ICFallDegrees.GrowthCircuitInterface
import Mathlib.Tactic

/-!
# Actual fixed-E6 finite-box common-denominator clearing

For the E6 rational section (X, (X+8)/(X-1)), multiplying every bounded
monomial X^u ((X+8)/(X-1))^v by the *same* polynomial (X-1)^(N-1)
produces X^u (X+8)^v (X-1)^(N-1-v).

This is the concrete input for the finite coefficient matrix (and for
positive-characteristic bounded-degree polynomial Wronskians). It uses the
actual RatFunc model, not an assumed model or an uninterpreted representation.
It does not yet prove the all-minor rank preservation or BoxSpreadBound.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6BoxClearing

open Polynomial
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncDenominatorIndependence
open ICFallDegrees.GrowthPrimeBoundedRationalWronskian

variable {K : Type*} [Field K]

/-- The single denominator for every exponent pair (u,v) with v<N. -/
def boxDenominator (N : ℕ) : Polynomial K :=
  (Polynomial.X - 1) ^ (N - 1)

/-- The actual integral-form E6 coefficient-polynomial column. -/
def boxNumerator (N : ℕ) (e : ℕ × ℕ) : Polynomial K :=
  Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C (8 : K)) ^ e.2 *
    (Polynomial.X - 1) ^ (N - 1 - e.2)

theorem boxDenominator_ne_zero (N : ℕ) :
    boxDenominator (K := K) N ≠ 0 := by
  have ht : (Polynomial.X - 1 : Polynomial K) ≠ 0 := by
    simpa only [Polynomial.C_1] using
      (Polynomial.X_sub_C_ne_zero (1 : K))
  exact pow_ne_zero _ ht

/-- The same denominator clears each actual E6 rational monomial.
This identity is valid over every field, including finite fields. -/
theorem box_monomial_cleared (N : ℕ) (e : ℕ × ℕ)
    (he : e.2 < N) :
    algebraMap (Polynomial K) (RatFunc K) (boxDenominator (K := K) N) *
        monomial (8 : K) e =
      algebraMap (Polynomial K) (RatFunc K) (boxNumerator (K := K) N e) := by
  let T : Polynomial K := Polynomial.X - 1
  let U : Polynomial K :=
    Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C (8 : K)) ^ e.2
  have hle : e.2 ≤ N - 1 := by omega
  have ht : (algebraMap (Polynomial K) (RatFunc K) T) ≠ 0 := by
    have hT : T ≠ 0 := by
      change (Polynomial.X - 1 : Polynomial K) ≠ 0
      simpa only [Polynomial.C_1] using
        (Polynomial.X_sub_C_ne_zero (1 : K))
    exact RatFunc.algebraMap_ne_zero hT
  change algebraMap (Polynomial K) (RatFunc K) (T ^ (N - 1)) *
      RatFunc.mk U (T ^ e.2) =
    algebraMap (Polynomial K) (RatFunc K) (U * T ^ (N - 1 - e.2))
  calc
    algebraMap (Polynomial K) (RatFunc K) (T ^ (N - 1)) *
        RatFunc.mk U (T ^ e.2) =
      (algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1) *
        (algebraMap (Polynomial K) (RatFunc K) U /
          (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) := by
          rw [RatFunc.mk_eq_div, map_pow, map_pow]
    _ = ((algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1 - e.2) *
        (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) *
        (algebraMap (Polynomial K) (RatFunc K) U /
          (algebraMap (Polynomial K) (RatFunc K) T) ^ e.2) := by
            rw [← pow_add, Nat.sub_add_cancel hle]
    _ = algebraMap (Polynomial K) (RatFunc K) U *
        (algebraMap (Polynomial K) (RatFunc K) T) ^ (N - 1 - e.2) := by
          field_simp [pow_ne_zero _ ht]
          ring
    _ = algebraMap (Polynomial K) (RatFunc K)
        (U * T ^ (N - 1 - e.2)) := by
          simp only [map_mul, map_pow, U] <;> ring

/-- Concrete polynomial columns are independent whenever their actual
rational-function monomials are independent. -/
theorem boxNumerator_linearIndependent (r N : ℕ)
    (e : Fin r → ℕ × ℕ)
    (he : ∀ j, (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j))) :
    LinearIndependent K (fun j => boxNumerator (K := K) N (e j)) := by
  exact polynomial_linearIndependent_of_clear r
    (fun j => monomial (8 : K) (e j))
    (boxDenominator (K := K) N)
    (fun j => boxNumerator (K := K) N (e j))
    (boxDenominator_ne_zero N)
    (fun j => box_monomial_cleared N (e j) (he j)) hind

/-- For any independent E6 box columns satisfying the explicit
below-characteristic degree guard, their genuine RatFunc Wronskian is
nonzero. No Wronskian or independent polynomial family is postulated. -/
theorem box_wronskian_ne_zero_of_degree_bound
    (r N k : ℕ) (e : Fin r → ℕ × ℕ)
    (he : ∀ j, (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hdeg : ∀ j, (boxNumerator (K := K) N (e j)).natDegree < k)
    (hchar : ringChar K = 0 ∨ k ≤ ringChar K) :
    (derivativeMatrix r (fun j => monomial (8 : K) (e j))).det ≠ 0 := by
  exact rational_wronskian_ne_zero_of_bounded_polynomial_representatives
    r k (fun j => monomial (8 : K) (e j))
    (boxDenominator (K := K) N)
    (fun j => boxNumerator (K := K) N (e j))
    (boxDenominator_ne_zero N)
    (fun j => box_monomial_cleared N (e j) (he j))
    (boxNumerator_linearIndependent r N e he hind)
    hdeg hchar

/-- Every actual E6 polynomial box column has degree strictly below 2N
when both exponents are in [0,N). This is a uniform all-subsets degree
bound, independent of any rank, determinant, or characteristic assumption. -/
theorem boxNumerator_natDegree_lt_two_mul (N : ℕ) (e : ℕ × ℕ)
    (hu : e.1 < N) (hv : e.2 < N) :
    (boxNumerator (K := K) N e).natDegree < 2 * N := by
  have hx : (Polynomial.X : Polynomial K).natDegree ≤ 1 :=
    Polynomial.natDegree_X_le
  have ha : (Polynomial.X + Polynomial.C (8 : K)).natDegree ≤ 1 :=
    Polynomial.natDegree_add_le_of_degree_le hx (by simp)
  have ht : (Polynomial.X - 1 : Polynomial K).natDegree ≤ 1 := by
    compute_degree
  have h0 : (Polynomial.X ^ e.1 : Polynomial K).natDegree ≤ e.1 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le e.1 hx)
  have h1 : ((Polynomial.X + Polynomial.C (8 : K)) ^ e.2).natDegree ≤ e.2 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le e.2 ha)
  have h2 : ((Polynomial.X - 1 : Polynomial K) ^ (N - 1 - e.2)).natDegree ≤
      N - 1 - e.2 := by
    simpa only [mul_one] using
      (Polynomial.natDegree_pow_le_of_le (N - 1 - e.2) ht)
  have hleft :
      (Polynomial.X ^ e.1 * (Polynomial.X + Polynomial.C (8 : K)) ^ e.2).natDegree ≤
        e.1 + e.2 :=
    Polynomial.natDegree_mul_le_of_le h0 h1
  have htot :
      (boxNumerator (K := K) N e).natDegree ≤
        (e.1 + e.2) + (N - 1 - e.2) := by
    change (Polynomial.X ^ e.1 *
      (Polynomial.X + Polynomial.C (8 : K)) ^ e.2 *
      (Polynomial.X - 1) ^ (N - 1 - e.2)).natDegree ≤
        (e.1 + e.2) + (N - 1 - e.2)
    exact Polynomial.natDegree_mul_le_of_le hleft h2
  omega

/-- **Unconditional within the finite exponent box:** a linearly
independent family of actual E6 rational monomials has a nonzero genuine
Wronskian whenever the characteristic is zero or at least 2N. In
particular this applies over large finite prime fields and eliminates
the polynomial-representative/degree premises in the preceding theorem. -/
theorem box_wronskian_ne_zero_of_box_char_bound
    (r N : ℕ) (e : Fin r → ℕ × ℕ)
    (hbox : ∀ j, (e j).1 < N ∧ (e j).2 < N)
    (hind : LinearIndependent K (fun j => monomial (8 : K) (e j)))
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K) :
    (derivativeMatrix r (fun j => monomial (8 : K) (e j))).det ≠ 0 := by
  exact box_wronskian_ne_zero_of_degree_bound r N (2 * N) e
    (fun j => (hbox j).2) hind
    (fun j => boxNumerator_natDegree_lt_two_mul N (e j)
      (hbox j).1 (hbox j).2) hchar

#print axioms boxNumerator_natDegree_lt_two_mul
#print axioms box_wronskian_ne_zero_of_box_char_bound

#print axioms boxDenominator_ne_zero
#print axioms box_monomial_cleared
#print axioms boxNumerator_linearIndependent
#print axioms box_wronskian_ne_zero_of_degree_bound

end ICFallDegrees.GrowthPrimeE6BoxClearing
