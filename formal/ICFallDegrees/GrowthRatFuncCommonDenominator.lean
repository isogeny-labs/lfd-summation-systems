import ICFallDegrees.GrowthRatFuncPolynomialWronskianMap
import Mathlib.FieldTheory.RatFunc.Basic
import Mathlib.Tactic

/-!
# N8: common polynomial denominator for arbitrary rational Wronskians

For an arbitrary finite family in K(X), clear all denominators with
the product of the actual RatFunc denominators. This constructs a
single nonzero polynomial H and polynomial numerators P_j with

  algebraMap H * f_j = algebraMap P_j.

Combined with the previously proved genuine rational common-scalar
Wronskian identity and polynomial embedding identity, this produces
an unconditional exact polynomial determinant after multiplication
by (algebraMap H)^r.

The independent-polynomial Wronskian nonvanishing theorem, cofactor
ratio, local valuations, and SpreadBound remain open.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncCommonDenominator

open Polynomial
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthRatFuncWronskianCommonScalar
open ICFallDegrees.GrowthRatFuncPolynomialWronskianMap

variable {K : Type*} [Field K]

/-- The product of all actual reduced denominators in a rational family. -/
def commonDenominator (r : ℕ) (f : Fin r → RatFunc K) : Polynomial K :=
  ∏ j : Fin r, (f j).denom

/-- The product denominator is nonzero, even for an empty family. -/
theorem commonDenominator_ne_zero (r : ℕ) (f : Fin r → RatFunc K) :
    commonDenominator r f ≠ 0 := by
  classical
  unfold commonDenominator
  exact Finset.prod_ne_zero_iff.mpr (fun j _ => RatFunc.denom_ne_zero (f j))

/-- Every member's actual denominator divides the common polynomial. -/
theorem denom_dvd_commonDenominator (r : ℕ)
    (f : Fin r → RatFunc K) (j : Fin r) :
    (f j).denom ∣ commonDenominator r f := by
  classical
  unfold commonDenominator
  exact Finset.dvd_prod_of_mem (f := fun i : Fin r => (f i).denom) (Finset.mem_univ j)

/-- Clear all rational denominators simultaneously, without assuming
that the input functions are independent or nonzero. -/
theorem exists_polynomial_numerators (r : ℕ) (f : Fin r → RatFunc K) :
    ∃ p : Fin r → Polynomial K, ∀ j,
      algebraMap (Polynomial K) (RatFunc K) (commonDenominator r f) * f j =
      algebraMap (Polynomial K) (RatFunc K) (p j) := by
  classical
  have hdvd : ∀ j : Fin r, (f j).denom ∣ commonDenominator r f :=
    denom_dvd_commonDenominator r f
  choose q hq using hdvd
  refine ⟨fun j => (f j).num * q j, ?_⟩
  intro j
  have hden : algebraMap (Polynomial K) (RatFunc K) (f j).denom ≠ 0 :=
    RatFunc.algebraMap_ne_zero (RatFunc.denom_ne_zero (f j))
  have hnum :
      algebraMap (Polynomial K) (RatFunc K) (f j).num =
        f j * algebraMap (Polynomial K) (RatFunc K) (f j).denom :=
    (div_eq_iff hden).mp (RatFunc.num_div_denom (f j))
  rw [hq j, map_mul, map_mul, hnum]
  ring

/-- An all-size genuine rational Wronskian becomes exactly a
polynomial Wronskian after clearing one single common denominator.
There is no hidden Wronskian nonvanishing or SpreadBound premise. -/
theorem exists_polynomial_wronskian_scale (r : ℕ)
    (f : Fin r → RatFunc K) :
    ∃ (H : Polynomial K) (p : Fin r → Polynomial K),
      H ≠ 0 ∧
      (algebraMap (Polynomial K) (RatFunc K) H)^r *
        (derivativeMatrix r f).det =
        algebraMap (Polynomial K) (RatFunc K)
          (polynomialDerivativeMatrix r p).det := by
  classical
  obtain ⟨p, hp⟩ := exists_polynomial_numerators r f
  refine ⟨commonDenominator r f, p, commonDenominator_ne_zero r f, ?_⟩
  have heq :
      (fun j => algebraMap (Polynomial K) (RatFunc K)
          (commonDenominator r f) * f j) =
        (fun j => algebraMap (Polynomial K) (RatFunc K) (p j)) := by
    funext j
    exact hp j
  have hdet := determinant_mul_scalar r
    (algebraMap (Polynomial K) (RatFunc K) (commonDenominator r f)) f
  rw [heq, determinant_polynomialMap] at hdet
  exact hdet.symm

end ICFallDegrees.GrowthRatFuncCommonDenominator
