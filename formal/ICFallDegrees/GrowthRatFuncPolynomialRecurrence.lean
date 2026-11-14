import ICFallDegrees.GrowthRatFuncScaledRecurrence
import ICFallDegrees.GrowthRatFuncLogDerivative
import ICFallDegrees.GrowthWronskianTightDegree
import Mathlib.Tactic

/-!
# N8: the polynomial numerator recurrence equals the genuine derivative recurrence

Combines the kernel-proved arbitrary-order differential identity on RatFunc K
with the previously verified polynomial recurrence. Crucially, the identities
below refer to actual iterates of the quotient-rule rational derivative:
no symbolic differentiation convention or Wronskian assumption is introduced.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncPolynomialRecurrence

open Polynomial
open ICFallDegrees.GrowthComplexRatFuncDerivative
open ICFallDegrees.GrowthRatFuncIteratedLinear
open ICFallDegrees.GrowthRatFuncScaledRecurrence
open ICFallDegrees.GrowthRationalWronskianRec
open ICFallDegrees.GrowthRationalLogNumerator

variable {K : Type*} [Field K]

/-- The polynomial numerator recursion is preserved under the actual
coefficient embedding K[X] → K(X), because the actual derivation extends
the formal polynomial derivative. -/
theorem rationalNumeratorRec_algebraMap (D Q : Polynomial K) (n : ℕ) :
    rationalNumeratorRec
        (algebraMap (Polynomial K) (RatFunc K) D)
        (algebraMap (Polynomial K) (RatFunc K) Q) n =
      algebraMap (Polynomial K) (RatFunc K) (numeratorRec D Q n) := by
  let φ : Polynomial K →+* RatFunc K := algebraMap (Polynomial K) (RatFunc K)
  change rationalNumeratorRec (φ D) (φ Q) n = φ (numeratorRec D Q n)
  induction n with
  | zero =>
      simp [rationalNumeratorRec, numeratorRec]
  | succ n ih =>
      have hn : (n : RatFunc K) = φ (Polynomial.C (n : K)) := by
        simp [φ]
      calc
        rationalNumeratorRec (φ D) (φ Q) (n + 1) =
            φ D * rationalDerivative (φ (numeratorRec D Q n)) +
              (φ Q - (n : RatFunc K) * rationalDerivative (φ D)) *
                φ (numeratorRec D Q n) := by
                  rw [rationalNumeratorRec_succ, ih]
        _ = φ (numeratorRec D Q (n + 1)) := by
          rw [rationalDerivative_algebraMap, rationalDerivative_algebraMap, hn]
          simp only [numeratorRec, map_add, map_mul, map_sub]
          rfl

/-- Actual higher rational derivatives of each two-torsion monomial
satisfy the *same* polynomial recurrence whose degree was audited
in GrowthWronskianTightDegree. This is an equality in K(X), not an
assumed differential bridge. -/
theorem concrete_higher_derivative (c : K) (u v n : ℕ) :
    iterDerivative n
        (RatFunc.mk
          (X ^ u * (X + C c) ^ v)
          ((X - 1 : Polynomial K) ^ v)) *
        (algebraMap (Polynomial K) (RatFunc K) (denominator c)) ^ n =
      RatFunc.mk
          (X ^ u * (X + C c) ^ v)
          ((X - 1 : Polynomial K) ^ v) *
        algebraMap (Polynomial K) (RatFunc K)
          (numeratorRec (denominator c) (logNumerator c u v) n) := by
  let f : RatFunc K := RatFunc.mk
    (X ^ u * (X + C c) ^ v) ((X - 1 : Polynomial K) ^ v)
  let d : RatFunc K := algebraMap (Polynomial K) (RatFunc K) (denominator c)
  let q : RatFunc K := algebraMap (Polynomial K) (RatFunc K) (logNumerator c u v)
  have hlog : rationalDerivative f * d = f * q :=
    GrowthRatFuncLogDerivative.concrete_rational_derivative c u v
  have h := iterDerivative_scaled f d q hlog n
  rw [rationalNumeratorRec_algebraMap] at h
  exact h

/-- Quantitative polynomial degree of each *actual* derivative numerator. -/
theorem concrete_higher_derivative_numerator_degree (c : K) (u v n : ℕ) :
    (numeratorRec (denominator c) (logNumerator c u v) n).natDegree ≤ 2 * n :=
  ICFallDegrees.GrowthWronskianTightDegree.numeratorRec_natDegree_le_two
    (denominator c) (logNumerator c u v)
    (denominator_natDegree_le c) (logNumerator_natDegree_le c u v) n

end ICFallDegrees.GrowthRatFuncPolynomialRecurrence
