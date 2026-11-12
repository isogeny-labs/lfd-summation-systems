import ICFallDegrees.ChainQuotient
import ICFallDegrees.BinaryS3Bounds

noncomputable section

/-!
# Unified last-fall-degree bounds for the three curve families

One kernel-checked statement per presentation family, each an instance of
the finite-quotient theorem with an explicit reduced-ideal-trivial subsystem:

* prime curves, `m` summands, subgroup chains of depth `k`
  (`ChainQuotient.chain_lfd_bound`): `LFD ≤ max (m*k+1) delta`;
* prime curves, two summands, the normalized S3 presentation at every depth
  (`ChainQuotient.prime_field_m2_lfd_bound`): `LFD ≤ 2*k+1`;
* binary and Koblitz curves, descended left-chain S3 systems on `r` Boolean
  coordinates (`BinaryS3Bounds.normalized_split_lfd_bound`):
  `LFD ≤ max (r+1) 3`, hence `LFD ≤ m*k+1` for `m` summands on coordinate
  spaces of dimension `k` and `m*k ≥ 2`.

Koblitz curves are the binary curves with `a ∈ {0,1}` and `b = 1`; the bound
is the same because it depends only on the Boolean field equations and the
input degrees, not on the curve coefficients or on the Frobenius structure
of the coordinate space.  Sharpness of the Boolean `r+1` bound as an
abstract statement is `BooleanSharp.exact_lfd`; the exact values of actual
curve systems are measured, not proved, and are typically far below.
-/

namespace ICFallDegrees.CurveBounds

open MvPolynomial

variable {K : Type*} [Field K]

/-- Prime curves, `m` summands of depth `k`, any further inputs of degree
at most `delta`: `LFD ≤ max (m*k+1) delta`. -/
theorem prime_chain_lfd_bound (m k : ℕ) (two_ne : (2 : K)≠0)
    (F : Set (MvPolynomial (Fin (m*k)) K)) (included : ChainQuotient.chainInputs K m k⊆F)
    (delta : ℕ) (input_degree : ∀ f ∈ F, f.totalDegree≤delta) :
    HasLFDBound F (max (m*k+1) delta) :=
  ChainQuotient.chain_lfd_bound m k two_ne F included delta input_degree

/-- Prime curves, two summands, normalized S3 at every depth `k ≥ 2`, prime
`p ≥ 5`, every curve `y² = x³ + a x + b` and target abscissa `t`. -/
theorem prime_m2_lfd_bound (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (a b t : ZMod p) :
    HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (2*k+1) :=
  ChainQuotient.prime_field_m2_lfd_bound p large k depth a b t

/-- Binary curves: descended normalized left-chain S3 system on `r` Boolean
coordinates with `r ≥ 2` has `LFD ≤ r+1`. For `m` summands on coordinate
spaces of dimension `k` (and `(m-2)*n` auxiliary coordinates) `r` is the
total number of Boolean coordinates. -/
theorem binary_descent_lfd_bound [Algebra (ZMod 2) K] {n r s : ℕ} (basis : Basis (Fin n) (ZMod 2) K)
    (constants : Fin s → Fin 3 → K) (coefficients : Fin s → Fin 3 → Fin r → K) (b : K)
    (dimension : 2≤r) :
    HasLFDBound
      (BinaryS3Bounds.descendedInputs basis (BinaryS3Bounds.normalizedConstraint constants coefficients b))
      (r+1) := by
  have bound := BinaryS3Bounds.normalized_split_lfd_bound basis constants coefficients b
  have collapse : max (r+1) 3=r+1 := by omega
  rwa [collapse] at bound

/-- Koblitz curves are binary curves with `b = 1`; the bound is unchanged. -/
theorem koblitz_descent_lfd_bound [Algebra (ZMod 2) K] {n r s : ℕ} (basis : Basis (Fin n) (ZMod 2) K)
    (constants : Fin s → Fin 3 → K) (coefficients : Fin s → Fin 3 → Fin r → K)
    (dimension : 2≤r) :
    HasLFDBound
      (BinaryS3Bounds.descendedInputs basis (BinaryS3Bounds.normalizedConstraint constants coefficients 1))
      (r+1) :=
  binary_descent_lfd_bound basis constants coefficients 1 dimension

/-- Prime chains, any number of summands: chain equations plus any reduced inputs. -/
theorem prime_reduced_inputs_lfd_bound (m k : ℕ) (two_ne : (2 : K)≠0)
    (R : Set (MvPolynomial (Fin (m*k)) K))
    (reduced : R⊆FiniteQuotient.reducedSpace (ChainQuotient.chainInputs K m k)) :
    HasLFDBound (ChainQuotient.chainInputs K m k∪R) (m*k+1) :=
  ChainQuotient.reduced_inputs_lfd_bound m k two_ne R reduced

end ICFallDegrees.CurveBounds
