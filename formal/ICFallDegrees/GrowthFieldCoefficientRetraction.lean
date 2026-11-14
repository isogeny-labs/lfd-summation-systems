import ICFallDegrees.GrowthIntegralBezoutSpecialization
import Mathlib.LinearAlgebra.Basis.VectorSpace
import Mathlib.Tactic

/-!
# N8: coefficientwise rational descent from a field extension

For any extension of fields K → L, a K-linear retraction exists because
K → L is injective and every vector subspace has a complement. Applying
it coefficientwise gives an explicit K-linear projection of multivariate
polynomials, compatible with multiplication by lifted K-polynomials.

This is *not* a ring homomorphism and is not asserted to preserve products
of arbitrary L-polynomials. Its special compatibility is sufficient for
descending a finite Bézout expression in K-polynomial generators.
-/

noncomputable section
namespace ICFallDegrees.GrowthFieldCoefficientRetraction

open MvPolynomial

variable {K L σ : Type*} [Field K] [Field L] [Algebra K L]

/-- A K-linear retraction of a field extension exists by choosing a
nonzero coordinate of 1 in a vector-space basis over K. -/
private theorem exists_scalarRetraction :
    ∃ φ : L →ₗ[K] K, ∀ a : K, φ (algebraMap K L a) = a := by
  classical
  let b := Basis.ofVectorSpace K L
  have hi : ∃ i, b.coord i (1 : L) ≠ 0 := by
    by_contra hn
    push_neg at hn
    have hz : b.repr (1 : L) = 0 := by
      ext i
      exact hn i
    exact one_ne_zero (b.repr.injective (by simpa using hz))
  obtain ⟨i, hnonzero⟩ := hi
  let c : K := b.coord i (1 : L)
  let φ : L →ₗ[K] K :=
    { toFun := fun x => b.coord i x / c
      map_add' := by
        intro x y
        simp only [map_add, add_div]
      map_smul' := by
        intro a x
        change b.coord i (a • x) / c = a • (b.coord i x / c)
        rw [map_smul]
        simp only [smul_eq_mul, mul_div_assoc] }
  refine ⟨φ, ?_⟩
  intro a
  have ha : (algebraMap K L a) = a • (1 : L) := by
    rw [Algebra.smul_def, mul_one]
  change b.coord i (algebraMap K L a) / c = a
  rw [ha, map_smul]
  simp only [smul_eq_mul]
  exact mul_div_cancel_right₀ a hnonzero

/-- A constructed K-linear left inverse of the inclusion K → L. -/
def scalarRetraction : L →ₗ[K] K :=
  Classical.choose (exists_scalarRetraction (K := K) (L := L))

/-- The selected linear retraction fixes every base-field scalar. -/
theorem scalarRetraction_apply (a : K) :
    (scalarRetraction (K := K) (L := L)) ((algebraMap K L) a) = a :=
  (Classical.choose_spec (exists_scalarRetraction (K := K) (L := L))) a

/-- Apply an arbitrary K-linear functional to polynomial coefficients. -/
def coefficientProjection (φ : L →ₗ[K] K) :
    MvPolynomial σ L →ₗ[K] MvPolynomial σ K where
  toFun p := Finsupp.mapRange φ φ.map_zero p
  map_add' p q := by
    ext e
    change φ (coeff e (p + q)) = φ (coeff e p) + φ (coeff e q)
    rw [coeff_add, map_add]
  map_smul' c p := by
    ext e
    change φ (coeff e (c • p)) = c • φ (coeff e p)
    rw [coeff_smul, map_smul]

/-- A coefficient projection commutes with multiplying by a variable. -/
theorem coefficientProjection_mul_X (φ : L →ₗ[K] K)
    (p : MvPolynomial σ L) (v : σ) :
    coefficientProjection φ (p * X v) = coefficientProjection φ p * X v := by
  classical
  ext e
  change φ (coeff e (p * X v)) =
    coeff e (coefficientProjection φ p * X v)
  rw [coeff_mul_X', coeff_mul_X']
  split_ifs
  · rfl
  · exact map_zero _

/-- Multiplication by a lifted base-field polynomial commutes with
coefficient projection, even though the projection is not multiplicative. -/
theorem coefficientProjection_mul_lift (φ : L →ₗ[K] K)
    (a : MvPolynomial σ L) (p : MvPolynomial σ K) :
    coefficientProjection φ (a * map (algebraMap K L) p) =
      coefficientProjection φ a * p := by
  induction p using MvPolynomial.induction_on with
  | C c =>
      have scalar : a * C ((algebraMap K L) c) = c • a := by
        ext e
        rw [mul_comm a (C _), coeff_C_mul, coeff_smul, Algebra.smul_def]
      simpa only [map_C, scalar, map_smul, smul_eq_C_mul, mul_comm] using
        (coefficientProjection φ).map_smul c a
  | add p q ihp ihq =>
      simp only [map_add, mul_add, map_add, ihp, ihq]
  | mul_X p i ihp =>
      simp only [map_mul, map_X, ← mul_assoc,
        coefficientProjection_mul_X, ihp]

/-- With a retraction fixing K, lifting and projecting a polynomial is
the identity. -/
theorem coefficientProjection_lift (φ : L →ₗ[K] K)
    (hφ : ∀ c : K, φ ((algebraMap K L) c) = c)
    (p : MvPolynomial σ K) :
    coefficientProjection φ (map (algebraMap K L) p) = p := by
  classical
  ext e
  change φ (coeff e (map (algebraMap K L) p)) = coeff e p
  rw [coeff_map]
  exact hφ _

/-- Fully unconditional scalar-extension/retraction identity. -/
theorem rationalRetraction_lift (p : MvPolynomial σ K) :
    coefficientProjection (scalarRetraction (K := K) (L := L))
      (map (algebraMap K L) p) = p := by
  exact coefficientProjection_lift
    (scalarRetraction (K := K) (L := L))
    (fun c => scalarRetraction_apply (K := K) (L := L) c) p

#print axioms scalarRetraction_apply
#print axioms coefficientProjection_mul_lift
#print axioms rationalRetraction_lift

end ICFallDegrees.GrowthFieldCoefficientRetraction
