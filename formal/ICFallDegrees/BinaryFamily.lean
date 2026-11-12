import ICFallDegrees.BooleanRaw
import Mathlib.Algebra.CharP.Defs

noncomputable section

/-!
# A parameterized characteristic-two summation-polynomial family

For E: Y²+XY=X³+aX²+u⁴, target (0,u²), coordinate base {0,u},
the raw S3 is u⁴(x²y²+1). Its explicit Boolean normalization is u⁴(xy+1).
Coefficient descent through any jointly injective family of linear coordinate
maps gives exact raw LFD 4 and normalized LFD 3. This includes the Koblitz
specialization u=1, a∈{0,1}. The group law and rational liftability are not
identified with the polynomial zero set.
-/

namespace ICFallDegrees.BinaryFamily

open MvPolynomial

variable {k K ι : Type*} [Field k] [Field K] [Algebra k K]

def coefficientMap (l : K →ₗ[k] k) (p : BooleanPair.P K) : BooleanPair.P k :=
  Finsupp.mapRange l l.map_zero p

theorem coefficientMap_scaled_monomial (l : K →ₗ[k] k) (a : K)
    (e : Fin 2 →₀ ℕ) :
    coefficientMap l (C a*(monomial e 1-1))=C (l a)*(monomial e 1-1) := by
  classical
  ext d
  change l (coeff d (C a*(monomial e 1-1)))=coeff d (C (l a)*(monomial e 1-1))
  simp only [coeff_C_mul,coeff_sub,coeff_monomial,coeff_one]
  split_ifs <;> simp

theorem coefficientMap_scaled_raw (l : K →ₗ[k] k) (a : K) :
    coefficientMap l (C a*BooleanPair.raw)=C (l a)*BooleanPair.raw := by
  simp only [BooleanPair.raw_monomial]
  exact coefficientMap_scaled_monomial l a _

theorem coefficientMap_scaled_normalized (l : K →ₗ[k] k) (a : K) :
    coefficientMap l (C a*BooleanPair.h)=C (l a)*BooleanPair.h := by
  simp only [BooleanPair.h,BooleanPair.x,BooleanPair.y,X,monomial_mul,mul_one]
  exact coefficientMap_scaled_monomial l a _

def s3 (b t : K) (x1 x2 : BooleanPair.P K) : BooleanPair.P K :=
  (x1*x2)^2+(x1*C t)^2+(x2*C t)^2+x1*x2*C t+C b

theorem target_on_curve (a u : K) :
    (u^2)^2+0*u^2=0^3+a*0^2+u^4 := by ring

theorem curve_b_nonzero (u : K) (hu : u≠0) : u^4≠0 := pow_ne_zero 4 hu

variable [CharP K 2]

theorem raw_s3_identity (u : K) :
    s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y)=
      C (u^4)*BooleanPair.raw := by
  have identity : s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y)-
      C (u^4)*BooleanPair.raw=C ((2 : K)*u^4) := by
    simp only [s3,BooleanPair.raw,map_zero,map_mul,map_ofNat,map_pow]
    ring
  have two : (2 : K)=0 := CharP.cast_eq_zero K 2
  have hz : s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y)-
      C (u^4)*BooleanPair.raw=0 := by rw [identity,two,zero_mul,map_zero]
  exact sub_eq_zero.mp hz

theorem normalization_identity (u : K) :
    s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y)-C (u^4)*BooleanPair.h=
      C (u^4)*(BooleanPair.y^2*BooleanPair.cx+BooleanPair.x*BooleanPair.cy) := by
  rw [raw_s3_identity]
  simp only [BooleanPair.raw,BooleanPair.h,BooleanPair.cx,BooleanPair.cy]
  ring

def rawDescendedInputs (l : ι → K →ₗ[k] k) (u : K) : Set (BooleanPair.P k) :=
  fun f => f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i,
    f=coefficientMap (l i) (s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y))

def normalizedDescendedInputs (l : ι → K →ₗ[k] k) (u : K) : Set (BooleanPair.P k) :=
  fun f => f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i,
    f=coefficientMap (l i) (C (u^4)*BooleanPair.h)

theorem raw_descent_identity (l : ι → K →ₗ[k] k) (u : K) :
    rawDescendedInputs l u=BooleanPair.rawInputs (fun i => l i (u^4)) := by
  funext f
  change (f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i,
    f=coefficientMap (l i) (s3 (u^4) 0 (C u*BooleanPair.x) (C u*BooleanPair.y))) =
      (f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i, f=C (l i (u^4))*BooleanPair.raw)
  rw [raw_s3_identity]
  simp only [coefficientMap_scaled_raw]

omit [CharP K 2] in
theorem normalized_descent_identity (l : ι → K →ₗ[k] k) (u : K) :
    normalizedDescendedInputs l u=BooleanPair.inputs (fun i => l i (u^4)) := by
  funext f
  change (f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i,
    f=coefficientMap (l i) (C (u^4)*BooleanPair.h)) =
      (f=BooleanPair.cx ∨ f=BooleanPair.cy ∨ ∃ i, f=C (l i (u^4))*BooleanPair.h)
  simp only [coefficientMap_scaled_normalized]

omit [CharP K 2] in
theorem nonzero_coordinate (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a))
    (u : K) (hu : u≠0) : ∃ i, l i (u^4)≠0 := by
  classical
  by_contra no
  have allzero : ∀ i, l i (u^4)=0 := by
    intro i
    by_contra hi
    exact no ⟨i,hi⟩
  have coordinates : (fun i => l i (u^4))=(fun i => l i 0) := by
    funext i
    simp [allzero i]
  exact pow_ne_zero 4 hu (injective coordinates)

theorem normalized_lfd_three (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a))
    (u : K) (hu : u≠0) : IsLastFallDegree (normalizedDescendedInputs l u) 3 := by
  rw [normalized_descent_identity]
  exact BooleanPair.exact_lfd_three _ (nonzero_coordinate l injective u hu)

theorem normalized_first_strict_fall_three (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a))
    (u : K) (hu : u≠0) : IsFirstStrictFall (normalizedDescendedInputs l u) 3 := by
  rw [normalized_descent_identity]
  exact BooleanPair.first_strict_fall_three _ (nonzero_coordinate l injective u hu)

theorem raw_lfd_four (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a))
    (u : K) (hu : u≠0) : IsLastFallDegree (rawDescendedInputs l u) 4 := by
  rw [raw_descent_identity]
  exact BooleanPair.raw_exact_lfd_four _ (nonzero_coordinate l injective u hu)

theorem presentations_same_ideal (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a))
    (u : K) (hu : u≠0) :
    Ideal.span (normalizedDescendedInputs l u)=Ideal.span (rawDescendedInputs l u) := by
  rw [normalized_descent_identity,raw_descent_identity]
  exact BooleanPair.same_ideal _ (nonzero_coordinate l injective u hu)

theorem normalized_ideal_proper (l : ι → K →ₗ[k] k) (u : K) :
    (1 : BooleanPair.P k)∉Ideal.span (normalizedDescendedInputs l u) := by
  rw [normalized_descent_identity]
  exact BooleanPair.proper_ideal _

omit [CharP K 2] in
/-- Squaring invariance is an extra condition on the factor base. -/
theorem frobenius_invariant_iff (u : K) (hu : u≠0) :
    (∀ z ∈ ({0,u} : Set K), z^2∈({0,u} : Set K)) ↔ u=1 := by
  constructor
  · intro invariant
    have hs := invariant u (by simp)
    simp only [Set.mem_insert_iff,Set.mem_singleton_iff] at hs
    rcases hs with zero | fixed
    · exact False.elim (pow_ne_zero 2 hu zero)
    · have product : u*(u-1)=0 := by calc
        _ = u^2-u := by ring
        _ = 0 := sub_eq_zero.mpr fixed
      rcases mul_eq_zero.mp product with zero | one
      · exact False.elim (hu zero)
      · exact sub_eq_zero.mp one
  · rintro rfl z hz
    simp only [Set.mem_insert_iff,Set.mem_singleton_iff] at hz
    rcases hz with rfl | rfl <;> simp

/-- The b=1 specialization includes both Koblitz choices a=0 and a=1. -/
theorem koblitz_lfd_three (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a)) :
    IsLastFallDegree (normalizedDescendedInputs l (1 : K)) 3 :=
  normalized_lfd_three l injective 1 one_ne_zero

theorem koblitz_raw_lfd_four (l : ι → K →ₗ[k] k)
    (injective : Function.Injective (fun a => fun i => l i a)) :
    IsLastFallDegree (rawDescendedInputs l (1 : K)) 4 :=
  raw_lfd_four l injective 1 one_ne_zero

end ICFallDegrees.BinaryFamily
