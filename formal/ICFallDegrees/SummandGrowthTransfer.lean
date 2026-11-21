import ICFallDegrees.AffineTransport
import ICFallDegrees.GrowthRationalUnit
import ICFallDegrees.GrowthRationalPrimeNative

/-!
Degree-controlled specialization for growth lower bounds.

These transport lemmas do not assert that a Semaev specialization satisfies
their input hypotheses, or that its source ideal is the unit ideal.
The preprocessing budget is explicit and may exceed the source budget.
-/
noncomputable section
namespace ICFallDegrees.SummandGrowthTransfer
open MvPolynomial
variable {σ τ K : Type*} [CommRing K]

theorem map_closure_budget
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D E : ℕ) (hDE : D ≤ E)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G E (φ f))
    {f} (hf : InClosure F D f) : InClosure G E (φ f) := by
  induction hf with
  | input member degree => exact inputs _ member degree
  | zero => simpa using (InClosure.zero : InClosure G E 0)
  | add _ _ ih ih' => simpa only [map_add] using ih.add ih'
  | smul c _ ih => simpa only [map_smul] using ih.smul c
  | mul a _ budget ih =>
      simpa only [map_mul] using ih.mul (φ a)
        (((Nat.add_le_add (AffineTransport.degree_le φ linear a)
          (AffineTransport.degree_le φ linear _)).trans budget).trans hDE)

/-- A fixed preprocessing cost does not become a cost at every derivation step. -/
theorem not_one_of_specialization
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D C : ℕ)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G (max D C) (φ f))
    (excluded : ¬ InClosure G (max D C) 1) : ¬ InClosure F D 1 := by
  intro h
  apply excluded
  simpa using map_closure_budget φ linear F G D (max D C)
    (Nat.le_max_left D C) inputs h

/-- Source unit-ideal status is a separate hypothesis, not implied by specialization. -/
theorem no_lfd_bound_of_specialization
    (φ : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K))
    (D C : ℕ)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → InClosure G (max D C) (φ f))
    (excluded : ¬ InClosure G (max D C) 1)
    (unit : (1 : MvPolynomial σ K) ∈ Ideal.span F) : ¬ HasLFDBound F D := by
  intro bound
  exact not_one_of_specialization φ linear F G D C inputs excluded
    (one_of_hasLFDBound unit bound)

/-- Explicit prime-field transfer at arbitrary depth, conditional on the actual
specialized inputs. This makes no higher-summand Semaev identity assumption implicit. -/
theorem prime_growth_transfer {σ : Type*}
    (p k : ℕ) [Fact p.Prime] (depth : 2 ≤ k) (a : ZMod p)
    (F : Set (MvPolynomial σ (ZMod p)))
    (φ : MvPolynomial σ (ZMod p) →ₐ[ZMod p]
      MvPolynomial (Fin (2*k)) (ZMod p))
    (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (D C : ℕ) (positive : 1 ≤ max D C)
    (hk : 14 * max D C + 10 ≤ k) (hp : 2 * 2^k + 1 ≤ p)
    (ha : a + 3 ≠ 0)
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D →
      InClosure (Set.range (PrimeM2General.inputs k depth a (-1-a) 1))
        (max D C) (φ f)) : ¬ InClosure F D 1 := by
  apply not_one_of_specialization φ linear F
    (Set.range (PrimeM2General.inputs k depth a (-1-a) 1)) D C inputs
  exact GrowthRationalPrimeNative.not_one_of_large_prime p k a depth
    (max D C) positive hk hp ha

/-- The determinant formula for the formal-degree-(2,2) Sylvester resultant. -/
def quadraticSylvester (a b c d e f : K) : K :=
  (a*f-c*d)^2-(a*e-b*d)*(b*f-c*e)

theorem quadraticSylvester_constant (a b c f : K) :
    quadraticSylvester a b c 0 0 f = a^2*f^2 := by
  unfold quadraticSylvester
  ring

/-- Specialize the universal quadratic Sylvester formula, retaining degree two.
This is a formula identity; it does not redefine resultant at the reduced degree. -/
theorem order_two_pair_specialization (a x y : K) :
    quadraticSylvester ((x-y)^2)
      (-2*((x+y)*(x*y+a)+2*(-1-a)))
      ((x*y-a)^2-4*(-1-a)*(x+y))
      ((1-1)^2) (-2*((1+1)*(1*1+a)+2*(-1-a)))
      ((1*1-a)^2-4*(-1-a)*(1+1)) =
      (a+3)^4*(x-y)^4 := by
  unfold quadraticSylvester
  ring

#print axioms map_closure_budget
#print axioms not_one_of_specialization
#print axioms no_lfd_bound_of_specialization
#print axioms prime_growth_transfer
#print axioms order_two_pair_specialization
end ICFallDegrees.SummandGrowthTransfer
