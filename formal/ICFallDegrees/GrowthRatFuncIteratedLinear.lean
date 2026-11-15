import ICFallDegrees.GrowthRatFuncDerivativeRules
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

/-!
# N8: genuine iterated rational derivation preserves constant-linear relations

The iterates here are the *actual* quotient-rule derivative on `RatFunc K`
from `GrowthComplexRatFuncDerivative`.  Every lemma is true in any
coefficient field and is a prerequisite for the characteristic-zero
Wronskian cofactor construction.  No determinant nonvanishing,
pole identity, circuit spread, or CurveModel is assumed or claimed.
-/

noncomputable section
namespace ICFallDegrees.GrowthRatFuncIteratedLinear

open ICFallDegrees.GrowthComplexRatFuncDerivative

variable {K : Type*} [Field K]

/-- Iterating the genuine rational-function derivation. -/
def iterDerivative : ℕ → RatFunc K → RatFunc K
  | 0, f => f
  | n + 1, f => rationalDerivative (iterDerivative n f)

@[simp]
theorem iterDerivative_zero (f : RatFunc K) : iterDerivative 0 f = f := rfl

@[simp]
theorem iterDerivative_succ (n : ℕ) (f : RatFunc K) :
    iterDerivative (n + 1) f = rationalDerivative (iterDerivative n f) := rfl

/-- Iterated rational derivatives are additive. -/
theorem iterDerivative_add (n : ℕ) (f g : RatFunc K) :
    iterDerivative n (f + g) = iterDerivative n f + iterDerivative n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    change rationalDerivative (iterDerivative n (f + g)) =
      rationalDerivative (iterDerivative n f) + rationalDerivative (iterDerivative n g)
    rw [ih, rationalDerivative_add]

/-- Iterated rational derivatives vanish on zero. -/
@[simp]
theorem iterDerivative_zero_fun (n : ℕ) :
    iterDerivative n (0 : RatFunc K) = 0 := by
  induction n with
  | zero => rfl
  | succ n ih =>
    change rationalDerivative (iterDerivative n (0 : RatFunc K)) = 0
    rw [ih, rationalDerivative_zero]

/-- Iterated rational differentiation is K-linear; coefficients are constants. -/
theorem iterDerivative_smul (n : ℕ) (c : K) (f : RatFunc K) :
    iterDerivative n (c • f) = c • iterDerivative n f := by
  induction n with
  | zero => rfl
  | succ n ih =>
    change rationalDerivative (iterDerivative n (c • f)) =
      c • rationalDerivative (iterDerivative n f)
    rw [ih, rationalDerivative_smul]

/-- Every derivative row commutes with finite sums. -/
theorem iterDerivative_sum {ι : Type*} (n : ℕ) (s : Finset ι)
    (f : ι → RatFunc K) :
    iterDerivative n (∑ i ∈ s, f i) = ∑ i ∈ s, iterDerivative n (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert i s hi ih =>
    simp only [Finset.sum_insert hi, iterDerivative_add, ih]

/-- Every derivative row commutes with finite constant-coefficient sums. -/
theorem iterDerivative_const_sum {ι : Type*} (n : ℕ)
    (s : Finset ι) (c : ι → K) (f : ι → RatFunc K) :
    iterDerivative n (∑ i ∈ s, c i • f i) =
      ∑ i ∈ s, c i • iterDerivative n (f i) := by
  rw [iterDerivative_sum]
  apply Finset.sum_congr rfl
  intro i hi
  exact iterDerivative_smul n (c i) (f i)

/-- A genuine constant-linear dependence persists in every Wronskian row.
The converse (nonzero Wronskian for independent columns) is a separate
essential theorem, still unproved in Lean. -/
theorem relation_in_every_row {ι : Type*}
    (s : Finset ι) (c : ι → K) (f : ι → RatFunc K)
    (hrel : (∑ i ∈ s, c i • f i) = 0) (n : ℕ) :
    (∑ i ∈ s, c i • iterDerivative n (f i)) = 0 := by
  rw [← iterDerivative_const_sum, hrel, iterDerivative_zero_fun]

end ICFallDegrees.GrowthRatFuncIteratedLinear
