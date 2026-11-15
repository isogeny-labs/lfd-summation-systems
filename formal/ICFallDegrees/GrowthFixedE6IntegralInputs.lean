import ICFallDegrees.PrimeM2General
import Mathlib.Tactic

/-!
# An actual integral input family for the fixed E6 N8 normalized chain

The existing `PrimeM2General.inputs` is field-typed. Here we construct the
same finite family with integer coefficients, and prove that its map to Q
is exactly the original normalized family, not merely a related S3 system.
No unit ideal or characteristic-p growth conclusion is asserted here.
-/

noncomputable section
namespace ICFallDegrees.GrowthFixedE6IntegralInputs

open MvPolynomial
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GeneralQuadraticTop

/-- The exact integral chain-tail term, including block endpoints. -/
def chainTailInt (k : ℕ) (i : Fin (2 * k)) :
    MvPolynomial (Fin (2 * k)) ℤ :=
  if h : i.val + 1 = k ∨ i.val + 1 = 2 * k then -1
  else -X ⟨i.val + 1, by have := i.isLt; omega⟩

/-- The actual fixed E6 normalized quadratic:
    quadratic cross-term coefficients (-14,-2);
    linear coefficients (1,16); and constant 64. -/
def normalizedInt (k : ℕ) (depth : 2 ≤ k) :
    MvPolynomial (Fin (2 * k)) ℤ :=
  top (active k depth) (-14 : ℤ) (-2 : ℤ) +
    C (1 : ℤ) * (X (active k depth 1) + X (active k depth 3)) +
    C (16 : ℤ) * (X (active k depth 0) + X (active k depth 2)) +
    C (64 : ℤ)

/-- Every coefficient is integral; no rational polynomial lift is assumed. -/
def e6InputsInt (k : ℕ) (depth : 2 ≤ k) :
    Fin (2 * k + 1) → MvPolynomial (Fin (2 * k)) ℤ :=
  Fin.lastCases (normalizedInt k depth)
    (fun i => X i ^ 2 + chainTailInt k i)

/-- The integer block endpoint conventions are the same over Q. -/
theorem map_chain_tail (k : ℕ) (i : Fin (2 * k)) :
    (chainTailInt k i).map (Int.castRingHom ℚ) =
      PrimeM2General.chainTail (K := ℚ) k i := by
  unfold chainTailInt PrimeM2General.chainTail
  split_ifs <;> simp

/-- Mapping the integer E6 quadratic to Q recovers the actual fixed E6
normalized polynomial, with t=1, a=6 and b=-7. -/
theorem map_normalized (k : ℕ) (depth : 2 ≤ k) :
    (normalizedInt k depth).map (Int.castRingHom ℚ) =
      PrimeM2General.normalized k depth (6 : ℚ) (-7 : ℚ) 1 := by
  unfold normalizedInt PrimeM2General.normalized PrimeM2General.tail
    GeneralQuadraticTop.top
  simp only [map_add, map_mul,
    MvPolynomial.map_C, MvPolynomial.map_X,
    map_sub, map_neg]
  simp only [map_ofNat, map_natCast, map_intCast] 
  norm_num <;> simp only [map_ofNat] <;> ring

/-- The *exact* E6 N8 normalized input list has an integral model. -/
theorem map_e6_inputs (k : ℕ) (depth : 2 ≤ k) (i : Fin (2 * k + 1)) :
    (e6InputsInt k depth i).map (Int.castRingHom ℚ) =
      PrimeM2General.inputs k depth (6 : ℚ) (-7 : ℚ) 1 i := by
  refine Fin.lastCases ?_ (fun i => ?_) i
  · simpa only [e6InputsInt, PrimeM2General.inputs, Fin.lastCases_last]
      using map_normalized k depth
  · simpa only [e6InputsInt, PrimeM2General.inputs,
      Fin.lastCases_castSucc, map_add,
      map_pow, MvPolynomial.map_X] using
      congrArg (fun f : MvPolynomial (Fin (2 * k)) ℚ =>
        (X i : MvPolynomial (Fin (2 * k)) ℚ) ^ 2 + f)
        (map_chain_tail k i)

#print axioms map_chain_tail
#print axioms map_normalized
#print axioms map_e6_inputs

end ICFallDegrees.GrowthFixedE6IntegralInputs
