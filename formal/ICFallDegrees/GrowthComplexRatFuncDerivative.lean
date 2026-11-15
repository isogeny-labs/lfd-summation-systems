/-
Copyright (c) 2026 Alok Singh. All rights reserved.
Released under Apache 2.0 (https://www.apache.org/licenses/LICENSE-2.0).
Original author: Alok Singh, JacobianTwo/RatFuncDerivative.lean,
https://github.com/alok/jacobian-two/blob/ded8e67fb47c155f83ec9bd68af6014499fc2d61/JacobianTwo/RatFuncDerivative.lean
Modified for isogeny-labs N8 formalization: namespace and integration target only.
-/
import Mathlib.Algebra.Polynomial.Derivation
import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Mathlib.RingTheory.Polynomial.Wronskian
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# The formal derivative on a rational-function field

This file packages the quotient-rule derivative on `K(X)` as a `K`-derivation,
including a characteristic-zero proof that its kernel is the constants.
It is a genuine analytic-algebraic prerequisite of the N8 rational Wronskian
bridge, not a proof of Wronskian nonvanishing or SpreadBound.
-/

noncomputable section

open Polynomial

namespace ICFallDegrees.GrowthComplexRatFuncDerivative

variable {K : Type*} [Field K]

/-- The quotient-rule derivative on the rational-function field `K(X)`. -/
def rationalDerivative (x : RatFunc K) : RatFunc K :=
  x.liftOn'
    (fun p q => RatFunc.mk (derivative p * q - p * derivative q) (q ^ 2))
    (by
      intro p q a hq ha
      apply (RatFunc.mk_eq_mk
        (pow_ne_zero 2 (mul_ne_zero ha hq)) (pow_ne_zero 2 hq)).2
      simp only [derivative_mul]
      ring)

theorem rationalDerivative_mk (p q : K[X]) :
    rationalDerivative (RatFunc.mk p q) =
      RatFunc.mk (derivative p * q - p * derivative q) (q ^ 2) := by
  rw [rationalDerivative, RatFunc.liftOn'_mk]
  · intro r
    simp

@[simp]
theorem rationalDerivative_algebraMap (p : K[X]) :
    rationalDerivative (algebraMap K[X] (RatFunc K) p) =
      algebraMap K[X] (RatFunc K) (derivative p) := by
  rw [← RatFunc.mk_one, rationalDerivative_mk]
  simp

theorem mk_add_mk (p r : K[X]) {q s : K[X]} (hq : q ≠ 0) (hs : s ≠ 0) :
    RatFunc.mk p q + RatFunc.mk r s =
      RatFunc.mk (p * s + q * r) (q * s) := by
  simp only [RatFunc.mk_eq_div, map_add, map_mul]
  field_simp [RatFunc.algebraMap_ne_zero hq, RatFunc.algebraMap_ne_zero hs]
  ring

theorem mk_mul_mk (p r : K[X]) {q s : K[X]} (hq : q ≠ 0) (hs : s ≠ 0) :
    RatFunc.mk p q * RatFunc.mk r s = RatFunc.mk (p * r) (q * s) := by
  simp only [RatFunc.mk_eq_div, map_mul]
  field_simp [RatFunc.algebraMap_ne_zero hq, RatFunc.algebraMap_ne_zero hs]

theorem rationalDerivative_add (x y : RatFunc K) :
    rationalDerivative (x + y) = rationalDerivative x + rationalDerivative y := by
  induction x using RatFunc.induction_on' with
  | _ p q hq =>
    induction y using RatFunc.induction_on' with
    | _ r s hs =>
      rw [mk_add_mk p r hq hs,
        rationalDerivative_mk, rationalDerivative_mk, rationalDerivative_mk,
        mk_add_mk _ _ (pow_ne_zero 2 hq) (pow_ne_zero 2 hs)]
      apply (RatFunc.mk_eq_mk
        (pow_ne_zero 2 (mul_ne_zero hq hs))
        (mul_ne_zero (pow_ne_zero 2 hq) (pow_ne_zero 2 hs))).2
      simp only [derivative_add, derivative_mul]
      ring

@[simp]
theorem rationalDerivative_zero : rationalDerivative (0 : RatFunc K) = 0 := by
  rw [← map_zero (algebraMap K[X] (RatFunc K)), rationalDerivative_algebraMap,
    derivative_zero, map_zero]

theorem rationalDerivative_mul (x y : RatFunc K) :
    rationalDerivative (x * y) =
      x * rationalDerivative y + y * rationalDerivative x := by
  induction x using RatFunc.induction_on' with
  | _ p q hq =>
    induction y using RatFunc.induction_on' with
    | _ r s hs =>
      rw [mk_mul_mk p r hq hs,
        rationalDerivative_mk, rationalDerivative_mk, rationalDerivative_mk,
        mk_mul_mk p _ hq (pow_ne_zero 2 hs),
        mk_mul_mk r _ hs (pow_ne_zero 2 hq),
        mk_add_mk _ _
          (mul_ne_zero hq (pow_ne_zero 2 hs))
          (mul_ne_zero hs (pow_ne_zero 2 hq))]
      apply (RatFunc.mk_eq_mk
        (pow_ne_zero 2 (mul_ne_zero hq hs))
        (mul_ne_zero
          (mul_ne_zero hq (pow_ne_zero 2 hs))
          (mul_ne_zero hs (pow_ne_zero 2 hq)))).2
      simp only [derivative_mul]
      ring

@[simp]
theorem rationalDerivative_one : rationalDerivative (1 : RatFunc K) = 0 := by
  rw [← map_one (algebraMap K[X] (RatFunc K)), rationalDerivative_algebraMap,
    derivative_one, map_zero]

@[simp]
theorem rationalDerivative_algebraMap_base (c : K) :
    rationalDerivative (algebraMap K (RatFunc K) c) = 0 := by
  rw [IsScalarTower.algebraMap_apply K K[X] (RatFunc K) c,
    rationalDerivative_algebraMap]
  simp

theorem rationalDerivative_smul (c : K) (x : RatFunc K) :
    rationalDerivative (c • x) = c • rationalDerivative x := by
  rw [Algebra.smul_def, Algebra.smul_def, rationalDerivative_mul,
    rationalDerivative_algebraMap_base, mul_zero, add_zero]

/-- In characteristic zero, the constants of the quotient-rule derivative on
`K(X)` are exactly the scalar rational functions. -/
theorem rationalDerivative_eq_zero_iff [CharZero K] (x : RatFunc K) :
    rationalDerivative x = 0 ↔ ∃ c : K, x = RatFunc.C c := by
  constructor
  · intro hx
    have hrepr : x = RatFunc.mk x.num x.denom := by
      rw [RatFunc.mk_eq_div, RatFunc.num_div_denom]
    have hmk :
        RatFunc.mk
          (derivative x.num * x.denom - x.num * derivative x.denom)
          (x.denom ^ 2) = 0 := by
      rw [← rationalDerivative_mk, ← hrepr]
      exact hx
    have hnumerator :
        derivative x.num * x.denom - x.num * derivative x.denom = 0 := by
      have hdisj :
          algebraMap K[X] (RatFunc K)
              (derivative x.num * x.denom -
                x.num * derivative x.denom) = 0 ∨
            algebraMap K[X] (RatFunc K) (x.denom ^ 2) = 0 := by
        simpa only [RatFunc.mk_eq_div, _root_.div_eq_zero_iff] using hmk
      have hmapped := hdisj.resolve_right
        (RatFunc.algebraMap_ne_zero (pow_ne_zero 2 x.denom_ne_zero))
      apply RatFunc.algebraMap_injective K
      simpa using hmapped
    have hwronskian : wronskian x.num x.denom = 0 := by
      rw [wronskian]
      calc
        x.num * derivative x.denom - derivative x.num * x.denom =
            -(derivative x.num * x.denom -
              x.num * derivative x.denom) := by ring
        _ = 0 := by rw [hnumerator, neg_zero]
    have hderivatives :=
      x.isCoprime_num_denom.wronskian_eq_zero_iff.mp hwronskian
    refine ⟨x.num.coeff 0 / x.denom.coeff 0, ?_⟩
    calc
      x = RatFunc.mk x.num x.denom := hrepr
      _ = RatFunc.mk (Polynomial.C (x.num.coeff 0))
          (Polynomial.C (x.denom.coeff 0)) := by
            rw [Polynomial.eq_C_of_derivative_eq_zero hderivatives.1,
              Polynomial.eq_C_of_derivative_eq_zero hderivatives.2]
            simp
      _ = RatFunc.C (x.num.coeff 0 / x.denom.coeff 0) := by
            simpa only [RatFunc.mk_eq_div, RatFunc.algebraMap_C] using
              ((map_div₀ (RatFunc.C : K →+* RatFunc K)
                (x.num.coeff 0) (x.denom.coeff 0)).symm)
  · rintro ⟨c, rfl⟩
    rw [← RatFunc.algebraMap_eq_C]
    exact rationalDerivative_algebraMap_base c

/-- The formal derivative on `K(X)`, bundled as a `K`-derivation. -/
def rationalDerivativeDerivation : Derivation K (RatFunc K) (RatFunc K) where
  toFun := rationalDerivative
  map_add' := rationalDerivative_add
  map_smul' := rationalDerivative_smul
  leibniz' := rationalDerivative_mul
  map_one_eq_zero' := rationalDerivative_one

@[simp]
theorem rationalDerivativeDerivation_apply (x : RatFunc K) :
    rationalDerivativeDerivation x = rationalDerivative x := rfl

@[simp]
theorem rationalDerivativeDerivation_algebraMap (p : K[X]) :
    rationalDerivativeDerivation (algebraMap K[X] (RatFunc K) p) =
      algebraMap K[X] (RatFunc K) (derivative p) :=
  rationalDerivative_algebraMap p

end ICFallDegrees.GrowthComplexRatFuncDerivative
