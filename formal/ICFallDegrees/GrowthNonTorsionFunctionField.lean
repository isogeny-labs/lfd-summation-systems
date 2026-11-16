import ICFallDegrees.GrowthNonTorsionRationalObstruction
import Mathlib.Algebra.Polynomial.SpecificDegree
import Mathlib.RingTheory.AdjoinRoot

noncomputable section

/-!
# The actual quadratic S3 function field

This constructs K(X)[y]/(y²-(X³+aX+b)) as a field, rather than assuming
an abstract curve model. It supplies the generic on-curve point and its
secant companion. Divisor valuations and circuit spread remain separate.
-/

namespace ICFallDegrees.GrowthNonTorsionFunctionField

open Polynomial GrowthSecant GrowthNonTorsionRationalObstruction

variable {K : Type*} [Field K]

def equation (a b : K) : Polynomial (RatFunc K) :=
  X ^ 2 - C (algebraMap (Polynomial K) (RatFunc K) (cubicPolynomial a b))

theorem equation_natDegree (a b : K) : (equation a b).natDegree = 2 := by
  exact Polynomial.natDegree_X_pow_sub_C

theorem equation_irreducible (a b : K) : Irreducible (equation a b) := by
  apply Polynomial.irreducible_of_degree_le_three_of_not_isRoot
  · simp [equation_natDegree]
  · intro y hy
    apply cubicPolynomial_nonsquare a b y
    simpa [Polynomial.IsRoot, equation, sub_eq_zero] using hy

instance equation_fact (a b : K) : Fact (Irreducible (equation a b)) :=
  ⟨equation_irreducible a b⟩

abbrev FunctionField (a b : K) := AdjoinRoot (equation a b)

theorem finrank_two (a b : K) :
    Module.finrank (RatFunc K) (FunctionField a b) = 2 := by
  calc
    _ = (AdjoinRoot.powerBasis (equation_irreducible a b).ne_zero).dim :=
      (AdjoinRoot.powerBasis (equation_irreducible a b).ne_zero).finrank
    _ = 2 := equation_natDegree a b

def constant (a b : K) : K →+* FunctionField a b :=
  (algebraMap (RatFunc K) (FunctionField a b)).comp RatFunc.C

def x (a b : K) : FunctionField a b :=
  algebraMap (RatFunc K) (FunctionField a b) RatFunc.X

def y (a b : K) : FunctionField a b := AdjoinRoot.root (equation a b)

theorem y_square (a b : K) :
    y a b ^ 2 = cubic (constant a b a) (constant a b b) (x a b) := by
  have hr := AdjoinRoot.eval₂_root (equation a b)
  simp only [equation, Polynomial.eval₂_sub, Polynomial.eval₂_pow,
    Polynomial.eval₂_X, Polynomial.eval₂_C, sub_eq_zero] at hr
  simpa [y, x, constant, cubicPolynomial, cubic, RingHom.comp_apply,
    RatFunc.algebraMap_X, RatFunc.algebraMap_C] using hr

theorem x_ne_zero (a b : K) : x a b ≠ 0 := by
  simpa only [x, map_zero] using
    (algebraMap (RatFunc K) (FunctionField a b)).injective.ne RatFunc.X_ne_zero

theorem x_ne_constant (a b t : K) : x a b ≠ constant a b t := by
  intro ht
  have he : (RatFunc.X : RatFunc K) = RatFunc.C t :=
    (algebraMap (RatFunc K) (FunctionField a b)).injective ht
  have hp : (Polynomial.X : Polynomial K) = Polynomial.C t := by
    apply RatFunc.algebraMap_injective K
    simpa only [RatFunc.algebraMap_X, RatFunc.algebraMap_C] using he
  have hc := congrArg (fun f : Polynomial K => f.coeff 1) hp
  simpa using hc

def z (a b t Y : K) : FunctionField a b :=
  secantX (constant a b a) (constant a b b) (x a b) (constant a b t)
    (y a b) (constant a b Y)

theorem s3_zero (a b t Y : K) (hR : Y ^ 2 = cubic a b t) :
    s3 (constant a b a) (constant a b b) (x a b) (z a b t Y)
      (constant a b t) = 0 := by
  apply secant_s3 _ _ _ _ _ _ (y_square a b) _ (x_ne_constant a b t)
  simpa only [cubic, map_add, map_mul, map_pow] using congrArg (constant a b) hR

theorem companionConstant_ne_zero (a b t : K) (ht : t ≠ 0) :
    companionConstant (constant a b a) (constant a b b) (x a b)
      (constant a b t) ≠ 0 := by
  let q : Polynomial K := C (t ^ 2) * X ^ 2 + C (-2*t*a-4*b) * X +
    C (a ^ 2 - 4*b*t)
  have hq2 : q.coeff 2 = t ^ 2 := by
    dsimp [q]
    simp only [coeff_add, coeff_C_mul_X_pow, coeff_C_mul_X, coeff_C,
      if_pos rfl, if_neg (by decide : (2:ℕ) ≠ 1), if_neg (by decide : (2:ℕ) ≠ 0),
      add_zero, if_true]
  have hq : q ≠ 0 := by
    intro hz
    have hc := congrArg (fun f : Polynomial K => f.coeff 2) hz
    have : t ^ 2 = 0 := by simpa only [hq2, coeff_zero] using hc
    exact (pow_ne_zero 2 ht) this
  have hm := (algebraMap (RatFunc K) (FunctionField a b)).injective.ne
    (RatFunc.algebraMap_ne_zero hq)
  have he : companionConstant (constant a b a) (constant a b b) (x a b)
      (constant a b t) = algebraMap (RatFunc K) (FunctionField a b)
        (algebraMap (Polynomial K) (RatFunc K) q) := by
    simp only [q, companionConstant, x, constant, RingHom.comp_apply,
      map_add, map_sub, map_mul, map_pow, map_neg, map_ofNat,
      RatFunc.algebraMap_X, RatFunc.algebraMap_C]
    ring
  rw [he]
  simpa only [map_zero] using hm

theorem z_ne_zero (a b t Y : K) (ht : t ≠ 0)
    (hR : Y ^ 2 = cubic a b t) : z a b t Y ≠ 0 := by
  apply GrowthNonTorsionBridge.secantX_ne_zero _ _ _ _ _ _ (y_square a b) _
    (x_ne_constant a b t) (companionConstant_ne_zero a b t ht)
  simpa only [cubic, map_add, map_mul, map_pow] using congrArg (constant a b) hR

#print axioms equation_irreducible
#print axioms finrank_two
#print axioms y_square
#print axioms x_ne_zero
#print axioms x_ne_constant
#print axioms s3_zero
#print axioms companionConstant_ne_zero
#print axioms z_ne_zero

end ICFallDegrees.GrowthNonTorsionFunctionField
