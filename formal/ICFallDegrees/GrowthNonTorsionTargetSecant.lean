import ICFallDegrees.GrowthNonTorsionTargetValuation

noncomputable section

/-!
# The second pole identity for the S3 secant companion

For a short Weierstrass curve and a target with Y nonzero in characteristic
other than two, the secant numerator evaluates to 4Y^2. Its denominator
(x-t)^2 has order two at the target, hence the companion has pole order
two. If t is nonzero, x is a unit there, giving the second monomial-ratio
identity. The infinity and target valuations are distinct. No global
nonsingularity, genus-one height bound or circuit-spread theorem is proved
or assumed by this module.
-/

namespace ICFallDegrees.GrowthNonTorsionTargetSecant

open Polynomial WeierstrassCurve.Affine
open GrowthSecant GrowthNonTorsionInfinity GrowthNonTorsionInfinitySecant
open GrowthNonTorsionTargetLocal GrowthNonTorsionTargetValuation
open scoped Multiplicative

variable {K : Type*} [Field K]

theorem target_equation (a b t Y : K) (hR : Y^2 = cubic a b t) :
    (shortCurve a b).Equation t Y := by
  simpa [equation_iff,shortCurve,cubic] using hR

theorem target_not_neg (a b t Y : K) (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) :
    Y ≠ (shortCurve a b).negY t Y := by
  apply sub_ne_zero.mp
  simpa only [negY,shortCurve,zero_mul,sub_zero,sub_neg_eq_add,← two_mul] using
    mul_ne_zero h2 hY

theorem target_nonsingular (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) : (shortCurve a b).Nonsingular t Y := by
  exact (nonsingular_iff t Y).mpr
    ⟨target_equation a b t Y hR, Or.inr (target_not_neg a b t Y h2 hY)⟩

def secantTargetValuation (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) : Valuation (shortCurve a b).FunctionField ℤₘ₀ :=
  targetValuation (shortCurve a b) t Y (target_equation a b t Y hR)
    (target_nonsingular a b t Y hR h2 hY) (target_not_neg a b t Y h2 hY)

def coordinateNumerator (a b t Y : K) : (shortCurve a b).CoordinateRing :=
  let W := shortCurve a b
  let c := algebraMap K W.CoordinateRing
  (coordinateX W + c t) * (coordinateX W * c t + c a) + 2*c b +
    2*coordinateY W*c Y

theorem evaluate_numerator (a b t Y : K) (hR : Y^2 = cubic a b t) :
    pointEval (shortCurve a b) t Y (target_equation a b t Y hR)
      (coordinateNumerator a b t Y) = 4*Y^2 := by
  simp only [coordinateNumerator,map_add,map_mul,map_ofNat,
    pointEval_constant,pointEval_X,pointEval_Y]
  calc
    _ = 2 * cubic a b t + 2*Y^2 := by dsimp [cubic]; ring
    _ = _ := by rw [← hR]; ring

theorem map_coordinateNumerator (a b t Y : K) :
    algebraMap (shortCurve a b).CoordinateRing (shortCurve a b).FunctionField
      (coordinateNumerator a b t Y) =
    secantNumerator (constants (shortCurve a b) a) (constants (shortCurve a b) b)
      (fieldX (shortCurve a b)) (constants (shortCurve a b) t)
      (fieldY (shortCurve a b)) (constants (shortCurve a b) Y) := by
  simp only [coordinateNumerator,secantNumerator,companion,map_add,map_mul,map_ofNat,
    fieldX,fieldY,constants,RingHom.comp_apply,← IsScalarTower.algebraMap_apply]
  rfl

theorem targetValuation_companion (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) :
    secantTargetValuation a b t Y hR h2 hY (companionX (shortCurve a b) a b t Y) =
      ↑(Multiplicative.ofAdd (2:ℤ)) := by
  let W := shortCurve a b
  let v := secantTargetValuation a b t Y hR h2 hY
  have hnum : v (algebraMap W.CoordinateRing W.FunctionField
      (coordinateNumerator a b t Y)) = 1 := by
    exact targetValuation_unit W t Y (target_equation a b t Y hR)
      (target_nonsingular a b t Y hR h2 hY) (target_not_neg a b t Y h2 hY) _ (by
        rw [evaluate_numerator a b t Y hR]
        have h4 : (4:K) ≠ 0 := by
          simpa only [show (2:K)*2=4 by norm_num] using mul_ne_zero h2 h2
        exact mul_ne_zero h4 (pow_ne_zero 2 hY))
  rw [map_coordinateNumerator] at hnum
  have hp : v (fieldX W - constants W t) = ↑(Multiplicative.ofAdd (-1:ℤ)) :=
    targetValuation_parameter W t Y _ _ _
  change v (secantX _ _ _ _ _ _) = _
  rw [secantX,v.map_div,hnum,map_pow,hp]
  norm_num [← WithZero.coe_pow,← ofAdd_nsmul,← WithZero.coe_inv,← ofAdd_neg]

theorem targetValuation_monomial (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) (ht : t ≠ 0) (i j : ℕ) :
    secantTargetValuation a b t Y hR h2 hY
      (monomial (shortCurve a b) a b t Y i j) =
        ↑(Multiplicative.ofAdd (2*(j:ℤ))) := by
  have hx : secantTargetValuation a b t Y hR h2 hY (fieldX (shortCurve a b)) = 1 :=
    targetValuation_X (shortCurve a b) t Y (target_equation a b t Y hR)
      (target_nonsingular a b t Y hR h2 hY) (target_not_neg a b t Y h2 hY) ht
  rw [GrowthNonTorsionInfinitySecant.monomial,map_mul,map_pow,hx,one_pow,one_mul,map_pow,
    targetValuation_companion a b t Y hR h2 hY,← WithZero.coe_pow,← ofAdd_nsmul]
  simp only [nsmul_eq_mul]
  congr 2
  ring

theorem targetValuation_monomial_ratio (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) (ht : t ≠ 0) (i j i' j' : ℕ) :
    secantTargetValuation a b t Y hR h2 hY
      (monomial (shortCurve a b) a b t Y i j /
        monomial (shortCurve a b) a b t Y i' j') =
      ↑(Multiplicative.ofAdd (2*((j:ℤ)-j'))) := by
  rw [Valuation.map_div,targetValuation_monomial a b t Y hR h2 hY ht,
    targetValuation_monomial a b t Y hR h2 hY ht,
    ← WithZero.coe_div,← ofAdd_sub]
  congr 2
  ring

theorem targetValuation_ne_infinity (a b t Y : K) (hR : Y^2 = cubic a b t)
    (h2 : (2:K) ≠ 0) (hY : Y ≠ 0) (ht : t ≠ 0) :
    secantTargetValuation a b t Y hR h2 hY ≠ infinityValuation (shortCurve a b) := by
  intro he
  have hx : secantTargetValuation a b t Y hR h2 hY (fieldX (shortCurve a b)) = 1 :=
    targetValuation_X (shortCurve a b) t Y (target_equation a b t Y hR)
      (target_nonsingular a b t Y hR h2 hY) (target_not_neg a b t Y h2 hY) ht
  rw [he,infinityValuation_X] at hx
  have hgt : (1:ℤₘ₀) < ↑(Multiplicative.ofAdd (2:ℤ)) := by
    rw [WithZero.one_lt_coe]
    change (0:ℤ) < 2
    norm_num
  exact (ne_of_gt hgt) hx

end ICFallDegrees.GrowthNonTorsionTargetSecant
