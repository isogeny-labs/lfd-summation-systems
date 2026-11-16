import ICFallDegrees.GrowthNonTorsionInfinity
import ICFallDegrees.GrowthSecant

noncomputable section

namespace ICFallDegrees.GrowthNonTorsionInfinitySecant

open GrowthSecant GrowthNonTorsionInfinity WeierstrassCurve.Affine
open scoped Multiplicative

variable {K : Type*} [Field K] (W : WeierstrassCurve.Affine K)

private def exp (n : ℤ) : ℤₘ₀ := ↑(Multiplicative.ofAdd n)

def companionX (a b t Y : K) : W.FunctionField :=
  secantX (constants W a) (constants W b) (fieldX W) (constants W t)
    (fieldY W) (constants W Y)

/-- The secant companion has order zero at infinity when t is nonzero.
The coefficient parameters may be arbitrary; S3 needs the curve equations
separately. This result uses the actual constructed infinity valuation. -/
theorem infinityValuation_companion (a b t Y : K) (ht : t ≠ 0) :
    infinityValuation W (companionX W a b t Y) = 1 := by
  let v := infinityValuation W
  let c := constants W
  let x := fieldX W
  let y := fieldY W
  have hx : v x = exp 2 := infinityValuation_X W
  have hy : v y = exp 3 := infinityValuation_Y W
  have hct : v (c t) = 1 := infinityValuation_constant W t ht
  have h23 : exp 2 ≤ exp 3 := by norm_num [exp]
  have h13 : (1:ℤₘ₀) ≤ exp 3 := by
    rw [exp,WithZero.one_le_coe]
    change (0:ℤ) ≤ 3
    norm_num
  have h34 : exp 3 < exp 4 := by norm_num [exp]
  have hx2 : v (x^2) = exp 4 := by
    rw [map_pow,hx]
    norm_num [exp,pow_two,← WithZero.coe_mul,← ofAdd_add]
  have hlead : v (c t * x^2) = exp 4 := by
    rw [map_mul,hct,hx2,one_mul]
  have hlowx : v (c (t^2+a) * x) ≤ exp 3 := by
    rw [map_mul]
    calc
      _ ≤ 1 * v x := mul_le_mul_right' (infinityValuation_constant_le_one W _) _
      _ ≤ exp 3 := by simpa only [one_mul,hx] using h23
  have hlowc : v (c (a*t+2*b)) ≤ exp 3 :=
    (infinityValuation_constant_le_one W _).trans h13
  have hlowy : v (c (2*Y) * y) ≤ exp 3 := by
    rw [map_mul]
    calc
      _ ≤ 1 * v y := mul_le_mul_right' (infinityValuation_constant_le_one W _) _
      _ ≤ exp 3 := by rw [one_mul,hy]
  let rest := c (t^2+a) * x + c (a*t+2*b) + c (2*Y) * y
  have hrest : v rest ≤ exp 3 := by
    apply (v.map_add _ _).trans
    exact max_le ((v.map_add _ _).trans (max_le hlowx hlowc)) hlowy
  have he : secantNumerator (c a) (c b) x (c t) y (c Y) = c t*x^2+rest := by
    dsimp [rest]
    simp only [secantNumerator,companion,map_add,map_mul,map_pow,map_ofNat]
    ring
  have hn : v (secantNumerator (c a) (c b) x (c t) y (c Y)) = exp 4 := by
    rw [he]
    apply (v.map_add_eq_of_lt_left _).trans hlead
    rw [hlead]
    exact hrest.trans_lt h34
  have htx : v (c t) < v x := by
    rw [hct,hx,exp,WithZero.one_lt_coe]
    change (0:ℤ) < 2
    norm_num
  have hd : v ((x-c t)^2) = exp 4 := by
    rw [map_pow,v.map_sub_eq_of_lt_left htx,hx]
    norm_num [exp,pow_two,← WithZero.coe_mul,← ofAdd_add]
  change v (secantX (c a) (c b) x (c t) y (c Y)) = 1
  rw [secantX,v.map_div,hn,hd]
  exact div_self (by simp [exp])

theorem companion_ne_zero (a b t Y : K) (ht : t ≠ 0) :
    companionX W a b t Y ≠ 0 := by
  intro hz
  have h := infinityValuation_companion W a b t Y ht
  rw [hz,map_zero] at h
  exact zero_ne_one h

def monomial (a b t Y : K) (i j : ℕ) : W.FunctionField :=
  fieldX W ^ i * companionX W a b t Y ^ j

theorem infinityValuation_monomial (a b t Y : K) (ht : t ≠ 0) (i j : ℕ) :
    infinityValuation W (monomial W a b t Y i j) =
      ↑(Multiplicative.ofAdd (2 * (i:ℤ))) := by
  rw [monomial,map_mul,map_pow,infinityValuation_X,map_pow,
    infinityValuation_companion W a b t Y ht,one_pow,mul_one,
    ← WithZero.coe_pow,← ofAdd_nsmul]
  simp only [nsmul_eq_mul]
  congr 2
  ring

theorem infinityValuation_monomial_ratio (a b t Y : K) (ht : t ≠ 0)
    (i j i' j' : ℕ) :
    infinityValuation W (monomial W a b t Y i j / monomial W a b t Y i' j') =
      ↑(Multiplicative.ofAdd (2 * ((i:ℤ) - i'))) := by
  rw [Valuation.map_div,infinityValuation_monomial W a b t Y ht,
    infinityValuation_monomial W a b t Y ht,← WithZero.coe_div,← ofAdd_sub]
  congr 2
  ring

theorem fieldX_ne_constant (t : K) : fieldX W ≠ constants W t := by
  intro he
  have h := infinityValuation_constant_le_one W t
  rw [← he,infinityValuation_X] at h
  rw [WithZero.coe_le_one] at h
  change (2:ℤ) ≤ 0 at h
  omega

def shortCurve (a b : K) : WeierstrassCurve.Affine K := ⟨0,0,0,a,b⟩

theorem shortCurve_equation (a b : K) :
    fieldY (shortCurve a b) ^ 2 =
      cubic (constants (shortCurve a b) a) (constants (shortCurve a b) b)
        (fieldX (shortCurve a b)) := by
  let W := shortCurve a b
  have h : (algebraMap W.CoordinateRing W.FunctionField)
      (CoordinateRing.mk W W.polynomial) = 0 := by
    rw [AdjoinRoot.mk_self,map_zero]
  simp only [WeierstrassCurve.Affine.polynomial,W,shortCurve,map_sub,map_add,map_mul,
    map_pow,map_zero,zero_mul,mul_zero,zero_add,add_zero] at h
  exact sub_eq_zero.mp h

theorem shortCurve_s3 (a b t Y : K) (hR : Y ^ 2 = cubic a b t) :
    s3 (constants (shortCurve a b) a) (constants (shortCurve a b) b)
      (fieldX (shortCurve a b)) (companionX (shortCurve a b) a b t Y)
      (constants (shortCurve a b) t) = 0 := by
  let W := shortCurve a b
  apply secant_s3 _ _ _ _ _ _ (shortCurve_equation a b) _ (fieldX_ne_constant W t)
  simpa only [cubic,map_add,map_mul,map_pow] using congrArg (constants W) hR


end ICFallDegrees.GrowthNonTorsionInfinitySecant
