import ICFallDegrees.GrowthNonTorsionBridge
import Mathlib.FieldTheory.RatFunc.Degree
import Mathlib.Algebra.Polynomial.Degree.SmallDegree

noncomputable section

/-!
# The generic non-two-torsion S3 curve has no rational section in K(X)

The cubic is monic of degree three. A nonzero rational square has even
integer degree, so it cannot equal that cubic. Combined with recoverY,
this rules out the E6-style rational section for a target with Y nonzero.
It does not construct the quadratic extension or prove circuit spread.
-/

namespace ICFallDegrees.GrowthNonTorsionRationalObstruction

open Polynomial GrowthSecant GrowthNonTorsionBridge

variable {K : Type*} [Field K]

def cubicPolynomial (a b : K) : Polynomial K := X ^ 3 + C a * X + C b

theorem cubicPolynomial_natDegree (a b : K) :
    (cubicPolynomial a b).natDegree = 3 := by
  simpa [cubicPolynomial] using
    (Polynomial.natDegree_cubic (a := (1 : K)) (b := 0) (c := a) (d := b) one_ne_zero)

/-- The concrete generic cubic is not a square even as a rational function. -/
theorem cubicPolynomial_nonsquare (a b : K) (y : RatFunc K) :
    y ^ 2 ≠ algebraMap (Polynomial K) (RatFunc K) (cubicPolynomial a b) := by
  intro hs
  have hp : cubicPolynomial a b ≠ 0 := by
    intro hz
    have hd := cubicPolynomial_natDegree a b
    rw [hz, Polynomial.natDegree_zero] at hd
    omega
  have hy : y ≠ 0 := by
    intro hz
    have hmap := RatFunc.algebraMap_ne_zero hp
    apply hmap
    simpa [hz] using hs.symm
  have hd := congrArg RatFunc.intDegree hs
  rw [pow_two, RatFunc.intDegree_mul hy hy, RatFunc.intDegree_polynomial,
    cubicPolynomial_natDegree] at hd
  omega

/-- No generic rational S3 companion exists for a target with Y nonzero.
Thus the rational E6 shortcut cannot discharge the genus-one obligation. -/
theorem generic_s3_ne_zero (a b t Y : K)
    (h2 : (2 : K) ≠ 0) (hY : Y ≠ 0) (hR : Y ^ 2 = cubic a b t)
    (z : RatFunc K) :
    s3 (RatFunc.C a) (RatFunc.C b) RatFunc.X z (RatFunc.C t) ≠ 0 := by
  have h2' : (2 : RatFunc K) ≠ 0 := by
    simpa using (RatFunc.C : K →+* RatFunc K).injective.ne h2
  have hY' : RatFunc.C Y ≠ 0 := by
    simpa using (RatFunc.C : K →+* RatFunc K).injective.ne hY
  have hR' : (RatFunc.C Y) ^ 2 = cubic (RatFunc.C a) (RatFunc.C b) (RatFunc.C t) := by
    simpa only [cubic, map_add, map_mul, map_pow] using congrArg (RatFunc.C : K →+* RatFunc K) hR
  apply no_section_of_nonsquare _ _ _ _ _ (RatFunc.C Y) h2' hY' hR'
  intro y hy
  apply cubicPolynomial_nonsquare a b y
  simpa only [cubic, cubicPolynomial, map_add, map_mul, map_pow,
    RatFunc.algebraMap_X, RatFunc.algebraMap_C] using hy

#print axioms cubicPolynomial_natDegree
#print axioms cubicPolynomial_nonsquare
#print axioms generic_s3_ne_zero

end ICFallDegrees.GrowthNonTorsionRationalObstruction
