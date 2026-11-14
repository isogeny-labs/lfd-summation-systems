import ICFallDegrees.GrowthRationalCoordinates
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Tactic

/-!
# N8: characteristic-p boundary relation for the rational section

The relation at a rational two-torsion target can be multiplied by
(X-1)^(p-1) and written with coefficients X^p-1 fixed under Frobenius.
This is an explicit large minimal-dependence candidate with p+1
monomials over F^p. Only the polynomial identity is proved here;
minimality and the Wronskian r=p determinant remain to be formalized.
-/

namespace ICFallDegrees.GrowthRationalBoundary

variable {R : Type*} [CommRing R]

/-- A characteristic-p identity on the rational section (X-1)Z=X+A.
It is the algebraic starting point for testing the r=p Wronskian case. -/
theorem pth_boundary_relation (p : ℕ) [Fact p.Prime] [CharP R p]
    (X A Z : R) (hsection : (X - 1) * Z = X + A) :
    (X ^ p - 1) * Z - (X ^ p - 1) -
      (A + 1) * (X - 1) ^ (p - 1) = 0 := by
  have hp : p = p - 1 + 1 :=
    (Nat.sub_add_cancel (Fact.out : p.Prime).one_le).symm
  have hpow : (X - 1) ^ p = (X - 1) ^ (p - 1) * (X - 1) := by
    conv_lhs => rw [hp]
    rw [pow_succ]
  have hfrob : (X - 1) ^ p = X ^ p - 1 := by
    simpa only [one_pow] using (sub_pow_char X (1 : R))
  calc
    (X ^ p - 1) * Z - (X ^ p - 1) -
        (A + 1) * (X - 1) ^ (p - 1) =
      (X - 1) ^ p * Z - (X - 1) ^ p -
        (A + 1) * (X - 1) ^ (p - 1) := by rw [hfrob]
    _ = (X - 1) ^ (p - 1) *
        ((X - 1) * Z - (X - 1) - (A + 1)) := by rw [hpow]; ring
    _ = 0 := by rw [hsection]; ring

end ICFallDegrees.GrowthRationalBoundary
