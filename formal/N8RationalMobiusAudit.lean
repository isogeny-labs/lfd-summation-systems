import ICFallDegrees.GrowthRationalChainUnit
import ICFallDegrees.GrowthRationalGridExclusion
import ICFallDegrees.GrowthRationalMobius
import ICFallDegrees.GrowthSecant
import ICFallDegrees.GrowthPoleArithmetic

-- Pinned Lean 4.19.0 smoke audit for PR #133 inherited source
#print axioms ICFallDegrees.GrowthSecant.secant_on_curve
#print axioms ICFallDegrees.GrowthPoleArithmetic.spread_from_two_pole_orders

-- New rational field identities.
#print axioms ICFallDegrees.GrowthRationalMobius.mul_rationalZ
#print axioms ICFallDegrees.GrowthRationalMobius.rationalZ_sub_one
#print axioms ICFallDegrees.GrowthRationalMobius.rationalZ_ne_one
#print axioms ICFallDegrees.GrowthRationalMobius.rationalZ_involutive
#print axioms ICFallDegrees.GrowthRationalMobius.ratFunc_xi2_involutive

#print axioms ICFallDegrees.GrowthRationalGridExclusion.s3_ne_zero_on_root_grid

#print axioms ICFallDegrees.GrowthRationalChainUnit.eval_rawS3_rational
#print axioms ICFallDegrees.GrowthRationalChainUnit.unit_of_primitiveRoot_frobenius
