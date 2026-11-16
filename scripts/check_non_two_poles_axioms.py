"""Require all 42 exact two-place theorem axiom answers."""
import sys
from pathlib import Path

import check_non_torsion_axioms as audit

GROUPS = {
    "GrowthNonTorsionInfinity": [
        "norm_degree_add_le",
        "coordinateValue_mul",
        "coordinateValue_add_le",
        "coordinateValue_ne_zero",
        "infinityValuation_map",
        "infinityValuation_X",
        "infinityValuation_Y",
        "infinityValuation_constant",
        "infinityValuation_constant_le_one",
        "infinityValuation_uniformizer"
    ],
    "GrowthNonTorsionInfinitySecant": [
        "infinityValuation_companion",
        "companion_ne_zero",
        "infinityValuation_monomial",
        "infinityValuation_monomial_ratio",
        "fieldX_ne_constant",
        "shortCurve_equation",
        "shortCurve_s3"
    ],
    "GrowthNonTorsionTargetLocal": [
        "pointEval_constant",
        "pointEval_X",
        "pointEval_Y",
        "pointEval_surjective",
        "pointEval_YClass",
        "pointEval_XClass",
        "pointIdeal_maximal",
        "pointIdeal_eq_ker",
        "pointLocal_maximal_principal",
        "localParameter_ne_zero",
        "pointLocal_not_field",
        "pointLocal_dvr"
    ],
    "GrowthNonTorsionTargetValuation": [
        "targetValuation_parameter",
        "targetValuation_unit",
        "targetValuation_X",
        "targetValuation_constant"
    ],
    "GrowthNonTorsionTargetSecant": [
        "target_equation",
        "target_not_neg",
        "target_nonsingular",
        "evaluate_numerator",
        "map_coordinateNumerator",
        "targetValuation_companion",
        "targetValuation_monomial",
        "targetValuation_monomial_ratio",
        "targetValuation_ne_infinity"
    ]
}
audit.EXPECTED = {
    "ICFallDegrees." + module + "." + name
    for module, names in GROUPS.items() for name in names
}

if __name__ == "__main__":
    print("PASS", audit.check("\n".join(Path(p).read_text() for p in sys.argv[1:])),
          "two-place theorem axiom answers")
