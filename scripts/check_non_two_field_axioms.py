"""Require the eight exact axiom answers for the constructed function field."""
import sys
from pathlib import Path

import check_non_torsion_axioms as audit

audit.EXPECTED = {'ICFallDegrees.GrowthNonTorsionFunctionField.' + name for name in (
    'equation_irreducible', 'finrank_two', 'y_square', 'x_ne_zero',
    'x_ne_constant', 's3_zero', 'companionConstant_ne_zero', 'z_ne_zero')}

if __name__ == '__main__':
    print('PASS', audit.check('\n'.join(Path(p).read_text() for p in sys.argv[1:])),
          'explicit function-field theorem axiom answers')
