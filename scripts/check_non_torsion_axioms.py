"""Check all eight explicit axiom answers from the two new algebraic modules."""
import re
import sys
from pathlib import Path

PREFIX = 'ICFallDegrees.'
EXPECTED = {PREFIX + 'GrowthNonTorsionBridge.' + x for x in
            ('recoverY_square', 'recoverY_secantX', 'secantX_recoverY',
             'no_section_of_nonsquare', 'secantX_ne_zero')}
EXPECTED |= {PREFIX + 'GrowthNonTorsionRationalObstruction.' + x for x in
             ('cubicPolynomial_natDegree', 'cubicPolynomial_nonsquare', 'generic_s3_ne_zero')}
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}


def check(text):
    if 'sorryAx' in text or 'Lean.ofReduceBool' in text or 'native_decide' in text:
        raise ValueError('placeholder or native evaluation in audit')
    answers = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text, re.S)
    answers += [(n, '') for n in re.findall(r"'([^']+)' does not depend on any axioms", text)]
    if len(answers) != len(EXPECTED) or {n for n, _ in answers} != EXPECTED:
        raise ValueError('missing, duplicate or unexpected theorem audit answer')
    for name, body in answers:
        extra = set(re.findall(r'[A-Za-z_][A-Za-z0-9_.]*', body)) - ALLOWED
        if extra:
            raise ValueError(f'{name}: unexpected transitive axioms {sorted(extra)}')
    return len(answers)


if __name__ == '__main__':
    print('PASS', check('\n'.join(Path(p).read_text() for p in sys.argv[1:])), 'explicit theorem axiom answers')
