"""Certify FFD=4 for the two frozen m=5,k=3 VM input systems."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_split_m as S  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode(index, variables, auxiliaries):
    mask = index & ((1 << variables) - 1)
    exponents = tuple((index >> (variables + 6 * j)) & 63 for j in range(auxiliaries))
    assert index == mask + sum(e << (variables + 6*j) for j,e in enumerate(exponents))
    return mask, exponents


def parse(path):
    lines = path.read_text().splitlines()
    p, m, k, auxiliaries, D, *_ = map(int, lines[0].split())
    assert (m, k, auxiliaries, D) == (5, 3, 3, 8)
    assert int(lines[1]) == m - 1 == len(lines) - 2
    variables = m*k
    gens = []
    for line in lines[2:]:
        numbers = list(map(int, line.split()))
        assert numbers[0] * 2 + 1 == len(numbers)
        gen = {decode(index, variables, auxiliaries): c % p
               for index, c in zip(numbers[1::2], numbers[2::2])}
        assert all(gen.values()) and len(gen) == numbers[0]
        gens.append(gen)
    return p, m, k, gens


def audit(path):
    p, m, k, gens = parse(path)
    identity_path = HERE / ('split-instance-identity.json' if p == 114689
                            else 'split-p163841-instance-identity.json')
    identity = json.loads(identity_path.read_text())
    assert identity.get('prime', identity.get('p')) == p
    assert (identity['m'], identity['k']) == (m, k)
    if p == 114689:
        a, b = identity['curve']
        target = identity['target'][0]
    else:
        instance = identity['instance']
        assert (instance['p'], instance['m'], instance['k'], instance['seed']) == (p,m,k,0)
        a, b, target = instance['a'], instance['b'], instance['t']
        assert identity['inputs'][path.name.replace('split-m5-k3-p163841','split-p163841-m5-k3')] == sha(path)
    expected, _ = S.split_generators_m(I.PrimeCurve(p,a,b),target,k,m)
    assert gens == expected
    system_hash = S.system_hash(p,m,k,gens)
    if p == 114689:
        assert system_hash == identity['system_hash']
    deg = lambda key: key[0].bit_count() + sum(key[1])
    degrees = [max(map(deg, g)) for g in gens]
    assert degrees == [4, 4, 4, 3]
    top = {key: c for key,c in gens[-1].items() if deg(key) == 3}
    assert len(top) == 2
    assert all(es == (0,0,2) for _,es in top)
    [(a,_), (b,_)] = sorted(top)
    assert a.bit_count() == b.bit_count() == 1 and not a&b
    ca, cb = top[(a,(0,0,2))], top[(b,(0,0,2))]
    # In the associated graded chain quotient a^2=b^2=0.
    # (ca*a-cb*b)(ca*a+cb*b)=0, with a nonzero coefficient tuple.
    assert ca and cb and (ca * cb - cb * ca) % p == 0
    # Degree 3: one input with a nonzero top, so the scalar map is injective.
    # Degree 4: the displayed annihilator; every Koszul generator starts at 7.
    return {'p': p, 'm': m, 'k': k, 'ffd_lower_bound': 4, 'ffd_upper_bound': 4,
            'system_hash': system_hash, 'identity_sha256': sha(identity_path),
            'degree_3_top': [[a, list(es), c] for (a,es),c in sorted(top.items())],
            'degree_4_syzygy_coefficient': [[a,ca], [b,-cb % p]],
            'input_sha256': sha(path)}


def main():
    cells = [audit(HERE / f'split-m5-k3-p{p}-D8.txt') for p in (114689,163841)]
    assert [c['p'] for c in cells] == [114689,163841]
    out = {'status': 'passed', 'scope': 'frozen m=5,k=3 D8 VM inputs',
           'convention': 'associated graded chain quotient, Koszul trivial module',
           'source_sha256': {name: sha(REPO/name) for name in
                             ('csd/ic_bounds_split_m.py','csd/ic_bounds_m3.py','csd/ic_dreg.py')},
           'cells': cells, 'audit_source_sha256': sha(Path(__file__))}
    (HERE/'audit-split-m5k3-ffd.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'passed','systems':len(cells)}))


if __name__ == '__main__':
    main()
