"""Recompute frozen product-space ranks and check the prime-degree ceiling."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import _f2_rank  # noqa: E402
from csd.ic_rootfind import product_space  # noqa: E402

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'subspace-remainder-cost-20261125' / 'manifest.json'
PINNED_MANIFEST_SHA256 = 'd30bf96488e1a382ca85f575ce9134453485a4fbf636bb6c49af33aa91c7f68a'
OUT = HERE / 'verification.json'


def main() -> None:
    parent_raw = PARENT.read_bytes()
    assert hashlib.sha256(parent_raw).hexdigest() == PINNED_MANIFEST_SHA256
    parent = json.loads(parent_raw)
    rows = []
    for c in parent['cells']:
        n, k = c['n'], c['k']
        assert n in (23, 29) and k > 0 and 2*k-1 <= n
        F = GF2n(n, c['modulus'])
        basis = [1 << i for i in range(k)]
        products = [F.mul(a, b) for i, a in enumerate(basis)
                    for b in basis[i:]]
        rank = _f2_rank(products)
        assert rank == len(product_space(F, basis))
        lower_bound = 2*k-1
        assert rank == lower_bound
        excess = n-k-rank
        assert excess == n-3*k+1
        rows.append({'n':n, 'k':k, 'modulus':c['modulus'],
                     'product_space_dimension':rank,
                     'linearization_excess':excess,
                     'max_excess_at_fixed_k':n-k-lower_bound})
    assert [(r['n'],r['k']) for r in rows] == [(23,8),(23,10),(23,11),
                                               (29,12),(29,13),(29,14)]
    out = {'parent_manifest_sha256':PINNED_MANIFEST_SHA256,
           'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'status':'passed', 'cells':rows}
    raw = json.dumps(out,sort_keys=True,indent=2)+'\n'
    if OUT.exists():
        assert OUT.read_text() == raw, 'frozen verification changed'
    else:
        OUT.write_text(raw)
    print('verified',len(rows),'prime-degree cells')


if __name__ == '__main__':
    main()
