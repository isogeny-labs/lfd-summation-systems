#!/usr/bin/env python3
"""Exact check: a proper coordinate space with small |Z| and LFD 5 (not preregistered).

V = F_(2^8) inside F_(2^16) (a proper, Frobenius-invariant, k = 8 < n = 16 subspace),
curve y^2 + xy = x^3 + a x^2 + 1, order-two target (0, 1).  Then
S_3(x1, x2, 0) = (x1 x2 + 1)^2 with x1 x2 + 1 in F_(2^8), so the span of the Boolean
equations is the span of the 8 coordinates of x1 x2 + 1 over F_(2^8): the inverse-family
ideal of n = 8 (docs/IC_INVERSE_PARITY_PROOF.md: written LFD = 8/2 + 1 = 5).  Here the
LFD is certified directly with the runner's certificate (csd.ic_bounds_binary.certify),
on |Z| = 255 solutions.  Also the same for V = F_(2^6) in F_(2^12) (predicted 4).
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_binary as B  # noqa: E402
from csd.ic_research_campaign import field  # noqa: E402


def subfield_basis(F, d):
    n = F.N
    imgs = [F.frob(1 << i, d) ^ (1 << i) for i in range(n)]
    # kernel of the F_2-linear map x -> x^(2^d) + x
    rows = []  # (image, preimage)
    basis = {}
    kern = []
    for i in range(n):
        v, pre = imgs[i], 1 << i
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p][0]; pre ^= basis[p][1]
            else:
                basis[p] = (v, pre); break
        if not v:
            kern.append(pre)
    assert len(kern) == d
    return kern


def main():
    out = []
    for n, d in ((12, 6), (16, 8)):
        F = field(n)
        V = subfield_basis(F, d)
        for a in (0, 1):
            E = I.BinaryCurve(F, a, 1)
            system = B.m2_system(E, V, 'frobenius-invariant', (0, 1), 'order-two', None)
            cert = B.certify(system, 7, 20000, 600_000_000)
            row = dict(n=n, k=d, a=a, b=1, target=[0, 1], basis=V, system_hash=system['system_hash'],
                       Z=len(system['points']), subfield=system['subfield'],
                       predicted_from_inverse_theorem=max(3, d // 2 + 1) if d % 2 == 0 else 3, **cert)
            out.append(row)
            print(json.dumps({k: row[k] for k in ('n', 'k', 'a', 'Z', 'subfield', 'status', 'exact_lfd',
                                                  'strict_falls', 'predicted_from_inverse_theorem',
                                                  'closure_seconds', 'peak_rss_mb')}))
    with open(os.path.join(HERE, 'subfield-check.jsonl'), 'w') as f:
        for r in out:
            f.write(json.dumps(r) + '\n')


if __name__ == '__main__':
    main()
