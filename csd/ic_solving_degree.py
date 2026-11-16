"""Reduced-Groebner-basis degree of the vanishing ideal of a point set on a cube.

For the chain systems (prime) and Boolean systems (binary) of the LFD paper the
ideal I is the vanishing ideal of its solution set Z (docs/proofs/N10, Lemma 3)
and every standard monomial is squarefree (N10, Lemma 2).  The Buchberger--Moeller
algorithm therefore only has to visit squarefree monomials: a monomial is
standard iff its evaluation vector on Z is independent of those of the smaller
standard monomials (any degree-compatible order, increasing order).  The
minimal generators of in(I) are the nonstandard squarefree monomials whose
proper divisors are all standard, and u_v^2 for each standard variable u_v.

Monomials are bit masks over r variables; bit v is u_v and u_0 > u_1 > ... .
Orders: 'drl' (degree reverse lexicographic) and 'grlex'.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

from csd.elim import fp_echelon


def order_key(mask: int, r: int, order: str):
    """Sort key, increasing in the term order, for squarefree masks over r variables."""
    deg = bin(mask).count('1')
    bits = [(mask >> v) & 1 for v in range(r)]
    if order == 'drl':
        # same degree: x^a < x^b iff the last nonzero entry of a - b is positive,
        # i.e. at the highest index where they differ, the smaller one has the bit.
        return (deg, tuple(-bits[v] for v in range(r - 1, -1, -1)))
    if order == 'grlex':
        return (deg, tuple(bits))
    raise ValueError(order)


class PrimeChainPoints:
    """Points of a chain system: x_i = c_i * g^{E_i}; chain monomial mask with block i
    in bits [i*k, (i+1)*k) evaluates to g^{sum_i a_i E_i} (a_i = block i of the mask)."""

    def __init__(self, exps: list[tuple[int, ...]], m: int, k: int, g: int, p: int):
        self.E = np.array(exps, dtype=np.int64).reshape(len(exps), m)
        self.m, self.k, self.g, self.p = m, k, g, p
        self.N = 1 << k
        self.r = m * k
        self.powers = np.array([pow(g, e, p) for e in range(self.N)], dtype=np.int64)

    def evaluate(self, masks: list[int]) -> np.ndarray:
        """len(masks) x |Z| matrix of monomial values."""
        out = np.zeros((len(masks), len(self.E)), dtype=np.int64)
        lowmask = self.N - 1
        for j, mask in enumerate(masks):
            e = np.zeros(len(self.E), dtype=np.int64)
            for i in range(self.m):
                a = (mask >> (i * self.k)) & lowmask
                if a:
                    e = e + a * self.E[:, i]
            out[j] = self.powers[e % self.N]
        return out


class BooleanPoints:
    """Points of {0,1}^r given as bit masks; a squarefree monomial evaluates to 1 iff mask is a subset."""

    def __init__(self, points: list[int], r: int):
        self.P = np.array(points, dtype=object)
        self.pts = [int(z) for z in points]
        self.r = r
        self.p = 2

    def evaluate(self, masks: list[int]) -> np.ndarray:
        assert self.r <= 63
        P = np.array(self.pts, dtype=np.uint64)
        M = np.array(masks, dtype=np.uint64)[:, None]
        return ((P[None, :] & M) == M).astype(np.int64)


def staircase(points, npoints: int, r: int, p: int, order: str = 'drl', chunk: int = 4000) -> dict:
    """Buchberger--Moeller on squarefree monomials.  Returns the standard monomials per
    degree (Hilbert function of A/I by degree), the minimal generators of in(I) by
    degree, and max.GB.deg."""
    if npoints == 0:
        return dict(points=0, standard_by_degree=[], generators_by_degree={0: 1}, max_gb_degree=0,
                    standard_count=0, squarefree_generators=[0], standard_variables=[])   # in(I) = (1)
    std = [0]
    std_set = {0}
    E = points.evaluate([0])                      # rows: evaluation vectors of standard monomials
    by_deg = {0: [0]}
    gens: dict[int, int] = {}
    gen_masks: list[int] = []
    d = 0
    while True:
        d += 1
        if d > r:
            break
        prev = by_deg.get(d - 1, [])
        cands = set()
        for s in prev:
            for v in range(r):
                if not (s >> v) & 1:
                    m = s | (1 << v)
                    if m in cands:
                        continue
                    if all((m & ~(1 << w)) in std_set for w in range(r) if (m >> w) & 1):
                        cands.add(m)
        if not cands:
            break
        cands = sorted(cands, key=lambda m: order_key(m, r, order))
        new_std = []
        if len(std) == npoints:
            gens[d] = gens.get(d, 0) + len(cands)
            gen_masks.extend(cands)
        else:
            for start in range(0, len(cands), chunk):
                block = cands[start:start + chunk]
                C = points.evaluate(block)
                M = np.vstack([E, C]) % p
                rank, piv, _ = fp_echelon(M.T.copy(), p)
                base = E.shape[0]
                assert piv[:base] == list(range(base)), 'standard rows lost independence'
                chosen = [c - base for c in piv[base:]]
                chosen_set = set(chosen)
                for i, m in enumerate(block):
                    if i in chosen_set:
                        new_std.append(m)
                    else:
                        gens[d] = gens.get(d, 0) + 1
                        gen_masks.append(m)
                if chosen:
                    E = np.vstack([E, C[chosen]])
            for m in new_std:
                std.append(m)
                std_set.add(m)
        by_deg[d] = new_std
    # squares of standard variables
    nvar_std = len(by_deg.get(1, []))
    if nvar_std:
        gens[2] = gens.get(2, 0) + nvar_std
    assert len(std) == npoints, (len(std), npoints)
    hf = [len(by_deg.get(i, [])) for i in range(max(by_deg) + 1)]
    while hf and hf[-1] == 0:
        hf.pop()
    return dict(points=npoints, standard_by_degree=hf, generators_by_degree=dict(sorted(gens.items())),
                max_gb_degree=max(gens) if gens else 0, standard_count=len(std),
                squarefree_generators=sorted(gen_masks), standard_variables=sorted(by_deg.get(1, [])))


def resolve(lfd: int, npoints: int | None, r: int | None) -> tuple[str | None, int | None]:
    """Tag by N10 Cor. 4: (tag, settled sd) or (None, None)."""
    if npoints == 0:
        return 'unit', lfd
    if npoints is not None and lfd >= max(2, (npoints + 1).bit_length() - 1):
        return 'count', lfd
    if r is not None and lfd >= max(2, r):
        return 'cube', lfd
    return None, None


def selftest():
    # 1. vanishing ideal of the whole cube {0,1}^3: in(I) = (Z^2): max.GB.deg 2, all squarefree standard
    pts = BooleanPoints(list(range(8)), 3)
    s = staircase(pts, 8, 3, 2)
    assert s['standard_by_degree'] == [1, 3, 3, 1] and s['max_gb_degree'] == 2, s
    # 2. single point: in(I) = (Z_0, Z_1, Z_2), degree 1
    s = staircase(BooleanPoints([5], 3), 1, 3, 2)
    assert s['max_gb_degree'] == 1 and s['generators_by_degree'] == {1: 3}, s
    # 3. N points of the cube in general position: 2^d - 1 <= N for every squarefree generator (N10 Lemma 2)
    import random
    rng = random.Random(1)
    for _ in range(200):
        r = rng.randint(2, 7)
        N = rng.randint(1, 2 ** r)
        Z = rng.sample(range(2 ** r), N)
        s = staircase(BooleanPoints(Z, r), N, r, 2, order=rng.choice(['drl', 'grlex']))
        sq = [d for d in s['generators_by_degree'] if d != 2 or True]
        assert s['max_gb_degree'] <= max(2, (N + 1).bit_length() - 1), (r, N, s)
    # 4. order keys: drl vs grlex on degree 2 in 3 variables
    drl = sorted([0b011, 0b101, 0b110], key=lambda m: order_key(m, 3, 'drl'))
    gl = sorted([0b011, 0b101, 0b110], key=lambda m: order_key(m, 3, 'grlex'))
    # u0 > u1 > u2. DRL: u0u1 > u0u2 > u1u2 ; grlex (= lex within degree): u0u1 > u0u2 > u1u2 too.
    assert drl == [0b110, 0b101, 0b011] and gl == [0b110, 0b101, 0b011], (drl, gl)
    # degree 3 in 4 variables separates them: DRL: u0u1u2 > u0u1u3 > u0u2u3 > u1u2u3 equals lex here;
    # degree 2 in 4 variables: grlex u0u3 > u1u2 ; DRL u1u2 > u0u3
    assert order_key(0b1001, 4, 'grlex') > order_key(0b0110, 4, 'grlex')
    assert order_key(0b0110, 4, 'drl') > order_key(0b1001, 4, 'drl')
    print('SELFTEST PASSED')


if __name__ == '__main__':
    selftest()
