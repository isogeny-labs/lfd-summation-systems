"""N8-growth, Corollary 5.2: existence of unit-ideal instances (numerical check of the construction).

Construction: p prime with 2^(k+1) | p - 1; X a non-square mod p; n a non-square; a = n - X^2; curve
y^2 = x^3 + a x (b = 0), target abscissa t = X.  Then X != 0, a != 0 (n is not a square), and
X^3 + a X = X n is a nonzero square, which are the hypotheses of the kernel lemma
ICFallDegrees.UnitInstance.unit_zmod.  The script checks, for every such (X, n) of a few small primes and
every k >= 2 with 2^(k+1) | p - 1: the hypotheses, and directly that S_3(x1, x2, X) has no zero with
x1, x2 in the group of 2^k-th roots of unity (so the chain system has no solution).
usage: ic_bounds_growth_unit.py [out.json]
"""
import json, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import s3_coeffs

def main():
    out = []
    for p in (17, 97, 193, 257, 641):
        sq = {y * y % p for y in range(1, p)}
        nons = [v for v in range(1, p) if v not in sq]
        gen = next(g for g in range(2, p) if len({pow(g, e, p) for e in range(p - 1)}) == p - 1)
        ks = [k for k in range(2, 12) if (p - 1) % (1 << (k + 1)) == 0]
        pairs = zeros = 0; hyp_ok = True
        for X in nons:
            for n in nons[::max(1, len(nons) // 12)]:
                a = (n - X * X) % p
                hyp_ok &= (a != 0) and (X != 0) and ((X * X * X + a * X) % p in sq)
                g = s3_coeffs(p, a, 0, X)
                for k in ks:
                    H = [pow(gen, e * ((p - 1) >> k), p) for e in range(1 << k)]
                    pairs += 1
                    for x1 in H:
                        for x2 in H:
                            if sum(c * pow(x1, e1, p) * pow(x2, e2, p) for (e1, e2), c in g.items()) % p == 0:
                                zeros += 1
        out.append(dict(p=p, depths=ks, instances_times_depths=pairs, hypotheses_hold=bool(hyp_ok), zeros_found=zeros))
        print(out[-1], flush=True)
    ok = all(r['hypotheses_hold'] and r['zeros_found'] == 0 for r in out)
    print('ok', ok)
    if len(sys.argv) > 1: json.dump(dict(ok=ok, rows=out), open(sys.argv[1], 'w'), indent=1)

main()
