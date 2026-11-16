# Check C1b: reduced deglex Groebner bases of the binary systems with sd > LFD, computed by PolyBoRi
# (BooleanPolynomialRing: the field equations are implicit).  usage: sage polybori_gb.sage in.json out.jsonl
# input: list of {id, r, equations: [[mask, ...], ...]} (a polynomial is the sum of its monomial masks)
import json, sys, time
src, dst = sys.argv[1], sys.argv[2]
with open(dst, 'w') as out:
    for s in json.load(open(src)):
        r = int(s['r'])
        B = BooleanPolynomialRing(r, 'u', order='deglex')
        u = B.gens()
        polys = []
        for eq in s['equations']:
            f = B(0)
            for m in eq:
                t = B(1)
                for v in range(r):
                    if (int(m) >> v) & 1:
                        t *= u[v]
                f += t
            polys.append(f)
        t0 = time.time()
        G = B.ideal(polys).groebner_basis()
        lead = sorted(int(sum(1 << v for v in range(r) if u[v] in g.lm().variables())) for g in G)
        out.write(json.dumps(dict(id=s['id'], r=r, max_degree=int(max(g.degree() for g in G)), leading_masks=lead,
                                  seconds=float(time.time() - t0))) + '\n')
        out.flush()
print('polybori done')
