# Reduced DRL Groebner bases of the input systems (check C1 of PROTOCOL.md), with Sage/Singular.
# usage: sage sage_gb.sage in.json out.json
# input: list of {id, p, r, polys: [[[e_0..e_{r-1}], c], ...] per polynomial}
# output: list of {id, max_degree, leading: sorted exponent lists, dimension, basis (if r <= 8)}
import json, sys
set_random_seed(20261114)
src, dst = sys.argv[1], sys.argv[2]
out = []
for sysm in json.load(open(src)):
    p, r = sysm['p'], sysm['r']
    R = PolynomialRing(GF(p), 'u', r, order='degrevlex')
    u = R.gens()
    polys = []
    for terms in sysm['polys']:
        f = R(0)
        for e, c in terms:
            f += R(c) * prod(u[i] ** e[i] for i in range(r))
        polys.append(f)
    I = R.ideal(polys)
    G = I.groebner_basis()
    lead = sorted([list(g.lm().exponents()[0]) for g in G])
    rec = dict(id=sysm['id'], max_degree=int(max(g.degree() for g in G)), leading=lead,
               dimension=int(0) if I.is_one() else int(I.vector_space_dimension()))
    if r <= 8:
        rec['basis'] = [[[list(e), int(c)] for e, c in g.dict().items()] for g in G]
    out.append(rec)
json.dump(out, open(dst, 'w'))
print('sage done', len(out))
