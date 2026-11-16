# Streaming variant of sage_gb.sage: one JSON line per system, flushed, with a per-system time limit.
# usage: sage sage_gb_stream.sage in.json out.jsonl seconds [maxr]
# Systems run in increasing number of variables; a timed-out system is written with "timeout": true.
import json, sys, time
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
set_random_seed(20261114)
src, dst, budget = sys.argv[1], sys.argv[2], int(sys.argv[3])
maxr = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
done = set()
try:
    done = {json.loads(l)['id'] for l in open(dst)}
except IOError:
    pass
systems = sorted(json.load(open(src)), key=lambda s: (s['r'], len(s['polys'])))
with open(dst, 'a') as out:
    for sysm in systems:
        if sysm['id'] in done or sysm['r'] > maxr:
            continue
        p, r = sysm['p'], sysm['r']
        R = PolynomialRing(GF(p), 'u', r, order='degrevlex')
        u = R.gens()
        polys = [sum(R(c) * prod(u[i] ** e[i] for i in range(r)) for e, c in terms) for terms in sysm['polys']]
        I = R.ideal(polys)
        t0 = time.time()
        try:
            alarm(budget)
            G = I.groebner_basis()
            dim = int(0) if I.is_one() else int(I.vector_space_dimension())
            cancel_alarm()
        except AlarmInterrupt:
            out.write(json.dumps(dict(id=sysm['id'], r=r, timeout=True, seconds=float(time.time() - t0))) + '\n')
            out.flush()
            continue
        rec = dict(id=sysm['id'], r=r, max_degree=int(max(g.degree() for g in G)),
                   leading=sorted([list(g.lm().exponents()[0]) for g in G]), dimension=dim,
                   seconds=float(time.time() - t0))
        if r <= 8:
            rec['basis'] = [[[list(e), int(c)] for e, c in g.dict().items()] for g in G]
        out.write(json.dumps(rec) + '\n')
        out.flush()
print('stream done')
