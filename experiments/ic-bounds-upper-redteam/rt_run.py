"""Red-team exact check of T1 (LFD <= min(k1,k2)+2) in the ORDINARY polynomial ring.
Engine: experiments/ic-curve-validation/code/lfd_engine.py (exact_lfd: iterated closure in the
full monomial space, ideal part from evaluation ranks, certificate of stabilisation).
Usage: python3 rt_run.py <job> ; writes rt-<job>.jsonl and rt-<job>-summary.json"""
import sys, json, time, random, math
from rt_lib import *

def plane(p):
    return [(a, b, t) for a in range(p) for b in range(p) if nonsingular(p, a, b) for t in range(p)]

def special(p, k1, k2, rng, ncurves=6, nrand=0):
    """Adversarial instances for the branches of the proof: t = 0, order-two targets (every fibre
    polynomial is a square), tangent targets 3t^2 + a = 0, t in H1 (a fibre drops degree), t in H1
    with f(t) = 0 (a constant fibre polynomial), j = 0 and j = 1728, plus nrand random ones."""
    H1 = mu(p, 1 << min(k1, k2)); Hbig = mu(p, 1 << max(k1, k2))
    out = set()
    curves = [(0, 1), (0, rng.randrange(1, p)), (1, 0), (rng.randrange(1, p), 0)]
    while len(curves) < 4 + ncurves:
        a, b = rng.randrange(p), rng.randrange(p)
        if nonsingular(p, a, b):
            curves.append((a, b))
    for a, b in curves:
        if not nonsingular(p, a, b):
            continue
        ts = {0}
        ts |= {t for t in range(p) if (t ** 3 + a * t + b) % p == 0}          # order two
        ts |= {t for t in range(p) if (3 * t * t + a) % p == 0}               # tangent case
        ts |= set(rng.sample(H1, min(3, len(H1)))) | set(rng.sample(Hbig, min(3, len(Hbig))))
        for t in ts:
            out.add((a, b, t))
    # constant fibre polynomial: t in H and f(t) = 0
    for _ in range(ncurves):
        t = rng.choice(Hbig); a = rng.randrange(p); b = (-t ** 3 - a * t) % p
        if nonsingular(p, a, b):
            out.add((a, b, t))
        t = rng.choice(H1); a = rng.randrange(p); b = (-t ** 3 - a * t) % p
        if nonsingular(p, a, b):
            out.add((a, b, t))
    got = 0
    while got < nrand:
        a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
        if nonsingular(p, a, b):
            out.add((a, b, t)); got += 1
    return sorted(out)

def jobs(name):
    rng = random.Random(987654321)
    J = []   # (family, p, k1, k2, s1, s2, raw, list of (a,b,t), Dcap)
    if name == 'eq12':
        for p in (3, 5, 7, 11, 13):
            for c in coset_reps(p, 1):
                J.append(('k1-plane', p, 1, 1, 1, c, False, plane(p), None))
        J.append(('k2-plane', 5, 2, 2, 1, 1, False, plane(5), None))
        for p in (13, 17):
            for c in coset_reps(p, 2):
                J.append(('k2-plane', p, 2, 2, 1, c, False, plane(p), None))
        for p in (29, 37):
            J.append(('k2-plane', p, 2, 2, 1, 1, False, plane(p), None))
        J.append(('k2-plane', 29, 2, 2, 1, coset_reps(29, 2)[1], False, plane(29), None))
        for p in (41, 53, 61, 73, 97, 101, 109, 113, 12289, 40961):
            for c in coset_reps(p, 2)[:3]:
                J.append(('k2-special', p, 2, 2, 1, c, False, special(p, 2, 2, rng, nrand=300), None))
    elif name == 'eq3':
        for c in coset_reps(17, 3):
            J.append(('k3-plane', 17, 3, 3, 1, c, False, plane(17), None))
        for p in (41, 73, 89, 97, 113, 12289):
            for c in coset_reps(p, 3)[:2]:
                J.append(('k3-special', p, 3, 3, 1, c, False, special(p, 3, 3, rng, nrand=400), None))
    elif name == 'eq45':
        def cut(insts, n):
            return insts if len(insts) <= n else sorted(rng.sample(insts, n))
        J.append(('k4-special', 17, 4, 4, 1, 1, False, cut(special(17, 4, 4, rng, ncurves=40, nrand=300), 500), None))
        for p in (97, 113, 193, 241, 257, 12289):
            for c in coset_reps(p, 4)[:2]:
                J.append(('k4-special', p, 4, 4, 1, c, False, cut(special(p, 4, 4, rng, nrand=40), 90), None))
        for p in (97, 193, 257):
            for c in coset_reps(p, 5)[:2]:
                J.append(('k5-special', p, 5, 5, 1, c, False, cut(special(p, 5, 5, rng, ncurves=2, nrand=8), 24), 6))
    elif name == 'uneq17':
        for (k1, k2) in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)):
            for c in coset_reps(17, max(k1, k2)):
                J.append(('uneq-plane', 17, k1, k2, 1, c, False, plane(17), None))
        J.append(('uneq-plane', 17, 2, 1, 1, 1, False, plane(17), None))      # builder symmetry control
        J.append(('uneq-plane', 13, 1, 2, 1, 1, False, plane(13), None))
        J.append(('uneq-plane', 5, 1, 2, 1, 1, False, plane(5), None))
        J.append(('uneq-special', 17, 3, 4, 1, 1, False, special(17, 3, 4, rng, ncurves=30, nrand=600), None))
    elif name == 'uneqdeep':
        # short block of depth 1, 2, 3 against a long block: the bound min(k1,k2)+2 is far below the
        # equal-depth values there, so this is where a defect of the selector argument would show.
        plan = []
        for p, kmax in ((97, 5), (193, 6), (257, 8), (769, 8), (12289, 12), (40961, 12)):
            for k2 in range(5, kmax + 1):
                plan.append((p, 1, k2, 80 if k2 <= 9 else 16))
                if k2 <= 10:
                    plan.append((p, 2, k2, 40 if k2 <= 6 else (16 if k2 <= 8 else 4)))
                if k2 <= 6:
                    plan.append((p, 3, k2, 8 if k2 == 5 else 4))
                if k2 == 7 and p in (257, 12289):
                    plan.append((p, 3, k2, 1))
        for p, k1, k2, n in plan:
            r, bound = k1 + k2, k1 + 2
            cap = max(D for D in range(2, bound + 2) if math.comb(r + D, D) <= 5200)
            if cap < bound:
                continue
            for c in coset_reps(p, k2)[:2]:
                insts = special(p, k1, k2, rng, ncurves=2, nrand=n)
                if len(insts) > n:
                    insts = sorted(rng.sample(insts, n))
                J.append(('uneq-deep', p, k1, k2, 1, c, False, insts, cap))
    elif name == 'uneqmid':
        # unequal depths with both blocks >= 4: the scan is capped below the bound when the
        # matrices get large; an instance not certified within the cap is reported as uncertified.
        for p, k1, k2, n in ((97, 4, 5, 12), (193, 4, 5, 12), (257, 4, 5, 12), (193, 4, 6, 8), (257, 4, 6, 8),
                             (193, 5, 6, 4), (257, 5, 6, 4), (257, 4, 7, 4), (257, 3, 8, 4), (12289, 3, 8, 2)):
            r, bound = k1 + k2, min(k1, k2) + 2
            cap = max(D for D in range(2, bound + 2) if math.comb(r + D, D) <= 5200)
            for c in coset_reps(p, max(k1, k2))[:2]:
                insts = special(p, k1, k2, rng, ncurves=2, nrand=n)
                if len(insts) > n:
                    insts = sorted(rng.sample(insts, n))
                J.append(('uneq-mid', p, k1, k2, 1, c, False, insts, cap))
    elif name == 'raw':
        for p in (5, 13):
            J.append(('raw-plane', p, 1, 1, 1, 1, True, plane(p), None))
            J.append(('raw-plane', p, 2, 2, 1, 1, True, plane(p), None))
        J.append(('raw-plane', 17, 1, 2, 1, 1, True, plane(17), None))
        J.append(('raw-special', 17, 3, 3, 1, 1, True, special(17, 3, 3, rng, nrand=150), None))
    else:
        raise SystemExit('unknown job')
    return J

def main():
    name = sys.argv[1]
    t0 = time.time()
    summ = {}
    viol, uncert = [], []
    total = 0
    with open(f'rt-{name}.jsonl', 'w') as fh:
        for fam, p, k1, k2, s1, s2, raw, insts, cap in jobs(name):
            key = f'{fam} p={p} k=({k1},{k2}) s=({s1},{s2}) raw={raw}'
            hist = {}
            tj = time.time()
            for (a, b, t) in insts:
                bound = min(k1, k2) + (4 if raw else 2)
                Dmax = bound + 1 if cap is None else min(bound + 1, cap)
                res, bound, z = exact(p, a, b, t, k1, k2, s1, s2, raw, Dmax=Dmax)
                lfd = res['lfd']
                rec = dict(fam=fam, p=p, a=a, b=b, t=t, k1=k1, k2=k2, s1=s1, s2=s2, raw=raw, Z=z, lfd=lfd,
                           cert=res['certified_at'], bound=bound, Dmax=Dmax)
                fh.write(json.dumps(rec) + '\n')
                total += 1
                hk = 'uncertified' if lfd is None else str(lfd)
                hist[hk] = hist.get(hk, 0) + 1
                if lfd is None:
                    uncert.append(rec)
                elif lfd > bound:
                    viol.append(rec)
                    print('VIOLATION', rec, flush=True)
            summ[key] = dict(instances=len(insts), bound=bound, lfd_hist=hist, seconds=round(time.time() - tj, 1))
            print(key, summ[key], flush=True)
            fh.flush()
    out = dict(job=name, total=total, violations=viol, uncertified=uncert, seconds=round(time.time() - t0, 1), families=summ)
    json.dump(out, open(f'rt-{name}-summary.json', 'w'), indent=1)
    print('DONE', name, 'total', total, 'violations', len(viol), 'uncertified', len(uncert), 'seconds', out['seconds'], flush=True)

if __name__ == '__main__':
    main()
