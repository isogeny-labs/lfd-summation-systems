#!/usr/bin/env python3
"""Refutation attempt for Corollary 2 (direct m=3, LFD >= 7 at k >= 3 when 5t^3+3at+6b != 0) and
for the 'LFD = 7' pattern, on primes/curves/targets OUTSIDE the producer's grid.  Seed 424242.
  python3 scripts/ic_bounds_m3_verify_newprimes.py [--k4]   -> verify/newprimes.jsonl
"""
import json, random, time, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_bounds_m3_verify_engine as V
from ic_bounds_m3_verify_symbolic import sqrt_mod

OUT = Path(__file__).resolve().parents[1] / 'experiments' / 'ic-bounds-m3' / 'verify' / 'newprimes.jsonl'


def main():
    k4 = '--k4' in sys.argv
    rng = random.Random(424242 + (1 if k4 else 0))
    terms = V.s4_terms()
    plan = [(113, 4), (193, 4), (241, 4)] if k4 else [(41, 3), (73, 3), (89, 3), (137, 3), (233, 3), (409, 3)]
    per = 1 if k4 else 4
    done = set()
    if OUT.exists():
        for l in open(OUT): d = json.loads(l); done.add((d['p'], d['a'], d['b'], d['t'], d['k']))
    with open(OUT, 'a') as f:
        for p, k in plan:
            n = 0
            while n < per:
                a, b = rng.randrange(p), rng.randrange(p)
                if (4 * a**3 + 27 * b * b) % p == 0: continue
                xs = [x for x in range(p) if pow((x**3 + a * x + b) % p, (p - 1) // 2, p) == 1]
                if not xs: continue
                t = rng.choice(xs)
                if (5 * t**3 + 3 * a * t + 6 * b) % p == 0: continue
                if (p, a, b, t, k) in done: n += 1; continue
                t0 = time.time()
                r = V.exact_profile('direct', p, a, b, t, k, terms, Dmax=3 * k + 1)
                r['seconds'] = round(time.time() - t0, 2); r['seed'] = 424242 + (1 if k4 else 0)
                r['y'] = sqrt_mod((t**3 + a * t + b) % p, p)
                f.write(json.dumps(r) + '\n'); f.flush()
                print(p, k, a, b, t, r['status'], r.get('lfd'), r['input_degrees'], r['dimA'], r['seconds'], flush=True)
                n += 1


if __name__ == '__main__':
    main()
