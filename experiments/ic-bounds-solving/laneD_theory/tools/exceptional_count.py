#!/usr/bin/env python3
"""Count exceptional targets t (unit instance with top coefficient of g^{-1} zero)
for fixed (a, b), m = 2, using the character sum
    c_top(t) = N^{-2} sum_{alpha in H^2} alpha_1 alpha_2 / S3(alpha_1, alpha_2, t),
and compare with the proved bound  #exceptional <= N^2 - 1  (Section 4 of SD1)
when c_top is not identically zero. Pure Python, single thread."""
import random, sys

def s3(x1, x2, t, a, b, p):
    return ((x1 - x2) ** 2 * t * t - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t
            + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p

def run(k, p, ncurves, rng):
    N = 2 ** k
    z = next(z for z in range(2, p) if pow(z, N, p) == 1 and pow(z, N // 2, p) != 1)
    H = [pow(z, j, p) for j in range(N)]
    worst = 0; tot_unit = tot_exc = 0
    for _ in range(ncurves):
        while True:
            a, b = rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p: break
        exc = unit = 0
        for t in range(p):
            s = 0; ok = True
            for x1 in H:
                for x2 in H:
                    v = s3(x1, x2, t, a, b, p)
                    if v == 0: ok = False; break
                    s = (s + x1 * x2 * pow(v, p - 2, p)) % p
                if not ok: break
            if not ok: continue
            unit += 1
            if s == 0: exc += 1
        tot_unit += unit; tot_exc += exc
        worst = max(worst, exc)
        print(f"k={k} p={p} a={a} b={b} unit_t={unit} exceptional_t={exc} bound={N*N-1}")
    print(f"SUMMARY k={k} p={p} curves={ncurves} unit={tot_unit} exceptional={tot_exc} "
          f"frac={tot_exc/max(tot_unit,1):.4f} 1/p={1/p:.4f} worst_per_curve={worst} bound={N*N-1}")
    return worst <= N * N - 1

rng = random.Random(15112026)
ok = all([run(1, 101, 6, rng), run(2, 97, 6, rng), run(2, 113, 4, rng), run(3, 97, 3, rng)])
sys.exit(0 if ok else 1)
