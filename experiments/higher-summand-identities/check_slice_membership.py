#!/usr/bin/env python3
"""Check the one new ingredient of the all-m growth corollary: after fixing the
abscissae of summands 3..m to q_3..q_m, the specialized summation polynomial
    s(x1,x2) = S_{m+1}(x1, x2, q_3, ..., q_m, t)
vanishes at (x(P), x(R'-P)) for every point P, where R' = R - Q_3 - ... - Q_m
and x(Q_i) = q_i. Hence s lies in the prime ideal p_{R'} of the two-summand
curve model with target R'.

Cases (random curves y^2 = x^3 + a x + b over primes p, rational points only):
  PAIR   m even, Q's in cancelling pairs (Q,-Q): R' = R; q arbitrary.
  ODD    m odd, pairs plus one rational Q: R' = R - Q.
  GEN    arbitrary rational Q_3..Q_m: R' = R - sum Q_i.
  CTRL   negative control: R'' = R' + W for a random point W (wrong target)
         must give s != 0 for most P (rate required > 0.9). (R' + 2Q_3 is NOT a
         valid control: it is the target of another sign pattern, and s
         vanishes there too, which a first version of this script showed.)
Evaluator: check_identities.S (formal-degree Sylvester recursion).
"""
import random, sys, json
from check_identities import S, ec_add, inv

def is_prime(n):
    if n < 2: return False
    for d in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % d == 0: return n == d
    dd, s = n - 1, 0
    while dd % 2 == 0: dd //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, dd, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def sqrt_mod(r, p):
    if r == 0: return 0
    if pow(r, (p - 1) // 2, p) != 1: return None
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, x = s, pow(z, q, p), pow(r, q, p), pow(r, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, x = i, b * b % p, t * b * b % p, x * b % p
    return x

def rand_point(a, b, p, rng):
    while True:
        x = rng.randrange(p); y = sqrt_mod((x ** 3 + a * x + b) % p, p)
        if y is not None: return (x, y if rng.random() < .5 else (p - y) % p)

def neg(P, p): return None if P is None else (P[0], (-P[1]) % p)

def main(seed=7):
    rng = random.Random(seed); res = {}
    def rec(k, ok): r = res.setdefault(k, [0, 0]); r[0] += ok; r[1] += 1
    for p in (10007, 65537, 1000003, 2**31 - 1):
        for _ in range(4):
            a, b = rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b ** 2) % p == 0: continue
            R = rand_point(a, b, p, rng); t = R[0]
            for m in range(3, 8):
                for case in ("PAIR", "ODD", "GEN"):
                    if case == "PAIR" and m % 2: continue
                    if case == "ODD" and m % 2 == 0: continue
                    if case == "GEN":
                        Qs = [rand_point(a, b, p, rng) for _ in range(m - 2)]
                    else:
                        Qs = []
                        if case == "ODD": Qs.append(rand_point(a, b, p, rng))
                        while len(Qs) < m - 2:
                            Q = rand_point(a, b, p, rng); Qs += [Q, neg(Q, p)]
                    Rp = R
                    for Q in Qs: Rp = ec_add(Rp, neg(Q, p), a, p)
                    if Rp is None: continue
                    qs = [Q[0] for Q in Qs]
                    for _k in range(3):
                        P = rand_point(a, b, p, rng)
                        D = ec_add(Rp, neg(P, p), a, p)
                        if D is None: continue
                        rec(case, int(S([P[0], D[0]] + qs + [t], a, b, p) == 0))
                    # negative control: wrong slice target
                    Rw = ec_add(Rp, rand_point(a, b, p, rng), a, p)
                    if Rw is None: continue
                    P = rand_point(a, b, p, rng); D = ec_add(Rw, neg(P, p), a, p)
                    if D is None: continue
                    rec("CTRL_nonzero", int(S([P[0], D[0]] + qs + [t], a, b, p) != 0))
    out = {k: {"passed": v[0], "total": v[1]} for k, v in res.items()}
    print(json.dumps(out, indent=1))
    ok = all(res[k][0] == res[k][1] for k in ("PAIR", "ODD", "GEN")) and \
        res["CTRL_nonzero"][0] > 0.9 * res["CTRL_nonzero"][1]
    print("ALL PASS" if ok else "FAILURE"); return 0 if ok else 1

if __name__ == "__main__": sys.exit(main())
