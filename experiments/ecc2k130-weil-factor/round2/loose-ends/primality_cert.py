#!/usr/bin/env python3
"""Independent primality proofs (pure Python, no PARI) for p2 = 146505763881528721 and 263.

Two independent proofs for p2:
  (A) Pocklington / Lehmer N-1 (Pratt-style, recursive): p2-1 fully factored.
  (B) Morrison N+1 Lucas test: p2+1 fully factored.
Every prime factor used is itself proved: by trial division if < 10^12 (sqrt < 10^6),
otherwise recursively by (A).  Factorizations are checked by multiplication, and the
factor lists are *derived here* by trial division / Pollard rho, not copied from PARI.
Controls: composites (Carmichael 561, strong pseudoprime 3215031751, p2+2, p2*263, 263^2)
must be rejected; a tampered factorization must be rejected.
"""
import math, random, sys, json

def is_prime_td(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    r = math.isqrt(n)
    f = 3
    while f <= r:
        if n % f == 0: return False
        f += 2
    return True

def pollard_rho(n):
    if n % 2 == 0: return 2
    while True:
        c = random.randrange(1, n); x = y = random.randrange(2, n); d = 1
        while d == 1:
            x = (x*x + c) % n; y = (y*y + c) % n; y = (y*y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n: return d

def factor(n, out=None):
    """Return dict prime->exp; primality of each factor PROVED by prove_prime."""
    if out is None: out = {}
    for p in (2, 3, 5, 7, 11, 13):
        while n % p == 0: out[p] = out.get(p, 0) + 1; n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1: continue
        if prove_prime(m)[0]:
            out[m] = out.get(m, 0) + 1; continue
        d = pollard_rho(m); stack += [d, m // d]
    return out

def prove_prime(n, depth=0, log=None):
    """Return (True, certificate) or (False, reason). Never returns True for composites."""
    if n < 2: return (False, "n<2")
    if n < 10**12:
        return (is_prime_td(n), {"n": n, "method": "trial division"})
    # Pocklington with N-1 fully factored (Lehmer): for each q | n-1 need a_q with
    # a_q^(n-1) = 1 and gcd(a_q^((n-1)/q) - 1, n) = 1.
    F = factor(n - 1)
    assert math.prod(q**e for q, e in F.items()) == n - 1
    wit = {}
    for q in F:
        for a in range(2, 1000):
            if pow(a, n - 1, n) != 1: return (False, f"Fermat fails base {a}")
            if math.gcd(pow(a, (n - 1)//q, n) - 1, n) == 1:
                wit[q] = a; break
        else:
            return (False, f"no witness for q={q}")
    sub = {q: prove_prime(q, depth+1)[1] for q in F}
    return (True, {"n": n, "method": "Lehmer N-1", "factors": F, "witnesses": wit, "sub": sub})

def check_lehmer(n, F, wit):
    """Re-check a given N-1 certificate (factor primality checked recursively)."""
    if math.prod(q**e for q, e in F.items()) != n - 1: return False
    for q, e in F.items():
        if not prove_prime(q)[0]: return False
        a = wit[q]
        if pow(a, n - 1, n) != 1 or math.gcd(pow(a, (n - 1)//q, n) - 1, n) != 1: return False
    return True

# ---------------- Morrison N+1 ----------------
def jacobi(a, n):
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0

def lucas_U(P, Q, k, n):
    """U_k(P,Q) mod n via the matrix [[P,-Q],[1,0]]^k."""
    def mul(A, B):
        return [[(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % n, (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % n],
                [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % n, (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % n]]
    R = [[1, 0], [0, 1]]; M = [[P % n, (-Q) % n], [1, 0]]
    while k:
        if k & 1: R = mul(R, M)
        M = mul(M, M); k >>= 1
    return R[1][0]  # M^k = [[U_{k+1}, -Q U_k],[U_k, -Q U_{k-1}]]

def morrison(n, F):
    """N+1 proof: N+1 = prod F fully factored; for each q find (P,Q) with (D/n) = -1,
    gcd(n, 2QD) = 1, U_{n+1} = 0 and gcd(U_{(n+1)/q}, n) = 1 (Morrison 1975)."""
    if math.prod(q**e for q, e in F.items()) != n + 1: return (False, "bad factorization")
    for q in F:
        if not prove_prime(q)[0]: return (False, f"factor {q} not prime")
    # One common discriminant D for all q (needed: the sign (D/p) must be the same for every q).
    for P in range(1, 60):
        for Q in range(-30, 30):
            D = P*P - 4*Q
            if Q == 0 or D == 0 or math.gcd(n, 2*Q*D) != 1 or jacobi(D, n) != -1: continue
            if lucas_U(P, Q, n + 1, n) != 0: return (False, f"U_(n+1) != 0 for P={P},Q={Q}")
            if all(math.gcd(lucas_U(P, Q, (n + 1)//q, n), n) == 1 for q in F):
                return (True, {"P": P, "Q": Q, "D": D})
    return (False, "no common Lucas witness")

if __name__ == "__main__":
    random.seed(1)
    p2 = 146505763881528721
    res = {}
    # sanity of the Lucas implementation: U_n(1,-1) = Fibonacci
    assert [lucas_U(1, -1, k, 10**9) for k in range(10)] == [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
    ok, cert = prove_prime(p2); res["p2_lehmer"] = ok
    print("p2 N-1 factorization:", cert["factors"], " witnesses:", cert["witnesses"])
    print("p2 Lehmer N-1 proof:", ok, " recheck:", check_lehmer(p2, cert["factors"], cert["witnesses"]))
    Fp = factor(p2 + 1)
    ok2, w2 = morrison(p2, Fp); res["p2_morrison"] = ok2
    print("p2 N+1 factorization:", Fp, " Lucas witnesses:", w2)
    print("p2 Morrison N+1 proof:", ok2)
    ok3, c3 = prove_prime(263); res["263_td"] = ok3
    F263 = {2: 1, 131: 1}
    res["263_lehmer"] = check_lehmer(263, F263, {2: 5, 131: 5})
    print("263 trial division:", ok3, " Lehmer N-1 (262=2*131, base 5):", res["263_lehmer"])
    # negative controls
    neg = {}
    for c in (561, 3215031751, p2 + 2, p2 * 263, 263**2, 4294967297, 2**61 - 3):
        if c < 10**12:
            neg[c] = prove_prime(c)[0]
        else:
            try: neg[c] = prove_prime(c)[0]
            except AssertionError: neg[c] = False
    print("negative controls (must all be False):", neg)
    # tampered certificate: replace a factor of p2-1 by a wrong one
    badF = dict(cert["factors"]); badF.pop(289510489); badF[289510487] = 1
    tamper = check_lehmer(p2, badF, {**cert["witnesses"], 289510487: 2})
    print("tampered N-1 cert accepted? (must be False):", tamper)
    tamperL = morrison(p2, {2: 1, 11: 1, 131: 1, 1777: 1, 28607065271: 1})[0]
    print("tampered N+1 cert accepted? (must be False):", tamperL)
    # positive control: a known 64-bit prime 2^61-1 must be proved
    pc = prove_prime(2**61 - 1)[0]; pc2 = morrison(2**61 - 1, factor(2**61))[0]
    print("positive control 2^61-1 (N-1, N+1):", pc, pc2)
    allgood = (res["p2_lehmer"] and res["p2_morrison"] and res["263_td"] and res["263_lehmer"]
               and not any(neg.values()) and not tamper and not tamperL and pc and pc2)
    print("ALL CHECKS PASS:", allgood)
    json.dump({"p2": p2, "N-1": {str(k): v for k, v in cert["factors"].items()},
               "N-1_witnesses": {str(k): v for k, v in cert["witnesses"].items()},
               "N+1": {str(k): v for k, v in Fp.items()},
               "N+1_witnesses": {str(k): v for k, v in w2.items()},
               "results": res, "negative_controls": {str(k): v for k, v in neg.items()},
               "all_pass": allgood}, open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w"), indent=1)
