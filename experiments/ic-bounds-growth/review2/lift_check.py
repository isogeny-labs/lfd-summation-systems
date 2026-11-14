#!/usr/bin/env python3
"""Independent brute-force checks of the lift identity (1.1) and the digit lemmas of
docs/proofs/N8-growth.md (review 2, 2026-11-12).  Imports nothing from the repository.

(1.1)  (v f)~ = v f~ - (x_i^N - 1) x_i^(2^j - N) f~_hi,   f~_hi = sum_{a_i >= N - 2^j} c_a x^a,
and its Lean form  x_i^(N-2^j) (v f)~ = x_i^N (f~ - f~_hi) + f~_hi  (GrowthReduction.lift_identity),
with A_k = K[x1,x2]/(x1^N - 1, x2^N - 1), N = 2^k, v = x_i^(2^j), 0 <= j < k.
Laurent polynomials are dicts {(e1, e2): c} with integer (possibly negative) exponents,
coefficients in GF(q).

Usage: python3 lift_check.py [KMAX_IDENTITY] [KMAX_DIGITS] [OUT.json]
"""
import json
import random
import sys

Q = 1000003  # coefficient field GF(Q)


def wt(n):
    return bin(n).count("1")


def add(f, g, s=1):
    h = dict(f)
    for e, c in g.items():
        h[e] = (h.get(e, 0) + s * c) % Q
        if h[e] == 0:
            del h[e]
    return h


def shift(f, i, s):
    """multiply by x_i^s (Laurent)"""
    out = {}
    for e, c in f.items():
        e2 = list(e)
        e2[i] += s
        out[tuple(e2)] = c
    return out


def mul_Ak(f, i, j, N):
    """u_{i,j} * f in A_k, returned as its lift (exponents in [0, N))"""
    out = {}
    for e, c in f.items():
        e2 = list(e)
        e2[i] = (e2[i] + 2 ** j) % N
        e2 = tuple(e2)
        out[e2] = (out.get(e2, 0) + c) % Q
        if out[e2] == 0:
            del out[e2]
    return out


def top_part(f, i, j, N):
    return {e: c for e, c in f.items() if e[i] >= N - 2 ** j}


def check_identity(f, i, j, N):
    lhs = mul_Ak(f, i, j, N)
    hi = top_part(f, i, j, N)
    # v f~ - (x_i^N - 1) x_i^(2^j - N) f~_hi = v f~ - x_i^(2^j) hi + x_i^(2^j - N) hi
    rhs = add(add(shift(f, i, 2 ** j), shift(hi, i, 2 ** j), -1), shift(hi, i, 2 ** j - N))
    ok_written = lhs == rhs
    # Lean form: x_i^(N - 2^j) * lhs == x_i^N (f - hi) + hi, as genuine polynomials
    lean_l = shift(lhs, i, N - 2 ** j)
    lean_r = add(shift(add(f, hi, -1), i, N), hi)
    ok_lean = lean_l == lean_r and all(min(e) >= 0 for e in lean_l)
    # divisibility claim: x_i^(N - 2^j) divides f~_hi
    ok_div = all(e[i] >= N - 2 ** j for e in hi)
    return ok_written, ok_lean, ok_div


def main():
    kmax_id = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    kmax_dig = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    out = sys.argv[3] if len(sys.argv) > 3 else None
    rng = random.Random(20261112)
    rep = {"identity": {}, "digits": {}, "lemma42": {}, "lemma44": {}, "lemma43_count": {}}

    # (1.1): every single monomial (the identity is linear), all i, j, k <= kmax_id
    for k in range(2, kmax_id + 1):
        N = 2 ** k
        cnt = bad = 0
        for a1 in range(N):
            for a2 in range(N):
                f = {(a1, a2): 1 + (a1 * 7 + a2) % (Q - 1)}
                for i in (0, 1):
                    for j in range(k):
                        cnt += 1
                        if not all(check_identity(f, i, j, N)):
                            bad += 1
        # random dense combinations, and boundary-heavy ones (a_i in {N-2^j-1, N-2^j, N-2^j+1, N-1, 0})
        rnd = 0
        for _ in range(200):
            f = {}
            for _ in range(rng.randint(1, 3 * N)):
                f[(rng.randrange(N), rng.randrange(N))] = rng.randrange(1, Q)
            for i in (0, 1):
                for j in range(k):
                    rnd += 1
                    if not all(check_identity(f, i, j, N)):
                        bad += 1
        for i in (0, 1):
            for j in range(k):
                vals = sorted({N - 2 ** j - 1, N - 2 ** j, N - 2 ** j + 1, N - 1, 0} & set(range(N)))
                f = {}
                for v in vals:
                    for o in range(N):
                        e = (v, o) if i == 0 else (o, v)
                        f[e] = rng.randrange(1, Q)
                rnd += 1
                if not all(check_identity(f, i, j, N)):
                    bad += 1
        rep["identity"][k] = {"monomial_cases": cnt, "combination_cases": rnd, "failures": bad}

    # digits: top run <=> a >= N - 2^j <=> a + 2^j >= N ; no double wrap ; weight of a shift
    for k in range(1, kmax_dig + 1):
        N = 2 ** k
        bad_eq = bad_wrap = bad_wt = 0
        for a in range(N):
            for j in range(k):
                toprun = all((a >> b) & 1 for b in range(j, k))
                if not (toprun == (a >= N - 2 ** j) == (a + 2 ** j >= N)):
                    bad_eq += 1
                if a + 2 ** j >= 2 * N:
                    bad_wrap += 1
                if wt((a + 2 ** j) % N) > wt(a) + 1:
                    bad_wt += 1
        rep["digits"][k] = {"toprun_equiv_failures": bad_eq, "double_wrap": bad_wrap,
                            "shift_weight_failures": bad_wt}

    # Lemma 4.2 (window): {a in [0,N) : ||2(a-a0)||_p <= 2 sigma} lies in <= 6 aligned blocks of size 2^l
    def dist(n, p):
        if p == 0:
            return abs(n)
        r = n % p
        return min(r, p - r)

    def primes_from(n, count):
        ps, m = [], n
        while len(ps) < count:
            if m > 1 and all(m % d for d in range(2, int(m ** 0.5) + 1)):
                ps.append(m)
            m += 1
        return ps

    for k in range(2, 8):
        N = 2 ** k
        worst = 0
        cases = 0
        for p in [0] + primes_from(2 * N + 1, 4):
            for sigma in range(0, N + 2):
                l = (2 * sigma).bit_length()  # least l with 2^l >= 2 sigma + 1
                assert 2 ** l >= 2 * sigma + 1 and (l == 0 or 2 ** (l - 1) < 2 * sigma + 1)
                for a0 in range(N):
                    S = [a for a in range(N) if dist(2 * (a - a0), p) <= 2 * sigma]
                    nb = len({a >> l for a in S})
                    worst = max(worst, nb)
                    cases += 1
        rep["lemma42"][k] = {"cases": cases, "max_blocks": worst}

    # Lemma 4.4 (no straddling), exhaustive over a, a' for every lambda, j, sigma (k <= 9)
    for k in range(1, 10):
        N = 2 ** k
        bad1 = bad2 = c1 = c2 = 0
        for lam in range(1, k + 1):
            for sigma in range(0, 2 ** (lam - 1)):
                pmin = 2 * N + 1  # the case p > 0: smallest admissible p makes the condition weakest
                for j in range(lam - 1, k):
                    cut = N - 2 ** j
                    for a in range(cut, N):
                        # first case: a - a' <= sigma with a' < cut
                        for ap in range(max(0, a - sigma), cut):
                            c1 += 1
                            if wt(ap) < k - lam:
                                bad1 += 1
                        # second case: a - a' >= (p - 2 sigma)/2, smallest p
                        for ap in range(0, cut):
                            if 2 * (a - ap) >= pmin - 2 * sigma:
                                c2 += 1
                                if wt(a) < k - lam + 1:
                                    bad2 += 1
        rep["lemma44"][k] = {"first_case_pairs": c1, "first_case_failures": bad1,
                             "second_case_pairs": c2, "second_case_failures": bad2}

    # Lemma 4.3 count: pairs of weight <= w inside one pair of aligned blocks <= C(2l, <= w)
    from math import comb
    for k in range(2, 7):
        N = 2 ** k
        worst_ratio = 0.0
        bad = 0
        for l in range(0, k + 1):
            B = 2 ** l
            for w in range(0, 2 * k + 1):
                bound = sum(comb(2 * l, i) for i in range(w + 1))
                for c1 in range(N // B):
                    for c2 in range(N // B):
                        n = sum(1 for a in range(c1 * B, (c1 + 1) * B) for b in range(c2 * B, (c2 + 1) * B)
                                if wt(a) + wt(b) <= w)
                        if n > bound:
                            bad += 1
        rep["lemma43_count"][k] = {"failures": bad}

    print(json.dumps(rep, indent=1))
    ok = (all(v["failures"] == 0 for v in rep["identity"].values())
          and all(v["toprun_equiv_failures"] == 0 and v["double_wrap"] == 0 and v["shift_weight_failures"] == 0
                  for v in rep["digits"].values())
          and all(v["max_blocks"] <= 6 for v in rep["lemma42"].values())
          and all(v["first_case_failures"] == 0 and v["second_case_failures"] == 0 for v in rep["lemma44"].values())
          and all(v["failures"] == 0 for v in rep["lemma43_count"].values()))
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
    if out:
        with open(out, "w") as fh:
            json.dump(rep, fh, indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
