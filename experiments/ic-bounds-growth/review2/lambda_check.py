#!/usr/bin/env python3
"""Independent checker for Lemma 4.1 of docs/proofs/N8-growth.md (review 2, 2026-11-12).

Written from the definitions only; it imports nothing from the repository.

Definitions checked:
  C(n, <=w)  = sum_{i<=w} binom(n, i)            (computed by Pascal rows, cross-checked by math.comb)
  Lambda(w)  = largest l >= 1 with 2^(l-1) < 7776 * C(2l, <=w)^2
  k_D        = D + Lambda(D-1)

Tail argument (independent of the ratio argument of the written proof and of the
Lean Chernoff argument): C(n, <=w) <= (n+1)^w.  Let L0 be an integer with
  (a) 7776 * (2*L0+1)^(2w) <= 2^(L0-1)  and
  (b) (2*L0+3)^(2w) <= 2 * (2*L0+1)^(2w).
Then g(l) = 7776*(2l+1)^(2w)/2^(l-1) is non-increasing for l >= L0 (the ratio
((2l+3)/(2l+1))^(2w)/2 decreases in l and is <= 1 at L0), so g(l) <= 1 and the
defining inequality of Lambda fails for every l >= L0.  Every l < L0 is tested
exactly with integers.

Usage: python3 lambda_check.py [WMAX] [OUT.json]
"""
import json
import math
import sys


def pascal_csum_table(nmax, wmax):
    """csum[n][w] = C(n, <=w) for n <= nmax, w <= wmax, by Pascal's rule (exact ints)."""
    row = [1] + [0] * wmax          # binom(0, i)
    table = []
    for n in range(nmax + 1):
        if n > 0:
            new = [1] + [row[i] + row[i - 1] for i in range(1, wmax + 1)]
            row = new
        acc, cs = 0, []
        for i in range(wmax + 1):
            acc += row[i]
            cs.append(acc)
        table.append(cs)
    return table


def tail_start(w):
    """Smallest L0 >= 1 satisfying (a) and (b) of the module docstring."""
    L0 = 1
    while True:
        a = 7776 * (2 * L0 + 1) ** (2 * w) <= 2 ** (L0 - 1)
        b = (2 * L0 + 3) ** (2 * w) <= 2 * (2 * L0 + 1) ** (2 * w)
        if a and b:
            return L0
        L0 += 1


def main():
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    out = sys.argv[2] if len(sys.argv) > 2 else None
    starts = {w: tail_start(w) for w in range(0, wmax + 1)}
    nmax = 2 * max(starts.values()) + 4
    cs = pascal_csum_table(nmax, wmax)
    # cross-check the Pascal table against math.comb on a sample
    for n in range(0, nmax + 1, 7):
        for w in range(0, wmax + 1, 5):
            assert cs[n][w] == sum(math.comb(n, i) for i in range(w + 1)), (n, w)

    def holds(w, l):
        return 2 ** (l - 1) < 7776 * cs[2 * l][w] ** 2

    lam = {}
    report = {"wmax": wmax, "rows": [], "failures": []}
    for w in range(0, wmax + 1):
        L0 = starts[w]
        sat = [l for l in range(1, L0) if holds(w, l)]
        assert sat and sat[0] == 1
        L = max(sat)
        lam[w] = L
        # the satisfying set restricted to [4w, L0) is an initial segment [4w, L]?
        seg = [l for l in range(max(4 * w, 1), L0) if holds(w, l)]
        initial = seg == list(range(max(4 * w, 1), L + 1))
        # rho strictly increasing on [4w, L0]: 2*C(2l)^2 > C(2l+2)^2
        incr = all(2 * cs[2 * l][w] ** 2 > cs[2 * l + 2][w] ** 2 for l in range(max(4 * w, 1), L0 + 1))
        all_below = [l for l in range(1, L + 1) if not holds(w, l)]
        row = {"w": w, "Lambda": L, "tail_start_L0": L0, "le_13w_23": L <= 13 * w + 23,
               "ge_4w": L >= 4 * w, "initial_segment_from_4w": initial,
               "rho_increasing_from_4w": incr, "gaps_below_Lambda": all_below}
        report["rows"].append(row)
        for key in ("le_13w_23", "initial_segment_from_4w", "rho_increasing_from_4w"):
            if not row[key]:
                report["failures"].append((w, key))
        if w >= 1 and not row["ge_4w"]:
            report["failures"].append((w, "ge_4w"))

    # k_D and the linear bound
    kD = {D: D + lam[D - 1] for D in range(1, wmax + 2)}
    report["k_D"] = {D: kD[D] for D in range(2, wmax + 2)}
    bad = [D for D in range(1, wmax + 2) if kD[D] > 14 * D + 10]
    report["k_D_le_14D_10_failures"] = bad
    doc_table = {2: (25, 27), 3: (36, 39), 4: (48, 52), 5: (59, 64), 6: (71, 77), 7: (83, 90),
                 8: (95, 103), 9: (106, 115), 10: (118, 128), 12: (142, 154), 16: (189, 205),
                 20: (237, 257)}
    report["doc_table_mismatches"] = [(D, v, (lam[D - 1], kD[D])) for D, v in doc_table.items()
                                      if (lam[D - 1], kD[D]) != v]
    report["Lambda_1_to_5"] = [lam[w] for w in range(1, 6)]
    assert report["Lambda_1_to_5"] == [25, 36, 48, 59, 71], report["Lambda_1_to_5"]

    # Corollary 5.2 (b): the linear form never exceeds the exact form; LFD >= 8 from k = 90
    lin_bad, exact_lb = [], {}
    for k in range(27, 14 * (wmax + 1) + 10):
        Ds = [D for D in range(2, wmax + 2) if kD[D] <= k]
        ex = 1 + max(Ds) if Ds else None
        exact_lb[k] = ex
        if k >= 38 and (k - 10) // 14 + 1 > ex:
            lin_bad.append(k)
    report["linear_form_exceeds_exact_form_at"] = lin_bad
    report["first_k_with_LFD_lb_8"] = min(k for k, v in exact_lb.items() if v is not None and v >= 8)

    # numerical constants quoted in the written proof of Lemma 4.1 and Remark (iii)
    l7776 = math.log2(7776)
    diffs_quoted = [13 * w + 10.07 - 2 * w * math.log2(5.437 * (13 + 24 / w)) for w in range(1, 8)]
    diffs_exact = [13 * w + 23 - l7776 - 2 * w * math.log2(2 * math.e * (13 + 24 / w)) for w in range(1, 8)]
    report["lemma41_constants"] = {
        "log2_7776": l7776, "two_e": 2 * math.e,
        "diffs_quoted_form_w1_7": [round(x, 3) for x in diffs_quoted],
        "diffs_exact_form_w1_7": [round(x, 3) for x in diffs_exact],
        "two_log2_5.437x16": 2 * math.log2(5.437 * 16),
        "two_log2_2e_x16": 2 * math.log2(2 * math.e * 16),
        "(8/7)^4": (8 / 7) ** 4, "2*(7/8)^4": 2 * (7 / 8) ** 4,
    }
    # direct exact check that 7776*C(2l,<=w)^2 <= 2^(l-1) at l = 13w+24 for w <= wmax
    report["rho_at_13w_24_ge_1"] = all(7776 * cs[2 * (13 * w + 24)][w] ** 2 <= 2 ** (13 * w + 23)
                                       for w in range(0, wmax + 1) if 2 * (13 * w + 24) <= nmax)

    # asymptotic constant c: 4*H2(1/(2c)) = 1
    def H2(x):
        return -x * math.log2(x) - (1 - x) * math.log2(1 - x)
    lo, hi = 2.0, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 4 * H2(1 / (2 * mid)) > 1:
            lo = mid
        else:
            hi = mid
    report["c_root"] = lo
    report["slope_1_over_1_plus_c"] = 1 / (1 + lo)
    report["kD_over_D_near_60"] = {D: round(kD[D] / D, 3) for D in range(55, min(62, wmax + 2))}

    # Lean numeric steps (GrowthCounting.numeric_base / numeric_all) and the Chernoff estimate
    report["lean_steps"] = {
        "17^52*16^2 <= 2^13*16^52": 17 ** 52 * 16 ** 2 <= 2 ** 13 * 16 ** 52,
        "17^4 <= 2*16^4": 17 ** 4 <= 2 * 16 ** 4,
        "base w=0: 15552*17^96 <= 2^24*16^96": 15552 * 17 ** 96 <= 2 ** 24 * 16 ** 96,
        "chernoff C(n,<=w)*16^n <= 17^n*16^w (n<=nmax, w<=wmax)":
            all(cs[n][w] * 16 ** n <= 17 ** n * 16 ** w for n in range(nmax + 1) for w in range(wmax + 1)),
    }

    print("Lambda(1..5) =", report["Lambda_1_to_5"])
    print("Lambda(w), w=1..%d:" % wmax, [lam[w] for w in range(1, wmax + 1)])
    print("k_D, D=2..%d:" % (wmax + 1), [kD[D] for D in range(2, wmax + 2)])
    print("max over w of Lambda(w) - (13w+23):", max(lam[w] - 13 * w - 23 for w in range(1, wmax + 1)))
    print("tight cases Lambda(w) = 13w+23:", [w for w in range(0, wmax + 1) if lam[w] == 13 * w + 23])
    print("failures:", report["failures"])
    print("k_D <= 14D+10 failures:", bad)
    print("doc table mismatches:", report["doc_table_mismatches"])
    print("linear form exceeds exact form at k:", lin_bad[:10])
    print("first k with exact lower bound LFD >= 8:", report["first_k_with_LFD_lb_8"])
    print("Lemma 4.1 constants:", json.dumps(report["lemma41_constants"]))
    print("rho(13w+24) >= 1 exactly for all w checked:", report["rho_at_13w_24_ge_1"])
    print("c =", report["c_root"], "slope =", report["slope_1_over_1_plus_c"])
    print("k_D/D near 60:", report["kD_over_D_near_60"])
    print("Lean numeric steps:", report["lean_steps"])
    gaps = {r["w"]: r["gaps_below_Lambda"] for r in report["rows"] if r["gaps_below_Lambda"]}
    print("w with some l < Lambda(w) failing the inequality:", gaps)
    if out:
        with open(out, "w") as fh:
            json.dump(report, fh, indent=1, default=str)
    ok = (not report["failures"] and not bad and not report["doc_table_mismatches"] and not lin_bad
          and report["rho_at_13w_24_ge_1"] and all(report["lean_steps"].values()))
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
