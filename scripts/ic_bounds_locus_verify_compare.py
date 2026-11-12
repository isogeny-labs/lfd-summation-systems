#!/usr/bin/env python3
"""Verifier step 6: compare the verifier's recomputations (experiments/ic-bounds-locus/verify/*.jsonl,
written by ic_bounds_locus_verify.py) with the producer's rows, and test the (A, B) invariance
with verifier-computed values on both sides where available."""
import glob
import json
import os
import collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXP = os.path.join(ROOT, "experiments", "ic-bounds-locus")
V = os.path.join(EXP, "verify")


def load_jsonl(pattern):
    out = []
    for fn in sorted(glob.glob(pattern)):
        out += [json.loads(l) for l in open(fn)]
    return out


def sig_mine(r):
    per = r["per"]
    s = {}
    for D in ("3", "4"):
        if D in per:
            s[f"V{D}"] = per[D]["dim"]
            s[f"low{D}"] = per[D]["low"]
            s[f"one{D}"] = per[D]["one"]
            s[f"sat{D}"] = per[D]["dim"] == per[D]["dimI"]
    return s


def main():
    prod = {}
    for sub in ("data", "heldout"):
        for r in load_jsonl(os.path.join(EXP, sub, "p*-k*-*.jsonl")):
            prod[(r["prime"], r["k"], r["a"], r["b"], r["t"])] = r
    plane = {}
    for r in load_jsonl(os.path.join(EXP, "plane", "p193-k6-plane-c*.jsonl")):
        plane[(r["A"], r["B"])] = r
    res = dict(targets={}, plane={}, invariance={})
    mine_t = load_jsonl(os.path.join(V, "p193-k6-*.jsonl")) + load_jsonl(os.path.join(V, "p257-k7.jsonl"))
    mism = []
    lfd_hist = collections.Counter()
    closed_all = True
    n = 0
    certs = collections.Counter()
    for m in mine_t:
        key = (m["prime"], m["k"], m["a"], m["b"], m["t"])
        if key not in prod:
            continue
        n += 1
        pr = prod[key]
        d = {x["D"]: x for x in pr["rows"]}
        for D, v in m["per"].items():
            D = int(D)
            closed_all &= v["closed"]
            if D in d and (d[D]["dim"], d[D]["low"], d[D]["one"]) != (v["dim"], v["low"], v["one"]):
                mism.append([list(key), D, [d[D]["dim"], d[D]["low"], d[D]["one"]], [v["dim"], v["low"], v["one"]]])
        if m["Z"] != pr["Z"]:
            mism.append([list(key), "Z", pr["Z"], m["Z"]])
        if m["k"] == 6:
            certs[str(m["cert"])[:10]] += 1
            lfd_hist[(m["lfd"], pr.get("exact_lfd"))] += 1
            if m["lfd"] != pr.get("exact_lfd"):
                mism.append([list(key), "lfd", pr.get("exact_lfd"), m["lfd"]])
            falls = m["falls"]
            if falls != pr["falls"]:
                mism.append([list(key), "falls", pr["falls"], falls])
        else:
            if m["lfd"] is not None or pr.get("exact_lfd") is not None:
                lfd_hist[(m["lfd"], pr.get("exact_lfd"))] += 1
    res["targets"] = dict(compared=n, mismatches=mism, all_closed=closed_all,
                          lfd_mine_vs_producer={str(k): v for k, v in lfd_hist.items()}, k6_certificates=dict(certs))
    # plane sample
    mine_p = load_jsonl(os.path.join(V, "plane-*.jsonl"))
    pm = []
    cls = collections.Counter()
    myplane = {}
    for m in mine_p:
        if m["t"] != 1:
            continue
        s = sig_mine(m)
        myplane[(m["prime"], m["k"], m["a"], m["b"])] = (s, m["Z"])
        if m["prime"] != 193 or m["k"] != 6:
            continue
        pr = plane[(m["a"], m["b"])]
        ok = (pr["V3"], pr["low3"], pr["V4"], pr["low4"], pr["one4"]) == (s["V3"], s["low3"], s["V4"], s["low4"], s["one4"])
        ok &= bool(pr["saturated4"]) == s["sat4"]
        ok &= pr["Z"] == m["Z"]
        if not ok:
            pm.append([m["a"], m["b"], pr, s])
        cls[("sat",) if s["sat4"] else (s["V3"], s["low3"], s["V4"], s["low4"])] += 1
    res["plane"] = dict(compared=len([m for m in mine_p if m["prime"] == 193 and m["k"] == 6]), mismatches=pm,
                        classes_recomputed=[[list(k), v] for k, v in cls.most_common()])
    # invariance: verifier target vs verifier plane point (same p, k), else vs producer plane (k = 6)
    inv = collections.Counter()
    bad = []
    for m in mine_t:
        p, k, a, b, t = m["prime"], m["k"], m["a"], m["b"], m["t"]
        if t == 0:
            continue
        A = a * pow(t, -2, p) % p
        B = b * pow(t, -3, p) % p
        H = pow(t, 1 << k, p) == 1
        s = sig_mine(m)
        if (p, k, A, B) in myplane:
            s2, z2 = myplane[(p, k, A, B)]
            src = "verifier-plane"
        elif k == 6:
            pr = plane[(A, B)]
            s2 = dict(V3=pr["V3"], low3=pr["low3"], V4=pr["V4"], low4=pr["low4"], one4=pr["one4"], sat4=bool(pr["saturated4"]))
            src = "producer-plane"
        else:
            continue
        same = (s["sat4"] and s2["sat4"] and s["V3"] == s2["V3"]) or all(s[x] == s2[x] for x in ("V3", "low3", "V4", "low4", "one4"))
        inv[(k, src, "t in H" if H else "t not in H", same)] += 1
        if not same:
            bad.append([p, k, a, b, t, A, B, s, s2])
    res["invariance"] = dict(counts={str(k): v for k, v in inv.items()}, failures=bad)
    json.dump(res, open(os.path.join(V, "compare.json"), "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str)[:5000])


if __name__ == "__main__":
    main()
