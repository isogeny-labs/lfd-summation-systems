#!/usr/bin/env python3
"""Verifier step 4: re-derive the held-out prediction scores with own code.
For each held-out curve: (i) re-apply the frozen rules (rules-k6-plane.json, rules-k7.json)
with own code to get the predicted target sets, compare with the producer's frozen
prediction files; (ii) score them against the held-out per-target rows: deviating set =
signature (V3, low3, V4, low4) != generic; value match per rule key.  Also records the
prediction files' made_at_unix against the held-out files' first/last log lines' mtimes."""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXP = os.path.join(ROOT, "experiments", "ic-bounds-locus")
GENERIC = {6: (13, 1, 605, 111), 7: (15, 1, 829, 115)}


def ev(coeffs, A, B, p):
    s = 0
    for key, c in coeffs.items():
        i, j = map(int, key.split(","))
        s += c * pow(A, i, p) * pow(B, j, p)
    return s % p


def main():
    plane = {}
    for fn in glob.glob(os.path.join(EXP, "plane", "p193-k6-plane-c*.jsonl")):
        for l in open(fn):
            r = json.loads(l)
            plane[(r["A"], r["B"])] = r
    out = []
    for p, k, rules_fn, pat in ((193, 6, "rules-k6-plane.json", "p193-k6-*-plane.json"), (257, 7, "rules-k7.json", "p257-k7-a*[0-9].json")):
        rules = json.load(open(os.path.join(EXP, rules_fn)))["rules"]
        cub = rules["cubic_C1"]["polynomial_AB"]["coefficients"]
        for pf in sorted(glob.glob(os.path.join(EXP, "predictions", pat))):
            pred = json.load(open(pf))
            a, b = pred["a"], pred["b"]
            data = [json.loads(l) for l in open(os.path.join(EXP, "heldout", f"p{p}-k{k}-a{a}-b{b}.jsonl"))]
            byt = {r["t"]: r for r in data}
            # own application of the rules
            mine = {"t0": set(), "tangent": set(), "order_two": set(), "C1": set(), "C2": set()}
            for t in byt:
                rhs = (t ** 3 + a * t + b) % p
                if t == 0:
                    mine["t0"].add(t)
                    continue
                if (3 * t * t + a) % p == 0:
                    mine["tangent"].add(t)
                elif rhs == 0:
                    mine["order_two"].add(t)
                A = a * pow(t, -2, p) % p
                B = b * pow(t, -3, p) % p
                if ev(cub, A, B, p) == 0:
                    mine["C1"].add(t)
                if p == 193:
                    r = plane[(A, B)]
                    if r["V3"] == 13 and r["V4"] > 700:
                        mine["C2"].add(t)
            mine_all = set().union(*mine.values())
            prod_all = set()
            for name, pr in pred["predictions"].items():
                prod_all |= set(pr["predicted_t"])
            dev = {}
            for t, r in byt.items():
                d = {x["D"]: x for x in r["rows"]}
                s = (d[3]["dim"], d[3]["low"], d[4]["dim"], d[4]["low"])
                if s != GENERIC[k]:
                    dev[t] = dict(V3=d[3]["dim"], V4=d[4]["dim"], lfd=r.get("exact_lfd"))
            # value check: predicted values per own class
            want = {6: dict(tangent=("V3", 34), order_two=("V4", 607), C1=("V4", 602), C2=("lfd", 4)),
                    7: dict(tangent=("V3", 40), order_two=("V4", 835), C1=("V4", 826))}[k]
            val_mis = []
            for cls, ts in mine.items():
                if cls == "t0" or cls not in want:
                    continue
                key, v = want[cls]
                for t in ts:
                    got = dev.get(t, {}).get(key)
                    if got != v:
                        val_mis.append([cls, t, v, got])
            out.append(dict(prime=p, k=k, curve=[a, b], targets=len(byt), deviating=sorted(dev),
                            own_predicted=sorted(mine_all), producer_predicted=sorted(prod_all),
                            same_prediction=mine_all == prod_all, missed=sorted(set(dev) - mine_all),
                            false_alarms=sorted(mine_all - set(dev)), value_mismatches=val_mis,
                            missed_values={t: dev[t] for t in sorted(set(dev) - mine_all)},
                            made_at_unix=pred.get("made_at_unix"),
                            heldout_mtime=int(os.path.getmtime(os.path.join(EXP, "heldout", f"p{p}-k{k}-a{a}-b{b}.jsonl")))))
    tot = {}
    for k in (6, 7):
        rs = [r for r in out if r["k"] == k]
        tot[k] = dict(curves=len(rs), targets=sum(r["targets"] for r in rs), deviating=sum(len(r["deviating"]) for r in rs),
                      predicted=sum(len(r["own_predicted"]) for r in rs), missed=sum(len(r["missed"]) for r in rs),
                      false_alarms=sum(len(r["false_alarms"]) for r in rs),
                      value_mismatches=sum(len(r["value_mismatches"]) for r in rs),
                      all_same_as_producer=all(r["same_prediction"] for r in rs),
                      made_before_data_mtime=all(r["made_at_unix"] < r["heldout_mtime"] for r in rs))
    json.dump(dict(per_curve=out, totals=tot), open(os.path.join(EXP, "verify", "predict-checks.json"), "w"), indent=1)
    for r in out:
        print(r["k"], r["curve"], r["targets"], "dev", r["deviating"], "missed", r["missed_values"], "FA", r["false_alarms"],
              "valmis", r["value_mismatches"], "same", r["same_prediction"], r["made_at_unix"], r["heldout_mtime"])
    print(json.dumps(tot, indent=1))


if __name__ == "__main__":
    main()
