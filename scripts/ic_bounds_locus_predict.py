"""Predict, then test, the exceptional locus on fresh (held-out) curves (bounds item 2).

predict: evaluate frozen descriptions (unions of conjunctions of feature literals
         from csd.ic_bounds_locus.features) on every rational target of a curve and
         write the predicted deviating targets, with the sha256 of the frozen rules,
         BEFORE the closures of that curve are computed.
check:   compare a prediction file with the computed closures
         (scripts/ic_bounds_locus.py output) and report hits, misses, false alarms.

Usage
  python3 scripts/ic_bounds_locus_predict.py predict --rules RULES.json --prime P --k K --curve a,b --out PRED.json
  python3 scripts/ic_bounds_locus_predict.py check --pred PRED.json --data DATA.jsonl --out CHECK.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from csd import ic_dreg as I  # noqa: E402
from csd.ic_bounds_locus import features  # noqa: E402
from ic_bounds_locus import rational_targets  # noqa: E402
from ic_bounds_locus_classify import summarize, deviations  # noqa: E402


def poly_hit(rule, p, a, b, t):
    """Rule on the normalized plane: t != 0 and P(a/t^2, b/t^3) = 0 over F_p."""
    spec = rule['polynomial_AB']
    assert spec['prime'] == p
    if t % p == 0:
        return False
    it = pow(t, -1, p)
    A, B = a * it * it % p, b * it ** 3 % p
    return sum(int(c) * pow(A, int(key.split(',')[0]), p) * pow(B, int(key.split(',')[1]), p)
               for key, c in spec['coefficients'].items()) % p == 0


_PLANES = {}


def plane_hit(rule, p, a, b, t):
    """Rule by lookup in a frozen plane sweep: t != 0 and the normalized point (a/t^2, b/t^3)
    has V_3 = 2k+1 and dim V_4 above `min_V4` in the plane files (t = 1 sheet)."""
    spec = rule['plane_lookup']
    assert spec['prime'] == p
    key = spec['glob']
    if key not in _PLANES:
        import glob as _glob
        table = {}
        for f in sorted(_glob.glob(str(ROOT / key))):
            for line in open(f):
                r = json.loads(line)
                table[(r['A'], r['B'])] = r
        _PLANES[key] = table
    if t % p == 0:
        return False
    it = pow(t, -1, p)
    r = _PLANES[key].get((a * it * it % p, b * it ** 3 % p))
    if r is None:
        raise KeyError('normalized point not in the plane sweep')
    return r['V3'] == spec['V3'] and r['V4'] > spec['min_V4']


def holds(clause, feats):
    for lit in clause:
        if lit.startswith('not(') and lit.endswith(')'):
            if feats[lit[4:-1]]:
                return False
        elif not feats[lit]:
            return False
    return True


def predict(args):
    rules_text = Path(args.rules).read_text()
    rules = json.loads(rules_text)
    p, k = args.prime, args.k
    a, b = [int(v) for v in args.curve.split(',')]
    E = I.PrimeCurve(p, a, b)
    out = dict(prime=p, k=k, a=a, b=b, rules_file=args.rules,
               rules_sha256=hashlib.sha256(rules_text.encode()).hexdigest(),
               made_at_unix=int(time.time()), predictions={})
    targets = rational_targets(E)
    for name, rule in rules['rules'].items():
        if rule.get('k') not in (None, k):
            continue
        if 'polynomial_AB' in rule:
            hit = [t for t in targets if poly_hit(rule, p, a, b, t)]
        elif 'plane_lookup' in rule:
            hit = [t for t in targets if plane_hit(rule, p, a, b, t)]
        else:
            hit = [t for t in targets if any(holds(c, features(p, a, b, t, k)) for c in rule['clauses'])]
        out['predictions'][name] = dict(key=rule['key'], clauses=rule.get('clauses'),
                                        polynomial_AB=rule.get('polynomial_AB'),
                                        plane_lookup=rule.get('plane_lookup'), predicted_t=hit,
                                        expected_value=rule.get('expected_value'))
    out['targets'] = len(targets)
    Path(args.out).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({n: v['predicted_t'] for n, v in out['predictions'].items()}))


def check(args):
    pred = json.loads(Path(args.pred).read_text())
    k = pred['k']
    recs = [json.loads(line) for line in open(args.data)]
    rows = {}
    for r in recs:
        s = summarize(r)
        rows[r['t']] = (s, deviations(s, k))
    report = dict(prediction_file=args.pred,
                  prediction_sha256=hashlib.sha256(Path(args.pred).read_bytes()).hexdigest(),
                  data_file=args.data, targets_computed=len(rows), targets_predicted_over=pred['targets'], rules={})
    for name, v in pred['predictions'].items():
        key = v['key']
        actual = sorted(t for t, (s, dev) in rows.items() if (bool(dev) if key == 'any' else key in dev))
        predicted = sorted(v['predicted_t'])
        report['rules'][name] = dict(key=key, predicted=predicted, actual=actual,
                                     misses=sorted(set(actual) - set(predicted)),
                                     false_alarms=sorted(set(predicted) - set(actual)),
                                     exact_match=predicted == actual)
        if v.get('expected_value') is not None:
            vals = {t: rows[t][1].get(key) for t in predicted if t in rows}
            report['rules'][name]['values_on_predicted'] = vals
            report['rules'][name]['values_match'] = all(x == v['expected_value'] for x in vals.values())
    Path(args.out).write_text(json.dumps(report, indent=1) + '\n')
    for n, r in report['rules'].items():
        print(n, 'exact' if r['exact_match'] else 'MISMATCH', 'pred', r['predicted'], 'actual', r['actual'],
              r.get('values_match', ''))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('predict')
    p.add_argument('--rules', required=True)
    p.add_argument('--prime', type=int, required=True)
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--curve', required=True)
    p.add_argument('--out', required=True)
    c = sub.add_parser('check')
    c.add_argument('--pred', required=True)
    c.add_argument('--data', required=True)
    c.add_argument('--out', required=True)
    args = ap.parse_args()
    predict(args) if args.cmd == 'predict' else check(args)


if __name__ == '__main__':
    main()
