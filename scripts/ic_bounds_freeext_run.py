"""Free-extension certificates for the two-summand chain closures (D = 3, 4, 5).

Uses csd/ic_bounds_freeext.py (W_{k0+1}(S), closedness sweep, equality, witnesses) on top of the
independent closure engine csd/ic_bounds_d5_verify.py.  Instances are drawn exactly as in
scripts/ic_bounds_d5_verify.py (random.Random(seed), ic_curve_sweep.random_rational_target,
redrawn with --want until the scalar solution set at depth kmax is nonempty), or given by an
explicit target abscissa X (--X, for the exceptional targets).

Subcommands (JSON lines on stdout; records in experiments/ic-bounds-freeext/):
  chain   --p P --a A --b B (--seed S [--want] | --X X) --D D --kmin K0 --kmax K1
          [--ext-top] [--witness] [--nosplit] --name NAME
      computes V_D(k) for k = kmin..kmax; for each k < kmax (and k = kmax with --ext-top) builds
      W_{k+1}(V_D(k)), runs the deterministic closedness sweep on it at depth k+1, checks the
      generator rows lie in it and 1 notin V_D(k), and compares it with V_D(k+1) when computed.
      --witness searches a shifted-generator witness e (tau_e(g) of weight <= D, notin V_D(kmax)).
  exceptional --p P --a A --b B --kind order2|tangent|zero
      prints an abscissa X of that kind (and whether it lifts to a rational point).
"""
import argparse
import hashlib
import importlib.util
import json
import resource
import sys
import time
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402
from csd import ic_bounds_d5_verify as V  # noqa: E402
from csd import ic_bounds_freeext as FX  # noqa: E402

OUT = ROOT / 'experiments/ic-bounds-freeext'
T0 = time.time()
SOURCES = ('csd/ic_bounds_freeext.py', 'csd/ic_bounds_d5_verify.py', 'scripts/ic_bounds_freeext_run.py',
           'scripts/ic_bounds_d5_verify.py')


def log(**kw):
    kw['t'] = round(time.time() - T0, 1)
    kw['rss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(json.dumps(kw), flush=True)


def load_verify_script():
    spec = importlib.util.spec_from_file_location('ic_bounds_d5_verify_script', ROOT / 'scripts/ic_bounds_d5_verify.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def target_X(args):
    if args.X is not None:
        return args.X % args.p, dict(X=args.X % args.p, drawn=False)
    VS = load_verify_script()
    inst = VS.instance(args.p, args.a, args.b, args.seed, args.draw_k or args.kmax, args.want)
    return inst['target'][0], dict(target=inst['target'], attempts=inst['attempts'], drawn=True, draw_k=args.draw_k or args.kmax)


def count(k, D):
    return sum(comb(2 * k, w) for w in range(D + 1))


def ideal_dim(p, a, b, X, k, D):
    """dim(I_k cap M_<=D) = |M_<=D| - rank(evaluation of M_<=D on Z) when 2^k | p-1, else None."""
    if (p - 1) % (1 << k):
        return None, None
    Z = V.solutions(p, a, b, X, k)
    L = V.Layout(k, D, False)
    return count(k, D) - V.eval_rank(L, Z, p, k, D), len(Z)


def space_record(sp):
    C = sp.C
    prof, parts = C.profile()
    return dict(rank=C.rank(), profile=prof, part_profiles=parts, contains_one=C.contains_one(),
                rref_ok=all(P.check_rref() for P in C.parts), basis_sha256=[P.sha() for P in C.parts])


def reference_sha(args, k):
    """recorded V_D(k) of scripts/ic_bounds_d5_verify.py (same engine source) for this instance."""
    if args.ref_name is None:
        return None
    path = ROOT / 'experiments/ic-bounds-d5/verify' / args.ref_name.format(k=k)
    if not path.exists():
        return None
    r = json.loads(path.read_text())
    eng = hashlib.sha256((ROOT / 'csd/ic_bounds_d5_verify.py').read_bytes()).hexdigest()
    if r.get('source_sha256', {}).get('csd/ic_bounds_d5_verify.py') != eng or not r.get('split'):
        return None
    return dict(file=str(path.relative_to(ROOT)), file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                rank=r['rank'], profile=r['profile'], basis_sha256=r['basis_sha256'],
                sweep_full_passed=r.get('sweep_full', {}).get('passed'), target=r.get('target'))


def chain(args):
    p, a, b, D = args.p, args.a % args.p, args.b % args.p, args.D
    X, tinfo = target_X(args)
    split = not args.nosplit
    rec = dict(p=p, a=args.a, b=args.b, seed=args.seed, want=args.want, D=D, split=split, target_info=tinfo,
               exceptional=dict(order_two=(X ** 3 + a * X + b) % p == 0, t_zero=X % p == 0,
                                tangent=(3 * X * X + a) % p == 0),
               closures={}, extensions={})
    log(stage='instance', X=X, **rec['exceptional'])
    prev = None
    for k in range(args.kmin, args.kmax + 1):
        g = V.chain_generator(p, a, b, X, k)
        t = time.time()
        C = V.Closure(g, k, p, D, split=split, gens='all', log=log)
        C.run()
        sp = FX.as_space(C)
        r = space_record(sp)
        r['sweep_full'] = C.sweep_full()
        r['seconds'] = round(time.time() - t, 1)
        r['generator_in_space'] = sp.residual_rows(FX.generator_matrix(g, C.L, D, k, p, 'g')) == 0
        idim, npts = ideal_dim(p, a, b, X, k, D)
        r.update(ideal_dim=idim, points=npts, columns=C.L.c,
                 deficient=None if idim is None else r['rank'] < idim)
        log(stage='closure', k=k, **{x: y for x, y in r.items() if x not in ('basis_sha256', 'part_profiles')})
        rec['closures'][str(k)] = r
        if prev is not None:
            kp, Wp = prev
            e = rec['extensions'][str(kp)]
            if Wp is not None:
                e['equals_V_next'] = FX.same_space(Wp, sp)
                e['equals_V_next_by'] = 'canonical RREF identical (same_space)'
            else:
                e['equals_V_next'] = e['W']['basis_sha256'] == r['basis_sha256']
                e['equals_V_next_by'] = 'SHA-256 of the canonical RREFs (W freed before V_D(k+1))'
            e['W_profile_equals_V_next_profile'] = e['W']['profile'] == r['profile']
            e['recursion_holds'] = FX.profile_recursion(rec['closures'][str(kp)]['profile']) == r['profile']
            log(stage='compare', k0=kp, equals=e['equals_V_next'], recursion=e['recursion_holds'])
            prev = None
            del Wp
        if k < args.kmax or args.ext_top:
            t = time.time()
            W = FX.free_extension(C) if (args.generic_ext or not split) else FX.free_extension_direct(C)
            wr0 = 'generic elimination' if (args.generic_ext or not split) else 'direct assembly (free_extension_direct)'
            wr = space_record(W)
            wr['expected_profile'] = FX.profile_recursion(r['profile'])
            wr['profile_formula_ok'] = wr['expected_profile'] == wr['profile']
            wr['sweep_full'] = W.C.sweep_full()
            wr['sweep_random'] = W.C.sweep_random(t=4)
            g1 = V.chain_generator(p, a, b, X, k + 1)
            wr['generator_rows'] = int(FX.generator_matrix(g1, W.C.L, D, k + 1, p, 'all').shape[0])
            wr['generator_rows_outside'] = W.residual_rows(FX.generator_matrix(g1, W.C.L, D, k + 1, p, 'all'))
            wr['construction'] = wr0
            wr['seconds'] = round(time.time() - t, 1)
            ref = reference_sha(args, k + 1)
            if ref is not None:
                wr['reference_V_next'] = ref
                wr['equals_reference_V_next'] = ref['basis_sha256'] == wr['basis_sha256']
            e = dict(k0=k, W=wr, one_notin_S=not r['contains_one'], g_in_S=r['generator_in_space'],
                     W_closed=bool(wr['sweep_full'].get('passed')))
            e['kernel_hypotheses_hold'] = bool(e['W_closed'] and e['one_notin_S'] and e['g_in_S'])
            rec['extensions'][str(k)] = e
            log(stage='extension', k0=k, closed=e['W_closed'], rank=wr['rank'], profile=wr['profile'],
                formula=wr['profile_formula_ok'], gen_outside=wr['generator_rows_outside'],
                hyp=e['kernel_hypotheses_hold'], sec=wr['seconds'])
            prev = (k, None if args.drop_w else W)
            if args.drop_w:
                del W
        if args.witness:
            t = time.time()
            wres = FX.shifted_generator_witnesses(C, g, D, limit=args.witness_limit, want=1)
            wres['seconds'] = round(time.time() - t, 1)
            wres['base_depth'] = k
            rec.setdefault('witness', {})[str(k)] = wres
            log(stage='witness', **{x: y for x, y in wres.items()})
        del C, sp
    rec['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    rec['seconds'] = round(time.time() - T0, 1)
    rec['command'] = ' '.join(sys.argv)
    rec['source_sha256'] = {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in SOURCES}
    rec['label'] = 'exact'
    rec['date'] = time.strftime('%Y-%m-%d')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'records').mkdir(exist_ok=True)
    path = OUT / 'records' / f'{args.name}.json'
    path.write_text(json.dumps(rec, indent=1) + '\n')
    log(stage='saved', file=str(path))


def exceptional(args):
    p, a, b = args.p, args.a % args.p, args.b % args.p
    out = []
    for X in range(p):
        ok = dict(order2=(X ** 3 + a * X + b) % p == 0, tangent=(3 * X * X + a) % p == 0, zero=X == 0)[args.kind]
        if ok:
            rhs = (X ** 3 + a * X + b) % p
            rational = rhs == 0 or pow(rhs, (p - 1) // 2, p) == 1
            out.append(dict(X=X, rational_point=rational))
    print(json.dumps(dict(p=p, a=a, b=b, kind=args.kind, abscissas=out)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['chain', 'exceptional'])
    ap.add_argument('--p', type=int, required=True)
    ap.add_argument('--a', type=int, default=1)
    ap.add_argument('--b', type=int, default=0)
    ap.add_argument('--seed', type=int, default=None)
    ap.add_argument('--want', action='store_true')
    ap.add_argument('--X', type=int, default=None)
    ap.add_argument('--D', type=int, default=5)
    ap.add_argument('--kmin', type=int, default=None)
    ap.add_argument('--kmax', type=int, default=None)
    ap.add_argument('--draw-k', type=int, default=None)
    ap.add_argument('--ext-top', action='store_true')
    ap.add_argument('--witness', action='store_true')
    ap.add_argument('--witness-limit', type=int, default=20000)
    ap.add_argument('--generic-ext', action='store_true')
    ap.add_argument('--drop-w', action='store_true')
    ap.add_argument('--ref-name', default=None, help='e.g. k{k}-p12289-s10-D5.json')
    ap.add_argument('--nosplit', action='store_true')
    ap.add_argument('--kind', default='order2')
    ap.add_argument('--name', default=None)
    args = ap.parse_args()
    if args.cmd == 'chain':
        assert args.kmin is not None and args.kmax is not None and args.name
        assert args.X is not None or args.seed is not None
        chain(args)
    else:
        exceptional(args)


if __name__ == '__main__':
    main()
