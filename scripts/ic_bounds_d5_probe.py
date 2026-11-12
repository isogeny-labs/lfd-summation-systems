"""Degree-D chain closure probe at any depth k: resumable, with a full closure verification pass.

Generalizes scripts/run_ic_chain_depth10_probe.py (fixed k = 10, D = 5) to --k and --degree.
Two summands, subgroup chains of depth k over F_p (2^k | p - 1), curve y^2 = x^3 + a x + b, target
drawn from random.Random(seed) exactly as in the depth-10 probe.  V_D (modulo the chains) is built
by one of two engines:

  --engine shift   the depth-10 probe's algorithm unchanged (csd/ic_chain_shift_closure.py rows,
                   fp_echelon rounds over the full int64 basis, every low pivot processed once,
                   then the 2k shifts of every final low row re-reduced against a reduced copy);
  --engine split   csd/ic_bounds_d5.py: summand-swap symmetry split, RREF stored on free columns
                   only (uint16), exact tiled float64 products, k block-1 shifts; a verification
                   sweep re-reduces all k shifts of all final low rows of both parts.

Both certify the same statements:
  * reached:          rank V_D = |M_<=D| - |Z| and M_<D evaluates onto Z, so V_D = I cap M_<=D and
                      LFD <= D (exact LFD = D when dim V_{D-1} < dim(I cap M_<=D-1));
  * deficient-closed: V_D verified closed (so it is exactly V_D) with rank < |M_<=D| - |Z|, so an
                      element of I of degree <= D is not in V_D and LFD >= D + 1 (unit ideal: 1 not
                      in V_D).
  * closed-full:      (only with --no-early-stop) V_D verified closed and equal to I cap M_<=D.
A paused run (--max-minutes reached) exits with status 3 after a checkpoint; rerun the same command.

V_{D-1} is computed exactly for the record and the strict fall at D: by the dense engine with its
exclusion certificate (--prev dense, the depth-10 probe's choice; k <= 10) or by the split engine run
to verified closure (--prev split, default for k >= 11).

The output experiments/ic-bounds-d5/runs/<tag>.json records the instance, the dimensions, the
filtration profile dim(V_D cap M_<=w), w = 0..D, timings, peak memory, and source hashes.
"""
import argparse, hashlib, json, math, os, resource, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from csd import ic_bounds_d5 as B
from csd import ic_chain_shift_closure as SC
from csd import ic_fall_certificates as C
from csd import ic_chain_block_closure as BC
from csd import elim

SOURCES = ['scripts/ic_bounds_d5_probe.py', 'csd/ic_bounds_d5.py', 'csd/ic_chain_shift_closure.py',
           'csd/ic_chain_block_closure.py', 'csd/ic_fall_certificates.py', 'csd/ic_curve_sweep.py',
           'csd/ic_dreg.py', 'csd/elim.py', 'csd/_elim.c']

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--k', type=int, required=True)
ap.add_argument('--degree', type=int, default=5)
ap.add_argument('--prime', type=int, default=12289)
ap.add_argument('--a', type=int, default=1)
ap.add_argument('--b', type=int, default=0)
ap.add_argument('--seed', type=int, default=10)
ap.add_argument('--want-solutions', action='store_true',
                help='redraw targets until the system has a scalar solution (a proper ideal)')
ap.add_argument('--engine', choices=['split', 'shift'], default='split')
ap.add_argument('--prev', choices=['dense', 'split', 'none'], default=None,
                help='how to compute V_{D-1} exactly (default: dense for k <= 10, split otherwise)')
ap.add_argument('--no-early-stop', action='store_true',
                help='split engine: continue to a verified closure even after the rank target is reached')
ap.add_argument('--batch', type=int, default=600, help='split engine: low rows processed per round')
ap.add_argument('--product-cap', type=int, default=6000, help='shift engine: rows of products per round')
ap.add_argument('--max-minutes', type=float, default=None, help='pause (checkpoint, exit 3) after this long')
ap.add_argument('--save-every', type=float, default=0.0,
                help='split engine: seconds between checkpoints (0 = after every round; pauses always save)')
ap.add_argument('--tag', default=None, help='output name (default k<k>-p<prime>-s<seed>-D<degree>-<engine>)')
ap.add_argument('--out-dir', type=Path, default=ROOT / 'experiments/ic-bounds-d5/runs')
ap.add_argument('--ckpt', type=Path, default=None, help='checkpoint directory (default ~/.cache/ic-bounds-d5/<tag>)')
args = ap.parse_args()
k, D, p = args.k, args.degree, args.prime
tag = args.tag or f'k{k}-p{p}-s{args.seed}-D{D}-{args.engine}' + ('-proper' if args.want_solutions else '')
CKPT = args.ckpt or Path.home() / '.cache' / 'ic-bounds-d5' / tag
OUT = args.out_dir / f'{tag}.json'
prev_mode = args.prev or ('dense' if k <= 10 else 'split')
started = time.time()
deadline = started + 60 * args.max_minutes if args.max_minutes else None


def log(**kw):
    kw['elapsed'] = round(time.time() - started)
    kw['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(json.dumps(kw), flush=True)


inst = B.make_instance(p, args.a, args.b, args.seed, k, args.want_solutions)
E, g, Z, g0 = inst['E'], inst['generator'], inst['points'], inst['g0']
log(stage='instance', k=k, degree=D, prime=p, a=args.a, b=args.b, seed=args.seed, target=inst['target'],
    points=len(Z), attempts=inst['attempts'], generator_terms=len(g))
count = lambda d: B.column_count(k, d)
target = count(D) - len(Z)
onto = True
if Z:
    masks = B.monomials(k, D)
    onto = B.evaluation_rank_fast(Z, [m for m in masks if m.bit_count() < D], k, g0, p) == len(Z)
log(stage='target', columns=count(D), target=target, evaluation_onto_below_D=onto)

# ---------------------------------------------------------------- exact V_{D-1}
prev = None
if prev_mode == 'dense':
    Vp = C.chain_degree_closure([g], 2, k, p, D - 1, max_columns=60000, max_cells=1_500_000_000)
    assert Vp.contains_one or C.check_closed_exclusion_certificate(Vp, [g])
    prev = dict(engine='dense+exclusion-certificate', dim=Vp.quotient_rank, contains_one=Vp.contains_one)
elif prev_mode == 'split':
    r = B.split_closure(g, k, p, D - 1, count(D - 1) - len(Z), early_stop=False, batch=args.batch)
    assert r['status'] == 'closed' and r['echelon_ok']
    prev = dict(engine='split-verified-closed', dim=r['rank'], contains_one=r['contains_one'], profile=r['profile'])
if prev:
    prev['ideal_dim_lower'] = count(D - 1) - len(Z)
    log(stage='prev', degree=D - 1, **prev)


def finish(record):
    args.out_dir.mkdir(parents=True, exist_ok=True)
    record.update(prime=p, a=args.a, b=args.b, seed=args.seed, target=inst['target'], k=k, degree=D,
                  points=len(Z), want_solutions=args.want_solutions, engine=args.engine, prev=prev,
                  columns=count(D), evaluation_onto_below_D=onto,
                  maxrss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024,
                  seconds_this_run=round(time.time() - started, 1),
                  source_sha256={n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in SOURCES})
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    log(stage='done', out=str(OUT), status=record['result']['status'])


def classify(rank, closed, contains_one, low_dim, reached):
    if reached:
        fall = prev is not None and prev['dim'] < count(D - 1) - len(Z)
        if not onto:
            return dict(status='reached-not-onto', rank=rank, target=target,
                        note='rank target reached but M_<D does not evaluate onto Z: no certificate')
        return dict(status='reached', lfd_upper_bound=D, exact_lfd=(D if fall else None) if prev else None,
                    rank=rank, target=target, contains_one=contains_one,
                    note=f'V_{D} = I cap M_<={D} (rank equals target, evaluation onto): LFD <= {D}'
                         + ((f'; strict fall at {D} since dim V_{D-1} < dim(I cap M_<={D-1}), so LFD = {D}' if fall
                             else f'; no fall at {D}') if prev else ''))
    if closed:
        return dict(status='deficient-closed', lfd_lower_bound=D + 1, rank=rank, target=target,
                    deficiency=target - rank, low_dimension=low_dim, contains_one=contains_one,
                    fall_at_D=(low_dim > prev['dim']) if prev else None,
                    note=f'V_{D} is exactly the degree-{D} closure and has dimension below dim(I cap M_<={D}), '
                         f'so an element of I of degree <= {D} is not in V_{D}: LFD >= {D + 1}'
                         + ('; 1 is not in V_D' if not contains_one else ''))
    return dict(status='not-closed', rank=rank)


# ---------------------------------------------------------------- split engine
if args.engine == 'split':
    L = B.Layout(k, D)
    log(stage='layout', columns=L.cols, part_columns=list(L.part_cols), low_columns=len(L.low))
    sc = B.SplitClosure(L, g, p, target, ckpt=CKPT, log=log, save_every=args.save_every)
    if sc.load():
        log(stage='resumed', phase=sc.state['phase'], rank=sc.rank, rounds=sc.state['rounds'])
    status = sc.run(batch=args.batch, deadline=deadline, early_stop=not args.no_early_stop)
    if status == 'paused':
        log(stage='paused', rank=sc.rank, phase=sc.state['phase'], checkpoint=str(CKPT))
        sys.exit(3)
    assert all(P.check_echelon() for P in sc.parts)
    rank = sc.rank
    res = classify(rank, status == 'closed', sc.contains_one(), sc.low_dimension(),
                   reached=(rank >= target and (status == 'reached' or args.no_early_stop)))
    if status == 'closed' and rank >= target:
        res['status'] = 'closed-full' if onto else res['status']
    h = hashlib.sha256()
    for P in sc.parts:
        h.update(P.piv.astype(np.int64).tobytes()); h.update(np.ascontiguousarray(P.B).tobytes())
    finish(dict(result=res, ranks=[sc.parts[0].rank, sc.parts[1].rank], part_columns=list(L.part_cols),
                profile=sc.profile(), sym_profile=sc.part_profile(0), anti_profile=sc.part_profile(1),
                low_dimension=sc.low_dimension(), rounds=sc.state['rounds'], sweeps=sc.state['sweeps'],
                rank_history=sc.state['history'], batch=args.batch,
                total_seconds=round(sc.state['seconds'] + time.time() - sc.started, 1),
                free_storage_entries=int(sum(P.B.size for P in sc.parts)),
                peak_free_entries=sc.state.get('peak_free_entries'), work=sc.state.get('work'),
                basis_sha256=h.hexdigest()))
    sys.exit(0)

# ---------------------------------------------------------------- shift engine (depth-10 probe algorithm)
CKPT.mkdir(parents=True, exist_ok=True)
PRODUCT_CAP, MAX_CELLS = args.product_cap, 1_300_000_000


def save_state(state):
    (CKPT / 'state.tmp.json').write_text(json.dumps(state))
    os.replace(CKPT / 'state.tmp.json', CKPT / 'state.json')


def save(basis, state):
    tmp = CKPT / 'basis.tmp.npy'
    np.save(tmp, basis.astype(np.uint16) if p < 65536 else basis.astype(np.uint32))
    os.replace(tmp, CKPT / 'basis.npy')
    save_state(state)


masks = SC.monomials(k, D)
index = {m: i for i, m in enumerate(masks)}
perms = [SC.shift_permutation(k, masks, index, v) for v in range(2 * k)]
weights = np.asarray([m.bit_count() for m in masks])
budget = max(1, PRODUCT_CAP // (2 * k))
state_path = CKPT / 'state.json'
if state_path.exists():
    state = json.loads(state_path.read_text())
    basis = np.load(CKPT / 'basis.npy').astype(np.int64)
    log(stage='resumed', phase=state['phase'], rounds=state['rounds'], rows=int(basis.shape[0]))
else:
    basis = SC.generator_rows(g, k, D, masks, index, p)
    state = dict(phase='closure', rounds=0, processed=[], history=[], verified_chunks=0)

if state['phase'] == 'closure':
    processed = set(state['processed'])
    while True:
        state['rounds'] += 1
        if basis.size > MAX_CELLS:
            raise ValueError('Shift-closure matrix budget exceeded.')
        out = []
        rank, pivots, _ = elim.fp_echelon(basis, p, out)
        basis = out[0]
        state['history'].append(rank)
        pending = [i for i, col in enumerate(pivots) if weights[col] < D and col not in processed]
        log(stage='round', round=state['rounds'], rank=rank, target=target, processed=len(processed), pending=len(pending))
        if rank >= target or not pending:
            state.update(phase='verify' if rank < target else 'reached', rank=rank, processed=sorted(processed))
            save(basis, state)
            break
        chosen = pending[:budget]
        state['processed'] = sorted(processed)
        save(basis, state)
        if deadline and time.time() > deadline:
            state['rounds'] -= 1
            state['history'].pop()
            save(basis, state)
            log(stage='paused', rank=rank, checkpoint=str(CKPT))
            sys.exit(3)
        processed.update(pivots[i] for i in chosen)
        basis = np.vstack([basis] + [SC.shift_rows(basis[chosen], perm) for perm in perms])

out = []
rank, pivots, _ = elim.fp_echelon(basis, p, out)
basis = out[0]
assert rank == state['rank']
low = [i for i, col in enumerate(pivots) if weights[col] < D]
one_col = len(masks) - 1
contains_one = one_col in pivots
profile = [int(sum(1 for c in pivots if weights[c] <= w)) for w in range(D + 1)]
if state['phase'] == 'reached':
    res = classify(rank, False, contains_one, len(low), reached=True)
else:
    chunks = list(range(0, len(low), budget))
    log(stage='verify-start', rank=rank, low_rows=len(low), chunks=len(chunks), done=state['verified_chunks'])
    closed = True
    low_rows = basis[low].copy()
    ref = BC.EchelonPart(basis.shape[1], p)
    for s0 in range(0, len(basis), 2000):
        ref.add(basis[s0:s0 + 2000])
    assert ref.rank == rank, (ref.rank, rank)
    del basis
    for ci, start in enumerate(chunks):
        if ci < state['verified_chunks']:
            continue
        if deadline and time.time() > deadline:
            log(stage='paused', phase='verify', checkpoint=str(CKPT))
            sys.exit(3)
        chunk_rows = low_rows[start:start + budget]
        missing = 0
        for perm in perms:
            shifted = SC.shift_rows(chunk_rows, perm)
            for s0 in range(0, len(shifted), 2000):
                missing += int(ref.reduce(shifted[s0:s0 + 2000]).any(axis=1).sum())
        if missing:
            closed = False
            log(stage='verify-chunk', chunk=ci, missing_rows=missing, closed=False)
            break
        state['verified_chunks'] = ci + 1
        save_state(state)
        log(stage='verify-chunk', chunk=ci, of=len(chunks), closed_so_far=True)
    res = classify(rank, closed, contains_one, len(low), reached=False)
finish(dict(result=res, rank_history=state['history'], rounds=state['rounds'], product_cap=PRODUCT_CAP,
            profile=profile, low_dimension=len(low),
            basis_sha256=hashlib.sha256((CKPT / 'basis.npy').read_bytes()).hexdigest()))
