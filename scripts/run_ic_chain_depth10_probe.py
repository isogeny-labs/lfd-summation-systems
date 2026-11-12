"""Depth-10 prime chain probe: resumable degree-5 shift closure with full closure verification.

Two summands, subgroup chains of depth k = 10 over F_p (2^10 | p-1), curve y^2 = x^3 + a x + b,
target drawn from random.Random(seed).  The degree-5 closure V_5 (modulo the chains) is built by
the shift-structured engine of csd/ic_chain_shift_closure.py, every low pivot processed once, and
then verified closed: the 2k variable shifts of every echelon row of weight < 5 are re-reduced
against the final basis.  V_5 lies in I ∩ M_{<=5}; the chain ideal is radical with all points
rational, so I is the vanishing ideal of the scalar solution set Z and, once the squarefree
monomials of weight < 5 evaluate onto Z, dim(I ∩ M_{<=5}) = |M_{<=5}| - |Z| (the target).

Outcomes:
  * reached:   rank V_5 = target, so V_5 = I ∩ M_{<=5} and LFD <= 5 (exact value from the falls);
  * deficient: V_5 verified closed with rank < target, so some h in I of degree <= 5 is not in
               V_5 and LFD >= 6 (for the unit ideal, Z empty: 1 is not in V_5).

Defaults reproduce the first instance (p = 12289, a = 1, b = 0, seed 10, unit ideal) and write
experiments/ic-chain-depth10/probe.json.  The run checkpoints the echelon basis (uint16,
about 770 MB) and its progress after every closure round and every verification chunk to
--ckpt, so a restart resumes.  One run takes about two hours and 6-8 GB of memory.
"""
import argparse, sys, time, random, json, os, hashlib, math
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from csd import ic_dreg as I
from csd import ic_curve_sweep as S
from csd import ic_chain_shift_closure as SC
from csd import ic_fall_certificates as C
from csd import elim
from csd import ic_chain_block_closure as BC

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--prime', type=int, default=12289)
ap.add_argument('--a', type=int, default=1)
ap.add_argument('--b', type=int, default=0)
ap.add_argument('--seed', type=int, default=10)
ap.add_argument('--want-solutions', action='store_true',
                help='redraw targets until the system has a scalar solution (a proper ideal)')
ap.add_argument('--tag', default='probe', help='output name: experiments/ic-chain-depth10/<tag>.json')
ap.add_argument('--ckpt', type=Path, default=None, help='checkpoint directory (default ~/.cache/ic-chain-depth10/<tag>)')
args = ap.parse_args()
p, a, b, k, D = args.prime, args.a, args.b, 10, 5
assert (p - 1) % (1 << k) == 0 and (4 * a ** 3 + 27 * b * b) % p
CKPT = args.ckpt or Path.home() / '.cache' / 'ic-chain-depth10' / args.tag
CKPT.mkdir(parents=True, exist_ok=True)
OUT = ROOT / 'experiments/ic-chain-depth10' / f'{args.tag}.json'
PRODUCT_CAP, MAX_CELLS = 6000, 1_300_000_000
started = time.time()


def log(**kw):
    kw['elapsed'] = round(time.time() - started)
    print(json.dumps(kw), flush=True)


def save_state(state):
    (CKPT / 'state.tmp.json').write_text(json.dumps(state))
    os.replace(CKPT / 'state.tmp.json', CKPT / 'state.json')


def save(basis, state):
    tmp = CKPT / 'basis.tmp.npy'
    np.save(tmp, basis.astype(np.uint16) if p < 65536 else basis.astype(np.uint32))
    os.replace(tmp, CKPT / 'basis.npy')
    save_state(state)


rng = random.Random(args.seed)
E = I.PrimeCurve(p, a, b)
g0 = S.subgroup_generator_value(p, k)
for attempt in range(200):
    R = S.random_rational_target(E, rng, p)
    Z = S.prime_points(E, k, g0, 1, 1, R[0])
    if Z or not args.want_solutions:
        break
else:
    raise SystemExit('No target with scalar solutions found.')
g = I.pkm_generator(E, R[0], k)
log(stage='instance', prime=p, a=a, b=b, target=list(R), points=len(Z), attempts=attempt + 1)

masks = SC.monomials(k, D)
index = {m: i for i, m in enumerate(masks)}
perms = [SC.shift_permutation(k, masks, index, v) for v in range(2 * k)]
weights = np.asarray([m.bit_count() for m in masks])
target = len(masks) - len(Z)
budget = max(1, PRODUCT_CAP // (2 * k))
onto = True
if Z:
    onto = S.prime_evaluation_rank(Z, [m for m in masks if m.bit_count() < D], k, g0, p, 10 ** 9) == len(Z)
log(stage='target', target=target, evaluation_onto_below_5=onto)
assert onto, 'Evaluation of M_<5 is not onto Z; the ideal dimension is not |M_<=5| - |Z|.'

# exact dense V_4 for the record and for the fall at 5
V4 = C.chain_degree_closure([g], 2, k, p, 4, max_columns=16384, max_cells=400_000_000)
assert not V4.contains_one and C.check_closed_exclusion_certificate(V4, [g])
count4 = sum(math.comb(2 * k, d) for d in range(5))
log(stage='dense-V4', dim_V4=V4.quotient_rank, fit=4 * k * k + 172 * k - 571, ideal_le4=count4 - len(Z))

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
        state['rounds'] += 1  # on resume the interrupted round is recomputed from the saved echelon basis
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
        save(basis, state)  # saved before marking `chosen`: a resumed round re-derives the same choice
        processed.update(pivots[i] for i in chosen)
        basis = np.vstack([basis] + [SC.shift_rows(basis[chosen], perm) for perm in perms])

out = []
rank, pivots, _ = elim.fp_echelon(basis, p, out)
basis = out[0]
assert rank == state['rank']
low = [i for i, col in enumerate(pivots) if weights[col] < D]
one_col = len(masks) - 1  # masks sorted by weight descending: the constant monomial is last
if state['phase'] == 'reached':
    fall_at_5 = V4.quotient_rank < count4 - len(Z)
    result = dict(status='reached', lfd_upper_bound=5, exact_lfd=5 if fall_at_5 else None,
                  note='V_5 = I ∩ M_<=5 (rank equals target, evaluation onto): LFD <= 5; '
                       + ('strict fall at 5 since dim V_4 < dim(I ∩ M_<=4), so LFD = 5' if fall_at_5 else 'no fall at 5'),
                  rank=rank, target=target)
else:
    chunks = list(range(0, len(low), budget))
    log(stage='verify-start', rank=rank, low_rows=len(low), chunks=len(chunks), done=state['verified_chunks'])
    closed = True
    # Memory-lean verification: a reduced echelon copy of the basis (uint16) and the low rows are
    # kept; every shift of every low row is reduced against it in tiles.  The span is closed under
    # the shifts iff every reduced row is zero, which is the same test as the rank of the stacked
    # matrix equalling `rank`, without materialising the stack.
    low_rows = basis[low].copy()
    ref = BC.EchelonPart(basis.shape[1], p)
    for s0 in range(0, len(basis), 2000):
        ref.add(basis[s0:s0 + 2000])
    assert ref.rank == rank, (ref.rank, rank)
    del basis
    log(stage='verify-reference', rank=ref.rank)
    for ci, start in enumerate(chunks):
        if ci < state['verified_chunks']:
            continue
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
    if closed:
        contains_one = one_col in pivots
        result = dict(status='deficient-closed', lfd_lower_bound=6,
                      note=('V_5 is exactly the degree-5 closure and has dimension below dim(I ∩ M_<=5), '
                            'so an element of I of degree <= 5 is not in V_5: LFD >= 6'
                            + ('; 1 is not in V_5' if not contains_one else '')),
                      rank=rank, target=target, deficiency=target - rank, low_dimension=len(low),
                      contains_one=contains_one)
    else:
        result = dict(status='not-closed', note='verification found a missing product; continue the closure loop', rank=rank)

record = dict(prime=p, a=a, b=b, seed=args.seed, target=list(R), k=k, degree=D, points=len(Z),
              dim_V4=V4.quotient_rank, rounds=state['rounds'], rank_history=state['history'],
              product_cap=PRODUCT_CAP, result=result,
              basis_sha256=hashlib.sha256((CKPT / 'basis.npy').read_bytes()).hexdigest())
OUT.write_text(json.dumps(record, indent=1) + '\n')
log(stage='done', **result)
