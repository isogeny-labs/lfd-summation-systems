"""Exact LFD of the depth-10 unit-ideal instance: decide whether 1 lies in the degree-6 closure.

Input: the verified degree-5 closure basis written by scripts/run_ic_chain_depth10_probe.py
(checkpoint directory --v5, files basis.npy and state.json with phase 'verify' and every
verification chunk passed).  The instance is the default probe instance (p = 12289,
y^2 = x^3 + x, target abscissa 9362, k = 10), a unit ideal with LFD >= 6 (V_5 misses 1).

Method: csd/ic_chain_block_closure.py (symmetry-split blocked upper certificate).  Every row added
lies in V_6, so a full-rank symmetric part proves 1 in V_6 and, with LFD >= 6, exact LFD = 6.
If the rounds stop below full rank the result is reported as undecided (no lower bound is claimed).

Checkpoints both parts after every round to --ckpt; a restart resumes.  Needs about 10 GB of
memory and several hours; run it alone.  Writes experiments/ic-chain-depth10/degree6.json.
"""
import argparse, hashlib, json, os, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from csd import ic_chain_block_closure as B
from csd.ic_chain_shift_closure import monomials

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--v5', type=Path, default=Path.home() / '.cache' / 'ic-chain-depth10' / 'rerun')
ap.add_argument('--ckpt', type=Path, default=Path.home() / '.cache' / 'ic-chain-depth10' / 'degree6')
ap.add_argument('--combos', type=int, default=6000)
ap.add_argument('--seed', type=int, default=20261105)
args = ap.parse_args()
p, k, D = 12289, 10, 6
started = time.time()


def log(**kw):
    kw['elapsed'] = round(time.time() - started)
    print(json.dumps(kw), flush=True)


v5_state = json.loads((args.v5 / 'state.json').read_text())
assert v5_state['phase'] == 'verify' and v5_state['rank'] == 18569, v5_state.get('phase')
v5 = np.load(args.v5 / 'basis.npy')  # uint16; embedded in batches by the engine
assert v5.shape == (18569, 21700)
v5_hash = hashlib.sha256((args.v5 / 'basis.npy').read_bytes()).hexdigest()
log(stage='input', v5_rank=len(v5), v5_sha256=v5_hash)

args.ckpt.mkdir(parents=True, exist_ok=True)
state = None
if (args.ckpt / 'state.json').exists():
    meta = json.loads((args.ckpt / 'state.json').read_text())
    try:
        sym = B.EchelonPart(meta['sym_cols'], p)
        anti = B.EchelonPart(meta['anti_cols'], p)
        sym.basis = np.load(args.ckpt / 'sym_basis.npy'); sym.pivots = np.load(args.ckpt / 'sym_pivots.npy')
        anti.basis = np.load(args.ckpt / 'anti_basis.npy'); anti.pivots = np.load(args.ckpt / 'anti_pivots.npy')
        assert sym.basis.shape == (len(sym.pivots), sym.cols) and anti.basis.shape == (len(anti.pivots), anti.cols)
        state = dict(round=meta['round'], sym=sym, anti=anti, embedded=meta.get('embedded', len(v5)))
        log(stage='resumed', round=meta['round'], embedded=state['embedded'], sym_rank=sym.rank, anti_rank=anti.rank)
    except Exception as error:  # a checkpoint torn by a restart: start over
        log(stage='checkpoint-unusable', error=str(error))
        state = None


def checkpoint(st):
    # every file is written to a temporary name and renamed, so a kill leaves either the old or the
    # new file; the shape assertion on resume rejects a mixed pair and the run starts over
    for name in ('sym', 'anti'):
        for part, arr in (('basis', st[name].basis), ('pivots', st[name].pivots)):
            tmp = args.ckpt / f'{name}_{part}.tmp.npy'
            np.save(tmp, arr)
            os.replace(tmp, args.ckpt / f'{name}_{part}.npy')
    meta = dict(round=st['round'], sym_cols=st['sym'].cols, anti_cols=st['anti'].cols,
                sym_rank=st['sym'].rank, anti_rank=st['anti'].rank, embedded=st.get('embedded', len(v5)))
    (args.ckpt / 'state.tmp.json').write_text(json.dumps(meta))
    os.replace(args.ckpt / 'state.tmp.json', args.ckpt / 'state.json')


# Resumable embedding of V_5 into the symmetric and antisymmetric parts at degree 6, checkpointed
# after every batch (the engine skips its own embedding when it is handed a state).
masks6 = monomials(k, D)
index6 = {m: i for i, m in enumerate(masks6)}
split = B.SymmetricSplit(k, masks6, index6)
pidx = np.asarray([index6[m] for m in monomials(k, 5)])
if state is None:
    state = dict(round=0, sym=B.EchelonPart(len(split.sym_reps), p), anti=B.EchelonPart(len(split.anti_reps), p),
                 embedded=0)
EMBED = 1000
while state['embedded'] < len(v5):
    s0 = state['embedded']
    prev = np.zeros((len(v5[s0:s0 + EMBED]), len(masks6)), dtype=np.int64)
    prev[:, pidx] = v5[s0:s0 + EMBED]
    state['sym'].add(split.sym(prev, p), chunk=1000)
    state['anti'].add(split.anti(prev, p), chunk=1000)
    del prev
    state['embedded'] = s0 + EMBED if s0 + EMBED < len(v5) else len(v5)
    checkpoint(state)
    log(stage='embed', embedded=state['embedded'], of=len(v5), sym_rank=state['sym'].rank, anti_rank=state['anti'].rank)


history = []
while True:
    res = B.block_closure_contains_one(k, p, D, v5, monomials(k, 5), seed=args.seed + (state['round'] if state else 0),
                                       combos_per_round=args.combos, max_rounds=(state['round'] if state else 0) + 1,
                                       log=lambda d: (history.append(d), log(stage='round', **d)), state=state)
    state = res['state']
    checkpoint(state)
    if res['contains_one'] or (history and history[-1]['added'] == [0, 0]) or state['round'] >= 40:
        break

result = dict(contains_one=res['contains_one'], sym_rank=res['sym_rank'], sym_cols=res['sym_cols'],
              anti_rank=res['anti_rank'], anti_cols=res['anti_cols'], rounds=res['rounds'])
if res['contains_one']:
    result.update(status='exact', exact_lfd=6,
                  note='every row lies in V_6; the symmetric part reached full rank, so 1 is in V_6; '
                       'with the verified LFD >= 6 (1 not in V_5) the exact LFD is 6')
else:
    result.update(status='undecided', note='the rounds stopped below full rank; no lower bound is claimed')
record = dict(prime=p, a=1, b=0, target_x=9362, k=k, degree=D, v5_sha256=v5_hash, combos_per_round=args.combos,
              seed=args.seed, round_history=history, result=result, seconds=round(time.time() - started))
out = ROOT / 'experiments/ic-chain-depth10/degree6.json'
out.write_text(json.dumps(record, indent=1) + '\n')
log(stage='done', **result)
