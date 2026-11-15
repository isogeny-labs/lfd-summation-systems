"""Rebuild the solution sets of the binary m-summand systems whose records keep only |Z| and a points hash
(S3 rows of rows.jsonl with no computed GB degree), check the hash, and compute max.GB.deg (DRL, grlex)."""
import hashlib, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
from csd import ic_bounds_msum as S, ic_dreg as I, ic_solving_degree as SD
from csd.ic_research_campaign import field
import ic_bounds_msum_run as RUN
EXP = ROOT / 'experiments'
out = []
for name in sys.argv[1:]:
    rec_path = EXP / 'ic-bounds-msummands/records' / (name + '.json')
    if rec_path.exists():
        rec = json.load(open(rec_path))
        F = field(int(rec['n'])); E = I.BinaryCurve(F, int(rec['a']), int(rec['b']))
        basis = json.loads(rec['basis']) if isinstance(rec['basis'], str) else rec['basis']
        target = json.loads(rec['target']) if isinstance(rec['target'], str) else rec['target']
        _, truth = S.binary_generator(E, basis, int(rec['m']), target[0])
        pts = S.binary_zero_set(truth); sha = rec['points_sha256']; lfd = int(rec['exact_lfd']); m, k = int(rec['m']), int(rec['k'])
    else:
        spec = json.load(open(EXP / 'ic-bounds-msummands/gpu/specs' / (name + '.json')))['spec']
        c = RUN.build(spec); pts = c['pts']
        rec = json.load(open(next((EXP / 'ic-bounds-msummands/pod-records').glob('*/out/%s.json' % name))))
        sha = rec['points_sha256']; lfd = int(rec['exact_lfd']); m, k = spec['m'], spec['k']
    ok = hashlib.sha256(json.dumps(pts).encode()).hexdigest() == sha
    r = m * k
    row = dict(id=name, m=m, k=k, r=r, lfd=lfd, points=len(pts), points_sha256_ok=ok)
    if ok:
        for order in ('drl', 'grlex'):
            st = SD.staircase(SD.BooleanPoints(pts, r), len(pts), r, 2, order)
            row['gb_' + order] = st['max_gb_degree']; row['hf_' + order] = st['standard_by_degree']
    print(json.dumps(row), flush=True)
    out.append(row)
