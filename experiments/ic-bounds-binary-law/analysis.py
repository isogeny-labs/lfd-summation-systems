#!/usr/bin/env python3
"""Candidate laws for the two-summand binary/Koblitz LFD, with counts.

Input: table.jsonl (build_table.py) and structure.jsonl (structure.py).
Output: laws.json and laws.md (counts only; no summary statistics).

Unit of counting: distinct systems (system_hash).  Proper-subspace systems
(k < n) from ic-curve-sweep and ic-bounds-binary; the full-field inverse family
(k = n) is reported separately.

System classes (from structure.jsonl, exact):
  const      1 lies in the span of the input equations (unit, LFD <= 2 trivially)
  order-two  target x_R = 0: x1 x2 = sqrt(b), homogeneous bilinear + constant
  affine     x_R != 0: affine bilinear, with the trace equation Tr(x1)+Tr(x2)=Tr(b/x_R^2)
             in the span unless V lies in the trace-zero hyperplane
  e = r_eff - m_eff  (log2 of the expected number of solutions)
"""
import json, math, os, collections

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    T = [json.loads(l) for l in open(os.path.join(HERE, 'table.jsonl'))]
    St = {}
    for l in open(os.path.join(HERE, 'structure.jsonl')):
        x = json.loads(l); St[x['system_hash']] = x
    S, inv = {}, []
    for t in T:
        if t['source'] == 'ic-binary-inverse':
            inv.append(t); continue
        if t['system_hash'] not in S:
            s = dict(t, **St[t['system_hash']])
            s['cls'] = 'const' if s['const'] else ('order-two' if s['target_x'] == 0 else 'affine')
            s['e'] = s['r_eff'] - s['m_eff']
            s['curves'] = set()
            s['bases'] = set()
            S[t['system_hash']] = s
        S[t['system_hash']]['curves'].add(t['curve_class'])
    return list(S.values()), inv


def coeffs(r, m, N):
    return [sum(math.comb(r, i - 2 * j) * (-1) ** j * math.comb(m + j - 1, j) for j in range(i // 2 + 1))
            for i in range(N + 1)]


def semi_regular_dreg(r, m):
    for i, c in enumerate(coeffs(r, m, r + 2)):
        if c <= 0:
            return i
    return r + 1


def crossover(r, m, Z):
    """Least D with sum_{i<=D} c_i(r, m) <= |Z| (signed Boolean semi-regular series)."""
    c = coeffs(r, m, r + 2)
    acc = 0
    for D, x in enumerate(c):
        acc += x
        if acc <= Z:
            return D
    return r + 1


def lookup_law(systems, key):
    """Fit LFD = mode over each cell of `key`; count satisfied/violated and cell sizes."""
    cells = collections.defaultdict(collections.Counter)
    for s in systems:
        cells[key(s)][s['lfd']] += 1
    sat = sum(c.most_common(1)[0][1] for c in cells.values())
    single = sum(1 for c in cells.values() if sum(c.values()) == 1)
    pure = sum(1 for c in cells.values() if len(c) == 1)
    return dict(systems=len(systems), cells=len(cells), singleton_cells=single, pure_cells=pure,
                satisfy=sat, violate=len(systems) - sat)


def closed_law(systems, pred):
    c = collections.Counter()
    for s in systems:
        p = pred(s)
        c['equal' if s['lfd'] == p else ('lfd_below' if s['lfd'] < p else 'lfd_above')] += 1
    return dict(systems=len(systems), satisfy=c['equal'], violate=c['lfd_below'] + c['lfd_above'],
                lfd_below_prediction=c['lfd_below'], lfd_above_prediction=c['lfd_above'])


def bound_law(systems, bound):
    v = [s for s in systems if s['lfd'] > bound(s)]
    return dict(systems=len(systems), satisfy=len(systems) - len(v), violate=len(v),
                examples=[dict(n=s['n'], k=s['k'], Z=s['Z'], lfd=s['lfd'], cls=s['cls'], base=s['base_kind'])
                          for s in v[:5]])


def B_threshold(x):
    """Frozen k+e law for the affine class (count-optimal thresholds on the existing data)."""
    return 2 + (x >= 3) + (x >= 9) + (x >= 19)


def A_threshold(e):
    """Frozen e-only law for the affine class (mode of LFD per e on the existing data)."""
    return 2 if e <= -1 else 3 if e <= 3 else 4 if e <= 8 else 5


def main():
    S, inv = load()
    nonconst = [s for s in S if s['cls'] != 'const']
    affine = [s for s in S if s['cls'] == 'affine']
    two = [s for s in S if s['cls'] == 'order-two']
    proper = [s for s in nonconst if s['Z'] > 0]
    out = dict(counts=dict(systems=len(S), const=len(S) - len(nonconst), affine=len(affine), order_two=len(two),
                           proper_ideal=len(proper),
                           unit_nonconst=len([s for s in nonconst if s['Z'] == 0]),
                           all_bilinear_symmetric=all(s['bilinear'] for s in S)))
    # structure facts (exact on every system)
    out['structure'] = dict(
        affine_m_eff_eq_n_minus_1=collections.Counter((s['m_eff'] == s['n'] - 1) for s in affine),
        affine_r_eff=collections.Counter(
            ('2k-1' if s['r_eff'] == 2 * s['k'] - 1 else '2k' if s['r_eff'] == 2 * s['k'] else 'other')
            for s in affine),
        affine_r_eff_2k_bases=collections.Counter(s['base_kind'] for s in affine if s['r_eff'] == 2 * s['k']),
        order_two_r_eff=collections.Counter(
            ('2k' if s['r_eff'] == 2 * s['k'] else 'other') for s in two),
        order_two_m_eff_eq_n=collections.Counter((s['m_eff'] == s['n']) for s in two))
    out['structure'] = {k: {str(a): b for a, b in v.items()} for k, v in out['structure'].items()}
    # lookup laws
    L = {}
    for name, key in [
        ('|Z|', lambda s: s['Z']),
        ('(k,|Z|)', lambda s: (s['k'], s['Z'])),
        ('n-k', lambda s: s['n'] - s['k']),
        ('2k-n', lambda s: 2 * s['k'] - s['n']),
        ('(n,k)', lambda s: (s['n'], s['k'])),
        ('(n,k,class)', lambda s: (s['n'], s['k'], s['cls'])),
        ('(n,k,class,base_kind,curve)', lambda s: (s['n'], s['k'], s['cls'], s['base_kind'], tuple(sorted(s['curves'])))),
        ('e,class', lambda s: (s['e'], s['cls'])),
        ('k+e,class', lambda s: (s['k'] + s['e'], s['cls'])),
        ('(k,e,class)', lambda s: (s['k'], s['e'], s['cls'])),
    ]:
        L[name] = dict(all_nonconst=lookup_law(nonconst, key), affine=lookup_law(affine, key),
                       order_two=lookup_law(two, key), proper_ideals=lookup_law(proper, key))
    out['lookup_laws'] = L
    # closed-form laws
    C = {}
    C['B: affine, LFD = 2+[k+e>=3]+[k+e>=9]+[k+e>=19]'] = closed_law(affine, lambda s: B_threshold(s['k'] + s['e']))
    C['A: affine, LFD = 2/3/4/5 for e<=-1 / <=3 / <=8 / >=9'] = closed_law(affine, lambda s: A_threshold(s['e']))
    C['B on order-two (same thresholds)'] = closed_law(two, lambda s: B_threshold(s['k'] + s['e']))
    C['semi-regular crossover min{D: sum_{i<=D} c_i(r_eff,m_eff) <= |Z|}'] = closed_law(
        nonconst, lambda s: crossover(s['r_eff'], s['m_eff'], s['Z']))
    C['semi-regular d_reg(r_eff, m_eff)'] = closed_law(nonconst, lambda s: semi_regular_dreg(s['r_eff'], s['m_eff']))
    C['k=10 grid rule: LFD 5 iff |Z|>470 else <=4 (proper, k>=8; equality to 5 / <=4 checked)'] = (lambda ss: dict(
        systems=len(ss),
        satisfy=sum((s['lfd'] == 5) == (s['Z'] > 470) for s in ss),
        violate=sum((s['lfd'] == 5) != (s['Z'] > 470) for s in ss)))(
        [s for s in proper if s['k'] >= 8])
    out['closed_laws'] = C
    # bounds
    out['bounds'] = {
        'LFD <= min(k,l)+2 = k+2': bound_law(nonconst, lambda s: s['k'] + 2),
        'LFD <= semi-regular d_reg(r_eff, m_eff)': bound_law(nonconst, lambda s: semi_regular_dreg(s['r_eff'], s['m_eff'])),
        'LFD <= max(2, floor(log2|Z|)+1) (proper)': bound_law(proper, lambda s: max(2, int(math.log2(s['Z'])) + 1)),
        'LFD <= max(2, floor(log2|Z|)+2) (proper)': bound_law(proper, lambda s: max(2, int(math.log2(s['Z'])) + 2)),
        'LFD <= 4 when |Z| <= 470 (proper)': bound_law([s for s in proper if s['Z'] <= 470], lambda s: 4),
        'LFD <= 4 when k+e <= 18 (nonconst)': bound_law([s for s in nonconst if s['k'] + s['e'] <= 18], lambda s: 4),
        'LFD <= 3 when k+e <= 7 (affine)': bound_law([s for s in affine if s['k'] + s['e'] <= 7], lambda s: 3),
        'certificate degree D* - LFD': dict(collections.Counter(s['cert_degree'] - s['lfd'] for s in proper)),
    }
    # pre-certificate closure profiles: same within (n, k, r_eff, m_eff, class)?
    groups = collections.defaultdict(list)
    for s in nonconst:
        groups[(s['n'], s['k'], s['r_eff'], s['m_eff'], s['cls'])].append(s)
    prof = collections.Counter(); prof_by = collections.defaultdict(collections.Counter)
    for g, ss in groups.items():
        # modal value of dim V_D among members not yet certified at D
        maxD = max(len(s['dims']) for s in ss)
        modal = []
        for D in range(maxD):
            c = collections.Counter(s['dims'][D] for s in ss if D < len(s['dims']) - 1)
            modal.append(c.most_common(1)[0][0] if c else None)
        for s in ss:
            pre = s['dims'][:-1]
            ok = all(pre[D] == modal[D] for D in range(len(pre)))
            prof['equal' if ok else 'differ'] += 1
            prof_by[(s['cls'], s['base_kind'])]['equal' if ok else 'differ'] += 1
    out['precertificate_profile_vs_group_mode'] = dict(
        all=dict(prof), by_class_base={'%s/%s' % k: dict(v) for k, v in prof_by.items()},
        groups=len(groups), singleton_groups=sum(1 for v in groups.values() if len(v) == 1))
    # does the same group mode hold across curves and bases?  (groups where >= 2 base kinds occur)
    mixed = [g for g, ss in groups.items() if len({s['base_kind'] for s in ss}) > 1]
    out['precertificate_profile_vs_group_mode']['groups_with_both_base_kinds'] = len(mixed)
    # inverse family
    out['inverse_family'] = [dict(n=t['n'], Z=t['Z'], lfd=t['lfd'], k_plus_e=t['n'] + t['n'],
                                  B_pred_order_two_thresholds=B_threshold(2 * t['n'])) for t in inv]
    json.dump(out, open(os.path.join(HERE, 'laws.json'), 'w'), indent=1, default=str)
    # markdown table
    lines = ['| law | systems | cells (singleton) | satisfy | violate |', '|---|---|---|---|---|']
    for name, d in L.items():
        for sub in ('all_nonconst', 'affine', 'order_two'):
            x = d[sub]
            lines.append('| `%s` lookup, %s | %d | %d (%d) | %d | %d |' % (
                name, sub, x['systems'], x['cells'], x['singleton_cells'], x['satisfy'], x['violate']))
    lines += ['', '| closed-form law | systems | satisfy | violate | LFD below / above prediction |',
              '|---|---|---|---|---|']
    for name, x in C.items():
        lines.append('| %s | %d | %d | %d | %s / %s |' % (name, x['systems'], x['satisfy'], x['violate'],
                                                       x.get('lfd_below_prediction', '-'), x.get('lfd_above_prediction', '-')))
    lines += ['', '| bound | systems | satisfy | violate |', '|---|---|---|---|']
    for name, x in out['bounds'].items():
        if 'satisfy' in x:
            lines.append('| %s | %d | %d | %d |' % (name, x['systems'], x['satisfy'], x['violate']))
    open(os.path.join(HERE, 'laws.md'), 'w').write('\n'.join(lines) + '\n')
    print(json.dumps(out['counts']))
    print(json.dumps(out['structure']))
    print('\n'.join(lines))
    print(json.dumps(out['bounds'], default=str)[:3000])
    print(json.dumps(out['precertificate_profile_vs_group_mode']))


if __name__ == '__main__':
    main()
