#!/usr/bin/env python3
"""Sensitivity control: the independent engine must NOT return 4 for every prime k=3 system.
Random generators of chosen degree g (dense over squarefree chain monomials of degree <= g), p=97, k=3."""
import os, sys, random, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lfd_engine as L
p, k = 97, 3
rng = random.Random(5)
masks_by_deg = {}
for m in range(1 << (2*k)):
    masks_by_deg.setdefault(bin(m).count('1'), []).append(m)
rows = []
for g in (1, 2, 3, 4):
    for t in range(4):
        gen = [(m, rng.randrange(1, p)) for d in range(g+1) for m in masks_by_deg[d] if rng.random() < 0.7 or d == g]
        inputs, G = L.prime_inputs(gen, k, p)
        pts = L.prime_points(G, k, p)
        res = L.exact_lfd(inputs, 2*k, p, pts, 8)
        rows.append((g, len(pts), res['lfd'], res['falls'], res['certified_at']))
        print('deg', g, 'points', len(pts), 'LFD', res['lfd'], 'falls', res['falls'], 'cert', res['certified_at'], flush=True)
