#!/usr/bin/env python3
"""Control for index_check.py's modular HNF: index of the smaller ring Z[pi] (pi generators only),
compared with PARI on the same lattice (control_zpi.gp). A wrong HNF would disagree here too."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from index_check import gens, hnf
for ell in [5, 11, 13, 19, 31]:
    n = 2 * (ell - 1); G, _ = gens(ell)
    piv = hnf(G[:n] + G[:n], n)   # Z[pi] = span of pi^0..pi^(2d-1)
    idx = 1
    for c in piv: idx *= piv[c][c]
    print(ell, idx)
