#!/usr/bin/env python3
"""Pack the code the m-summand pod runs need into one archive: csd/*.py, the three driver scripts, gpu/ (chain,
wrapper, cells, chains plan, cost with frozen hashes, specs).  usage: make_bundle.py OUT.tgz"""
import glob
import os
import sys
import tarfile
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
G = os.path.join(ROOT, 'experiments', 'ic-bounds-msummands', 'gpu')
with tarfile.open(sys.argv[1], 'w:gz') as tar:
    # the C sources are compiled on first use (csd/elim.py reads csd/_elim.c); the first pod run failed without them
    for f in sorted(glob.glob(os.path.join(ROOT, 'csd', '*.py')) + glob.glob(os.path.join(ROOT, 'csd', '*.c'))):
        tar.add(f, arcname='csd/' + os.path.basename(f))
    for f in ('ic_bounds_msum_run.py', 'ic_bounds_bin_big_run.py', 'ic_bounds_m3_depth.py'):
        tar.add(os.path.join(ROOT, 'scripts', f), arcname='scripts/' + f)
    for f in sorted(glob.glob(os.path.join(G, '*'))):
        if os.path.isfile(f):
            tar.add(f, arcname='gpu/' + os.path.basename(f))
    for f in sorted(glob.glob(os.path.join(G, 'specs', '*.json'))):
        tar.add(f, arcname='gpu/specs/' + os.path.basename(f))
print('bundle', sys.argv[1], os.path.getsize(sys.argv[1]))
