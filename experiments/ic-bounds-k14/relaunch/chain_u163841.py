#!/usr/bin/env python3
"""Pod-side chain for job B of ../RELAUNCH.md only: the unit instance over p = 163841 at k = 14
(closure D = 6 with the basis kept, lift D = 6 -> 7, independent checker). Uses chain_k14.py, with one change:
the closure writes a single checkpoint, at the end of the closure phase, instead of one every 30 minutes
(two checkpoints near the peak size are about 240 GB; on 2026-11-11 that filled a 150 GB disk and the closure
exited with code 1 at about 4,100 s)."""
import chain_k14 as C

C.INSTANCES = [I for I in C.INSTANCES if I['tag'] == 'unit-p163841-k14']
assert len(C.INSTANCES) == 1
_run = C.run


def run(name, argv, nopool=True):
    if 'closure' in argv:
        i = argv.index('--ckpt-minutes')
        argv = argv[:i + 1] + ['100000'] + argv[i + 2:]
    return _run(name, argv, nopool)


C.run = run
C.main()
