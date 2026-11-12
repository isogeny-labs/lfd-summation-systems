#!/usr/bin/env python3
"""Run a script with the CuPy DEVICE memory pool disabled and the pinned-memory pool left at its default.

usage: python3 nopool_dev_only.py SCRIPT [args...]
This is the setting of the recorded runs of 2026-11-10 (scratch wrapper `nopool.py`). `nopool_run.py` also
disables the pinned-memory pool; with it the self-test cases ran about five times slower on 2026-11-11.
On the pod this file is uploaded under the name `nopool_run.py`, which `chain_k14.py` calls.
"""
import runpy
import sys

try:
    import cupy
except Exception as e:   # noqa: BLE001
    print('nopool_dev_only: cannot import cupy: %r' % (e,), flush=True)
    sys.exit(3)
cupy.cuda.set_allocator(None)
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
