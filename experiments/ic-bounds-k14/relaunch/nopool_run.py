#!/usr/bin/env python3
"""Run a script with the CuPy device and pinned memory pools disabled.

usage: python3 nopool_run.py SCRIPT [args...]
The recorded k = 14 runs of 2026-11-10 were made with the pools disabled; the switch was not in the
committed drivers, so this wrapper supplies it. Exit 3 if CuPy cannot be imported.
"""
import runpy
import sys

try:
    import cupy
except Exception as e:   # noqa: BLE001
    print('nopool_run: cannot import cupy: %r' % (e,), flush=True)
    sys.exit(3)
cupy.cuda.set_allocator(None)
cupy.cuda.set_pinned_memory_allocator(None)
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
