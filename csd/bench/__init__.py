"""Benchmark harness for relation finders of index calculus on binary curves.

    spec.py     instance specifications and their deterministic construction
    solvers.py  solver registry (LRF, enumeration, meet in the middle, plain Python, MutantXL)
    runner.py   isolated runs: forked child, timeout, peak RSS, A/A arms, randomized order, cross-checks
    stats.py    medians, bootstrap intervals, paired ratios, log2-linear fits, noise floor
    model.py    end-to-end index-calculus cost model against Pollard rho (assumptions are explicit parameters)
    grids.py    named parameter grids

Driver: scripts/bench.py.  Protocol and records: experiments/ic-bench/.
fixtures  frozen instance sets with expected solutions and recorded reference timings
"""
