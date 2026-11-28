# Quadratic root lift: frozen residual-query pilot

## Result

The lifted Boolean equations represent the **same two-point relations**
as native LRF on the eight selected residual queries. They have input
degree two, but SageMath/Singular did not make a complete query cheaper
on this frozen panel. At `k=11`, both positive cases hit the 30-second
solver cap. This is a backend and workload result, not a lower bound on
all auxiliary-variable solvers, nor a measurement of complete ECDLP
cost. The maximal intermediate solver degree was not exposed by the
backend and is **unmeasured**.

## Frozen fixture and independent audit

The [protocol](PROTOCOL.md) fixes source hashes, seed 1, curve recipe,
`n=23`, `(k,delta)=(10,6),(11,9)`, and the first two positive and first
two empty solvable residual systems per cell. The curve coefficient is
`b=1127129` for both cells. Selected attempt indices are in the table
below. Selection by outcome makes this a diagnostic panel rather than
an unbiased relation-yield estimate.

`audit.py` enumerated **all 4,456,448 Boolean assignments** across the
eight systems and found exactly the saved solution pairs. It compared
every signed point row with the Python LRF reference and native C LRF,
and verified point addition and factor-base membership. The audit
passed even for both timed-out solver cases. On every completed solver
case, Sage's extracted variety exactly matched the direct Boolean
solution set. The timeout cases have no extracted solver result.

## Timings

Seconds on macOS 27.2 arm64, Python 3.13.1, SageMath 10.10.rc0.
Each row is one deterministic execution, with wall and CPU retained in
`results.json`; no uncertainty interval is inferred. The native column
times the C `Solver.lrf(2,T)` query including relation recovery. The
Sage query column starts after interpreter import and includes equation
construction, Gröbner basis, variety extraction, and answer packaging.
The worker-total column also charges per-query Sage process startup;
it is a conservative cold cost. Shared native `Solver` setup was
0.0115 s at `k=10` and 0.00143 s at `k=11`, recorded separately.

| `k,delta` | Outcome | Attempt | Boolean solutions | Native C wall | Sage query wall | Worker total wall | Sage status |
|---|---|---:|---:|---:|---:|---:|---|
| `10,6` | empty | 7 | 0 | 0.000056 | 0.0902 | 2.408 | exact match |
| `10,6` | empty | 8 | 0 | 0.002176 | 0.0595 | 1.848 | exact match |
| `10,6` | positive | 19 | 2 | 0.004412 | 1.936 | 3.779 | exact match |
| `10,6` | positive | 80 | 2 | 0.002571 | 1.114 | 2.879 | exact match |
| `11,9` | empty | 2 | 0 | 0.002649 | 0.3785 | 2.097 | exact match |
| `11,9` | empty | 3 | 0 | 0.002509 | 0.4201 | 2.180 | exact match |
| `11,9` | positive | 29 | 2 | 0.002756 | >30 | >30 | censored at 30 s |
| `11,9` | positive | 34 | 4 | 0.002966 | >30 | >30 | censored at 30 s |

For `k=10` positives, Gröbner basis construction itself took about
0.047 s; extracting the solutions took 1.05–1.87 s. The two `k=11`
positive timeout records do not locate the slowdown within the worker.
The solver's 30-second cap makes their time a lower bound, not a
completed solve. `results.json` retains all CPU, setup, native stats,
solution rows, and status fields.

## Scientific interpretation and next gate

The high-degree projected good-set ANF does **not** prove that kernel
scanning is necessary. The auxiliary-root presentation is an explicit
quadratic counterexample to that inference. This pilot instead shows
that this particular Sage/Singular execution is slower than the native
scan on its matched residual systems. A different solver or an
output-sensitive extractor would need a new frozen comparison with
identical relation sets and full setup, failure, recovery and
verification cost. The lift leaves the outer factor-base scan and
relation yield unchanged, so even a residual-query gain would still
need a cold same-target collection and linear-algebra comparison with
generic solving.

## Replay

From the repository root:

```sh
python3 experiments/lrf-quadratic-lift-gate-20261127/run.py
python3 experiments/lrf-quadratic-lift-gate-20261127/audit.py
```

The first command rewrites `results.json` and timings vary by host. The
second command compares every Boolean assignment against the saved
solution sets and native references without invoking the Gröbner
solver.
