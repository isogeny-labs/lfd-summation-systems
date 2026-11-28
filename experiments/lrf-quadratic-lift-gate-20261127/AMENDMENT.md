# Timing and memory audit amendment

The frozen protocol specified a 30-second solver-call cap and memory
cost. The first implementation applied the 30-second cap to the
**entire Sage worker process**, including interpreter startup. This is
a stricter cap than the wording implied. The two `k=11` positive
timeouts are therefore censored at 30 seconds of worker wall time;
their internal solver time is not known and must not be reported as
exceeding 30 seconds.

Before re-running the same deterministic fixture, the worker was
extended to record `resource.RUSAGE_SELF.ru_maxrss` on macOS, in bytes,
for completed cases. A killed timeout worker cannot return this value;
its peak memory is explicitly unknown. All fixture selection rules,
equations, backend, and comparison arms remain frozen. The re-run
updates the raw timing record rather than pooling two host runs.
