# A1: macOS memory-limit preflight

28 November 2026, before any msolve process was launched or result recorded.
The first invocation of `measure.py` failed in Python's `preexec_fn`: this
macOS host rejects setting `RLIMIT_AS` to 2 GiB with `ValueError: current
limit exceeds maximum limit`, despite reporting the current soft and hard
limits as `9223372036854775807`. A direct `setrlimit` preflight reproduced
the failure. The process exit occurred before `exec`, so no arm was timed.

Amend the measurement rule for this host: retain the 20-second process cap,
one thread, fixed solver seed, and three repetitions, but do not claim a hard
memory cap. Record the cumulative child `ru_maxrss` only as a diagnostic;
it is not a per-arm peak. If a process causes memory pressure, stop the panel,
record which cells were not run, and do not assign a timing win to censored
cells. The original protocol remains intact for custody.
