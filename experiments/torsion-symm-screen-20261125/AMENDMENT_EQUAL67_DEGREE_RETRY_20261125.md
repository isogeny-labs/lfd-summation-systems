# Equal-67 degree-run retry

The first `equal67-degree.jsonl` run failed before invoking any closure:
macOS rejected `setrlimit(RLIMIT_AS, 512 MiB)` with `ValueError: current
limit exceeds maximum limit`. All eight environment-error cells remain in
that append-only store. The source now passes `max_bytes=512 MiB` directly to
`local_profile`, which bounds each closure and evaluation matrix as frozen in
`PROTOCOL_EQUAL67_20261125.md`; the 60-second process wall cap and 8192-column
cap remain unchanged. The retry uses a separate append-only store,
`equal67-degree-retry.jsonl`.
