# Note on verification.json (2026-11-12)

For three source files, the SHA-256 recorded in `verification.json` no longer matches the current
tree: those files were edited after this run. The hashes are left as recorded, because they
identify the sources the run actually used.

## Replay with the current sources (2026-11-14)

`replay-20261114.json`: `python scripts/replay_with_current_sources.py experiments/ic-chain-depth8 -- --workers 3`
recomputed all 22 systems with the current sources and found every result equal to `results.jsonl` (all fields
except the timing). Two of the source hashes in `verification.json` differ from the current files
(`csd/ic_chain_shift_closure.py`, `scripts/run_ic_chain_depth8.py`); the data hashes all match.
