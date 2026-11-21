# ic-bench: benchmark harness for relation finders

Code: `csd/bench/` (instances, solver registry, isolated runner, statistics, cost model), driver `scripts/bench.py`,
tests `tests/test_bench.py`. Protocol and predictions: `PROTOCOL.md`. Outcomes: `RESULTS.md`.

    python3 scripts/bench.py solvers
    python3 scripts/bench.py run --grid smoke --solvers py_ref,enum,lrf --aa --pin 0 --out experiments/ic-bench/x.jsonl
    python3 scripts/bench.py summarize experiments/ic-bench/x.jsonl --baseline lrf
    python3 scripts/bench.py model --preset ecc2k130 --m 4 --from-records experiments/ic-rootfind/n131.jsonl

Records: `c1-calib*.jsonl` (first calibration), `c1b-calib*.jsonl` (after the output-buffer fix), `c2-planted.jsonl` (planted
targets, all solvers), `c2-summary.md`, `model-ecc2k130-m{3,4,5}.md`, `model-sensitivity.txt`; one `*.meta.json` per
record (commit, CPU, load, source hashes). Each JSONL row is one (instance, solver, arm).

Adding a solver: subclass `csd.bench.solvers.Solver` (`available`, `applicable`, `prepare`, `run`), decorate with
`@register`; return decompositions as tuples of `(x, y)` so the runner can re-add them to `R` and compare sets.
