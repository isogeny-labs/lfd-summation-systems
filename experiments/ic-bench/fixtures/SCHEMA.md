# Frozen benchmark fixtures

Instances written out once, so that a benchmark does not depend on Python's `random`, on `csd.ic_dreg`'s choice of
field polynomial, or on the generator code. Each `*.jsonl` file is one set, one instance per line; `MANIFEST.json`
holds the SHA-256 and row count of every file (`bench.py fixtures verify` checks it).

## Row (format 1)

| field | content |
|---|---|
| `id` | `<key>@<set>/<label>`, the id a `Spec.fixture` carries |
| `spec` | the `Spec` fields (n, k, m, seed, basis, j, curve, plant) the instance was made from |
| `modulus` | hex field polynomial; `instance()` refuses to build if `csd.ic_dreg` now picks another |
| `a`, `b` | curve coefficients y^2 + xy = x^3 + a x^2 + b (`b` hex) |
| `basis` | factor-base basis, hex field elements (V = their F_2-span) |
| `R` | target point `[x, y]`, hex |
| `expected` | `complete` (an independent solver agreed with LRF: py_ref, enum or mitm4), `verified_by`, `count`, `sha` (of the canonical set), `decs` (the set itself when at most 64), `must_contain` (planted decomposition) |
| `reference` | recorded timings of the same instance (`t_lrf`, `t_brute`, `t_mitm`, `ns_per_target`) from `experiments/ic-rootfind/`, or `{}` |
| `provenance` | where the instance came from |

`rust-records` has no expected set (`complete` false): the instances are for the Rust kernel up to n = 163, with
the recorded `ns_per_target`.

## Use

    python3 scripts/bench.py fixtures list
    python3 scripts/bench.py fixtures verify [--deep]       # --deep re-runs every verifying solver
    python3 scripts/bench.py run --fixture m3-natural --filter k10 --solvers enum,lrf --reps 3 --pin 0

A fixture run compares each solver's output with the frozen expected set (`fixture_check`; a difference is status
`mismatch`) and reports `ratio_to_reference` against the recorded timing, which is machine-dependent: read it as
"this machine versus the one that recorded it", not as a result.

## Regenerating

`bench.py fixtures make [--set NAME]` rewrites sets (`m3-large-planted` runs enumeration, about 3 minutes). A
regenerated set that differs from the committed one means the generator or the solvers changed: look at the diff
before committing it.
