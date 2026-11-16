# Stage B readings of an independent instrument (three to six summands, split presentation)

This directory holds a copy of scored records from another repository and a
comparison with the certified values of `experiments/ic-bounds-m3/`. Nothing
here was computed in this repository.

## Source

`isogeny-labs/crypto`, directory `research/pkm_falldegree_unblock_20261107`,
commit `cb4dbc2` (`SOURCE.json` has the paths and SHA-256 of each copied
file). The study is preregistered there
(`research/prime_fall_degree_20261027/PREREGISTRATION.md`); "stage B" is its
ladder in the depth for three to six summands. Its records are labelled
unverified in that repository: no certificate of the kind used in this
paper exists for them.

## Systems

The "Kummer tower" arm of that study is the split presentation of
`experiments/ic-bounds-m3/README.md`, generalized to `m` summands:

- `m` squaring chains of length `t` (the depth `k` of the paper) for a coset
  `g·μ_{2^t}` of abscissae;
- `m − 2` free auxiliary variables, the abscissae of the partial sums;
- the left chain `S_3(x_1, x_2, u_1)`, `S_3(u_j, x_{j+2}, u_{j+1})`,
  `S_3(u_{m−2}, x_m, t)`.

Source: `examples/pkm_tower_pilot.rs`, function `summation`, of that
repository. Targets are random, so the systems are inconsistent (unit
ideals). Primes `786433` and `2013265921`, two draws per cell and prime.

## Instrument

A Rust degree-closure engine (`hky_fp_tower`) for the last fall degree, an
F4 implementation for the solving degree `D`, and a gate that requires
`d_F ≤ D` on every cell. It shares no code with the Python engines here.

## Comparison (`python3 scripts/check_crossrepo_stage_b.py`)

| three summands, depth | certified here (unit instances) | stage B draws |
|---|---|---|
| 2 | 5 on 63 of 63 | 5, 5, 5, 5 |
| 3 | 6 on 51 of 51 | 6, 6, 6, 6 |
| 4 | 6 on 27 of 27 | 6, 6, 6, 6 |

The script checks the digests of the copied files, rebuilds both columns
from the records and fails on any difference; the result is
`comparison.json`. It runs no closure.

## Readings that exist only in stage B

Last fall degree, both draws at both primes unless noted:

| summands \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| 3 | 5 | 6 | 6 | 7 | ≥ 7 |
| 4 | 6 | 7 | 7 | ≥ 7 | |
| 5 | 6 | 8 | ≥ 8 | | |
| 6 | 7 | | | | |

- Four summands, depth 4: from `M4-t4-k{0,1}-c7.jsonl` (prime `2013265921`):
  the closures of degree 3 to 6 do not contain 1 and the closure of degree
  7 does.
- Five summands at depth 3 and six at depth 2: the closure is complete
  without 1 through degree 7 (respectively 6) and F4 ends at degree 8
  (respectively 7), so the value follows from `d_F ≤ D`.
- Three summands, depth 5, prime `786433`: one draw stopped on the time
  budget (`≥ 7`), the other reads 7.
- Five and six summands at depth 2 use a base of 4 points, fewer than the
  number of summands.
- First fall degree: 4 on every cell.

These are readings of one instrument without certificates. The paper
reports them as such (Section "Three to eight summands").
