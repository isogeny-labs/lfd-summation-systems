# Equal-size Koblitz torsion-frame screen — 2026-11-25

**Decision:** this frozen equal-size `x`/`u` comparison does not pass either
part of the proposed joint gate. The symmetrized base admits fewer of the
same 24 targets, and exact iterated LFD is equal on the first four targets.
The earlier binary Edwards colored-base screen, frozen n=9 Boolean screen,
ECC2K-130 recipe screen, and prime-field subgroup question remain separate.

## Frozen instance and constructions

The protocol is `PROTOCOL_EQUAL67_20261125.md`. On
`K_1/F_(2^15): y²+xy=x³+x²+1`, field modulus `0x8003`, the independently
recounted curve order is `32,494=154·211`. The fixed point
`Q=(15115,1293)` has order 211. Targets are exactly `Q,...,24Q`.

The two bases use the same frozen dimension-seven Frobenius-invariant
subspace `V=span_F2[1,892,1694,3678,5682,8344,22848]` but different
coordinate frames. The direct base uses `x∈V`; the symmetrized base uses
`u=1/(x+1)∈V`. Each has 67 rational abscissae, 133 signed rational points,
and 60 nonzero columns after projection by 154. Their abscissa sets overlap
in seven values. The verifier checked `x(P+(0,1))=1/x(P)` and
`u(P+(0,1))=u(P)+1` at every finite rational `u`-base abscissa; `x=0`
is the identity boundary.

The direct Boolean system is the repository's complete Weil descent of
`S3(x1,x2,xR)` (14 unknown bits). The symmetrized system uses the native
`w1*w2*wR+w1+w2+wR+s=0` formula (13 unknown bits), with `w=u²+u` and
one torsion-parity bit. Both have input degree 2 on these target-fixed
systems. All 16,384 direct and 8,192 symmetrized assignments were checked
against independent finite-field evaluations for each of the first four
targets. Every symmetrized polynomial root across all 24 targets had a
finite rational group-law lift; no boundary or nonrational root choice
occurred in this panel.

## Independent useful relations

Both arms attempted all `133²=17,689` ordered signed base pairs and retained
all 24 target outcomes, including failures. Relation vectors are deduplicated
modulo 211; the independent verifier re-enumerated unordered pairs, checked
every saved vector against discrete logs of the cofactor-projected points,
and recomputed rank modulo 211.

| Arm | Admitted targets / 24 | Failed targets | Distinct checked rows | Final rank |
|:--|--:|--:|--:|--:|
| Direct `x` | 12 | 12 | 12 | 12 |
| Symmetrized `u` | 6 | 18 | 12 | 12 |

The direct arm admits scalars `1,5,7,10,13,14,15,16,19,21,23,24`;
the symmetrized arm admits `5,7,8,12,17,18`. Each arm has 24 checked
ordered pairs in total, but the candidate concentrates two independent rows
on each of six targets. Neither admitted-target count nor rank increases.
These are exact finite-group counts for the frozen target stream, not a
relation probability estimate on larger fields.

## Exact LFD on the frozen degree panel

The first degree-run store records eight environment errors before closure:
macOS rejected the requested process `RLIMIT_AS` setting. The
`AMENDMENT_EQUAL67_DEGREE_RETRY_20261125.md` retry uses a distinct store,
retains the same 60-second per-cell wall cap, 8,192-column cap, and the
512 MiB cap applied by `local_profile` to each closure/evaluation matrix.
Every retry cell completed without censoring.

| Target scalar | Direct roots | Symmetrized roots | Direct exact LFD | Symmetrized exact LFD |
|--:|--:|--:|--:|--:|
| 1 | 2 | 0 | 3 | 3 |
| 2 | 0 | 0 | 3 | 3 |
| 3 | 0 | 0 | 3 | 3 |
| 4 | 0 | 0 | 3 | 3 |

The first target has a proper direct ideal and a unit symmetrized ideal;
targets 2–4 are failed decomposition cells in both arms. The exact LFD
values are equal on those matched failed cells. These are iterated Boolean
closure measurements for the saved equations, not the native F4 solver's
total solving degree.

## Charged components and limits

The saved shared base/target preflight took 0.0061 wall / 0.0060 process-CPU
seconds. Constructing and exhaustively checking both systems on the first
four targets took 0.8062 / 0.8050 seconds. The full 24-target pair enumeration,
group checks, row deduplication and rank updates took 0.1028 / 0.1027 seconds
for the direct arm and 0.1032 / 0.1030 for the symmetrized arm. The four
bounded closure calls took 0.3723 wall / 0.3532 child-CPU seconds for direct
and 0.3550 / 0.3411 for symmetrized. These timings charge failed pair
attempts and setup, but were collected on a shared host and are component
measurements. Summing the saved preflight, both equation/root checks, both
24-target relation inventories and all eight successful closure calls gives
1.7458 wall and 1.7110 process/child-CPU seconds for the **combined panel**;
the source records overlap no timed phases. This is an audit total for the
screen, not an end-to-end logarithm cost or an arm-level speed ratio. The
first degree-run errors are retained separately.

A complete factor-base linear solve, descent and cold same-target rho were
not run. The candidate already fails the strict yield and LFD gates, so the
saved component timings do not support a favorable total-cost claim. This
panel does not imply a general obstruction for other torsion-equivariant
bases or for the prime-field subgroup systems.

## Evidence and replay

- `equal67-koblitz.json`: bases, equations, complete first-four root sets,
  all target relation rows, timings; SHA-256
  `75291f8313ccc6acde200efd098cea9f8104287d26b273f46798c08be5a40d23`.
- `equal67-verification.json`: curve order, torsion action, subgroup-log row
  checks, rank, and all 24 symmetrized root counts; SHA-256
  `80bc68d15e2fbce24e4578cdc04097607480e581cef79d6735791980e10855c8`.
- `equal67-degree.jsonl`: retained macOS limit errors; SHA-256
  `7a8a679b95d3cae710cae369ecf519f462f8e1304e79780ac5e8f7706e303fad`.
- `equal67-degree-retry.jsonl`: exact LFD profiles and per-cell wall/CPU;
  SHA-256 `fccb8e03e794da06d9d1d4c1cc5cb740ec00b50c363eb440daa355baa4fa254e`.

```sh
python3 experiments/torsion-symm-screen-20261125/equal67_koblitz.py
python3 experiments/torsion-symm-screen-20261125/verify_equal67.py
python3 experiments/torsion-symm-screen-20261125/equal67_degree.py \
  --store experiments/torsion-symm-screen-20261125/equal67-degree-retry.jsonl
```
