# Equal-size three-summand Koblitz torsion-frame screen — 2026-11-25

**Decision:** the native symmetrized `S4` presentation reduces exact Boolean
LFD from 7 to 5 on the frozen first four targets, including a target with
rational decompositions in both arms. The equal-size symmetrized base does
not improve independent useful-relation yield on the same 24-target stream:
both arms reach factor-base rank 10, while target admission falls from 24
to 16. This candidate does not pass the joint gate.

## Frozen construction

`PROTOCOL_EQUAL11_M3_20261125.md` fixes `K_1/F_(2^9)`, modulus 515, group
order `518=14·37`, subgroup generator `Q=(135,123)`, and targets
`Q,...,24Q`. The base-size-only first-match rule fixes
`V=span_F2[1,2,4,56]`. The direct `x∈V` and symmetrized
`u=1/(x+1)∈V` bases each have 11 rational abscissae, 21 signed points,
and 10 nonzero cofactor-projected columns. They share one abscissa.
The symmetrized system uses the native 18-term `S4` exponent list from
`koblitz_symmetrised.rs` and one torsion-parity bit.

## Checked relations and degree

Both arms attempted all `21³=9,261` ordered signed triples. An independent
unordered-triple enumerator checked the saved rows against the group law and
the order-37 subgroup logs. Sage separately confirmed rank 10 for each
factor-base matrix over `GF(37)`.

| Arm | Admitted / 24 | Failed | Checked ordered triples | Distinct rows | Base-column rank |
|:--|--:|--:|--:|--:|--:|
| Direct `x` | 24 | 0 | 400 | 67 | 10 |
| Symmetrized `u` | 16 | 8 | 348 | 67 | 10 |

All 4,096 direct and 1,024 symmetrized Boolean assignments were checked
against independent finite-field polynomial evaluation on each of the first
four targets. Exact iterated closure was bounded by degree 8, 4,096
normal-form columns, 512 MiB per matrix and 60 seconds per cell; all eight
cells completed.

| Target scalar | Direct input degree / exact LFD | Symmetrized input degree / exact LFD | Direct polynomial roots / group lifts | Symmetrized roots / group lifts |
|--:|:--|:--|:--|:--|
| 1 | 6 / 7 | 4 / 5 | 12 / 12 | 0 / 0 |
| 2 | 6 / 7 | 4 / 5 | 9 / 9 | 0 / 0 |
| 3 | **6 / 7** | **4 / 5** | 61 / 46 | 9 / 9 |
| 4 | 6 / 7 | 4 / 5 | 15 / 15 | 0 / 0 |

Target 3 is a same-target positive-class comparison; the other three
targets are positive for direct and failed for the symmetrized base.
Fifteen of the direct `S4` coordinate roots on target 3 have no rational
three-point lift, so polynomial roots were never counted as useful
relations without group verification. The exact LFD is a property of the
saved Boolean presentations; the native F4 solver's total solving degree
was not measured.

## Charged components and limits

Shared base/target preflight took 0.0011 wall/CPU seconds. Constructing and
exhaustively checking both systems on four targets took 0.8708 wall and
0.8689 CPU seconds. The complete 24-target triple enumeration, group
checks, row deduplication and rank updates took 0.0369 wall / 0.0367 CPU
seconds for direct and 0.0358 / 0.0358 for symmetrized. The four exact
closure calls took 14.9242 wall / 14.1900 child-CPU seconds for direct and
0.8716 / 0.6475 for symmetrized. These are shared-host component timings;
the different positive/failed mixes prevent a same-target solver-speed
conclusion from the aggregate closure times. Failed relation attempts and
setup are charged in the saved records. A complete factor-base solve,
descent and cold same-target rho were not run because the candidate does
not pass the useful-yield gate.

## Evidence and replay

- `equal11-m3.json`: frozen bases, equations, complete first-four roots,
  24-target relation rows and costs; SHA-256
  `c81944225f8fb2a4ece9030875942cc9fd01dd45ecf8ba6094aa9580ef60b8a2`.
- `equal11-m3-degree.jsonl`: eight exact closure profiles and wall/CPU;
  SHA-256 `0a589beade156149f4322b98c6afe2a7e0fce6b8689bd1bb32947ef5e1972b34`.
- `equal11-m3-verification.json`: independent curve order, rows,
  subgroup-log checks and polynomial-root group lifts; SHA-256
  `bca2f0eea0bb9759f81f99ab8c40315e3adcf57590c051f7b4ffc5c94ee8a51b`.

```sh
python3 experiments/torsion-symm-screen-20261125/equal11_m3.py
python3 experiments/torsion-symm-screen-20261125/verify_equal11_m3.py
python3 experiments/torsion-symm-screen-20261125/equal11_m3_degree.py
```
