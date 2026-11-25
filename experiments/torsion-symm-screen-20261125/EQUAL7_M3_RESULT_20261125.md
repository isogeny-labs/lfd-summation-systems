# Rank-unsaturated equal-size Koblitz `S4` screen — 2026-11-25

**Decision:** the native torsion-symmetrized three-summand presentation
lowers exact Boolean LFD on the frozen targets, including a target where
both systems have no roots. It does not increase useful relation yield at
equal base size: independent factor-base rank is 6 in both arms, and
target admission falls from 12 to 6. The joint gate fails again.

This is an adaptive exploratory follow-up to the equal-11 panel, whose
control saturated rank and admitted all 24 targets. The
`PROTOCOL_EQUAL7_M3_20261125.md` base-size-only rule was committed before
target outcomes were measured. An independent verifier replayed the rule
and confirmed `V=span_F2[1,2,4]` is its first match. This panel uses the
same `K_1/F_(2^9)` curve, `Q=(135,123)`, and first 24 targets as the
equal-11 panel. Direct `x` and symmetrized `u=1/(x+1)` bases each have
seven rational abscissae, 13 signed points and six nonzero projected
columns; they share one abscissa.

| Arm | Admitted / 24 | Failed | Checked ordered triples | Distinct rows | Base-column rank |
|:--|--:|--:|--:|--:|--:|
| Direct `x` | 12 | 12 | 96 | 20 | 6 |
| Symmetrized `u` | 6 | 18 | 96 | 20 | 6 |

Both arms attempted all `13³=2,197` ordered signed triples. The verifier
independently enumerated unordered triples, checked every saved row by the
group law and subgroup logs modulo 37, and recomputed the ranks. The
distinct rows and rank are for the factor-base columns with known target
multiples on the right-hand side. No repeated sign, torsion or permutation
image was counted as a new row.

All 512 direct and 128 symmetrized Boolean assignments were checked against
finite-field polynomial evaluation on each of the first four targets.
Every exact iterated LFD cell finished under the degree-8, 512-column,
512 MiB per-matrix and 60-second caps.

| Target scalar | Direct roots / group lifts | Symmetrized roots / group lifts | Direct input degree / exact LFD | Symmetrized input degree / exact LFD |
|--:|:--|:--|:--|:--|
| 1 | 6 / 6 | 0 / 0 | 6 / 7 | 4 / 4 |
| 2 | 6 / 6 | 0 / 0 | 6 / 7 | 4 / 4 |
| 3 | **0 / 0** | **0 / 0** | **6 / 7** | **4 / 4** |
| 4 | 9 / 9 | 0 / 0 | 6 / 7 | 4 / 4 |

Target 3 supplies a same-target failed-class comparison. The measured
advantage is exact LFD of these Boolean presentations, not the native F4
solver's total solving degree. The absence of a rank or admission gain
rules out the requested *joint* improvement for this frozen base and
target stream, regardless of an unmeasured total solver degree.

The saved shared preflight took 0.0011 wall/CPU seconds; construction and
exhaustive checking of both first-four systems took 0.0498 wall seconds.
Full 24-target group checks, failed triple attempts, row deduplication and
rank updates took 0.0086 wall / 0.0086 CPU seconds per arm. Four bounded
closure calls took 0.6326 wall seconds for direct and 0.3187 for
symmetrized. These are component costs under shared host load and
different positive/failed mixes. Base construction and every failed
attempt are included; a complete factor-base solve, descent and matched
cold same-target rho were not run after the useful-yield gate failed.

## Raw evidence and replay

- `equal7-m3.json`: bases, equations, complete first-four roots, all
  24-target relation rows and component costs; SHA-256
  `2c22a9d138a0ea5902dbaddeed2b90156f990a8e62f04fac1c906c2ed124aa58`.
- `equal7-m3-degree.jsonl`: eight exact LFD profiles, caps and wall/CPU;
  SHA-256 `00125cdeeb4a84a5e25b65210b5ef052e85e5818966e0b7c8bcf7daef58a732a`.
- `equal7-m3-verification.json`: replayed base selection, curve order,
  subgroup-log row checks and root lifts; SHA-256
  `fdcbff714c383775375aad66ca93e3e4ab4afea1d8651cf0b66641213895596c`.

```sh
python3 experiments/torsion-symm-screen-20261125/equal7_m3.py
python3 experiments/torsion-symm-screen-20261125/verify_equal7_m3.py
python3 experiments/torsion-symm-screen-20261125/equal7_m3_degree.py
```
