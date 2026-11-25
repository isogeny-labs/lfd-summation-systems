# Cold same-target cost gate

The [frozen protocol](COLD_TARGET_PROTOCOL.md) and [target manifest](cold_target_manifest.json) define three nonzero upper-half subgroup scalars per cell. For each public point Q=sG, the generic arm recovered s and verified [s]G=Q. The image arm constructed the exact fixed factor base and solver, found all three-summand relations for Q, and checked every signed triple against the group law. Each arm ran in a fresh native ARM process three times per target, with order rotated. All 27 pairs completed.

| Field, base dimension | Matched pairs where BSGS was faster on wall and CPU | Median BSGS wall | Median image first-query wall | Median BSGS CPU | Median image CPU | Image relations per target |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 2^23, 11 | 9/9 | 10.2 ms | 22.6 ms | 10.2 ms | 22.4 ms | 162, 180, 204 |
| 2^29, 13 | 9/9 | 57.5 ms | 85.0 ms | 57.2 ms | 84.1 ms | 190, 177, 171 |
| 2^29, 14 | 9/9 | 56.9 ms | 1,024.1 ms | 56.8 ms | 1,020.5 ms | 1,404, 1,370, 1,406 |

These figures are complete worker costs from field/curve construction through verification. The parent also measured 6.4 ms wall and 6.0 ms CPU to regenerate and check the entire frozen target manifest, a common cost across both arms. The generic solver is written in Python; the image kernel uses native C. Native solver compilation happened before both worker timings and was not charged to the image arm, which favors the image method. Median peak process RSS was about 36 MB for generic BSGS and 37–43 MB for the image arm. [Raw paired receipts](cold_target_results.json.gz) include full relation sets, recovered scalars, source hashes, wall/CPU and RSS.

**Interpretation.** A full relation collection that uses this exact image implementation to query these public targets must do at least the measured image setup and first verified query, followed by more relations and linear algebra. That measured lower bound already exceeds the complete generic solve on every matched run. The earlier three-summand improvement versus native LRF remains a valid relation-query component result, but it does not cross the cold same-target cost gate here. The gate does not rule out a different relation family, factor base, field size, target transformation, or improved implementation. No large-field or asymptotic conclusion follows.
