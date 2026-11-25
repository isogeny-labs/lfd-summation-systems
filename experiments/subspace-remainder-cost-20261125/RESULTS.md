# Matched binary two-summand component result

The frozen [protocol](PROTOCOL.md), [manifest](manifest.json), and [raw receipts](results.json) define six identical-target, identical-factor-base cells. Each timed arm ran in a fresh process on the same ARM host. The conditional-image arm returned the same signed point pairs as native LRF on every completed target; censored prefixes also matched. Three repetitions were run per cell, including empty and planted targets.

| Field, base dimension | Factor-base points | Independent random-target rows | Image wins vs both arms | Median native LRF wall (ms) | Median image wall (ms) | Median image CPU (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2^23, 8 | 273 | 0 | 0/3 (image censored) | 16.9 | — | — |
| 2^23, 10 | 1,057 | 5 | 1/3 | 19.3 | 18.8 among completed | 18.7 among completed |
| 2^23, 11 | 2,087 | 18 | 3/3 | 22.7 | 19.8 | 19.8 |
| 2^29, 12 | 4,033 | 2 | 1/3 | 22.9 | 20.9 among completed | 20.9 among completed |
| 2^29, 13 | 8,193 | 3 | 3/3 | 28.8 | 27.1 | 27.0 |
| 2^29, 14 | 16,337 | 15 | 3/3 | 52.5 | 39.1 | 39.1 |

The image times above are **arm-specific** totals including solver construction, conditional-image table construction, queries, and verification. The common fixture/target-generation time is separately recorded in each raw cell and added equally to all arms. Median common-charged wall time is 121.9 versus 119.0 ms at (23,11), 445.9 versus 444.2 ms at (29,13), and 838.2 versus 824.8 ms at (29,14), native LRF versus image respectively. Image-table allocations are 180,224, 851,968, and 1,835,008 bytes for these cells.

This is a component crossover under the frozen gate. The conditional-image condition does not change the relation family or skip the LRF kernel enumeration. The (29,13) margin is especially close. The small number of random-target rows means this panel does not establish a useful linear algebra system. The original worker attributed lazy image-table construction to its first query rather than its setup phase; **total wall, CPU, and memory include it**. Common fixture generation also prepares experimental inputs that each solver rebuilds, so the common-charged column is intentionally conservative and should not be interpreted as algorithmic preprocessing alone. No solver-degree reduction or large-field extrapolation follows from these data. The next gate is a frozen three-summand comparison with complete relation sets, failed attempts, and linear algebra charged.
