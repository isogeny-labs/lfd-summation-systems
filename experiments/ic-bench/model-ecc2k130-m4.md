n=131 m=4  t_lin=1849.72 ns  t_cand=260 ns  t_add=150.05 ns  t_la=10 ns (assumed)  rho: 2^60.81 iterations

| k | log2 mu | log2 relations (s) | log2 algebra (s) | log2 total (s) | total core-years | algebra memory | gap to rho (log2) |
|---|---|---|---|---|---|---|---|
| 42 | 32.4 | 73.6 | 58.4 | 73.6 | 4.57e+14 | 158 TiB | 35.5 |
| 43 | 36.4 | 72.6 | 60.4 | 72.6 | 2.29e+14 | 316 TiB | 34.5 |
| 44 | 40.4 | 71.6 | 62.4 | 71.6 | 1.14e+14 | 632 TiB | 33.5 |
| 45 **←** | 44.4 | 71.6 | 64.4 | 71.6 | 1.14e+14 | 1.26e+03 TiB | 33.5 |
| 46 | 48.4 | 72.9 | 66.4 | 72.9 | 2.85e+14 | 2.53e+03 TiB | 34.8 |
| 47 | 52.4 | 74.8 | 68.4 | 74.8 | 1.05e+15 | 5.06e+03 TiB | 36.7 |
| 48 | 56.4 | 76.8 | 70.4 | 76.8 | 4.17e+15 | 1.01e+04 TiB | 38.7 |

rho: 3.03e+11 core-seconds = 9.6e+03 core-years
floor with a FREE relation finder (t_lin = 0): 2^35.7 core-seconds at k = 30 (-2.5 vs rho in log2; algebra memory 2^-4.9 TiB)
required speed-up of the relation finder to reach rho (k re-optimized): 2^47.0 = 1.4e+14x, i.e. t_lin <= 1.32e-11 ns, break-even at k = 31
per-target budget for ANY relation finder (all decompositions of one R, one core), best k = 31: 60.9 ms  (mu = 2^-11.6, 2^30.1 relations, 187 s per relation)
Assumptions: |FB| = 2^k, mu = |FB|^m/(m! 2^n), t_la is not measured, kernel limits (e.g. the Rust kernel needs k >= 35 at n = 131) are not modelled.
