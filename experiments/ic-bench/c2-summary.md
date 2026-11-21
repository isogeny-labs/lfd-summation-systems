## Row status

| solver | ok | skipped |
|---|---|---|
| enum | 24 | 0 |
| lrf | 24 | 0 |
| mitm4 | 6 | 18 |
| mutantxl | 6 | 18 |
| py_ref | 24 | 0 |

## Per cell: median of t_min over seeds [IQR]; ndec = decompositions summed over seeds

| n | k | m | basis | enum | lrf | mitm4 | mutantxl | py_ref | ndec |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 4 | 3 | standard+planted | 110 µs [106 µs–120 µs] | 94.1 µs [87.5 µs–106 µs] | — | 1.58 s [1.58 s–1.66 s] | 590 µs [536 µs–726 µs] | 14 |
| 13 | 5 | 3 | standard+planted | 139 µs [125 µs–182 µs] | 155 µs [126 µs–155 µs] | — | — | 3.4 ms [2.89 ms–4.21 ms] | 22 |
| 16 | 6 | 3 | standard+planted | 207 µs [202 µs–226 µs] | 141 µs [140 µs–149 µs] | — | — | 14.8 ms [13.8 ms–15.7 ms] | 7 |
| 13 | 4 | 4 | standard+planted | 278 µs [262 µs–421 µs] | 310 µs [241 µs–345 µs] | 266 µs [236 µs–356 µs] | — | 10.4 ms [9.72 ms–14 ms] | 21 |

## Growth: log2 of median t_min against k, per (solver, m, basis); slope = doublings per step of k

| solver | m | basis | slope | 95% CI | R² | points |
|---|---|---|---|---|---|---|
| enum | 3 | standard | 0.45 | 0.33–0.58 | 0.977 | 3 |
| lrf | 3 | standard | 0.29 | -0.14–0.72 | 0.583 | 3 |
| py_ref | 3 | standard | 2.32 | 2.12–2.53 | 0.997 | 3 |

## A/A noise floor (same solver, same instance, two arms): B/A median 0.998, 5–95% 0.763–1.369 over 42 pairs. Ratios inside this band are not resolved.

## Peak RSS of the run (MiB, max over instances): enum 27, lrf 27, mitm4 27, mutantxl 64, py_ref 27

skipped mitm4: defined for m = 4 only

skipped mutantxl: 15 Boolean variables > 12; 16 Boolean variables > 12; 18 Boolean variables > 12
