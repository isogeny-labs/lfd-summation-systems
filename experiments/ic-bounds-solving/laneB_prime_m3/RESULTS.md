# Lane B: three-summand prime chains, static solving degree sd vs d_ff, d_lf

Every result statement below carries [pcc:unverified] (no independent validator has run; all numbers are from this lane's own code).
Run 2026-11-15, host 4 cores shared, OMP_NUM_THREADS=1..2, `ulimit -v 6000000`, python3 + numpy 2.5 (no sympy/galois/flint available; none needed). Seeds are fixed in the scripts (`random.Random(str)`), so every scan is reproducible.

## What was built (all under `tools/`)

| file | role |
|---|---|
| `build_m3.py` | own S_4 = Res_X(S_3(x1,x2,X), S_3(x3,t,X)) as dict polynomial, exponents reduced mod N=2^k into an (N,N,N) array over F_p, delta = max weight of supp g, exact NTT over the 3 axes (int64, p < 2^31), instance samplers (random curve and target x-coordinate of an F_p point; the locus 5t^3+3at+6b=0), `--selftest` |
| `sd_inverse.py` | g^-1 = IDFT(1/DFT(g)); maxwt, zeros per weight class, top coefficient (N-1,N-1,N-1), predicted sd = delta + maxwt; optional check g * g^-1 = 1 by direct convolution (done k <= 3 on the first 3 units per scan) |
| `static_bruteforce.py` | ordinary-ring (3k variables) static Macaulay sd by exact sparse echelon elimination of all rows m*f (f in 3k chains + g, deg(mf) <= D), degree by degree; `sd_sparse` (tuple keys) and `sd_sparse_int` (packed int keys, 3x faster, identical ranks); `lfd_ff` calls `lfd_engine.exact_lfd` (independent engine) for d_lf and d_ff |
| `run_k2.py`, `bf_exceptions.py`, `run_k3_bf.py`, `run_k3_lfd.py` | drivers for the brute-force / LFD runs |
| `scan.py`, `run_scans.sh`, `run_extra.sh` | scans k = 2..8 |
| `summarize.py`, `plot.py` | `data/summary.json`, `data/summary_table.md`, `data/sd_vs_k_m3.svg` |

Commands (from this directory): `python3 tools/build_m3.py --selftest`; `python3 tools/run_k2.py 2` (12 min); `tools/run_scans.sh`; `tools/run_extra.sh`; `python3 tools/bf_exceptions.py`; `python3 tools/run_k3_bf.py 0` (40 min, 0.5 GB); `python3 tools/run_k3_lfd.py` (33 min, 2.7 GB); `python3 tools/summarize.py > data/summary_table.md; python3 tools/plot.py`.

## Task 1: construction and unit tests (`data/selftest_build_m3.json`)

* Cross-check against the repository (`csd.ic_bounds_m3.direct_generator`, which expands the same resultant formula, then `chain_nf`): on 4 instances ((p,k) = (97,2), (97,3), (257,4), (193,3)) the normalized coefficient dictionaries are identical and delta = 6 in all four [pcc:unverified].
* g vanishes exactly on Z: for the same 4 instances (k <= 4, 4096 points at k=4; zero counts 9, 7, 18, 3) the zero set of the NTT-evaluation of g agrees with the zero set of the repository's numeric `s4_value` (independent resultant formula) at every point of H^3: 0 mismatches over 64+512+4096+512 points [pcc:unverified].
* Unit instance has Z empty: p = 4013, k = 2, (a,b,t) = (1561, 2811, 887): the repository's `s4_value` has 0 zeros on H^3, and the NTT evaluation has none [pcc:unverified].
* The inverse is a true inverse: g * g^-1 = 1 in F_p[(Z/2^k)^3] verified by cyclic convolution on the first 3 units of every scan at k <= 3 (all `inverse_verified_by_convolution = true`) [pcc:unverified].
* delta = 6 on every one of the generic unit instances scanned (k = 2..8, 6 distinct primes ... see tables), delta = 5 on every instance of the locus 5t^3+3at+6b = 0, as stated in ic-bounds-m3 [pcc:unverified].

## Task 2 + 3: P4 first half, brute force vs delta + maxwt(supp g^-1)

Ordinary ring, k = 2 (6 variables, D up to 12, 18564 monomials). 27 instances in `data/bruteforce_k2.jsonl` (primes 17, 4013, 65537; 21 sampled generically, of which 1 at p = 17 fell on the locus by chance, plus 6 deliberately on the locus), plus 12 exceptional instances (`data/bruteforce_k2_exceptions.json`; maxwt(g^-1) = 5):

* sd by brute force equals delta + maxwt(supp g^-1) on **39 / 39** instances [pcc:unverified]. Values: sd = 12 = 6 + 6 on 19 generic instances (top coefficient nonzero); sd = 11 on 20 instances: 1 sampled instance with maxwt 5 (p = 17), 1 sampled instance that landed on the locus (delta 5 + maxwt 6), 6 deliberate locus instances (5 + 6) and 12 further exceptional instances (6 + 5). sd = 10 (delta 5, maxwt 5) occurs once in the scan but was not brute-forced.
* The profile (rank of Mac_D) is in each record (`sd_profile`); e.g. p = 4013: ranks D = 2..12: 6, 42, 153, 399, 861, 1659, 2961, 4983, 8001, 12375, 18564 (full), 1 first in Mac_12 [pcc:unverified].
* **k = 3 brute force was feasible and was run** (9 variables, D up to 15, 1,307,504 monomials, 2406 s, peak RSS 0.53 GB, packed-int sparse elimination), one unit instance (p = 65537, (a,b,t) = (17881, 16412, 56632), delta 6, maxwt(g^-1) = 9): ranks D = 10..15 = 92122, 167830, 293884, 497410, 817189, 1307504; 1 not in Mac_D for D <= 14, 1 in Mac_15. sd = 15 = delta + 3k = predicted (`data/bruteforce_k3_0.json/.log`) [pcc:unverified]. Only one instance at k = 3 (cost); not a statistical check.
* k = 4 brute force skipped: 12 variables, D up to 18 gives C(30,12) = 86.5 million monomials (about 100x the k=3 job in rows and, at the observed 4x per degree growth, many hours and far above 6 GB).

Verdict part 1 (P4 first half, brute force agreeing at k = 2): not falsified: 0 mismatches in 39 instances at k = 2 and 1 at k = 3 [pcc:unverified].

## Task 4: scan k = 2..8 (data/scan_*.jsonl, data/summary.json, data/summary_table.md)

Instances: random nonsingular y^2 = x^3 + ax + b, t the x-coordinate of an F_p point; draws until 20 unit instances (g has no zero on H^3), at least; delta+maxwt computed exactly from g^-1. "sd = delta+3k" means the coefficient of x^(N-1,N-1,N-1) in g^-1 is nonzero.

Primes: "small" tier = smallest p = 1 mod 2^k above 8^k/4 (k <= 5), 8^k/2 (k = 6), 4 * 8^k (k = 7, 8) and extra very small primes (13, 17, 41, 97, 113) to expose exceptions; "mid" tier 65537 (2^16+1, about 2^16 as asked). Deviation: at k = 7 and 8, p about 2^16 gives almost no unit instances (k = 7: 2 of 22, k = 8: 0 of 30, because |H^3| = 2^21 / 2^24 >> p), so the 20-unit rows at k = 7, 8 use p = 8390273 and 67109633 (about 4 * 8^k; 18 s per k = 8 instance, 368 s for 20 units). Small primes are sampled with replacement, so the table counts distinct (p,a,b,t) only.

| k | p | draws (distinct) | units | delta | sd = delta+3k | exceptions (sd = delta+3k-1) |
|---|---|---|---|---|---|---|
| 2 | 13 | 594 | 163 | 6 | 156 / 163 | 7 |
| 2 | 17 | 650 | 242 | 6 | 216 / 242 | 26 |
| 2 | 29 | 53 | 19 | 6 | 19 / 19 | 0 |
| 2 | 65537 | 20 | 20 | 6 | 20 / 20 | 0 |
| 3 | 17 | 2019 | 76 | 6 | 76 / 76 | 0 |
| 3 | 41 | 1851 | 284 | 6 | 275 / 284 | 9 |
| 3 | 137 | 36 | 20 | 6 | 20 / 20 | 0 |
| 3 | 65537 | 20 | 20 | 6 | 20 / 20 | 0 |
| 4 | 97 | 11261 | 296 | 6 | 296 / 296 | 0 |
| 4 | 113 | 8311 | 296 | 6 | 294 / 296 | 2 |
| 4 | 1153 | 32 | 20 | 6 | 20 / 20 | 0 |
| 4 | 65537 | 20 | 20 | 6 | 20 / 20 | 0 |
| 5 | 8353 | 39 | 20 | 6 | 20 / 20 | 0 |
| 5 | 65537 | 20 | 20 | 6 | 20 / 20 | 0 |
| 6 | 131713 | 23 | 20 | 6 | 20 / 20 | 0 |
| 6 | 65537 | 39 | 20 | 6 | 20 / 20 | 0 |
| 7 | 8390273 | 21 | 20 | 6 | 20 / 20 | 0 |
| 7 | 65537 | 22 | 2 | 6 | 2 / 2 | 0 |
| 8 | 67109633 | 20 | 20 | 6 | 20 / 20 | 0 |
| 8 | 65537 | 30 | 0 | - | - | - |

All claims in this table [pcc:unverified]. Findings:

* Total 45 distinct exceptions, all at primes <= 113. In every one maxwt(supp g^-1) = 3k - 1 exactly (top coefficient zero, the next weight class nonzero), so sd = delta + 3k - 1 = 3k + 5; never lower than that among the generic ones. Exceptional fraction: 4.3% (p = 13, k=2), 10.7% (p = 17, k=2), 3.2% (p = 41, k=3), 0.7% (p = 113, k=4), i.e. between about 0.6/p and 1.8/p; every row satisfies fraction >= 1 - 3/p (and hence >= 1 - 10/p) [pcc:unverified]. At primes >= 1153 no exception occurred in any of the k = 4..8 rows (20 units each, so this says only that the exception rate is below about 15% there; the only quantitative evidence about the rate is at the small primes).
* Exceptions found by this scan: none at p >= 1153, k = 2..8.
* The locus 5t^3 + 3at + 6b = 0 (delta = 5): at p = 65537, k = 2..6 and p = 8390273 (k = 7), 67109633 (k = 8): 6 units each (k = 8: 3), all with delta = 5 and maxwt = 3k, so sd = 3k + 5 (11, 14, 17, 20, 23, 26, 29), i.e. delta + 3k = 5 + 3k. Small-prime locus rows: one exception at p = 17 (k = 2, sd = 10 = 5 + 5). Locus brute force at k = 2: 7/7 agree with 11 [pcc:unverified].
* sd reaches delta + 3k = 3k + 6 at every depth k = 2..8 on generic unit instances, a linear growth with slope 3 = m, unbounded in k, unlike d_lf [pcc:unverified].

Wall times: k=2..5 scans < 1 s each; k=6: 3-4 s; k=7: 33 s (20 units), k=8: 368 s (20 units); extra small-prime scans <= 17 s each; locus k=8: 68 s. All within the 10-minute budget per (k, prime).

## Task 5: d_ff, d_lf, sd, gap

Measured (lfd_engine, ordinary ring, exact, unit instances): k = 2 on 27 instances (`bruteforce_k2.jsonl`): d_ff = d_lf = 7 on all generic instances (24), 6 on the 6 + 1 locus instances (delta 5). The ordinary-ring V_D profile for p = 4013: dim V_6 = 861 = rank Mac_6, dim V_7 = 1716 = all monomials (1 in V_7) while Mac_7 has rank 1659, so Mac_7 is a proper subspace of V_7 (sd = 12 > 7). k = 3 on one instance (p = 65537, above): d_ff = d_lf = 7, certified at D = 7 (`data/lfd_k3.json`, 1988 s, 2.7 GB); at D = 7 rank Mac_7 = 10948 < 11440 = dim V_7 [pcc:unverified].

Recorded (not re-measured, from `/home/user/ml-cryptanalysis/experiments/ic-bounds-m3/README.md` and `ic-bounds-m3-depth/README.md`): exact direct m = 3 LFD = 7 at k = 2, 3, 4 on all cases of that record, 7 at k = 5 on 12 instances; at k = 6 one non-unit instance (p = 7681, |Z| = 48) has V_7 behind the ideal, LFD(6) >= 8 there (so "LFD 7 for k >= 2" is exact only up to k = 5 and is NOT established for k = 6; it is not a unit instance). The carry lemma gives LFD >= 7. Not measured at k >= 4 by me.

| k | d_ff | d_lf | sd (generic unit) | sd - d_lf | 3k+1 (K1) | delta+3k |
|---|---|---|---|---|---|---|
| 2 | 7 (measured) | 7 (measured) | 12 (brute force + inverse) | 5 | 7 | 12 |
| 3 | 7 (measured, 1 instance) | 7 (measured, 1 instance) | 15 (brute force 1 instance + inverse on all scanned units) | 8 | 10 | 15 |
| 4 | n/m | 7 (recorded) | 18 (inverse) | 11 | 13 | 18 |
| 5 | n/m | 7 (recorded) | 21 (inverse) | 14 | 16 | 21 |
| 6 | n/m | >= 8 (recorded, non-unit instance) | 24 (inverse) | <= 16 | 19 | 24 |
| 7 | n/m | n/m | 27 (inverse) | - | 22 | 27 |
| 8 | n/m | n/m | 30 (inverse) | - | 25 | 30 |

Remarks [pcc:unverified]: (1) sd exceeds the kernel bound 3k+1 (K1) for every k >= 2 (12 > 7, 15 > 10, ... 30 > 25) and its excess over 3k+1 is the constant 5 = delta - 1. (2) sd - d_lf grows by 3 per level over k = 2..5 where d_lf is recorded; P3-like "sd - d_lf >= k" holds for m = 3 trivially (5, 8, 11, 14 >= k). (3) Always d_lf <= sd, as required: 7 <= 12, 15, ...; at k = 2 the Mac space is strictly inside V_7 (rank 1659 vs 1716), so the static space needs five more degrees than the iterated closure.

## Plot

`data/sd_vs_k_m3.svg`: sd (inverse prediction; max = 3k+6, exceptions as open diamonds at 3k+5; brute-force values as orange squares at k = 2, 3), d_lf (7 for k = 2..5 plus the k = 6 lower bound 8), 3k+1 (K1 bound) and delta + 3k.

## P4 verdict

* First half (brute force at k = 2 agrees with delta + maxwt(g^-1)): 39/39 at k = 2, 1/1 at k = 3: not falsified [pcc:unverified].
* Second half (sd = delta + 3k generically, delta = 6, and 5 on the locus): all generic unit instances at primes >= 1153 (k = 2..8, 20 or more units per k) have delta = 6 and sd = 6 + 3k; all locus instances at large primes have delta = 5 and sd = 5 + 3k; exceptions occur only at primes <= 113 and all have sd = delta + 3k - 1; the exceptional fraction is within 1 - 3/p everywhere. Not falsified (falsifier: a fraction below 1 - 10/p; minimum observed fraction 0.893 at p = 17 against 1 - 10/17 = 0.41 and 1 - 3/17 = 0.82) [pcc:unverified].
* P4 holds in the tested range; it is a finite test: p >= 1153 samples have only 20-296 units per row.

## Anomalies, caveats, skips

* delta = 6 never failed generically; delta = 5 exactly on the locus. At p = 17, k = 2 one sampled instance landed on the locus by chance (delta 5, sd 11).
* At k = 6 the recorded d_lf is only a lower bound (>= 8) on one non-unit instance; the "LFD 7 for all k >= 2" statement quoted in the task is exact only for k <= 5 in the repository records.
* k = 7, 8 at p about 2^16: no (k=8) or very few (k=7) unit instances since |H^3| >> p; larger primes (8.4e6, 6.7e7) used and listed.
* Skips: k = 4 ordinary-ring brute force (86.5 M monomials); d_ff / d_lf measured only at k = 2 (27 instances) and k = 3 (1 instance); k >= 4 d_lf cited from the repository record.
* Small primes sample (p,a,b,t) with replacement, so unit counts in `scan_k*_p{13,17,41,97,113}.jsonl` include repeats; `summary.json` and the table above deduplicate.
* Nothing was committed; all files are inside `laneB_prime_m3/`.
