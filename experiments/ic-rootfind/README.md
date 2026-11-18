# Point decomposition by root finding instead of a Gröbner basis (binary curves)

**Question.** The relation search of index calculus asks for every `R = P_1 + ... + P_m` with `x(P_i)` in a
subspace `V` of `F_(2^n)` (dimension `k`). The usual route is the Weil descent of `S_(m+1)` and F4/F5. Is there a
faster route that finds the roots directly?

**Answer.** Yes, when the factor base is `V = <1, w, ..., w^(k-1)>` (or any `V` whose product space `V2 = span(V·V)`
is small). For such `V` the two-summand problem is one F_2-linear system in `k + 1` unknowns and one quadratic, so
`m` summands cost `|FB|^(m-2)` small linear systems. Enumeration costs `|FB|^(m-1)` curve additions. For a random
`V` with `k(k+1)/2 >= n` the method gives nothing; the table below shows the cliff.

Code: `csd/ic_rootfind.py` (mathematics, tables, Python reference), `csd/_rootfind.c` (kernels, `n <= 63`),
`rust/ic_rootfind/` (Rust kernel, `n <= 191`, used for `n = 131`),
`tests/test_ic_rootfind.py`, driver `scripts/run_ic_rootfind.py`.

## The algorithm (LRF, "linearized root finding")

`P_1 + P_2 = T` iff `P_1, P_2, -T` lie on a line `y = λ(x + xT) + yT + xT`. Substituting into
`y^2 + xy = x^3 + a x^2 + b`:

    e1 = x_1 + x_2 = λ^2 + λ + a + xT          e2 = x_1 x_2 = xT λ^2 + ((yT + xT)^2 + b)/xT

Both are F_2-affine in `λ`. `x_1, x_2 ∈ V` forces `e1 ∈ V` and `e2 ∈ V2`. Writing `e1 = Σ t_i v_i` and
`λ = Q(e1 + a + xT) + ε`, with `Q` a fixed linear right inverse of `z ↦ z^2 + z`, the condition `e2 ∈ V2` plus one
trace condition is a linear system in `(t, ε)`: `k + 1` unknowns, `n - dim V2 + 1` equations. Each solution gives
`x^2 + e1 x + e2`; it is kept when its roots lie in `V`. The points come with it (`y_i = λ x_i + μ`), so only
rational decompositions with the right signs are produced.

For `m` summands, enumerate `m - 2` factor-base points and solve the two-summand problem for
`T = R - P_3 - ... - P_m`.

Per target the kernel does `k + 1` field multiplications, a `(k + 1)`-column elimination on 64-bit words, and about
six more multiplications per candidate; inversions are shared across blocks of 256 targets and 256 candidates
(Montgomery's trick). The enumeration baseline shares its inversions the same way.

The linear structure of `S_3` in the symmetric functions is not new (it is behind the symmetrized systems of
Faugère–Perret–Petit–Renault and Huang–Petit–Shinohara–Takagi, and Semaev's splitting into `S_3`'s). What is here
is the solver built on it and its measurement against the baselines.

## Checks

- `tests/test_ic_rootfind.py`: on 480 instances (`n` from 7 to 17, `m` = 2, 3, 4, standard and random `V`, random
  and Koblitz curves; 87 of them have decompositions) LRF, the C enumeration, the C meet in the middle and a
  plain Python enumeration on `csd.ic_dreg`'s group law return the same sets. The C field arithmetic and factor base
  are compared with `csd.ic_dreg`'s. The targets with `x(T) = 0` (the one case the line cannot take) are tested
  separately.
- Every row of every record below: all methods run on the instance return the same set, and every decomposition is
  re-added to `R` with `csd.ic_dreg`'s group law before the row is written.

## Results

Single thread, wall clock, minimum of three runs for the C kernels, median over three seeded instances per size
(random `a ∈ {0,1}`, random `b`, random `R`). Machine: 4-core Intel Xeon (PCLMULQDQ), gcc 13.3 `-O3 -march=native`.
Set-up of the tables is excluded (it is below 0.2 s everywhere).

### m = 3, V = <1, w, ..., w^(k-1)>, n = 3k - 2

| n | k | \|FB\| | decompositions (sum over seeds) | LRF | enumeration | speed-up | LRF per target |
|---|---|---|---|---|---|---|---|
| 13 | 5 | 31 | 2 | 50 µs | 59 µs | 1.3× | 1620 ns |
| 16 | 6 | 67 | 1 | 56 µs | 156 µs | 2.3× | 829 ns |
| 19 | 7 | 127 | 2 | 66 µs | 414 µs | 6.4× | 518 ns |
| 22 | 8 | 259 | 0 | 67 µs | 860 µs | 13× | 259 ns |
| 25 | 9 | 517 | 2 | 189 µs | 0.00693 s | 38× | 367 ns |
| 28 | 10 | 1015 | 0 | 185 µs | 0.0122 s | 66× | 182 ns |
| 31 | 11 | 2093 | 0 | 743 µs | 0.109 s | 145× | 355 ns |
| 34 | 12 | 4059 | 0 | 0.00106 s | 0.324 s | 293× | 261 ns |
| 37 | 13 | 8167 | 1 | 0.00301 s | — | — | 368 ns |
| 40 | 14 | 16565 | 1 | 0.0049 s | — | — | 296 ns |
| 43 | 15 | 32651 | 1 | 0.0134 s | — | — | 410 ns |
| 46 | 16 | 65461 | 1 | 0.0293 s | — | — | 447 ns |
| 49 | 17 | 130865 | 2 | 0.0579 s | — | — | 443 ns |
| 52 | 18 | 262397 | 2 | 0.166 s | — | — | 632 ns |
| 55 | 19 | 524657 | 0 | 0.264 s | — | — | 504 ns |
| 58 | 20 | 1048695 | 3 | 0.435 s | — | — | 415 ns |
| 61 | 21 | 2097669 | 0 | 1.18 s | — | — | 564 ns |

### m = 4, V = <1, w, ..., w^(k-1)>, n = 4k - 3

| n | k | decompositions | LRF | meet in the middle | its table | enumeration | LRF vs MITM | LRF vs enumeration |
|---|---|---|---|---|---|---|---|---|
| 13 | 4 | 3 | 49 µs | 112 µs | 0.0 MiB | 126 µs | 2.1× | 2.2× |
| 17 | 5 | 2 | 94 µs | 234 µs | 0.0 MiB | 287 µs | 3.4× | 4.1× |
| 21 | 6 | 4 | 340 µs | 0.00143 s | 0.0 MiB | 0.00263 s | 4.4× | 8.2× |
| 25 | 7 | 1 | 0.00129 s | 0.00612 s | 0.2 MiB | 0.0177 s | 4.7× | 14× |
| 29 | 8 | 0 | 0.00562 s | 0.0338 s | 0.8 MiB | 0.146 s | 6.1× | 26× |
| 33 | 9 | 1 | 0.0222 s | 0.153 s | 3.0 MiB | — | 7.0× | — |
| 37 | 10 | 0 | 0.0911 s | 0.66 s | 10.9 MiB | — | 7.3× | — |
| 41 | 11 | 3 | 0.46 s | 3.49 s | 47.1 MiB | — | 7.6× | — |

### The factor base decides: n = 31, k = 10, m = 3, V = <1, ..., w^(j-1)> + (10 - j) random vectors

| j | dim V2 | k + dim V2 - n | candidates per target | LRF | enumeration | speed-up |
|---|---|---|---|---|---|---|
| 10 | 19 | -2 | 0.25 | 234 µs | 0.0257 s | 110× |
| 9 | 27 | 6 | 64.18 | 0.00434 s | 0.0321 s | 6.3× |
| 8 | 31 | 10 | enumeration fallback | 0.0289 s | 0.0291 s | 1.0× |
| 7 | 31 | 10 | enumeration fallback | 0.0267 s | 0.0272 s | 1.0× |
| 6 | 31 | 10 | enumeration fallback | 0.0274 s | 0.0285 s | 1.0× |
| 5 | 31 | 10 | enumeration fallback | 0.0274 s | 0.0274 s | 1.0× |
| 4 | 31 | 10 | enumeration fallback | 0.0275 s | 0.0278 s | 1.0× |
| 3 | 31 | 10 | enumeration fallback | 0.026 s | 0.0264 s | 1.0× |
| 2 | 31 | 10 | enumeration fallback | 0.027 s | 0.0278 s | 1.0× |
| 0 | 31 | 10 | enumeration fallback | 0.0292 s | 0.029 s | 1.0× |

### Below and above n = 3k: k = 12, m = 3, V = <1, w, ..., w^11>

| n | unknowns | equations | targets with a solvable system | candidates per target | decompositions (sum) | LRF | enumeration | speed-up |
|---|---|---|---|---|---|---|---|---|
| 27 | 13 | 5 | 0.998 | 256 | 266 | 0.052 s | 0.358 s | 6.9× |
| 29 | 13 | 7 | 0.992 | 64.1 | 55 | 0.0135 s | 0.357 s | 26× |
| 31 | 13 | 9 | 0.968 | 16 | 11 | 0.00436 s | 0.385 s | 89× |
| 33 | 13 | 11 | 0.887 | 4.02 | 2 | 0.00175 s | 0.289 s | 165× |
| 35 | 13 | 13 | 0.615 | 1 | 0 | 0.00118 s | 0.358 s | 298× |
| 37 | 13 | 15 | 0.22 | 0.249 | 0 | 0.00101 s | 0.327 s | 322× |
| 39 | 13 | 17 | 0.0626 | 0.0648 | 0 | 984 µs | 0.364 s | 370× |
| 41 | 13 | 19 | 0.0167 | 0.0172 | 1 | 0.00105 s | 0.426 s | 409× |
| 43 | 13 | 21 | 0.00517 | 0.00517 | 0 | 916 µs | 0.32 s | 349× |
| 45 | 13 | 23 | 0.000763 | 0.000763 | 0 | 979 µs | 0.363 s | 371× |

### Against the repository's MutantXL on the Weil descent (m = 3)

| n | k | Boolean variables | Boolean solutions | decompositions of R | closing degree | MutantXL at that degree | MutantXL, all degrees | enumeration | LRF | LRF vs MutantXL |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 4 | 12 | [0, 0, 6] | [0, 0, 1] | [7] | 1.73 s | 1.73 s | 23 µs | 20 µs | 82,088× |
| 13 | 5 | 15 | [0, 3, 6] | [0, 1, 1] | [7] | 94.2 s | 94.2 s | 38 µs | 28 µs | 3,357,844× |

## Reading

- **Three summands.** LRF costs about 0.2 to 0.6 µs per factor-base point, flat in `k` up to `n = 61` (the word
  size of the kernel). Enumeration grows by 4 per step of `k`, LRF by 2: the speed-up doubles with each `k`
  (293× at `k = 12`, the last size where enumeration was run). At `n = 61, k = 21` one target `R` takes 1.18 s on one
  thread and 0.33 s on four; enumeration there is about 2^42 curve additions and was not run.
- **Four summands.** LRF beats the sort-based meet in the middle by 2× to 8× and stores nothing; both are
  `|FB|^2` up to the factor `O(k)`. The ratio grows with the table (cache misses), not with an exponent.
- **The factor base decides.** With one random vector in place of `w^9` the product space jumps from dimension
  19 to 27, the system has 6 free unknowns, and the speed-up falls from 110× to 6×. With two it is gone and the
  kernel falls back to enumeration. The quantity that matters is `k + dim V2 - n`: each unit doubles the candidates
  per target (the `k = 12` table shows 256, 64, 16, 4, 1 going from `n = 27` to `35`).
- **Against the algebraic route.** The repository's MutantXL on the Weil descent of `S_4` closes at degree 7 and
  takes 1.7 s at 12 Boolean variables and 94 s at 15. LRF and enumeration take tens of microseconds there (mostly
  the Python call). This baseline is the in-repo Python/C MutantXL, not Magma's F4 or an F5; no Gröbner
  package is installed in this environment. Its direction is not in doubt (at these sizes plain enumeration already
  wins by five orders of magnitude) but the ratio says nothing precise about F4. The Boolean system also has more
  solutions than `R` has decompositions (abscissas of points over the quadratic extension, sign patterns summing
  to `-R`); LRF returns only the decompositions.
- Sizes below `k = 8` are dominated by the Python call (about 20 to 50 µs); read ratios from `k >= 9`.

## n = 131: the ECC2K-130 curve (Rust kernel, multi-word field)

`rust/ic_rootfind/` is a second implementation of LRF for `m = 2, 3`, in Rust with no dependencies and no
Python in its execution path (field arithmetic, table construction, solver, self-test and benchmark are all in the
crate). Fields up to `n = 191` (three words, PCLMULQDQ). It differs from the C kernel in three ways:

- **No stored factor base.** The abscissas `x(P_3)` are walked by Gray-code index and lifted on the fly, so a range
  of a factor base far too large to list (`k = 44`, about 2^44 points) can be scanned and split across threads.
- **Row form of the linear system.** For `V = <1, w, ..., w^(k-1)>` and a modulus whose low terms have degree at
  most `2k - 1`, the annihilator of `V2` under the trace form is `δ·<1, w, ..., w^(n-2k)>` with `δ = 1/f'(w)`.
  Equation `j` is then the trace form against `z_j = w^j δ xT`: its whole coefficient row is one table lookup, and
  `z_(j+1) = w z_j` is a shift. The `k + 1` field multiplications per target of the column form become two. The
  crate checks the annihilator numerically at set-up and falls back to the column form when it does not hold.
- **Gauss–Jordan over all rows**, written so that the compiler vectorizes it.

Checks: `ic_rootfind selftest` (field arithmetic against a bit-by-bit product at 14 degrees from 7 to 191; 144
random instances against enumeration, 51 with decompositions, row and column forms compared on 42; planted
decompositions at `n` = 67, 71, 97, 131, 163 in both forms). `tests/test_ic_rootfind.py` runs the crate on the C
kernel's instances (`n` up to 61) and compares the decomposition sets.

Measurement (`rust/ic_rootfind/run_n131.sh`, record `n131.jsonl`): the curve `y^2 + xy = x^3 + 1` over
`F_2[x]/(x^131 + x^8 + x^3 + x^2 + 1)`, `k = 44`, so 45 unknowns and 45 equations. One decomposition
`R = P_1 + P_2 + P_3` is planted, a range of 2^22 abscissas containing a planted one is scanned, the planted
decomposition must be found and every row must re-add to `R`. Median of three seeds:

| form | threads | per target | targets with a solvable system | candidates per target |
|---|---|---|---|---|
| rows | 1 | 1850 ns | 0.61 | 1.00 |
| columns | 1 | 2757 ns | 0.61 | 1.00 |
| rows | 4 | 463 ns (wall) | 0.61 | 1.00 |

The enumeration inner loop on the same field (shared-inversion additions and membership tests, 2^24 of them)
takes 150 ns per addition.

Extrapolation to the whole search, **not measured**:

| per target `R` | work | single-thread time at the measured rate |
|---|---|---|
| LRF, rows | 2^44 abscissas, about 2^44 linear systems | about 1.0 core-years |
| enumeration | about 2^87 additions (unordered pairs) | about 7e+11 core-years |

The expected number of decompositions of a random `R` is about `|FB|^3 / (6 · 2^131) ≈ 1/3`, so one relation costs
about three such scans, and index calculus needs about 2^43 relations. LRF is about twelve orders of magnitude
below enumeration per target, and the relation stage would still be about 3·10^13 core-years. Pollard rho
on this curve needs about 2^61 group operations in total. Nothing here changes the security of ECC2K-130.

Cost per target of the Rust kernel at other sizes (`rust-scaling.jsonl`, one seed, ranges of at most 2^21 abscissas):

| n | k | per target | candidates per target |
|---|---|---|---|
| 34 | 12 | 1049 ns | 4.02 |
| 46 | 16 | 937 ns | 4.00 |
| 61 | 21 | 1096 ns | 2.00 |
| 97 | 33 | 1548 ns | 2.00 |
| 113 | 38 | 1581 ns | 1.00 |
| 131 | 44 | 1820 ns | 1.00 |
| 163 | 55 | 2417 ns | 2.00 |

At `n = 61, k = 21` the C kernel takes 564 ns per target (table above) and the Rust kernel about twice that: it
always computes on three words. The Rust kernel needs `n - dim V2 <= 62` (the system's rows are single words),
so at `n = 131` it needs `k >= 35`.

A stage ablation at `n = 131` (temporary early returns, not in the committed code) put the time per target at
about 320 ns for lifting and the target addition, 850 ns for building the 45 rows (45 lookups of 17 bytes each),
580 ns for the elimination and 260 ns for the candidates.

## Limits

- **No ECDLP speed claim.** The cost per target is `|FB|^(m-2)` linear systems, and the expected number of
  decompositions per target is about `|FB|^m / (m! 2^n)`. For `n = 131` (ECC2K-130) and `m = 3` this needs
  `k ≈ 44`: about 2^44 linear systems for about a third of a relation, and about 2^44 relations. That is far above
  Pollard rho (about 2^61 group operations in total). The per-target rate at `n = 131` is measured (section above);
  the whole search is extrapolated.
- **Koblitz structure is unused.** For prime `n` such as 131 no proper Frobenius-invariant subspace exists
  (`csd.ic_dreg.koblitz_invariant_dims`), so `V = <1, w, ..>` is not Frobenius-stable.
- **Binary curves only.** For the prime-field chain systems (`x^(2^k) = 1`) `S_3` has no such linear structure;
  nothing here applies.
- **Baseline quality.** Enumeration and meet in the middle share the field arithmetic and the shared-inversion
  trick with LRF. Meet in the middle uses `qsort` and binary search; a hash table would be faster by a constant.

## Open threads

- Three summands in fewer than `|FB|` linear systems. The Riemann–Roch form (`f = a(x) + y c(x)`, norm with
  roots in `V`) gives equations that are linear in squares plus bilinear terms `a_i c_j`; guessing less than a
  whole point and linearizing the rest has not been tried.
- Building the rows of the linear system by carry-less correlation of `δ xT` with fixed bit strings (about six
  `PCLMULQDQ` per column) instead of 45 table lookups of 17 bytes; the ablation says the rows are the largest stage.
- The related Rust work in `isogeny-labs/crypto` (pair-sum tables and enumeration: `src/cryptanalysis/ic_framework/rootfind.rs`,
  `RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md`) was found only after this solver was written. I did not find the
  linear two-summand step there (not an exhaustive search of that repository). Comparing LRF with its
  `rootfind-mitm` on the same instances, inside its harness and its matrix-F4, is the natural next step.
- F4 (Magma or msolve-style) on the same seeded instances, to replace the MutantXL stand-in.

## Commands

    python3 -m pytest -q tests/test_ic_rootfind.py
    python3 scripts/run_ic_rootfind.py all          # about 6 minutes, almost all of it the MutantXL baseline
    python3 experiments/ic-rootfind/summarize.py    # tables.md (the tables above)
    rust/ic_rootfind/run_n131.sh                    # Rust self-test, n131.jsonl, rust-scaling.jsonl (about 90 s)

Records: `scaling-m3.jsonl`, `m4.jsonl`, `basis-n31-k10.jsonl`, `regime-k12.jsonl`, `algebraic-m3.jsonl`, `n131.jsonl`, `rust-scaling.jsonl` (one row
per seeded instance, with the curve, the target, the timings and the kernel's counters).
