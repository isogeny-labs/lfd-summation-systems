# ic-bounds-upper-redteam: independent checks for the review of U1 (T1: LFD <= k + 2)

Scripts written by the Red Team for `docs/reviews-U1-upper-bound-redteam-20261110.md`
(2026-11-10). They import nothing from `experiments/ic-bounds-upper/` (the producer's code).
The closure and exact-LFD engine is the independent validator's
`experiments/ic-curve-validation/code/lfd_engine.py`. It works in the ORDINARY polynomial ring:
all monomials of degree <= D in the chain variables, the literal closure rule, the ideal part
from evaluation ranks, and a stabilisation certificate (V_D = I_D and evaluation of the
monomials of degree < D onto K^Z) that does not use Theorem 1.

Run from this directory:

    ulimit -v 1500000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
    python3 rt_sanity.py        > out-sanity.txt       # S_3 is the summation polynomial; Lemma 2.2 by brute force; engine = literal Lean rule
    python3 rt_sweep_cross.py   > out-sweep-cross.txt  # builder + engine reproduce the certified sweep values, k <= 4
    python3 rt_run.py eq12                             # k = 1, 2: exhaustive planes, every coset ratio, special targets
    python3 rt_run.py eq3                              # k = 3
    python3 rt_run.py eq45                             # k = 4, 5
    python3 rt_run.py uneq17                           # unequal depths, exhaustive planes at p = 17, 13, 5
    python3 rt_run.py uneqdeep                         # short block of depth 1, 2, 3 against depth 5..12
    python3 rt_run.py uneqmid                          # unequal depths, both blocks >= 3
    python3 rt_run.py raw                              # unreduced degree-4 generator
    python3 rt_control.py       > out-control.txt      # negative control: random degree-2 generators
    python3 rt_univariate.py    > out-univariate.txt   # Lemma 2.3 as a statement: every q != 0 of degree <= 2
    python3 rt_singular.py      > out-singular.txt     # singular (a, b): outside T1
    python3 rt_p5.py            > out-p5.txt           # Proposition 5 example (k = 8, p = 257), D*(k), E2 arithmetic
    python3 rt_e1.py            > out-e1.txt           # P3 / E1 spot check in the ordinary ring
    python3 -I sweep_compare.py > out-sweep-compare.txt  # all certified two-summand prime sweep values against k + 2

`rt_run.py <job>` writes `rt-<job>.jsonl` (one row per instance; about 45 MB in total, ignored
by the local `.gitignore`) and `rt-<job>-summary.json` (kept here). A row is a violation if the
exact LFD exceeds min(k1, k2) + 2 (+4 for the raw generator); it is "uncertified" if the scan
ended before the stabilisation certificate. Both lists are in the summary. The test target was
fixed in `rt_run.py` before any run: any violation or any uncertified instance.

Files: `rt_lib.py` (system builder: chains, generator from the S_3 formula, coset scales, unequal
depths, raw option, zero set), `rt_run.py` (families), the single-purpose scripts above, the
`rt-*-summary.json` outputs and the `out-*.txt` console outputs.

Total run time on one shared core: about 1.5 hours. All scripts stay under 1.5 GB.
