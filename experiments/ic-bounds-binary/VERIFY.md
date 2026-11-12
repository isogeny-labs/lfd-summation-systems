# Independent verification of item 6 (binary and Koblitz at larger size)

Date 2026-11-08. Verifier code is new and does not import the producer's modules or scripts, or any
closure engine in the repository (`csd.ic_bounds_binary`, `scripts/ic_bounds_binary_{run,validate,related,summarize}.py`,
`csd.ic_boolean_closure`, `csd.ic_local_certificate`, `csd/_elim.c`). From `csd.ic_dreg` it uses only
`binary_modulus`, for the single side check on the full-field family.

| file | what it does |
|---|---|
| `scripts/ic_bounds_binary_verify.py` | Rebuilds each system from the stored spec (modulus, b, basis, target). It uses its own GF(2^n) arithmetic, the value table of S_3 (or of S_4 = Res_X(S_3(x1,x2,X), S_3(x3,xR,X))) and a Moebius transform to the ANF. It recomputes the system hash in the stored format, reads Z off the value table and runs its own C closure kernel (compiled at run time). It reports dim V_D, dim V_D ∩ M_<D, the evaluation ranks of M_<D on Z, falls, D* and the exact LFD. A pure-Python saturating closure (`naive_dims`) cross-checks the C kernel. |
| `scripts/ic_bounds_binary_verify_static.py` | Static (non-iterated) Boolean Macaulay ranks and `d_full` |
| `scripts/ic_bounds_binary_verify_data.py` | Recomputes the headline counts from the stored jsonl |
| `scripts/ic_bounds_binary_verify_summary.py` | Summarizes the re-runs into `verify/summary.json` |
| `verify/` | `m2.jsonl`, `m3-direct.jsonl`, `m3-split.jsonl` (one row per case), `static.jsonl`, `data-checks.json`, `p4-check.json`, `binary-inverse-fullfield.json`, `summary.json`, `m2-k10.jsonl` (first timing probe) |

**Engine definition.** V_D is the least subspace of M_<=D (squarefree monomials) that contains the
inputs of Boolean degree <= D and is closed under h -> x_i h (normal form) for deg h <= D-1. Rows go
into a semi-echelon basis with columns in degree-descending order. Rows are never changed after
insertion, and every row whose leading monomial has degree <= D-1 is shifted by every variable once.
So on exit the low rows, which span V_D ∩ M_<D, have all their shifts in the span, and V_D is closed.
A fall at D means dim(V_D ∩ M_<D) > dim V_{D-1}. Certificates:
- proper ideal: the first D with codim V_D = |Z| and M_<D onto F_2^Z; LFD is then the last fall <= D*;
- unit ideal: the least D with 1 in V_D.

Every degree also asserts the soundness condition codim V_D >= rank(eval M_<=D on Z).

## Commands (all run with the resource prefix)

```sh
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VERIFY_BUILD=<scratch dir for the compiled kernel>
V=experiments/ic-bounds-binary/verify
python3 scripts/ic_bounds_binary_verify.py all --family m3-direct --k 2 3 --naive-max-r 9 --extra 1 --out $V/m3-direct.jsonl
python3 scripts/ic_bounds_binary_verify.py all --family m3-direct --k 4 5 --extra 1 --extra-max-r 12 --out $V/m3-direct.jsonl
python3 scripts/ic_bounds_binary_verify.py all --family m3-split --extra 1 --extra-max-r 15 --out $V/m3-split.jsonl
python3 scripts/ic_bounds_binary_verify.py all --family m2 --k 8 9 10 --out $V/m2.jsonl     # resumable, chunks of <= 24 min
python3 scripts/ic_bounds_binary_verify_static.py --family m3-direct --case-ids <14 ids in static.jsonl> --out $V/static.jsonl
python3 scripts/ic_bounds_binary_verify_static.py --family m2 --max-columns 15000 --case-ids <9 ids> --out $V/static.jsonl
python3 scripts/ic_bounds_binary_verify_data.py > $V/data-checks.json
python3 scripts/ic_bounds_binary_verify_summary.py > $V/summary.json
```

In `all` mode a case whose stored system hash was already certified in the same run is still
rebuilt and hash-checked, and its |Z| and point hash are recomputed. Only the closure is reused. Total
closure time was about 50 CPU-minutes. The peak RSS was 137 MB per process. The slowest case
(m2, k = 10, LFD 5) took 64 s.

## What was re-checked

**The full grid, not a subsample.** The engine was fast enough to run all 1871 cases (1726 distinct
systems, all certified). Compared with the stored values:

| family | cases | closures | hash | point set (sha) | stored dim V_D at every D <= D* | D* | LFD | falls | naive cross-check |
|---|---|---|---|---|---|---|---|---|---|
| m2 | 568 | 532 | 568/568 | 568/568 | 568/568 | 568/568 | 568/568 | 568/568 | n/a (r >= 16) |
| m3-direct | 695 | 634 | 695/695 | 695/695 | 695/695 | 695/695 | 695/695 | 695/695 | 393/393 (all r = 6, 9 closures) |
| m3-split | 608 | 560 | 608/608 | 608/608 | 608/608 | 608/608 | 608/608 | 608/608 | n/a |

Beyond D*, one extra degree was scanned for 880 closures (m3-direct r <= 12, m3-split r <= 15). No
fall appears after the certificate degree, which agrees with the local-to-global lemma.

Recomputed histograms (exact, from my own certificates; they agree with README and tables.md):

| family | k | LFD histogram | max LFD | D* - LFD |
|---|---|---|---|---|
| m2 | 8 | 2:3, 3:70, 4:119 | 4 | 0:190, 1:2 |
| m2 | 9 | 2:16, 3:30, 4:146 | 4 | 0:190, 1:2 |
| m2 | 10 | 2:17, 3:4, 4:133, 5:30 | 5 | 0:180, 1:4 |
| m3-direct | 2 | 3:43, 4:12, 5:42, 6:104 | 6 | 0:197, 1:4 |
| m3-direct | 3 | 4:54, 5:6, 6:30, 7:141 | 7 | 0:229, 1:2 |
| m3-direct | 4 | 4:23, 5:41, 6:11, 7:180 | 7 | 0:232, 1:23 |
| m3-direct | 5 | 5:2, 7:6 | 7 | 0:6, 1:2 |
| m3-split | 2 | 3:80, 4:121 | 4 | 0:194, 1:7 |
| m3-split | 3 | 3:14, 4:213, 5:4 | 5 | 0:230, 1:1 |
| m3-split | 4 | 4:166, 5:10 | 5 | 0:169, 1:7 |

D* - LFD totals 1817 zeros and 54 ones. The maximum D* is 5 (m2, m3-split) and 7 (m3-direct).

**Static Macaulay (P5 data).** 23 rows were recomputed with my own static ranks. They are the
largest-`d_full` and random rows of m3-direct k = 2, 3, 4, the six m2 rows with a finite `d_full`,
and three capped m2 k = 8 rows. `d_full` and the rank at every stored degree agree 23/23. The
m3-direct k = 4 rows reach `d_full = 12 = r` against D* = 7, and the capped m2 rows do not meet the
certificate up to the 15,000-column cap. For unit ideals the static criterion is "1 in the row
space" (with that reading `d_full = 2` for five m2 unit rows). For D > r the rows keep growing while
the columns stay fixed, which is why `d_full` = 9 > r = 6 occurs at m3-direct k = 2.

**README commands.** These were run in a scratch copy of the tree, so no repository file was rewritten:
- `ic_bounds_binary_summarize.py` regenerates `summary.json` and `tables.md` byte-identically.
- `ic_bounds_binary_run.py --family m3-direct --n 5 --k 2 --max-degree 14` gives 28/28 rows, and
  `--family m2 --n 14 --k 8 --limit 8` gives 8/8 rows. Both match the stored hash, LFD, D*, ranks,
  falls and point hash.
- `--family m2 --n 11 --k 10 --list` against the stored file reports "24 cases to run", which
  matches the stated gap.
- `ic_bounds_binary_validate.py --hash-all` on a 90-row subset rehashed 90 rows, and engine A
  agreed 4/4.

The stored validation files contain what the claim says: 1871 rehashed, engine A 20/20, and MXL
7 + 4 + 6 = 17/17. Note that 3 of the 11 "m3" MXL rows are m3-split rows, not m3-direct.

**Data counts.** The following were recomputed from the jsonl and all agree:
- HKYY >= r: 1041/1355 (m2 403/525, m3-direct 406/449, m3-split 232/381), including every proper m2 row at n = 11;
- P1 bound = D* in 144 proper rows, and never above D*;
- paired order-two direct vs split: 54/103/4 of 161;
- largest certificate matrices: 21,700 columns x 16,212 rows (m2), 16,384 columns (m3-direct);
- mean |Z| against 2^(2k-n);
- unordered/|Z| = 0.449 to 0.453.

**P4.** Computed in `verify/p4-check.json`. F_(2^11) and F_(2^13) each have exactly one
Frobenius-invariant hyperplane, and it equals {z^2 + z}. The order of 2 is 10 mod 11 and 12 mod 13,
and a normal element exists, so the invariant subspaces have dimensions 0, 1, n-1 and n only.

**P3.** All 121 Koblitz (0,1)/(1,1) order-two pairs have equal system hashes. On 8 sampled m2
order-two rows (random and Koblitz b), Z = {x1 x2 = sqrt b} holds exactly.

## Problems found

1. **"First certified binary m = 2 systems with LFD 5 in the repository" is false.**
   `experiments/ic-binary-inverse` and `docs/IC_BINARY_INVERSE_FAMILY.md` already contain a
   binary/Koblitz two-summand system with certified exact LFD 5: the full-field base k = n = 8, curve
   `y^2+xy=x^3+ax^2+1`, target (0,1), in original coordinates. That family also has model LFD 6 and
   7 at n = 10, 12. My engine independently gives LFD 5 for that n = 8 system (a = 0 and a = 1), with
   |Z| = 255, dims [0,0,8,472,2260,6630], falls [3,5] and D* = 5
   (`verify/binary-inverse-fullfield.json`). What holds is narrower: these are the first LFD-5 rows
   among proper-subspace bases (k < n). The k <= 7 sweep in `experiments/ic-curve-sweep` does have a
   maximum of 4.
2. **Systems and cases are conflated in the LFD-5 counts.** There are 30 LFD-5 *cases* but 28
   distinct *systems*. At n = 11, k = 10 there are 20 proper cases but 18 systems, because two
   Koblitz order-two pairs coincide. The LFD <= 4 side is 137 cases and 127 systems. The separation
   itself holds on systems: LFD 5 for |Z| in [472, 1056], LFD <= 4 for |Z| in [46, 352].
3. **The m3-direct maximum does not exceed k+2 at every k.** The claim says it does, but at k = 5
   the maximum LFD is 7 = k+2.
4. **The P3 explanation of the duplicates is incomplete.** P3 says the Koblitz coincidence "is why
   1871 cases are 1726 systems". It explains 117 of the 139 duplicate groups (121 of the coinciding
   pairs). The remaining groups come from targets with the same x-coordinate on one curve (random and
   planted, R and -R) and from three 4-way groups of m3-direct k = 2 F_4-subfield bases, whose
   systems coincide across n = 6/8 and n = 6/9.
5. Minor wording. P2's algorithmic remark needs B chosen without Z. That is possible: take the
   standard monomials of the degree-descending echelon form of V_D*, which lie in M_<D* under the
   certificate. The written proof picks B through Z but stays valid as an existence argument.
   The README's "11/11 m3-direct" MXL line includes 3 m3-split rows.

No numerical value labelled exact failed to reproduce. No stall is labelled as a bound: every value
is certified and the caps are stated as caps. The m = 3 k = 5 and n = 11, k = 10 coverage gaps are
disclosed.

## Verdicts

| # | claim | verdict |
|---|---|---|
| 1 | all 1871 exact (1355 local-to-global, 516 unit), none capped | confirmed (all 1871 recertified independently) |
| 2 | m2 max LFD 4/4/5 at k = 8/9/10, far below k+2 and 2k+1 | confirmed |
| 3 | first certified binary m = 2 LFD 5 in the repository; 30 systems | refuted (prior full-field n = 8 LFD 5; 28 systems / 30 cases); the LFD values themselves confirmed |
| 4 | k = 10 separation by |Z| (empirical) | confirmed as an empirical pattern in this grid (counts are cases, not systems; n = 11 grid half done) |
| 5 | mean |Z| tracks 2^(2k-n), unordered ~0.45|Z| (empirical) | confirmed |
| 6 | m3-direct max LFD 6, 7, 7, 7; exceeds k+2 | overstated (values exact and confirmed; 7 = k+2 at k = 5, so it does not exceed there) |
| 7 | m3-direct maximum stays at 7 to k = 5, thin k = 5 evidence (empirical) | confirmed (8/8 rows recertified; correctly labelled as a pattern) |
| 8 | m3-split max 4, 5, 5; paired 54/103/4 of 161 | confirmed (k = 4 coverage partial, disclosed) |
| 9 | P1 counting bound; equals D* in 144 rows | confirmed |
| 10 | P2: Z read off V_D* by |Z|-size multiplication matrices | confirmed (proof valid; minor remark on choosing B) |
| 11 | P3: S_3, S_4 depend only on b; order-two Z | confirmed (the duplicate-count explanation is incomplete) |
| 12 | P4 invariant subspaces of F_(2^11), F_(2^13) | confirmed |
| 13 | P5: D* <= d_full; d_full = 12 at m3-direct k = 4; m2 111/117 not met below cap | confirmed (proof valid; 23/23 static rows reproduced) |
| 14 | D* - LFD 0 in 1817, 1 in 54; D* maxima; Macaulay sizes | confirmed |
| 15 | HKYY degree >= r in 1041/1355 | confirmed |
| 16 | independent validation (1871 hashes, engine A 20/20, MXL 17/17) | confirmed (stored files match; my own builder and engine also agree on all 1871) |

Evidence labels: the recertification of every exact value is **exact**. It uses an independent
engine with the codim = |Z| plus evaluation-onto certificate, or 1 in V_D. P1 to P5 are **proof**
(the arguments were checked as written). The |Z| patterns are **empirical**. No Lean was run.
