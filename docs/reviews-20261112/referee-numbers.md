# Number audit: `paper/curve-lfd-bounds.tex` against the record files

Repository: `/home/user/lfd-summation-systems` (single commit `705279f`, 2026-11-12; shallow clone, no history).
Nothing in the repository was modified. Scripts and their output are in this directory:

- `c1.py` (claim 1, sweep), `verify_all.py` → `verify_all.out` (claims 2–14), `c6_rank.py` (evaluation rank at k = 14).
- All recomputations read the JSON/JSONL/CSV records directly; "paper line" is the line in `curve-lfd-bounds.tex`.

Verdicts: **MATCH**, **MISMATCH** (both values given), **NO RECORD** (where I looked), **PROSE ONLY** (the number
appears in a Markdown file but in no data record).

## 1. Result table

| # | Claim | Paper line | Record | Paper | Record | Verdict |
|---|---|---|---|---|---|---|
| 1 | cases / distinct systems / capped | `curve-sweep-data.tex` | `ic-curve-sweep/cases.jsonl`, `results.jsonl` | 3163 / 3005 / 0 | 3163 / 3005 / 0 | MATCH |
| 1 | prime systems at k = 2, 3 | 399 | same | 455 | 220 + 235 = 455 | MATCH |
| 1 | binary + Koblitz systems with ≤ 8 Boolean variables | 399 | `systems.jsonl` (`variables`) | 985 | 985 | MATCH |
| 1 | k = 2 systems with LFD 3 | 928–930 | same | 3, all 3t²+a = 0; two unit, one proper with 2 solutions | 3 (all random-coset; p = 97, 97, 193); 3t²+a = 0 in all; unit, unit, proper with 2 points | MATCH |
| 1 | value-4 systems at k = 6, 7 | 930–932 | same | 21 (12 + 9), all y² = x³+x, target (0,0) | 12 + 9; every case a = 1, b = 0, target (0,0); every such system has value 4 | MATCH |
| 1 | dim V_4 deviations at k = 6, 7 | 1032–1037 | `results.jsonl` (`quotient_ranks[4]`) | 86 of 333; values 607, 835, 730, 794, 1343, 1471, 603, 857 | 86 of 333; 607×35, 835×28, 730×4, 794×8, 1343×5, 1471×4, 603×1, 857×1 | MATCH |
| 1 | dim V_3 ≠ 2k+1 | 1024–1025 | `quotient_ranks[3]` | three systems at k = 2, one at k = 6 | 3 at k = 2 (15, 15, 13), 1 at k = 6 (34; p = 193, 3t²+a = 0) | MATCH |
| 2 | validator agreement | 399–402 | `ic-curve-validation/all_compare.csv` | 1440, no mismatch | 1440 `required`, all `match_all` = True (455 prime, 535 + 450 Boolean) | MATCH |
| 2 | further systems at k = 4, 5, 6 | 405–406 | same (`group = bonus`) | 14, two with LFD 5 and falls 4, 5 | 14 (k = 4: 6, k = 5: 4, k = 6: 4); two with LFD 5, falls [4, 5] | MATCH |
| 3 | depth 8 | 902–905, 937–938 | `ic-chain-depth8/cases.jsonl`, `results.jsonl` | 20 systems; 12 subgroup, 8 coset; 19 proper, 1 unit; LFD 5, falls 4, 5 | 20; 12 / 8; 19 / 1; all (5, [4, 5]); primes 257 (4), 769 (8), 7681 (8) | MATCH |
| 3 | k = 7 controls | 904, 941 | depth8 and depth9 `cases.jsonl` | two + two, y² = x³+x, p = 7681, LFD 5 | 2 + 2, all (7681, 1, 0), LFD 5, falls [4, 5] | MATCH |
| 3 | depth 9 | 904, 939–940 | `ic-chain-depth9/` | 8 systems, 4 + 4, 24 to 44 solutions | 8; 4 subgroup, 4 coset; points 24, 28, 32, 32, 40, 40, 44, 44 | MATCH |
| 3 | "same curves" at depth 9 | 904 | same | same curves as depth 8 | random curve is (6068, 5176) at k = 9, (6667, 4191) at k = 8 over 7681 | minor MISMATCH (wording) |
| 4 | first depth-10 instance | 1136–1142 | `ic-chain-depth10/probe.json`, `rerun.json` | abscissa 9362; dim V_4 1549; rank 18,569 after 17 rounds; 3163 low rows | 9362; 1549; 18569, 17 rounds; 3163 | MATCH |
| 4 | second / third instance | 19–25 | `proper-13313.json`, `random-15361.json`; `ic-bounds-d5/runs/k10-*` | 7412, 98 sol.; 7776, 60 sol.; rank 18,569; targets 21,602 and 21,640; parts 9311 + 9258 | 7412, 98; 7776, 60; 18569; 21602, 21640; [9311, 9258] on all three | MATCH |
| 4 | degree 6 | 36–41 | `ic-chain-depth10/degree6.json` | 8 rounds; 30,318 + 30,142; 59,371 after round 7; 18,569 + 2042 low rows | 8; 30318, 30142; 29776 + 29595 = 59371; 20611 = 18569 + 2042 | MATCH |
| 4 | independent engine, D = 6 | 42–45 | `ic-bounds-d5/verify/k10-p12289-s10-D6.json` | 453 rounds; 21,700 low rows, no residual | 453; rows_checked 21700, 0 residuals | MATCH |
| 5 | k = 11 | 59–61 | `ic-bounds-d5/verify/k11-*.json` | 25,033 on two instances; proper one has 348 solutions | 25033 on both; points 0 and 348 | MATCH |
| 5 | k = 12 | 65–68 | `verify/k12-*.json` | 32,057; 1352 and 1399; (0,0,2,146,3727,32057); 16,061 + 15,996 | identical | MATCH |
| 5 | k = 13 | 68–71 | `verify/k13-p40961-*.json` | 39,657; 1484 and 1582; 19,864 + 19,793 | identical; profile (0,0,2,150,4021,39657) | MATCH |
| 5 | dim V_4 = 4k²+172k−571, k = 6..13, two instances | 64 | `ic-bounds-d5/verify/d4.json`, `d4-check.json` | holds | 16 rows (2 per k), all equal the formula | MATCH |
| 5 | dim(V_4 ∩ M_≤3) = 4k+87 | 1031 | same (`profile[3]`) | holds | 111, 115, …, 139 at k = 6..13 | MATCH |
| 5 | reduction "checked exactly at k = 7..12 on three instances" | 1047 | `ic-bounds-fall/verify/free-p40961.json`, `k11/`, `k12/`, `cross.json`; README §8 | k = 7..12 on three instances | k = 8..12 on the primary instance only; the two cross-check instances have k = 8, 9 only | MISMATCH |
| 5 | depth-step curves (five curves, p ∈ {257, 7681, 12289}) | 1019–1022 | `experiments/ic-depth-step/growth.json` named by `scripts/run_ic_depth_step.py` | five curves, k = 3..10 | directory is not in the repository | NO RECORD |
| 6 | dim V_5 and profiles, unit instance, k = 10..15 | 173–179 | `ic-bounds-k14/runs/unit-p65537-k{10..15}-D5.json` | 18,569 … 56,649 with the six profiles | identical | MATCH |
| 6 | dim V_6, k = 10..13 | 182–183 | `runs/unit-p65537-k{10..13}-D6.json` | 60,460 / 110,056 / 190,051 / 313,912, full | identical, `contains_one` true | MATCH |
| 6 | random curve p = 65537 | 273–279 | `runs/proper-p65537-k13-D{5,6}.json`, `instances/proper-p65537-s12-k1{3,4}.json` | target (1424, 826); 1028 and 4068 solutions; V_6(13) = 312,884; 82,654 | (1424, 826); 1028, 4068; 312884; profile[5] = 82654 | MATCH |
| 6 | evaluation rank 4068 at k = 14 (hence 118,370, 495,110, 46,188) | 276, 284–286 | `RESULTS.md` line 54 (prose); the k = 14 instance JSON has no `eval_rank` field | rank 4068 | no stored field; recomputed here: rank 4068 already on 4500 random monomials of weight ≤ 4 | MATCH by recomputation; record is prose only |
| 6 | p = 114689 | 302–305 | `runs/proper-p114689-k13-D{5,6}.json`, instances | 594 and 2370 solutions; 313,318; 83,088 | 594, 2370; 313318; profile[5] = 83088 | MATCH |
| 6 | p = 163841 | 291–300 | `runs/unit-p163841-k1{3,4}-D5.json`, `-k13-D6.json` | 39,657; V_6(13) full; 47,849 | identical; target (73730, 28291), 0 points | MATCH |
| 6 | k = 15 and k = 16, D = 5 | 319–321 | `runs/proper-p65537-k1{5,6}-D5.json`, `unit-p65537-k15-D5.json` | 56,649; 66,073, (0,0,2,162,4951,66073) | identical | MATCH |
| 6 | "eighteen closures … four execution modes" | 137–138 | `validation/small18.json` | 18; k = 6, 7; D = 3, 4, 5; 4 modes | 18 rows (3 instances × 2 × 3), 4 modes, all match | MATCH |
| 6 | "ten cases over p = 65537 at k = 6, 7, 8" | 139 | `validation/p65537.json` | 10 | 10 (three modes each) | MATCH |
| 6 | GPU self-test "on eight cases" | 141–142 | `runs/selftest-gpu-h200*.json` | 8 | 8 in the first self-test; 10 in the three later ones (job B, relaunch, spill engine) | MATCH for the first; later ones differ |
| 6 | reference engine reproduces V_5(14), V_6(10), V_6(11) | 154–157 | `runs/reference-engine-unit-p65537-*.json` | hashes equal, 47,849 in 122,438 columns | hashes equal to the GPU records in all three; 47849 / 122438 | MATCH |
| 7 | predictions committed before the runs | 376–380 | `ic-bounds-freeext/predictions-2026-11-10.json`, `outcomes-2026-11-10.json` | committed before the runs | self-attested only (section 2, item B) | NOT VERIFIABLE in the repo |
| 7 | "three predictions"; value 7 "was not among them" | 376–384 | `predictions-2026-11-10.json` | three; 7 at k = 14 not predicted | six entries; `lfd-table` has `"14": 7` and `d6-falls-at-14` says "(LFD(14) = 7)" | MISMATCH |
| 7 | free extension of V_5(10) to depth 11 not closed | 339–341 | `ic-bounds-freeext/records/D5-*.json` | four instances | five records, `W_closed` false at k0 = 10 in all five | MISMATCH (4 vs 5) |
| 7 | free extension of V_5(11) to depth 12 | 342–344 | same | closed on four; same hash as V_5(12) on two | closed, g in, 1 not in on five; hash equal on three (two over 12289, one over 40961) | MISMATCH (4 / 2 vs 5 / 3) |
| 8 | lift validation | 231–233 | `ic-bounds-k14/validation/lift-*.json` (10 sets named in `LIFT.md` line 287) | 170 runs, 37 cases, D ≤ 5, k ≤ 11, four primes | 170 rows, 37 distinct cases, primes 12289, 15361, 65537, 16777213, max k 11, max D 5, all `ok` | MATCH |
| 8 | `\LIFTCONTROL` | 435–440 | `validation/lift-k10d5-proper.json` | 14 solutions; 60,446; no 1; K codimension 7 per part | 14; 60446; false; dim T − dim K = 7, 7; all comparisons true | MATCH |
| 9 | k = 3, 4 | 1157–1158 | `ic-bounds-m3/rows.jsonl` | 234 cases, all LFD 7 | 117 + 117, all 7 | MATCH |
| 9 | k = 2 | 1159–1162 | same | 105 of 116 are 7; exceptions have degree 5 exactly when 5t³+3at+6b = 0 | 105 of 116; the 11 others have LFD 6, degree 5 and vanishing coefficient; no other case does | MATCH |
| 9 | recomputed values | 1161 | `ic-bounds-m3/verify/recheck.jsonl` | 1092 | 1092 rows | MATCH |
| 9 | split presentation | 1185–1187 | `rows.jsonl` (`split_free_lfd`) | 6 at k = 4; 111 of 117 at k = 3, 5 in six | 117 of 117; 111 and 6 | MATCH |
| 9 | twelve further instances | 1189–1195 | `ic-bounds-m3-depth/closures.jsonl`, `summary.json` | LFD 7 at k = 2..5; eight primes; \|Z\| 12 to 544 at k = 5 | 12 instances × 4 depths, all 7; 8 primes; \|Z\| 12 … 544; V_6 rank 1 | MATCH |
| 9 | k = 6 instance | 1202–1208 | `closures.jsonl` (p7681-r, k = 6) | p = 7681, (1876, 2678), t = 6223, \|Z\| = 48, 19,958, (0,0,0,0,0,3,2375), 62,956 | identical | MATCH |
| 10 | Table `tab:msum` | 1251–1256 | `ic-bounds-msum/results.jsonl` | 13 cells | all 13 cells equal | MATCH |
| 10 | 409 of 412; three exceptions | 1276–1282 | same | 409 / 412; unit ideals on y² = x³+x, deg g = mk−1 | 409 / 412; (4,3), (7,2), (8,2); a = 1, b = 0, zeros 0, degree mk−1 | MATCH |
| 10 | "96 instances", "56 instances" | 1284–1286 | same | 96; 56 | 48 + 48; 48 + 8 | MATCH |
| 10 | "(eight with 2^20 or 2^21 grid points, two with 2^24)" | 1289–1290 | same (`points`) | 8; 2 | 2^20: 8, 2^21: 8 (16 together); 2^24: 4 | MISMATCH |
| 11 | k = 6, p = 193 targets | 1432–1437 | `ic-bounds-locus/data/`, `heldout/` | 1276 targets, 13 curves; 12 with value 4 (3 + 9) | 1276, 13; 12 (3 with t = 0, ab = 0; 9 others, all in the 190-point class) | MATCH |
| 11 | C_2 | 1435–1436 | `c-fit-p193.json` | 190 points, no curve of degree ≤ 18 | count 190, least degree 19 | MATCH |
| 11 | stratified sample | 1431–1432 | `verify/sample.json` | 96 | 44 deviating + 52 generic = 96 | MATCH |
| 11 | 730 against 794 | 1449 | `data/p193-k6-a1-b0.jsonl`, `a5-b0.jsonl` (t = 0) | 730, 794 | 730 (a = 1), 794 (a = 5) | MATCH |
| 11 | denominator degree ≥ 480; 52 CPU-hours | 1445–1446, 1457–1458 | `docs/proofs/N6-genericity.md` lines 40, 69, 409, 547; data directory `experiments/ic-bounds-genericity/` is not in the repository | 480; 52 | prose says 480 and "about 52 CPU-hours" (an estimate) | PROSE ONLY |
| 12 | selector-bound test | 851, 856, 859 | `docs/reviews-U1-upper-bound-redteam-20261110.md` (table at lines 276–291, line 77); `experiments/ic-bounds-upper-redteam/` and `ic-bounds-upper/` are not in the repository | 149,122; 132,303; 2570 | table rows sum to 132,303 and 149,122; 2,570 in line 77 | PROSE ONLY, no record file |
| 13 | value 5 on the first grid | 1508–1509 | `ic-bounds-binary/results-m2.jsonl` | 28 systems, all k = 10, exactly those with > 470 solutions | 28 distinct systems (30 cases), all k = 10; min 472; largest count below value 5 is 352 | MATCH |
| 13 | new systems | 1514 | `ic-bounds-binary-law/results-prereg.jsonl`, `score.json` | 80 new systems | 80 cases, 75 distinct systems | minor MISMATCH |
| 13 | law scores | 1516–1517 | `score.json` | 11 / 36 and 36 / 11 on 47 | A: 11 hit, 36 miss; B: 36 hit, 11 miss; 47 | MATCH |
| 13 | trace-zero hyperplane at k = 9 | 1518–1519 | `score.json` (cell P1, affine) | 11 of 11 with value 5, about 490 solutions | 11 of 11; \|Z\| = 482 (6), 488 (2), 499 (3) | MATCH |
| 13 | subfield F_{2^8} ⊂ F_{2^16} | 1520–1522 | `subfield-check.jsonl` | 255 solutions, value 5, k = 8 | 255, 5 (both curves) | MATCH |
| 13 | 1678 of 1775 | 1525 | `laws.json` (`precertificate_profile_vs_group_mode`) | depend only on (n, k) and target class | 1678 / 1775, but grouped by (n, k, r_eff, m_eff, class) | MATCH on the count; grouping differs |
| 14 | `\LEANDECLS` | 1585 | `ic-lean/verification.json`, `axioms.log` | 513 | 513 declarations, 513 log lines, axioms within the three allowed | MATCH |
| 14 | 3220 memberships | 1336 | `ic-bounds-cert/*/certificate.json` (`shift_checks`) | 3220 | 3220 (both instances); 829 rows, 1471 columns | MATCH |
| 14 | "43 `decide +kernel` chunks" | 1337–1338 | `formal/ICFallDegrees/N5CertKernelP40961.lean`, `…P11.lean` | 43 | 42 theorems closed by `decide +kernel` (29 + 9 + 4); `grep -c` gives 43 only because of the docstring; the cert README says 43 | MISMATCH (43 vs 42) |
| 14 | 83 rows, 39,126 terms | 1372 | `ic-bounds-exact-k6/…/summary.json` | 83; 39,126 | 83; sum of `terms_per_row` = 39126 | MATCH |
| 14 | 1332 shifts of 111 low rows | 1376–1377 | same | 1332; 111 | 111 low rows; 1333 member tests = 1332 + generator | MATCH |
| 14 | "125 `decide +kernel` theorems" | 1378–1379 | `ExactK6KernelLowerP1153.lean`, `ExactK6KernelUpperP1153.lean` | 125 | 40 + 85 = 125 | MATCH |
| 14 | 131 candidates; deg h = (46, 46); "about 150 monomials" | 1347–1349 | `ic-bounds-cert/p40961-s7/witness.json` | 131; (46, 46); about 150 | 131 admissible, 0 outside; [46, 46]; product has 161 monomials | MATCH (161 for "about 150") |
| 15 | dim V_6(14) = 448,922 "on three instances" | 1472 | lines 314, 355, Table `tab:profile` row 14; four `*-k14-D6.json` records | three | four everywhere else, and four records | MISMATCH (internal) |
| 15 | `tab:twosided`, binary lower bound | 801 | `paper/binary-lfd.tex` Theorem `thm:inverse`; `ic-bounds-binary-law/table.jsonl` (inverse rows) | ≥ ⌊k/2⌋+1 for k ≥ 6 | exact value 3 for odd n; rows k = 7, 9, 11 have LFD 3 | MISMATCH (internal) |
| 15 | `tab:twosided` and abstract, binary "≤ 5 through k = 10" | 477, 557, 801–802, 1507, 1613 | same | at most 5 for dimension up to 10 | inverse family at k = n = 10 has LFD 6 (theorem; `table.jsonl` row n = k = 10) | MISMATCH unless restricted to the seeded grids |
| 15 | `tab:profile` rows 2–15 against the body | 979–988 | body text | — | consistent, except row 13 "as k = 11" (which names p = 12289; the k = 13 instances are over 40961) | MATCH, one wording slip |

## 2. Mismatches and unbacked claims, by priority

**A. Selector-bound test has no record (claim 12).** 149,122 systems, 132,303 at k = 2 and 2570 singular systems
appear only in the prose of `docs/reviews-U1-upper-bound-redteam-20261110.md`. Its table rows do sum to the paper's
figures (121,060 + 11,243 = 132,303; adding 9,248 + 5,847 + 1,580 + 144 gives 149,122). The directory it names,
`experiments/ic-bounds-upper-redteam/`, is absent, and the review states that the per-instance rows "were not kept".
Nothing can be recomputed.

**B. Ordering of the predictions is self-attested (claim 7).** Evidence inside the repository:
- `predictions-2026-11-10.json` has `frozen_on: 2026-11-10`; its SHA-256 (`b612cbbc…`) equals the hash stored in the
  outcomes file.
- `outcomes-2026-11-10.json` says the predictions were "first committed in WIP snapshot 3f2c93d2 at 2026-11-10
  03:06:50 UTC" and the run records at 06:59 UTC and later, and adds "not an external timestamp".
- The clone has one commit (`705279f`, 2026-11-12, `.git/shallow` present), so `3f2c93d2` cannot be inspected. The run
  JSONs carry no timestamps. File mtimes are the clone time.
So the claim "committed to the repository before the runs" cannot be checked from this clone.

**C. "The value 7 … was not among them" contradicts the predictions file (claim 7).** The file has six entries,
not three. Entry `lfd-table` lists `"14": 7`, and `d6-falls-at-14` ends "(LFD(14) = 7)". The three predictions the
paper names are in the file and their outcomes are recorded as held.

**D. Free-extension counts (claim 7).** Paper: not closed on four instances; closed on four; hash equal on two.
Records: five D = 5 records (12289 unit, 12289 random, 13313, 15361, 40961 unit); not closed at k0 = 10 on five;
closed at k0 = 11 on five; hash equal to V_5(12) on three. The paper understates.

**E. Grid-point sentence (claim 10).** Paper: "eight with 2^20 or 2^21 grid points, two with 2^24". Records: eight
at 2^20 ((m, k) = (5, 4)) and eight at 2^21 ((7, 3)), so sixteen; four at 2^24 (two at (6, 4), two at (8, 3)).

**F. "Three instances" against "four" (claim 15).** Line 1472 says dim V_6(14) = 448,922 "on three instances".
Lines 314 and 355, Table `tab:profile` and the four records say four. (`\CERTTAB` "on three" for the checker and
"lift on three" are correct.)

**G. Binary row of Table `tab:twosided` (claim 15).**
- "≥ ⌊k/2⌋+1, inverse family, k ≥ 6" fails for odd k: the companion theorem gives 3 for odd n, and `table.jsonl`
  has LFD 3 at k = 7, 9, 11.
- "≤ 5 through k = 10" (also abstract, introduction, item 5, conclusion) is true on the seeded grids (k < n), but
  the inverse family at k = n = 10 has value 6. The same cell's "5 occurs from k = 8" relies on an off-grid system
  (the subfield check), so the cell mixes the two scopes.

**H. Missing record directories.**
- `experiments/ic-depth-step/` (five depth-step curves, line 1019): absent; `ic-curve-validation/REPORT.md` also
  cites it. The formula itself is backed on other instances (`ic-bounds-d5/verify/d4.json`, `ic-bounds-fall`).
- `experiments/ic-bounds-genericity/` (denominator degree 480): absent; prose only.
- `experiments/ic-binary-inverse/`: absent; only its rows inside `ic-bounds-binary-law/table.jsonl`.

**I. "Checked exactly at k = 7..12 on three instances" (line 1047).** Records give k = 8..12 on one instance and
k = 8, 9 on two more.

**J. "43 `decide +kernel` chunks" (claim 14).** The two Lean files have 42 such theorems each.

**K. Smaller items.**
- "80 new systems": 80 cases, 75 distinct systems.
- "1678 of 1775 depend only on (n, k) and the target class": the record groups by (n, k, r_eff, m_eff, class).
- Depth 9 "same curves": the random curve differs from the depth-8 one.
- GPU self-test "eight cases": the three later self-tests have ten.

**L. Evaluation rank 4068 at k = 14 on the random curve over p = 65537 (claim 6).** The only record is a sentence in `RESULTS.md` ("evaluation at Z is onto"); the k = 14 instance JSON has no
`eval_rank` field, unlike the k = 10 files. I recomputed it (`c6_rank.py`, 413 s): all 4068 stored points are zeros
of S_3 for target abscissa 1424, and the evaluation of 4500 randomly chosen monomials of weight ≤ 4 on them has rank
4068. So the rank is 4068 on M_<6 and M_≤6, and 118,370, 495,110 and 46,188 follow. I did not re-enumerate Z, so
completeness of the stored point list rests on the instance file.

## 3. Stale statements in the experiment READMEs

1. `ic-lean/README.md`: "The 339 audited declarations"; the record and the paper say 513. It also links to files
   that are not in the repository (`docs/IC_*.md`, `../ic-general-prime/verification.json`).
2. `ic-bounds-freeext/README.md` ends with the literal line `RESULTS_PLACEHOLDER`; it has no results section.
3. `ic-bounds-m3/README.md` line 226: "Depth k >= 5 is not run … nothing here bounds the m = 3 LFD for k >= 5".
   The paper reports k = 5 (twelve instances) and k = 6.
4. `ic-bounds-m3-depth/README.md` lines 161–187: the k = 6 run "was still waiting for the lock when this report was
   written", with the weak conjecture LFD(6) = 7. The outcome paragraph above it reports the run (LFD ≥ 8) and
   labels this text as the earlier plan.
5. `ic-bounds-d5/README.md` lines 86 and 171, `VERIFY.md` line 110: k = 12, 13 and degree 6 at k = 10 "not run
   here" / "beyond the limits", "(k = 10, 11 done)". The records `verify/k12-*`, `k13-*`, `k10-…-D6.json` exist and
   the paper uses them.
6. `ic-bounds-exact10/README.md` lines 235–244, `VERIFY.md` line 66: exact LFD at k = 10 and depths 11–13 listed as
   "Heavy runs needed (not run here)".
7. `ic-bounds-k14/RELAUNCH.md` lines 10–13 and 94: "Jobs B1 and B2 (p = 163841) were not run", "the sentence that
   says the k = 14 control has not been completed". A later section of the same file and `RESULTS.md` report job B.
8. `ic-bounds-k14/RESULTS.md`: "Limits: One instance …, one engine at k = 14, one run"; the lift label "one run, one
   implementation, one instance"; "checked independently on two instances". Later sections give four closures,
   three checks and three lifts. `LIFT.md` §9a lists two k = 14 lifts, not the third (p = 163841).
9. `ic-bounds-locus/README.md` line 196: "Exact LFD at k = 7 is not computed". `VERIFY.md` and
   `verify/p257-k7-d5-probe.jsonl` give LFD 5 for one C_2-class target, which the paper cites.
10. `ic-bounds-freeext/outcomes-2026-11-10.json`, `d5-14`: "on both instances" beside five records.
11. `ic-bounds-binary-law/README.md` line 7: "80 new preregistered cases (60 distinct systems)"; `score.json` has 75.
12. Stored source hashes that no longer match the tree: `ic-chain-depth8/verification.json`
    (`csd/ic_chain_shift_closure.py`, `scripts/run_ic_chain_depth8.py`), `ic-chain-depth10/verification.json`
    (`scripts/run_ic_chain_depth10_probe.py`); `ic-curve-sweep/verification.json` names `csd/gf2n.py`, which is
    absent. Data hashes all match.

## 4. Not checked

- Items the brief excluded (binomials, the k = 14 D = 6/7 numbers, the sweep tables, the msum self-test).
- I did not rerun any closure. Every "record" value above is what the stored file says, recomputed from its rows
  where the file has rows (counts, sums, formulas) and read from its fields otherwise.
- The membership of the sweep's 857 system in C_2 (line 1037) and the base-class sign tests beyond 22 of 30 and
  3 against 5 (which match `base-class-analysis.json`).
