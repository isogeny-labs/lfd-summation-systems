# Lane C: binary (Weil-descent / Koblitz) summation systems, static Macaulay degree d_full next to LFD, d_ff and d_reg

Frozen protocol: `../PROTOCOL.md` (prediction P5). Date 2026-11-15. Nothing committed. Every result statement below carries [pcc:unverified]
(no independent validator has run; the scan was cross-checked only against the repository's own code, see "Checks").
Setting: Boolean ring F_2[x]/(x_i^2+x_i); m = 2 rows have n degree-2 equations in r = 2k variables; m = 3 rows are read from the records only.
Definitions as in PROTOCOL.md: d_full = least D with Mac_D = I cap M_{<=D} and monomials of degree < D onto Z (proper ideal) / least D with 1 in Mac_D (unit ideal);
d_ff = least strict fall = min of the recorded `strict_falls`; d_reg = Boolean semi-regular degree d_reg(r_eff, m_eff) (laws.md; first non-positive coefficient of (1+t)^r/(1+t^2)^m, `analysis.semi_regular_dreg`);
kernel bound r+1.

## Verdict on P5 [pcc:unverified]

* **First half holds.** On every row with an exact d_full (m2: 400 cases / 359 distinct systems = 185 proper + 174 unit, k = 2..10; m3-direct 687; m3-split 53 cases) there is no row with d_full < LFD and no row with D* > d_full
  (the programme's assertion, rechecked on all of them, plus the systems rescanned here). No d_full lower bound that is <= LFD remains open on a scanned system (14 + 4 k = 10 LFD-5 systems were rescanned to close the repository's cap-5 = LFD cases).
* **Second half holds on the scanned sample, with k >= 8 only partly measured.** d_full > LFD for every proper-ideal system scanned with k >= 8 (6 exact, 124 by lower bound: k=8 40 of 40 proper, k=9 41 of 41, k=10 49 of 49 distinct systems);
  the only systems with d_full = LFD at k >= 8 are unit ideals (1 at k=8, 3 at k=9, 9 at k=10; all d_full = LFD = 2). Over all scanned k >= 8 distinct systems (unit included): 42 of 43 (k=8), 41 of 44 (k=9), 49 of 58 (k=10) have d_full > LFD, so "more than half" holds (of the 568 k >= 8 cases / 532 distinct systems only those with scan data are counted; the rest have no static record).
* **d_full is not constant in k, and the growth is measured, not extrapolated, up to k = 8.** On proper ideals d_full = k+1 or k+2 on every one of the 185 exact proper systems, k = 2..8 (143 at k+1, 42 at k+2; table "d_full - k" in `data/growth_table.md`). The split is by LFD and sharpens with k: LFD 3 gives k+2 on 17 of 17 systems at k = 7 and 2 of 2 at k = 8 (10 of 24 at k = 6, 5 of 14 at k = 5), LFD 4 gives k+1 on 20 of 20 at k = 7 and 4 of 4 at k = 8 (16 of 17 at k = 6), LFD 2 gives k+1 on all 29. Measured exactly at k = 8: LFD 3 -> d_full 10 (2 systems), LFD 4 -> d_full 9 (4 systems). At k = 9 and 10 only lower bounds exist: d_full >= 8 and >= 7 (consistent with k+1, not a measurement of it).
  LFD stays 2..4 for k <= 9 and reaches 5 only at k = 10 (law B): so d_full/LFD grows (median proper k=4: 5/3, k=8: 9/4 on the scanned rows), d_ff stays in {2,3}.
* Falsifiers: no row with d_full < LFD; d_full not constant in k. P5 not falsified. Caveat: d_full at k = 6..8 is a stratified sample (per (k, n, LFD, unit) group, hash order), not a random sample; at k = 9, 10 only a lower bound; the 18 scanned LFD-5 systems (k = 10) have d_full >= 7 only.
* Unit ideals do NOT show growth: on 174 unit systems scanned, k = 2..10, d_full <= LFD + 2 and at most 5 (k=7: LFD 4 -> d_full 5; k=8: LFD 4 -> 5, LFD 3 -> 5, LFD 2 -> 2; k = 9, 10: d_full = LFD = 2).

## Task 1: inventory (`data/binary_table.jsonl`, 3772 rows) [pcc:unverified]

Built by `tools/build_table.py` from `results-{m2,m3-direct,m3-split}.jsonl` (1871 cases, ic-bounds-binary) plus the 1901 m2 rows k = 2..7 of `ic-bounds-binary-law/table.jsonl` (source ic-curve-sweep; these have no related-*.jsonl, so d_full exists only where my scan ran). [pcc:unverified]
Columns: family (m2 / m3-direct / m3-split), n, k, r, e = r_eff - m_eff (m2 only; from `structure.jsonl`), Z, class (`target_class`, `unit`), LFD, D*, d_lin, d_full or d_full_lb (lower bound D with d_full >= D), d_ff, d_reg (m2 only), kernel r+1, provenance.
* d_ff: **computed from the recorded `strict_falls`** (= min of the list; null if no strict fall is recorded: 66 distinct m2 systems, k = 2: 41, k = 3: 19, k = 4: 3, k = 5: 1, k = 6: 2; all are unit ideals with LFD 0 except 2 proper k = 2 systems whose LFD is recorded as 0 = "no strict fall"). The closure dimension profile is not needed because the strict-fall list is stored.
  d_ff is 2 or 3 on every row (m2), median 2, independent of k. The equivalence of the Boolean-closure strict fall with PROTOCOL's affine d_ff (field equations as degree-2 generators) is assumed, not re-proved here.
* e, d_reg, class k+e are m2 only. For m3 rows e/d_reg are null (the law's r_eff, m_eff are defined for the m2 structure files); m3 rows are reported with d_lin, d_full, LFD, D*.
* Rows are cases (Koblitz (0,1) and (1,1) give identical systems: 1871 cases = 1726 systems); all statistics below use distinct system hashes unless stated "cases".
* P5 first half over all families (`tools/p5_check.py`, counts are cases): m2 2469 rows (exact d_full 400, lower-bound only 130), m3-direct 695 (687 exact, 8 lb at k=5), m3-split 608 (53 exact, 86 lb). D* > d_full: 0 everywhere. d_full < LFD: 0 everywhere.
  d_full > LFD by k (cases, exact + lb>LFD): m2 k=2: 52 (30 equal); k=3: 41 (22); k=4: 38 (22); k=5: 31 (11); k=6: 62 (14); k=7: 48 (8); k=8: 42 (1); k=9: 43 (3); k=10: 53 (9). m3-direct k=2: 173 (28 equal); k=3: 231 (0); k=4: 255 (0); k=5: 8 (lb). m3-split k=2: 46 (2 equal); k=3: 54 (0); k=4: 37 (lb). The "equal" cases are unit ideals (d_full = LFD = 2) or small-k rows.

## Task 2: static scan extended (`tools/static_scan.py`, `tools/gf2rank.c`) [pcc:unverified]

Rebuilds Mac_D for the stored equations (k = 2..7 from `ic-curve-sweep/systems.jsonl` by system hash; k = 8..10 with the repository builder `ic_bounds_binary_validate.rebuild`, system hash asserted equal), uses the repository's `csd.ic_dreg.boolean_macaulay`
for the matrix, and my own in-place bit-packed GF(2) echelon (8-pivot blocked elimination, OpenMP) for the rank. Z is recomputed by brute force over all 2^r assignments and asserted equal to the recorded |Z| on every system. Per degree it records rows, columns, rank, deficiency = cols - |Z| - rank, evaluation rank, 1-in-row-space, affine forms.
Memory model: matrix = rows x words x 8 bytes (rows = n * sum_{i<=D-2} C(r,i)); cap 4.3 GB at k = 8 (D = 10 is 4.0 GB), 2.5 GB at k = 9, 10. Per-degree wall: k = 8, D = 9: 8 to 10 min, D = 10: 25 to 29 min (1 thread, shared host).

| k | r | systems scanned | exact d_full (proper) | new bound where not exact | previous repo bound |
|---|---|---|---|---|---|
| 6 | 12 | 41 proper + 28 unit | 7 (16 LFD-4, 14 LFD-3 systems) or 8 (10 LFD-3, 1 LFD-4) | exact | none |
| 7 | 14 | 37 proper + 15 unit | 8 (20 LFD-4 systems) or 9 (17 LFD-3 systems) | exact | none (k <= 7 had no related file) |
| 8 | 16 | 6 proper + 3 unit | 10 (LFD 3, 2 systems), 9 (LFD 4, 4 systems); unit: 2, 5, 5 | exact | d_full > 6 (i.e. >= 7) |
| 9 | 18 | 41 proper (D <= 7) + 3 unit | not reached | **d_full >= 8** (deficiency > 0 at D = 7) | >= 6 |
| 10 | 20 | 49 proper (D <= 6) + 9 unit | not reached (D = 7 needs 5.2 GB, ~1e13 word ops) | **d_full >= 7** (deficiency > 0 at D = 6), incl. 18 of the 28 distinct LFD-5 systems (the other 10 have no static record) | >= 5 |

The k = 9 and 10 lower bounds are lower bounds on the number of the Macaulay degree needed, with deficiency (dim of (I cap M_{<=D}) / Mac_D) recorded in the scan files: at k = 8, the deficiency just below d_full is 9 = k+1 (LFD 4, D = 8) and 2 (LFD 3, D = 9);
the same pattern holds at k = 5..7. For systems with d_full = k+1 the deficiency at D = k equals k+1 in 92 of 143 (not universal). Unit rows: d_full in {2,...,5}.
**Not done:** exact d_full for k = 9 (D >= 8 needs ~5.8 GB for rows alone), k = 10 (D = 7: 5.2 GB, then D >= 11), m3 rows (README's caps were left as they are), the 2 k = 8 systems of the grid-order list not reached (n = 13 LFD 3; n = 14 LFD 3 proper), and any k = 8..10 system outside the first-row-per-group sample.
Scanned systems: `data/static_scan_k2_7.jsonl` (338), `static_scan_k8.jsonl` (5), `static_scan_k8b.jsonl` (4), `static_scan_k9_10.jsonl` (18), `static_scan_k10_lfd5.jsonl` (29: 14 + 4 chosen LFD-5 systems + 11 added by a stray default-selection run that I killed; all valid, none discarded). Peak RSS 3.4 GB (k = 8, D = 10).
d_reg yardstick column: in the table below and in `binary_table.jsonl` (`d_reg`); d_reg(r_eff, m_eff) = 5 at k = 8, 5..6 at k = 9, 6..7 at k = 10, i.e. d_full > d_reg on every exact proper system with k >= 4 (equality on 7 systems at k = 2, 3 and 1 more at k = 3; d_full - d_reg is 4..5 at k = 8).

## Task 3: growth table (m2, distinct systems) (`tools/growth.py`; full tables `data/growth_table.md`, `data/growth_summary.json`) [pcc:unverified]

Proper ideals (|Z| > 0), median [min..max] over systems with a value (d_full: scanned exact systems only):

| k | systems | LFD | D* | d_ff | d_reg | d_full exact | d_full lower bound (n systems, median) | kernel r+1 |
|---|---|---|---|---|---|---|---|---|
| 2 | 53 | 2 [0..3] | 2 [2..3] | 2 [2..3] | 2 [2..3] | 3 [3..4] (26) | | 5 |
| 3 | 122 | 3 [2..3] | 3 [2..3] | 2 [2..3] | 3 [2..4] | 4 [4..5] (27) | | 7 |
| 4 | 226 | 3 [2..4] | 3 [2..4] | 2 [2..3] | 3 [2..4] | 5 [5..6] (28) | | 9 |
| 5 | 192 | 3 [2..4] | 3 [2..4] | 2 [2..3] | 4 [3..4] | 6 [6..7] (20) | | 11 |
| 6 | 216 | 3 [3..4] | 3 [3..4] | 2 [2..3] | 4 [3..5] | 7 [7..8] (41) | | 13 |
| 7 | 188 | 4 [3..4] | 4 [3..4] | 2 [2..3] | 5 [4..5] | 8 [8..9] (37) | | 15 |
| 8 | 170 | 4 [3..4] | 4 [3..4] | 2 [2..3] | 5 [5..5] | 9 [9..10] (6) | >= 7 (31) | 17 |
| 9 | 164 | 4 [3..4] | 4 [3..4] | 2 [2..3] | 5 [5..6] | not reached | >= 6 (39), best >= 8 | 19 |
| 10 | 155 | 4 [3..5] | 4 [4..5] | 2 [2..3] | 6 [6..7] | not reached | >= 5 (46), best >= 7 | 21 |

By k+e class (law B classes; proper ideals; LFD / d_full, medians): k+e in 3..8: k=2..6 LFD 3 (3 [2..4]), d_full = 3..7 (= k+1 mostly); k+e in 9..18: LFD 4 for k = 5..10 (3..4), d_full 6.5 (k=5), 7 (k=6), 8 (k=7), 9 (k=8, 6 systems); k+e >= 19 (k = 10 only, 28 systems): LFD 5 (all), D* 5, d_ff 3 [2..3], d_full >= 7 (best lb, 18 systems rescanned). Per (k, class) rows are in `data/growth_table.md`. The class carries LFD (law B) but d_full depends on k and on LFD (k+1 at LFD 4, k+2 at LFD 3), not on the class by itself.
Figure `data/binary_sd_vs_k.svg`: LFD (median, min-max band), semi-regular d_reg, d_full exact (scanned systems), d_full lower bounds (open triangles, k = 9, 10), kernel bound r+1 = 2k+1. Hand-written SVG (no matplotlib on this host); not rendered to PNG/PDF.
Also from item05 (`../../leads_followup_20261114/item05_binary_balanced/RESULTS.md`, read only, not rerun): balanced unit systems k = 18, 19 (|Z| = 0) have 1 not in V_4, i.e. LFD >= 5 (not certified exactly); no static scan was run on them (degree-5 Macaulay at r = 36 has 443,704 columns). m3 numbers (item10) were not read beyond the repository's records.

## Task 4: theory check, one sentence [pcc:unverified]

In F_2[x]/(x_i^2+x_i) (the ring of functions on F_2^{2k}) the descended system has n generators f_1..f_n of degree 2, the ideal is the set of functions vanishing on Z, and Mac_D = sum_i f_i * M_{<=D-2} reduced mod the field equations, so there is no single-generator inverse and no group-algebra image: unit means a Bezout identity 1 = sum g_i f_i whose degree is a minimum over a coset of the syzygy module (not the support of one inverse), the quotient ring is a product of 2^r copies of F_2 rather than K[(Z/2^k)^m] with p odd, and the prime-case NF(Mac_D) = g * M_{<=D-delta} fails;
the data agree: binary unit rows have d_full <= 5 (<= LFD + 2) for k = 2..10 (no growth with k, unlike delta + mk in the prime unit case), while proper rows have d_full = k+1 or k+2 for k = 2..8 and the deficiency profile (k = 8, D = 6..9: 462, 64, 9, 0 on an LFD-4 system; D = 7..10: 120, 25, 2, 0 on an LFD-3 system) is set by the bilinear (a_i b_j) structure and by Z, so the prime theorem, which is about a unit instance, does not transfer.

## Task 5: commands, reruns vs records, checks, skips

Environment: `ulimit -v 6000000`, at most 2 threads in total (two single-thread processes, or one 2-thread process); `python3 -I`. Shared host (load 3 to 6): wall times are inflated. ml-cryptanalysis commit 4575505b (read only).
Inputs (sha256 prefix): results-m2.jsonl 3902fa9f4e2aae74, related-m2.jsonl 7c05334bf8f3cb61 (= README), table.jsonl c6f025aeb1a5ffc9, structure.jsonl b449f7e57914d110. Output `binary_table.jsonl` 97c8d840bf72bf6c (built after all scans).
```
export GF2RANK_SO=<scratch>/gf2rank.so OMP_NUM_THREADS=1;  ulimit -v 6000000
python3 -I tools/static_scan.py --k 2 3 4 5 6 7 --per-group 2 --threads 1 --out data/static_scan_k2_7.jsonl
python3 -I tools/static_scan.py --k 6 7 --per-group 4 --threads 1 --out data/static_scan_k2_7.jsonl      # adds 54
python3 -I tools/static_scan.py --k 8 --per-group 1 --threads 1 --max-gb 4.3 --max-cols 70000 --time-cap 2400 --out data/static_scan_k8.jsonl   # killed after 5 of 11
python3 -I tools/static_scan.py --k 8 --hashes <4 hashes> --threads 1 --max-gb 4.3 --out data/static_scan_k8b.jsonl
python3 -I tools/static_scan.py --k 9 10 --per-group 1 --threads 1 --max-gb 2.5 --max-cols 140000 --out data/static_scan_k9_10.jsonl
python3 -I tools/static_scan.py --k 10 --hashes <LFD-5 systems> --threads 1 --max-gb 2.5 --max-cols 140000 --out data/static_scan_k10_lfd5.jsonl
python3 -I tools/build_table.py; python3 -I tools/growth.py; python3 -I tools/p5_check.py; python3 -I tools/compare_related.py; python3 -I tools/crosscheck.py 8
```
* **Rerun here:** Mac_D rank and d_full on 338 + 5 + 4 + 18 + 29 = 394 scan records (k = 2..10, 386 distinct systems), Z by brute force. **Read from records:** all LFD, D*, strict falls, |Z| (as a check), d_lin and the capped d_full of the other m2 rows, all m3 values, structure (r_eff, m_eff). The LFD values were not recomputed.
* **Checks:** (a) new kernel vs `csd.elim.gf2_echelon` on random matrices (5 shapes, identical rank and pivots) and on the degree-4 Macaulay matrix of 8 systems (identical, `data/crosscheck.json`); (b) `tools/crosscheck.py`: my d_full equals `ic_bounds_binary_related.static_degrees` on 8 of 8 systems k = 5..7; (c) `tools/compare_related.py`: 3 exact-vs-exact agreements with `related-m2.jsonl` (k = 8 unit d_full 5; k = 9, 10 unit d_full 2) and 18 consistent caps (our value >= the repository's capped_at), no disagreement (`data/related_agreement.json`); (d) |Z| brute force equals the recorded value on every scanned system; (e) D* <= d_full asserted per record.
* **Skips and limits:** no k = 9, 10 exact d_full (memory/time); k = 8 only 9 systems (all LFD-3/4 combos at n = 11..14 not covered; no LFD 2 proper); m3 not rescanned; d_lin only from the repository's capped scan (my scan stores affine-form counts per degree, not used in the table); figure SVG not rendered or visually inspected; no KB/ecdb queries (the research-kb server was not used; this lane only reads repository records); no Conductor task / resume checkpoint / benchmark record was created (local development analysis, no experiment routed to EKS); nothing committed.
* Anomalies noted: three k = 2 proper rows record LFD 0 (no strict fall; D* = 3), outside the "D* - LFD in {0,1}" statement of the README for ic-bounds-binary rows (those 3 are ic-curve-sweep rows). A stray default-selection run of the scanner (empty `--hashes` list is falsy) added 11 extra valid k = 10 scans before I killed it.
Note: tools/gf2rank.so is not shipped (rebuilt from gf2rank.c on demand; runs used GF2RANK_SO=<scratch>/gf2rank.so built from the same .c).
