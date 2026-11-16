# Validator report: solving-degree study (SD1 Theorem, lanes A–D, Rust check)

Validator: independent session (Mythos 5.1), 2026-11-15, 05:18–06:05 UTC (about 47 minutes wall, 2 processes
at most, `ulimit -v 8000000`, Python 3.13.18 + numpy 2.5.3 only, cargo 1.97.0, gcc 13.3). Nothing was committed;
every file written is inside this directory. No producer tool was imported; the producers' Rust binary was run
only as a black box (check 3). Every number below that I did not recompute is tagged `[pcc:unverified]`.

## 0. Frozen artifact set and drift

| repo | dispatch commit | what I validated from | drift observed during the run |
|---|---|---|---|
| crypto `claude/cool-davinci-szm3om` | `ef468c885b89ec7f537abfc0194f743bdd385fbc` | `git archive ef468c885b research/solving_degree_20261115` extracted to the scratchpad; hashes of all 25 files in `out/producer_files_sha256_crypto_ef468c885b.txt` | branch moved to `2678160ab0` (commit `293278ad59` rewrote the d_lf k = 9..14 cells of RESULTS.md and moved `m2_p17_k3_many.jsonl`, `m2_p257.jsonl` to the artifact store) and then to `c162fdb151` (merge from main; no change under `research/solving_degree_20261115`). The working-tree RESULTS.md (sha256 `ad488b6f…`) differs from the committed one (`b081c844…`). |
| ml-cryptanalysis `claude/lfd-proofs-20261108` | `bc98cdc71ed6944f5caf7b6db0f552004d4bdacf` | files on disk, each verified byte-identical to the `bc98cdc7` blob (hashes in `out/producer_files_sha256_mlcrypt.txt`, committed hashes checked with `git show bc98cdc7:<path>`) | branch moved to `77447e8d` ("SD1: Lean formalization plan; remark on delta = 2 for k >= 2"): `docs/proofs/SD1-solving-degree.md` changed from `ec7d5476…` to `7c5d7597…` (one added remark in §4.1, correct as far as I checked: the normalized generator contains u_{1,1}u_{2,1} with coefficient 1, so δ = 2 for k ≥ 2); `laneD_theory/SD1-solving-degree.md` is still the dispatch version. Untracked `experiments/ic-bounds-solving/laneC_k9/` appeared (out of scope). |

Finding F1 (process): the producers kept editing both branches while the validation ran. The dispositions below
are for the dispatch commits only; the later commits are not validated.

Key producer file hashes (sha256, from the dispatch commits): SD1 `ec7d54767d2f21ee…`; PROTOCOL.md
`b266b00b6064b2d5…`; crypto RESULTS.md `b081c8448adb5f3d…` (full values in the hash files); rust_check `src/main.rs`
`62848b0c11bfb3f5…`; `data/CHECK_LIST.md` `a452b42b57948aa4…`; lane A RESULTS.md `4d30d59ab63a9f4a…`;
`all_instances.json` `889fe91e1c8281d5…`; lane B RESULTS.md `d2ecc854a083c3fa…`; `bruteforce_k2.jsonl`
`84df606995527138…`; lane C RESULTS.md `ae5b435ef749d2ce…`; `binary_table.jsonl` `97c8d840bf72bf6c…`;
`static_scan_k2_7.jsonl` `1c75fe2016a5ccc0…`; REDTEAM.md `416fca9606e275e2…`; `ic-curve-sweep/systems.jsonl`
`99248b4338e1dd13…`. Built Rust binary (release, from the frozen source): `3a65dc5bf31463ea…`.

## 1. Proof reading (SD1 §2–§4)

Read in full; each step re-derived by hand. Result: **no gap found in the Lemma, the Theorem or Corollaries
4.1–4.3.** Details:

| step | verdict |
|---|---|
| Lemma 2.1(a): chains are a Gröbner basis (pairwise coprime leading monomials u_{i,j}²), division is degree-non-increasing, reduced representatives are unique, NF injective on reduced polynomials | verified; also observed computationally at every degree D on all 344 brute-force instances (chain-row pivot count = #monomials≤D − #squarefree≤D, no squarefree lead) |
| Lemma 2.1(b): NF(μ) is a monomial of weight ≤ deg μ (weight inequality), every x^e with wt(e) ≤ d is reached, deg(μĝ) = deg μ + δ in a domain, hence NF(Mac_D) = g·M_{≤D−δ} | verified |
| Theorem 3.1 (i)⇔(ii)⇔(iii), and 1 ∈ Mac_D ⇔ D ≥ δ + maxwt(supp g^{-1}) (⇒ via uniqueness of the inverse; ⇐ via P = ι(g^{-1})ĝ and (a); edge case D < 2) | verified |
| Cor. 4.1 (bounds, c_top criterion) | verified |
| Cor. 4.2 (character sum for c_top; needs p ∤ N, H ⊂ K) | verified; checked numerically on all 344 + 30 + 17 + 12 + 6 instances (`c_top_charsum_agrees`) |
| Cor. 4.3 (dichotomy, |E| ≤ N² − 1) | conclusion verified; **wording slip**: the second case says "min_O deg q_O ≤ 1 (a diagonal orbit), so again deg P ≤ N² − 1"; the inequality needed is min_O deg q_O ≥ 1 (off-diagonal q_α have t-degree exactly 2 because (α₁−α₂)² ≠ 0, so once no diagonal q_O is constant every orbit has degree ≥ 1). Not a gap. |
| Cor. 4.4 and the d_lf side (U1: d_lf ≤ k + 2; CG Thm 3.1: sd^CG = d_lf; N8; HKYY) | cited, **not verified here** [pcc:unverified] |

Steps I could not verify: only the external citations above. Cosmetic: a stray `[pcc:unverified]` tag sits in
the proof text of Lemma 2.1(a) (line 101).

## 2. Own static-Macaulay brute force vs own g^{-1} (check 2)

Code: `vcommon.py`, `v_static.py`. The ordinary-ring matrix is literal (every chain multiple and every
generator multiple of degree ≤ D as a row over all monomials of degree ≤ D); the elimination inserts the chain
rows into an echelon table, fully reduces the generator rows by it, and tests e_1 against the reduced generator
rows (justified in the file header; nothing from the Lemma is assumed). g^{-1} by evaluation on H^m, pointwise
inversion and my own inverse DFT (16-bit-split matmul, exact int64), with an exact cyclic-convolution check
g·g^{-1} = 1 on every instance. Instances: my own RNG (seeds 20261115–20261123), random nonsingular (a, b),
random t, unit test on H^m. The fully dense elimination of the whole matrix was also run where it fits
(m = 2, k = 1, 2 and m = 3, k = 1: 11 instances, every degree) and agreed with the sparse answer at every D.

| set | instances | sd_static = δ + maxwt | values |
|---|---|---|---|
| m = 2, k = 3, p ∈ {97, 113, 193, 257} × 10 | 40 | **40/40** | sd = 8 (all), c_top ≠ 0 |
| m = 2, k = 3, p ∈ {17, 41, 73} × 40 | 120 | **120/120** | 117 × sd = 8; **3 exceptions at p = 17 with sd = 7** (c_top = 0) |
| m = 2, k = 4, p ∈ {97, 113, 193, 257} × 5 (+1 probe) | 21 | **21/21** | sd = 10 |
| m = 2, k = 4, p = 113 × 100 | 100 | **100/100** | 98 × sd = 10; **2 exceptions with sd = 9** |
| m = 2, k = 5, p ∈ {97, 193, 257} × 1 | 3 | **3/3** | sd = 12 (63–67 s each, 646,646 monomials) |
| m = 3, k = 2, p ∈ {97, 113, 193} × 3, 257 × 1 | 10 | **10/10** | δ = 6, sd = 12 |
| m = 3, k = 2, p ∈ {13, 17} × 25 | 50 | **50/50** | 44 × (δ = 6, sd = 12); 5 × (δ = 5, sd = 11, locus 5t³+3at+6b = 0); **1 exception (δ = 6, sd = 11) at p = 13** |

Total **344 instances, 0 mismatches**; in all 344, sd = δ + mk ⇔ c_top ≠ 0, every exception has
sd = δ + mk − 1 exactly, and the convolution and character-sum checks pass. Data: `out/static_*.jsonl`
(per-degree rank profile recorded; e.g. m = 2, k = 3: rank Mac_7 = 1715 of 1716 monomials, Mac_8 = 3003 = all).

## 3. Producers' Rust binary (check 3)

Built from the frozen source with `cargo build --release -j 2` (0.8 s, no crates). `--selftest` passes
(`out/rust_selftest.txt`). All 7 rows of `rust_check/data/CHECK_LIST.md` re-run with the listed commands
reproduce maxwt, δ, sd and the top-coefficient flag exactly, including row 6 (p = 17 exception, sd = 7) and row 7
(k = 1 coefficients 91, 81, 81, 95): `out/rust_checklist_rerun.jsonl`. On 7 of my own brute-forced instances
(m = 2: k = 3 p = 97, k = 4 p = 97, k = 3 p = 17 generic and exceptional; m = 3: k = 2 p = 97, p = 13 generic and
exceptional) the binary agrees with my maxwt, δ, sd_static and c_top flag **7/7** (`out/rust_my_instances.json`).
`main.rs` was grepped for the check-list constants: none are hard-coded (only the documented k = 1 self-test row).

## 4. Lane A `data/all_instances.json` (check 4)

Code: `v_spot.py`, output `out/spot_laneA.json`.

Recounts from the rows: **29,775 rows** (k = 3: 19,179; 4: 6,970; 5: 2,723; 6: 438; 7: 122; 8: 73; 9: 85;
10..12: 60 each; 13: 5). **159 exceptions** (sd_inverse ≠ 2k + 2), all with sd = 2k + 1, all with
`zero_top_coeff` = true, and `zero_top_coeff` ⇔ exception on every row. **k = 3, p = 73: 40 / 3060 = 1.307 %**
(1/p = 1.370 %). Brute-force rows 18,828 (k = 3: 18,492; 4: 282; 5: 48; 6: 6), 0 internal mismatches, 127
exceptional rows among them. `inverse_verified` true on all rows. All 63 (k, p) cells have generic fraction
≥ 1 − 3/p; none is below 1 − 10/p. Minor: **33 duplicate (p, a, b, t, k) rows** (the 29,775 are not all distinct
instances; the lane does not say so), and the cell (k = 3, p = 113) with 1 exception in 60 rows is missing from
the lane's exception table (it is in the data and consistent with the claims).

30 stratified random rows (k = 3: 6, 4: 6, 5: 4, 6: 3, 7: 3, 8: 3, 9: 2, 10: 3; one exception forced in at k = 3
and at k = 4) recomputed with my g^{-1}: **30/30 agree** on unit status, maxwt, sd, top-coefficient flag and the
number of zero coefficients, convolution check true. The 12 rows with k ≤ 4 also got my ordinary-ring brute
force: **12/12** equal to the row's sd (and to its `sd_bruteforce` where present), including the two exceptions
(k = 3 p = 97 sd = 7; k = 4 p = 401 sd = 9).

Six further lane A rows at k = 11, 12, 13 (two each, p up to 2^31 − 2^11) recomputed with my own radix-2 NTT
(`v_growth.py`, self-checked against the matmul DFT): **6/6 agree** (sd = 24, 26, 28; `out/growth.json`).

## 5. Lane C binary d_full (check 5)

Code: `v_binary.py` (own Boolean Macaulay rows, own GF(2) bitset elimination, Z by brute force over 2^r),
equations found as lane C did (`experiments/ic-curve-sweep/systems.jsonl` by system hash). Five proper systems:

| system hash | k | n | |Z| | lane C d_full (LFD, D*) | mine | rank profile | still full at D+1, D+2 |
|---|---|---|---|---|---|---|---|
| 00357f759d06… | 3 | 10 | 2 | 4 (2, 2) | **4** | identical | yes |
| 02669d0a927b… | 4 | 6 | 8 | 5 (3, 3) | **5** | identical | yes |
| 01979f83c236… | 5 | 8 | 4 | 6 (3, 3) | **6** | identical | yes |
| 04ec65496036… | 6 | 9 | 12 | 7 (3, 3) | **7** | identical | yes |
| 0898244a2d55… | 6 | 7 | 32 | 8 (3, 3) | **8** | identical | yes |

**5/5 agree**, |Z| agrees, ranks and columns agree at every degree lane C recorded (`out/binary_dfull.json`).
Recount of `binary_table.jsonl`: 2,469 m2 rows, 2,277 distinct systems, **185 proper systems with an exact
d_full, k = 2..8, 143 at k + 1 and 42 at k + 2** (k = 7: LFD 4 → k + 1 on 20, LFD 3 → k + 2 on 17; k = 8: 4 and 2),
359 exact m2 systems with **0** rows d_full < LFD and **0** rows D* > d_full, 174 unit systems with d_full ≤ 5 and
d_full − LFD ≤ 2, lower bounds ≥ 7 (k = 8, 10) and ≥ 8 (k = 9). Not recomputed by me: any k = 7, 8 system and the
LFD / D* values themselves [pcc:unverified]. Definitional note: lane C reports the first D where the containment
condition holds, SD1 §5 defines d_full^stat as the least D₀ from which it holds for all D ≥ D₀ (REDTEAM A12); on
the five systems the condition stayed true at D + 1 and D + 2, so no numerical difference was observed.

## 6. d_lf by my own iterated closure (check 6)

Code: `v_dlf.py` (ordinary ring, variable-wise closure with degree cap, echelon basis in a degree-compatible
order). On 15 lane A instances carrying a recorded d_lf (k = 3 p = 97 × 5, k = 4 p = 1009 × 5, k = 5 p = 3137 × 5):
**d_lf = 4 and d_ff = 4 on all 15, agreeing 15/15 with lane A** (dim V_3 = 49/81/121 < all, dim V_4 = 210/495/1001 =
all monomials). Gaps sd − d_lf = 4, 6, 8 at k = 3, 4, 5. On 3 of my m = 3, k = 2 instances: **d_lf = d_ff = 7**
(dim V_6 = 861, V_7 = 1716 = all), matching lane B and sd = 12 (gap 5). One lane A instance at k = 6
(p = 12289, a = 8806, b = 10016, t = 9788; 12 variables): **d_lf = 5, d_ff = 4, agreeing with lane A's recorded
values** (dim V_4 = 1631 < 1820 = all monomials of degree ≤ 4, so 1 ∉ V_4; dim V_5 = 6188 = all, so 1 ∈ V_5;
933 s for D = 5; `out/dlf_k6.jsonl`). Gap sd − d_lf = 14 − 5 = 9 at k = 6.

## 7. m = 3 growth beyond the brute-force range (through the Theorem)

`v_growth.py`: 17 lane B scan rows at k = 3..8 (p = 65537, 8390273, 67109633) recomputed: **17/17 agree**
(δ = 6, maxwt = 3k, sd = 3k + 6); lane B's single k = 3 brute-force instance agrees (maxwt 9, sd 15); 12 fresh
unit instances of my own at p = 65537, k = 3..6: all δ = 6, sd = 3k + 6, none on the δ = 5 locus.

## 8. Dispositions

| claim | disposition | basis and limits |
|---|---|---|
| (a) Theorem: sd^stat = δ + maxwt(supp g^{-1}) for unit instances | **validated** | proof re-derived with no gap (§1); 344/344 own literal brute-force agreements at m = 2, k ≤ 5 and m = 3, k = 2, including 6 exceptional and 5 locus instances; dense cross-check at k ≤ 2. Limits: computational depth k ≤ 5 (m = 2), k = 2 (m = 3); only the normalized generator ĝ = ι(g) was tested (δ = 2, 6, 5); raw S_3 and random generators were tested by lane D only [pcc:unverified]. |
| (b) growth: sd^stat = 2k + 2 (m = 2) and 3k + 6 (m = 3) on generic unit instances, exceptions with c_top = 0 at rate ≈ 1/p | **validated** for m = 2, k = 3..13 and m = 3, k = 2..8 as a Theorem-backed statement: exact iff c_top ≠ 0 (proved); every exception seen (my 6, lane A's 159, lane B's rows I touched) has c_top = 0 and sd = δ + mk − 1; all 63 lane A cells pass 1 − 3/p. Limits: for k ≥ 6 (m = 2) and k ≥ 3 (m = 3) the value rests on the Theorem plus my g^{-1} recomputation of 30 + 6 + 17 rows, not on brute force; "≈ 1/p" is empirical at k = 3, 4 (m = 2) and k ≤ 4 (m = 3) only, and SD1 itself says the N² − 1 bound is vacuous in the regime p ≈ c·2^k. |
| (c) d_lf = 4, 4, 4, 5, 5, 5 at k = 3..8 (m = 2); 7 (m = 3) | **validated for k = 3, 4, 5, 6 (m = 2; 15/15 at k = 3..5, 1/1 at k = 6) and k = 2 (m = 3; 3/3); insufficient_evidence for k = 7, 8** | the k = 7, 8 values are lane A's quotient-closure measurements on 2 + 1 instances [pcc:unverified]; my k = 6 check is a single instance (ordinary-ring closure, 933 s). m = 3 for k ≥ 3 is lane B's one instance and the repository record [pcc:unverified]. |
| (d) binary proper m = 2: static d_full ∈ {k + 1, k + 2} for k = 2..8 | **validated on the recomputed systems (k = 3..6, 5/5) and as a consistent reading of the recorded table (185/185); insufficient_evidence as an independent recomputation at k = 7, 8** | no k = 7, 8 system was rescanned by me; LFD and D* are read from the records [pcc:unverified]. |

## 9. Findings and disagreements

1. **F1 (process)**: both branches moved during the validation (§0). No disagreement in content, but the artifact
   set named by the dispatch is not what the working trees now show.
2. **F2 (data)**: 33 duplicate instance keys in `all_instances.json`; the (3, 113) exception cell is absent from
   lane A's table. Neither changes any count claimed in RESULTS.md (29,775 rows and 159 exceptions are row counts
   and reproduce).
3. **F3 (text)**: Cor. 4.3 proof wording "min_O deg q_O ≤ 1" should be "≥ 1" (§1). Conclusion unaffected.
4. **F4 (definition)**: lane C's d_full is "first D with the condition", SD1's is "least D₀ from which it holds
   for all larger D"; equal on the 5 systems I checked (condition persisted two degrees further).
5. No numerical disagreement was found anywhere: 344/344, 7/7 + 7/7, 30/30 + 12/12, 15/15 + 3/3, 5/5, 17/17,
   12/12, 6/6.
6. The PROTOCOL's control "k = 12 instance p = 12289" is infeasible (no unit instance exists there); lane A says so
   and substituted a k = 7 control. Not re-run by me [pcc:unverified].
7. k = 6 d_lf closure (`out/dlf_k6.jsonl`, `out/log_dlf_k6.txt`): finished 06:02 UTC, d_lf = 5, d_ff = 4, agreeing
   with lane A (one instance; §6). Final wall time of the whole validation: 05:18–06:05 UTC.

## 10. Files written (all under this directory)

`vcommon.py` (own library: primes, roots of unity, S_3/S_4, exact DFT by matmul and by radix-2 NTT, inverse,
weights, monomial table, sparse echelon, dense rank), `v_static.py`, `v_dlf.py`, `v_spot.py`, `v_binary.py`,
`v_rust_compare.py`, `v_growth.py`; outputs in `out/`: `static_m2_k3.jsonl`, `static_m2_k3_smallp.jsonl`,
`static_m2_k4.jsonl`, `static_m2_k4_p113.jsonl`, `static_m2_k4_probe.jsonl`, `static_m2_k5.jsonl`,
`static_m3_k2.jsonl`, `static_m3_k2_smallp.jsonl`, `smoke_*.jsonl` (dense cross-checks), `dlf.jsonl`,
`dlf_k6.jsonl`, `spot_laneA.json`, `binary_dfull.json`, `rust_checklist_rerun.jsonl`, `rust_selftest.txt`,
`rust_my_instances.json`, `growth.json`, `producer_files_sha256_crypto_ef468c885b.txt`,
`producer_files_sha256_mlcrypt.txt`, `log_*.txt`. The frozen crypto tree and the built binary live in the
session scratchpad (`frozen_crypto/`), not in either repository.
