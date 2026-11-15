# Growth in depth: numerical checks (`docs/proofs/N8-growth.md`)

Outputs of the tests listed in Section 6 of `docs/proofs/N8-growth.md`, the written proof of the growth
theorem. The proof is the evidence; these files check its finite lemmas and its reduction on small cases.
The engine is `scripts/ic_bounds_growth_lib.py` (exact linear algebra mod p, numpy only).

| Record | Script | What it checks |
|---|---|---|
| `semaev.json`, `lemmas.json` | `scripts/ic_bounds_growth_lemmas.py semaev`, `... lambda` (also `window`, `straddle`, `sparse`, which print only) | Semaev identity for the repository's formula; Lemma 4.1 (table of Λ); Lemmas 4.2 and 4.4 by brute force |
| `ideal-D3-*.json`, `ideal-D4-*.json` | `scripts/ic_bounds_growth_ideal.py p a b t D kmin kmax` | Propositions 2.1 and 2.2: dim V_D(k), dim S_D(k), V_D(k) ⊂ S_D(k), closedness of S_D(k) against FEI_k(D−1) |
| `lowpart-D3-*.json`, `lowpart-D4-*.json` | `scripts/ic_bounds_growth_lowpart.py` | Section 7: the low part V_D(k) ∩ M_<D against S_{D−1}(k) |
| `fei-w4-k*-p40961.json` | `scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 4 k` | the splitting property FEI_k(4) alone at k = 8, ..., 12 |
| `fei-w2-k*-p40961.json`, `fei-w3-k*-p40961.json` | `scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 w k` | FEI_k(w) for (w, k) = (2, 10..12) and (3, 10..11), run for the review of 12 November (`docs/proofs/N8-growth-review.md`); recomputed with equal records on 2026-11-14 |
| `nowrap-D3-k3.json`, `nowrap-D4-k7.json` | `scripts/ic_bounds_growth_nowrap.py` | Section 7: V_D(k) as the low part of a wrap-free closure |
| `circuits-w2-p12289.json`, `circuits-w3-p12289.json` | `scripts/ic_bounds_growth_circuits.py` | Theorem 3.3 over the constants, on computed circuits |
| `charp.json` | `scripts/ic_bounds_growth_charp.py` | Theorem 3.3 in characteristic p, dependence over the constants only, with a control |
| `fpdep.json` | `scripts/ic_bounds_growth_fpdep.py` | Theorem 3.3 over F^p itself (Hasse–Wronskian criterion), with the twisted control |
| `obstruction-D3-p12289-s12.json`, `obstruction-D4-p12289-s12.json` | `scripts/ic_bounds_growth_obstruction.py` | rank of the top-run map (a basis-independent measure of the failure of FEI) |
| `unit-existence.json` | `scripts/ic_bounds_growth_unit.py` | Corollary 5.2 (c): the construction of unit instances |
| `prime-grid-exceptions.json` | `scripts/check_n8_prime_grid_exceptions.py` (`--write` to regenerate; run in CI) | not about the written proof: for the fixed curve y² = x³ + 6x − 7, t = 1 of the kernel-checked theorem, every prime p ≡ 1 mod 2^(k+1) below 64·4^k at k = 2, …, 10, and whether the system is a proper ideal over F_p (405 of 8,567 primes; 124 of the 174 below 4^k) |
| `review2/` | scripts in that directory (independent of the repository code) | second red-team review, `docs/reviews-N8-growth-redteam2-20261112.md` |

Most scripts print their results and write a record only when an output path is given as the last argument
(for example `python scripts/ic_bounds_growth_lemmas.py semaev experiments/ic-bounds-growth/semaev.json`).
The exact arguments of each run and its result are given in the row of the table in Section 6 of
`docs/proofs/N8-growth.md` that names the record. Records marked "over the constants only" there test a weaker
statement than Theorem 3.3; the proof says which.
