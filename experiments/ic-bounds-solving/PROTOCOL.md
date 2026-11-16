# Solving degree alongside d_ff and d_lf for the ECDLP index-calculus systems — protocol

Frozen 2026-11-15T02:05Z, before any data. Goal: a rigorous statement of the (static, Macaulay)
solving degree and its growth in the depth k, side by side with the first fall degree d_ff and the
last fall degree d_lf, for the prime two-summand and three-summand chain systems and for the binary
(Weil-descent) systems. All results are [pcc:unverified] until an independent validator runs.

## Definitions (fixed for this study)

Ring R = K[u_1..u_r] (ordinary polynomial ring), generating list F (chains plus the normalized
summation generator for prime systems; the descended equations plus field equations for binary).

* V_{F,D}: iterated degree-D closure (Closure.lean; lfd_engine.py). d_lf = least d with
  I ∩ R_{≤D} = V_{F,D} for all D ≥ d (unit ideal: least D with 1 ∈ V_{F,D}).
* d_ff: least D with V_{F,D} ∩ R_{≤D−1} ≠ V_{F,D−1} (affine first strict fall).
* Mac_D = span{m·f : f ∈ F, m a monomial, deg(m f) ≤ D} (static Macaulay row space; no iteration).
  **sd** (static solving degree, Caminata–Gorla for a unit ideal) = least D with 1 ∈ Mac_D.
  For a proper ideal, record d_full = least D with Mac_D ⊇ I ∩ R_{≤D} and monomials of degree < D
  onto Z (the binary programme's definition), and d_lin.
* Always Mac_D ⊆ V_{F,D}, hence d_lf ≤ sd for unit ideals.

## Claim to test (prime chains, m summands, unit ideal)

Let A = K[G], G = (Z/2^k)^m, g = NF(generator) of degree δ. Lemma: NF(Mac_D) = g·M_{≤D−δ}
(M_{≤d} = span of chain monomials of Hamming weight ≤ d) and ker NF ∩ R_{≤D} ⊆ Mac_D.
Theorem: for a unit instance, **sd = δ + maxwt(supp g^{-1})**, g^{-1} the inverse in K[G].
Corollary: d_lf ≤ sd ≤ δ + mk, with equality sd = δ + mk iff the coefficient of x^{(N−1,…,N−1)}
in g^{-1} is nonzero.

## Preregistered predictions

* P1 (lane A). Brute-force sd in the ordinary ring (k = 3, 4; k = 5 if it fits) equals
  2 + maxwt(supp g^{-1}) on every unit two-summand instance tested. Falsifier: one mismatch.
* P2 (lane A). At k = 3..12, at least a fraction 1 − 3/p of random unit instances have
  sd = 2k + 2; the exceptions have sd ≤ 2k + 1. Falsifier: a fraction below 1 − 10/p at any k.
* P3 (lane A). On the same instances the exact d_lf (lfd_engine / recorded values) satisfies
  d_lf ≤ k + 2 < sd for k ≥ 1, with sd − d_lf ≥ k at every k ≥ 3.
* P4 (lane B). Three summands: sd = δ_3 + 3k generically (δ_3 = deg NF(S_4), expected 6; 5 on
  5t³ + 3at + 6b = 0), brute force agreeing at k = 2. Falsifier as in P1/P2.
* P5 (lane C). Binary m = 2 rows: d_full ≥ d_lf on every row (assertion of the programme, rechecked);
  d_full > d_lf on more than half of the rows with k ≥ 8 and sd/d_full grows with k while d_lf
  stays in {4, 5} (law B). Falsifier: a row with d_full < d_lf, or d_full constant in k.
* P6 (lane D). The written proof of the Lemma and Theorem survives an adversarial read: no step
  needs K to contain more than the 2^k-th roots of unity and p odd. Falsifier: a gap.

## Controls, budget, rules

k = 12 control instance p = 12289 reproduces the recorded dim V_4 and LFD values of the growth
record where available. Every lane: at most 2 threads, ulimit -v 6000000, writes only inside its
lane directory, RESULTS.md with commands and numbers, no commits. Statements of results carry
[pcc:unverified].

## Amendments (dated; predictions unchanged)

* 2026-11-15T02:40Z, before any lane reported. Naming: Caminata–Gorla's solving degree (Definition 1.1,
  arXiv:2112.05579, re-read) adds multiples of the lower-degree rows found by elimination (mutants), so
  by their Theorem 3.1 sd^CG = max(d_lf, max.GB.deg); for a unit ideal sd^CG = d_lf. The quantity
  defined above as "sd" (no iteration) is the static / XL Macaulay degree and is written **sd^stat**
  from here on. Exact check at m = 2, k = 1, p = 101 on 8 random unit instances: sd^stat = 4,
  d_lf = 3, 2 + maxwt(supp g^{-1}) = 4 (scratch script, main session) [pcc:unverified].
* 2026-11-15T02:45Z. Language and location: the crypto repository forbids Python for result
  generation (AGENTS.md "Implementation language"). The in-repository computation is therefore a
  standalone Rust program (`rust_check/`), and the Python lanes A–C (which use the ml-cryptanalysis
  engines, where Python is the convention) are moved into ml-cryptanalysis
  `experiments/ic-bounds-solving/` when the study lands; the crypto directory keeps the protocol,
  results summary, report, figures and the Rust program.
