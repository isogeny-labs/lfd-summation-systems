# Protocol: solving degree of the paper's exact-LFD systems (frozen 2026-11-14, before any GB degree was computed)

## Question

For every system whose exact last fall degree the paper reports, what is the solving degree, and does it equal
the last fall degree?

Here the solving degree is `sd_σ` in the sense of Caminata–Gorla (CG), Definition 1.1. That is the Lazard /
mutant-iterated Macaulay solving degree, not the step degree of an F4 run.

**Primary order.** `σ` = degree reverse lexicographic order (DRL), variables ordered by record index
(`u_0 > u_1 > …`, i.e. bit 0 of a monomial mask is the largest variable).

**Secondary order.** Graded lexicographic order (grlex), the order the paper uses for leading terms. It is reported
alongside DRL.

## Method

**The bridge.** By CG Thm 3.1, `sd_σ = max(LFD, max.GB.deg_σ)` (`docs/proofs/N10-solving-degree.md`, Theorem 1).
Each system is resolved in the first of these ways that applies:

| tag | rule | computation |
|---|---|---|
| `unit` | `Z = ∅`, so the reduced GB is `{1}` and `sd = LFD` (N10, Cor. 4(a)) | none |
| `count` | `LFD ≥ max(2, ⌊log2(|Z|+1)⌋)`, so `sd = LFD` (N10, Cor. 4(b)) | `|Z|` only |
| `cube` | `LFD ≥ max(2, r)`, so `sd = LFD` (N10, Cor. 4(c)) | none |
| `computed` | none of the above: `max.GB.deg_σ` computed exactly | as below |

**Where the `computed` rule is applied anyway.** It is applied to every proper system with at most `|Z| ≤ 3000`
points, regardless of tag, so that the exact `max.GB.deg` is known on all of them. The tag says which argument
settles `sd`; the computation is a check on that argument.

**How `max.GB.deg_σ` is computed.**
- `I` is the vanishing ideal of `Z`. N10 Lemma 3 proves this for every system here, because each input contains
  the chain or Boolean equations.
- The engine is the Buchberger–Möller algorithm on the squarefree staircase, applied to the record's solution
  set. It runs over `F_p` for prime chains and over `F_2` for Boolean systems (`csd/ic_solving_degree.py`).
- Monomials are processed degree by degree, in increasing `σ` order within a degree. A monomial is standard if
  and only if its evaluation vector on `Z` is independent of those of the smaller standard monomials. Rank
  profiles use `csd.elim.fp_echelon`.
- The minimal generators of `in_σ(I)` are:
  - the nonstandard squarefree monomials whose proper divisors are all standard;
  - `u_v^2` for each standard variable `u_v`.
- `max.GB.deg_σ` is their largest degree, or 0 if `Z = ∅`.

## Systems (all records with an exact LFD and a stored solution set or solution count)

| set | source | systems |
|---|---|---|
| S1 | `experiments/ic-curve-sweep/` (prime two-summand `k = 2..7`, binary/Koblitz) | all 3005 |
| S2 | `experiments/ic-chain-depth8/`, `experiments/ic-chain-depth9/` | all 32 |
| S3 | `experiments/ic-bounds-msummands/records/` (prime and binary, `m = 3..5`) with an exact LFD | all |
| S4 | `experiments/ic-bounds-msum/results.jsonl` with `lfd_exact` | all; tag `cube` or `count` from the stored `points`, else reported unresolved |
| S5 | `experiments/ic-bounds-binary/results-m2.jsonl`, `results-m3-*.jsonl` with an exact LFD | all |
| S6 | `experiments/ic-bounds-binary-deep/` instances with an exact LFD in `outcomes-2026-11-12.json` | all |
| S7 | the depth-10 unit instance and the depth-11 to 14 instances with an exact LFD (`ic-chain-depth10/`, `ic-bounds-k14/`) | tag only |

**Systems with no stored solution set.** They are tagged by `unit` / `count` / `cube` from their stored counts. If
none applies, they are reported as unresolved. They are never dropped.

**Systems with points outside `F_p`.** Prime systems with `2^k ∤ p − 1` have no `F_p`-rational solution grid. Any
such system with an exact LFD is listed separately: tagged if a rule applies, otherwise unresolved.

## Independent checks (all reported, pass or fail)

- **C1, Gröbner basis.** On a seeded sample (seed 20261114) of 150 S1 systems and every S2 system, compute the
  reduced DRL Gröbner basis of the input system (generators plus chain or Boolean equations) with Sage
  (Singular `std`, `degrevlex`, same variable order). Compare:
  - `max.GB.deg`;
  - the leading-monomial set, against the Buchberger–Möller staircase;
  - `vector_space_dimension`, against `|Z|`.

  This does not use the solution set.
- **C2, the bridge itself.** On 40 small S1 systems (`r ≤ 8`, prime and binary), compute `sd_DRL` directly from
  CG Definition 1.1, using an independent dense mutant-iterated Macaulay closure written for this check: the least
  `d` such that every Sage reduced-GB element lies in the closure. Compare it with `max(LFD_record, max.GB.deg)`.
- **C3, completeness of `Z`.** On 100 S1 prime systems, enumerate `Z` by brute force over `(cH)^2` and compare
  with the stored solutions.
- **C4, the counting lemma.** Check `|Z| ≤ 2^{k+1}` (N10 Lemma 5 / 5′) on every two-summand system of S1, S2 and
  S5.

## Reporting rule

- For each set, report the number of systems per tag.
- For the systems resolved by computation, report the distribution of `max.GB.deg` and of `max.GB.deg − LFD`.
- List every system where `max.GB.deg > LFD`, i.e. where `sd > LFD`. If there are none, say "none", with the
  counts.
- There is no success criterion. The result is reported whichever way it comes out.
- A failed check (C1–C4) is reported as a failure and blocks any statement in the paper that depends on it.

## Outputs

`experiments/ic-solving-degree/`:
- `rows.jsonl`: one row per system, with source, id, `r`, `|Z|`, LFD, tag, DRL and grlex `max.GB.deg` where
  computed, `sd`;
- `checks.json`;
- `summary.json`;
- `README.md`.

Driver: `scripts/run_ic_solving_degree.py`.
