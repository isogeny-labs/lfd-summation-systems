# Solving degree of the LFD paper's systems

The question: for every system in `paper/curve-lfd-bounds.tex` with an exact last fall degree, what is the
solving degree, and is it equal to the LFD?

`PROTOCOL.md` (frozen in commit 5a4c2307, before any Gröbner-basis degree was computed) fixes the systems, the
orders, the tags and the checks. The written proofs are in `docs/proofs/N10-solving-degree.md`. The paper section
is `\section{The solving degree}` (`sec:solving`).

## The argument

- **Solving degree.** This is Caminata–Gorla, Def. 1.1: the Lazard / mutant-iterated Macaulay solving degree. It
  is not the step degree of an F4 or F5 run. In particular it is not the Magma "degree of regularity" of Petit,
  Kosters and Messeng.
- **The bridge (CG Thm 3.1).** For every degree-compatible order σ, `sd_σ = max(LFD, max.GB.deg_σ)`.
- **The cube lemma (N10 Lemma 2–3).** Every system here contains the chain or Boolean equations. Hence:
  - the ideal is the vanishing ideal of its solution set `Z`;
  - standard monomials are squarefree;
  - `max.GB.deg ≤ max(2, min(r, ⌊log2(|Z|+1)⌋))`.
- **Two summands (N10 Lemma 5/5′).** `|Z| ≤ 2^{k+1}`, hence `sd ≤ max(LFD, k+1) ≤ k+2`.
- **m summands.** `sd ≤ max(LFD, mk) ≤ mk+1`.
- **Exactness and transfer.**
  - `sd = LFD` exactly when `Z = ∅` or `LFD ≥ max(2, ⌊log2(|Z|+1)⌋)`.
  - Every lower bound on the LFD transfers to the solving degree, including unboundedness for unit instances.
- **Prior work.**
  - The linear bound `max.GB.deg ≤ |Z|` is folklore.
  - Steiner (arXiv:2310.03637, Lemma 7.4) proves the DRL upper bound.
  - Mészáros (arXiv:1908.03045) treats Boolean standard-monomial down-sets.
  - We found no statement of the logarithmic form, and claim no novelty for it.

## Commands

All commands need Python ≥ 3.10 (they use `int.bit_count`).

```bash
python3.13 -m csd.ic_solving_degree                       # self-test
python3.13 scripts/run_ic_solving_degree.py               # rows.jsonl (cap 3000 solutions, as frozen)
python3.13 experiments/ic-solving-degree/rebuild_binary_msum.py   # rows-s3-binary.jsonl (9 binary m>=3 systems, sha-checked)
python3.13 scripts/run_ic_solving_degree.py --sets S6 --max-points 9000 \
    --out experiments/ic-solving-degree/rows-s6-cap9000.jsonl     # beyond the frozen cap; labelled
python3 experiments/ic-solving-degree/summarize.py        # rows-merged.jsonl, summary.json, paper/solving-degree-data.tex
python3.13 scripts/check_ic_solving_degree.py export      # sage-in.json
sage experiments/ic-solving-degree/sage_gb_stream.sage experiments/ic-solving-degree/sage-in.json \
     experiments/ic-solving-degree/sage-out.jsonl         # Singular std, degrevlex; one JSON line per system,
                                                          # per-system alarm, resumable (re-run skips finished ids)
python3.13 scripts/check_ic_solving_degree.py polybori-export     # polybori-all-in.json (not committed; regenerable)
sage experiments/ic-solving-degree/polybori_gb.sage experiments/ic-solving-degree/polybori-all-in.json \
     experiments/ic-solving-degree/polybori-all-out.jsonl # PolyBoRi, deglex
python3.13 scripts/check_ic_solving_degree.py polybori-compare    # checks-polybori-all.json (C1b)
python3.13 scripts/check_ic_solving_degree.py compare     # checks.json (C1-C4); reads sage-out.json[l] and sage-shard*-out.jsonl
```

The Singular run was split into shards (`sage-shard*-in.json` / `-out.jsonl`) to run in parallel; `compare` merges
them. Sage needs to run outside the sandbox with `TMPDIR`, `DOT_SAGE` and `MPLCONFIGDIR` pointed at a writable
directory and `PYTHONNOUSERSITE=1`.

## Results (DRL; see `summary.json`)

- **Settled.** The solving degree is settled on all 5,813 systems with an exact LFD (three of them at
  k = 14, whose upper bound 7 is lifted; `DEEP_LFD` in `scripts/run_ic_solving_degree.py`).
  - 2,496 settled by `unit`, 2,467 by `count`, 850 by computation (corrected 18 November 2026, see
    *Erratum* below; previously 1,879 and 3,084).
  - `sd = LFD` on 5,805.
  - Five of the 850 computed systems are labelled `beyond_cap` (see below).
- **`sd > LFD` on 8 systems.** The same 8 in both orders:
  - Two binary sweep systems have `LFD = 0` (the input is already closed), so `sd = max.GB.deg = 2`.
  - Six binary two-summand descents (`k = 8, 10`; 129–197 solutions) have `LFD = 3` and `max.GB.deg = 4`, so
    `sd = 4`.
- **`max.GB.deg − LFD` (DRL), over all 3,141 systems with a computed degree.** −6: 3, −5: 65, −4: 116, −3: 274,
  −2: 939, −1: 1447, 0: 289, +1: 6, +2: 2.
  DRL and grlex give different `max.GB.deg` on 325 systems, but `sd_DRL = sd_grlex` on every system.
- **Deep prime instances.**
  - `max.GB.deg_DRL = 3` at k=13 (594 solutions, LFD 6).
  - `max.GB.deg_DRL = 4` at k=14 (2370 solutions, LFD 7).
  - Unit instances have degree 0.
- **Beyond the frozen cap.** 5 binary-deep systems have 4004–8234 solutions, above the frozen cap of 3000.
  - They were computed in a separate run after the freeze (`rows-s6-cap9000.jsonl`, 19 minutes) and are labelled
    `beyond_cap` in `rows-merged.jsonl`.
  - Results: `max.GB.deg_DRL` is 5, 5, 5, 5, 6, against LFD 6, 6, 6, 6, 7.
  - Grlex: 5, 5, 6, 6, 7.
  - So `sd = LFD` on all five.
- **Tag checks.**
  - Tag violations (a `count`/`cube` tag contradicted by a computed degree): 0.
  - The Boolean zero sets from the zeta transform equal field enumeration on 1355/1355 systems.

## Checks

See `checks.json`; the results are summarized in the paper's Appendix `app:checks`.

- **C1.** Singular reduced DRL basis vs Buchberger–Möller: 182/182 agree in max degree, leading monomials and
  dimension = |Z|. The sample was 192 systems: 150 seeded sweep systems, the 32 at depth 8/9, the 8 with `sd > LFD`
  and the two proper k=13 instances. Not finished (`not_run` in `checks.json`):
  - 8 depth-9 prime systems (r=18), 20-minute limit each;
  - the 2 proper k=13 instances (r=26), 40-minute limit each;

  The 6 binary m=2 systems with `sd > LFD` (r=16, 20) were not reached by the first run. They were run separately
  (`sage-bin*-in.json` / `-out.jsonl`, at most 20 s each) and agree: `max.GB.deg = 4`.

  So the prime degrees at depth 9 and at k=13, 14 have no second engine.
- **C1b.** PolyBoRi (`BooleanPolynomialRing`, deglex) vs grlex Buchberger–Möller: 2,352/2,352 binary systems,
  including all 8 with `sd > LFD`. The Boolean ring hides the field equations, so `u_v^2` is added for each standard
  variable before comparing.
- **C2.** An independent dense mutant closure gives `sd_DRL` directly from CG Def. 1.1, compared with
  `max(LFD, max.GB.deg)`: 40/40 systems with r ≤ 8. The first version failed on 29/40 because it ordered the Macaulay
  columns by increasing degree; they must be in decreasing degree.
- **C3.** Brute-force `Z` on 100 prime sweep systems: 100/100.
- **C4.** `|Z| ≤ 2^{k+1}` on every two-summand system: 1,932 checked, 0 violations.

## Files

| File | Contents |
|---|---|
| `rows.jsonl` | one row per system: set, id, family, m, k, r, points, LFD, tag, `gb_drl`, `gb_grlex`, `hf_drl`, `sd_*` |
| `rows-s3-binary.jsonl` | the nine binary m-summand systems rebuilt from parameters (system hash asserted) |
| `rows-merged.jsonl` | the merged table that the paper's numbers come from |
| `summary.json`, `checks.json` | aggregate results, check results |
| `sage_gb_stream.sage`, `sage-in.json`, `sage-out.jsonl`, `sage-shard*` | the C1 Singular run (seed 20261114) |
| `polybori_gb.sage`, `polybori-all-out.jsonl`, `checks-polybori-all.json` | the C1b PolyBoRi run |

## Erratum, 18 November 2026

Set S4 (direct m-summand prime systems, `experiments/ic-bounds-msum/results.jsonl`) passed the record's
`points` field as the solution count. That field is the grid size 2^(mk); the solution count is `zeros`.
With the grid size every S4 system was tagged `count`, including the 617 with no solutions. The driver now
passes `zeros` (`msum()` in `scripts/run_ic_solving_degree.py`); only the S4 rows were regenerated
(`--sets S4`) and spliced into `rows.jsonl` in generation order, and `summarize.py` was rerun. Changed fields:
`points`, `count_bound` on all 733 S4 rows and `tag` on 617 (`count` -> `unit`). No Groebner degree, LFD or
solving degree changed: S4 has no stored solution set, so none of its systems was computed, and the bound
|Z| <= 2^(mk) used before was valid. The table row is relabelled "Direct m-summand, prime, m = 2..8" because
the set includes 47 systems with m = 2.

