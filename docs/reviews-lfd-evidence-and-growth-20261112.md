# LFD manuscript: evidence corrections and proof obligations

This review companion records changes required before upgrading claims. It does **not** change experimental data or assert new theorems.

## Evidence labels

- **Kernel-checked:** a named Lean declaration, exact project revision, and audited axioms.
- **Independently certified computation:** exact arithmetic plus a separately replayable witness/checker.
- **Written theorem:** complete mathematical proof with all hypotheses, but no kernel check.
- **Computed observation:** exact-arithmetic run with insufficient independent verification for the claimed certificate.

Do not call a lift-produced value independently certified solely because the implementation reports full rank.

## Two-summand upper bound

The manuscript's Theorem 9 asserts LFD <= k+2 for nonsingular prime-field two-summand chain systems. Before marking it kernel-checked, verify the exact declaration on the target commit and its hypotheses. The PDF's 513-declaration count is a snapshot, not a live project count. Regenerate the audit and PDF together.

Proof checklist:
1. Character idempotents of the first block have filtered degree at most k.
2. Fiber restrictions of the normalized S3 generator have degree at most two.
3. Every fiber ideal element admits a derivation within the asserted filtered budget, including constant, linear, repeated-root and no-root fibers.
4. Fiber reconstruction and the quotient-to-ordinary closure bridge preserve the budget.
5. The local-to-global criterion covers all ideal elements, not just the generator.

## Depth fourteen: separate lower and upper evidence

For the selected k=14 unit instances, a verified exclusion 1 not in V6 proves LFD >= 7. The written selector bound gives LFD <= 16. The reported equality LFD=7 additionally relies on the degree-seven lift, for which the manuscript discloses no independent certificate.

A **smaller target certificate** than the entire V7 basis is a derivation of 1 in V7. A verifier should check:
- the exact input generator and chain quotient conventions;
- every source row and degree-budgeted variable shift;
- each linear combination over the specified prime field;
- the final equality to the constant polynomial 1;
- compatibility with the ordinary-ring closure bridge.

The checker must not trust the lift's claimed rank, pivot partition, or fall counts. If a compact witness cannot be produced, retain the conditional status of the upper-side result.

## Multi-summand regimes

Distinguish **k <= m-1**, where degree mk of the normalized generator can make the general mk+1 bound sharp, from **k >= m**, where generator degree is bounded by m(m-1) and the same sharpness argument does not apply. The m=3, k=6 observation is LFD >= 8 on one instance, not an exact value.

## Genericity and growth

The proposed floor((k+14)/4) law concerns a generic-profile family, not all curves, targets, ideals or finite fields. Exceptional instances are documented in the manuscript.

At each fixed (k,D), symbolic matrix minors can define a rank-constant Zariski-open stratum; prove nonemptiness on the actual nonsingular curve parameter family and address finite-field specializations separately. Do not infer a depth-uniform exceptional locus from fixed-depth rank statements.

Define the filtered profile H_{D,k}(z) from dimensions of V_D(k) intersected with degree <= j. Under a **proved** free filtered extension, its profile transforms by truncated multiplication with (1+z)^2. Observed equality of dimensions does not establish that the actual closure is such an extension: the depth-10-to-11 counterexample in the manuscript must remain explicit.

The next conjectural transition to test on an appropriate unit-ideal family is 1 in V7(17) and 1 not in V7(18). No claim is made that these computations or proofs have been completed.

## Release checks

1. Verify exact Lean head, theorem statements, dependency audit and declaration count.
2. Compare manuscript tables and figures against immutable experiment records.
3. Regenerate paper and inspect equations, footnotes, captions and evidence labels.
4. Replay independent checks for any claim promoted to certified.
5. Do not modify existing experimental records to match the narrative.
