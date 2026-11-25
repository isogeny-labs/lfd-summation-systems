# Frozen binary subspace-remainder cost gate (2026-11-25)

## Question

On identical binary two-summand relation systems, can the exact
subspace-polynomial remainder condition locate all signed factor-base point
pairs for less total work than native LRF's kernel-candidate scan? This is a
component falsifier for the broader end-to-end relation-finding goal. A win
here is not an end-to-end index-calculus or three-summand speedup.

## Fixture and target panel

Use ordinary curves E/F_(2^n): y²+xy=x³+b, with n,k in
{(23,8),(23,10),(23,11),(29,12),(29,13),(29,14)} and the irreducible
moduli x^23+x^5+1 and x^29+x^2+1, respectively, as used by
`experiments/ic-rootfind-yield/goodset_structure.py`.
Freeze b=894692 for n=23 and b=326057084 for n=29. Their respective group
orders are 8390136=24·349589 and 536882624=64·8388791, checked by the
repository's group-order routine before measurement. The prime subgroup
orders q are 349589 and 8388791, and cofactors h are 24 and 64. Let G be
the first nonzero [h]P found by ascending x and `BinaryCurve.lift_x`; assert
[q]G=O, and set Q=[113]G. V=span_F2(1,w,...,w^(k-1)); the factor base is every
rational affine curve point with x in V, with both signs when distinct.
For each (n,k), generate 64 subgroup targets T=uG+vQ, with u,v drawn from a
local Python PRNG seeded 20261125+100n+k, and 16 planted sums of two factor-base points
using seed 20262125+100n+k. Preserve duplicate targets and every failed
attempt. Target generation and common factor-base construction are reported
and charged equally to all arms. Record the exact b, modulus, point panel,
factor-base size, and source/protocol hashes before timing. Only the random
subgroup targets enter relation-rank calculation. For their verified pairs,
cofactor-project the factor-base points, identify P and -P with opposite
signs, omit columns for [h]P=O, and form rows over F_q with known right-hand
side -h(u+113v). Report the resulting matrix rank, including zero-rank cases.

Execution environment gate: the repository C kernel uses x86 intrinsics and
cannot load into the local ARM Python process. Run every timed arm on the
same native x86_64 Linux host or pinned x86_64 container; record CPU model,
OS, compiler, Python, solver versions, and container digest where applicable
before the first timing. Local Rosetta compilation is an untimed diagnostic,
not a comparable benchmark receipt.

## Arms and correctness

1. Native LRF C kernel, including its table setup, linear solve, all
   kernel-candidate tests, root/sign recovery, and group validation.
2. Native C direct point enumeration on the same base and targets. This is a
   correctness and competing-time control; if a reusable direct pair-sum
   table is built, charge its setup and include it as an additional control.
3. A compact remainder-constraint solver: construct L_V and its recurrence
   circuit, solve the exact symmetric image equations from PR #63 (including
   e1=0), recover signed points, and group-validate every pair. Charge all
   circuit/solver setup, failed targets, retries, and recovery.

Each arm must produce the identical canonical signed point-pair set for
every target. Failed/empty targets count. Record setup and query wall time,
total process plus child CPU time, peak memory, result counts, solver status,
and raw per-target receipts. Run three repetitions with rotated direct-arm
order and the candidate after the direct arms. A candidate may stop as
`dominated_censored` only after both its elapsed wall and CPU time exceed the
faster completed direct-arm full-batch totals; retain the processed prefix
and lower-bound timings. Completed candidate runs require full equality.

The initial cost gate passes only if the candidate has identical relations
and lower total wall and CPU time than *both* direct arms in all three
repetitions for at least one (n,k). Report independent signed relation rank
and cost per rank only when rank is positive. A nonpassing or censored arm
cannot support a speedup claim. If a case passes, freeze a separate full
three-summand same-target campaign including factor-base enumeration and
linear algebra before claiming end-to-end gain. Solver degree alone is a
diagnostic. No extrapolation to n=131 or other curves from this panel.
