# Frozen quadratic root-lift solver screen

## Question and scope

Can the exact quadratic root-membership presentation beat native LRF on the
same two-summand relation queries? This is a bounded falsification screen,
not a relation-collection or scalar-recovery experiment. The construction and
its proof are in `../quadratic-root-lift-screen-20261125/SCREEN.md`.

The parent is `../subspace-remainder-cost-20261125/manifest.json`. Use its
curves, monomial factor bases and target indices 0, 1, 64 and 65 for cells
(23,11), (29,13) and (29,14). Add one deterministic diagonal target per cell
by doubling the first sorted factor-base point whose double has nonzero x.
The parent indices cover random and planted targets; the extra target covers
the repeated-root case. Freeze the exact points and hashes in `manifest.json`
before timing. Do not replace failed or empty queries.

## Arms and accounting

For every target use one fresh worker per arm. Charge field/curve and factor
base setup, affine LRF-system construction, CNF construction, Kissat invocations
and blocking clauses, root/sign recovery, exact group-law verification, and
canonicalization. Record wall, CPU, peak RSS, number of SAT models, result
status, and the full signed pair set. A hard 10 s wall cap per SAT target is
reported as censored; never count a partial pair set as complete. Native LRF
and conditional-image LRF use the repository C kernel with their own setup.
Direct enumeration is a correctness control when it finishes under the same
cap. Worker and solver process startup are charged. Include failed systems.

The SAT engine uses a Tseitin encoding of the n Boolean equations
`x(z)^2 + e1(g)*x(z) + e2(g) = 0`. It exposes no polynomial solving-degree
trace, so this screen **cannot** claim a reduction in observed maximal solving
degree. The input equations have degree at most two. A subsequent algebraic
solver gate is required if the SAT cost leaves a plausible advantage.

Exact relation-set parity and total charged query cost decide this screen.
Any end-to-end advantage requires a separate collection, linear algebra and
cold same-target scalar comparison; the already measured BSGS reference is
the current downstream gate. Keep all raw receipts and censored rows.
