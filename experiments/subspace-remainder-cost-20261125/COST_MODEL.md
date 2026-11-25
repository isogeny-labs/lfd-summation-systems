# Why the component win does not imply a single-target crossover

Let N=#E(F_(2^n)), let B be the number of **signed** points in the fixed factor base, and consider unordered three-point multisets. Under a uniform-sum heuristic, the expected number of triples for a random target is

    lambda = binom(B+2,3)/N ≈ B³/(6N).

This is a yield estimate, not an independence theorem. The measured random-target counts in the frozen three-summand pilot are consistent with it:

| (n,k) | B | Uniform-sum lambda | Pilot random-target counts |
| --- | ---: | ---: | --- |
| (23,11) | 2,087 | 180.8 | 178, 180 |
| (29,13) | 8,193 | 170.8 | 172, 187 |
| (29,14) | 16,337 | 1,353.8 | 1,372, 1,381 |

For fixed bounded cofactor and a sign-quotiented factor base with c=Theta(B) unknown subgroup logarithms, at least c independent relation rows are required to determine all of them. Even if **every** triple row is independent, the uniform heuristic therefore asks for at least Omega(B/lambda)=Omega(N/B²) known-scalar target queries. The current m=3 LRF and conditional-image kernels enumerate all B first factor-base points on each query, giving Omega(B) work per query before considering kernel candidates, root recovery, memory, sparse elimination, failed queries, or target resolution. Their relation-query work is consequently Omega(N/B) under these assumptions.

If B is at most on the order of sqrt(N), that expression is at least on the order of sqrt(N). If B exceeds sqrt(N), the **first** query's B-point scan is already above that scale. With subgroup order q comparable to N, a generic square-root method therefore sets a difficult cold-cost floor for this fixed-arity, full-scan architecture that solves the whole factor-base system. This conditional argument does not prove a lower bound for arbitrary relation methods: highly nonuniform yields, sublinear indexed queries, target-specific resolution with only a partial row span, a smaller effective set of logarithms, a different group/factor-base construction, or growing arity can change the premises. Setup reuse across targets also changes an amortized question, but not a cold single-target comparison.

The [cold matched receipts](COLD_TARGET_RESULTS.md) supply a finite-workload check of the cost issue: a Python generic solver finishes before the native image method's first verified query on every target tested. A future formulation intended to lower solver degree must also change at least one **total-cost** premise above. Lower degree alone is insufficient if it leaves the B-point scan and relation yield intact.
