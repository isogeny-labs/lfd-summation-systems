# Frozen cold same-target lower-bound gate (2026-11-25)

## Question

Does the conditional-image three-summand method's **first verified target query**, including factor-base and solver setup, already cost more than a generic scalar-recovery method on the identical target? If so, a complete relation collection through this oracle cannot beat that generic method in these exact cells. This is a bounded cost falsifier, not a proof about all curves, factor bases, or future relation algorithms.

## Inputs and arms

Use the exact ordinary binary curves, prime subgroup orders, generators, and factor bases in the existing frozen manifest for (n,k)=(23,11),(29,13),(29,14). A separate manifest frozen before measurement chooses three independent nonzero secret scalars from the upper half of each prime subgroup, using Python PRNG seed 20264125+100n+i for i=0,1,2. The n=29 cells share each target Q=sG. Record each Q and scalar for post-run validation; solvers receive only Q, G, q, and curve parameters. Run three fresh processes per target per arm, rotating arm order, on the same native ARM host. Exclude Python process startup equally; include all work from field/curve construction through verification. Record source hashes, host, compiler and Python versions.

1. Generic baby-step–giant-step using the same Python `BinaryCurve` group law: build the full baby table, search giant steps, recover s mod q, and independently verify [s]G=Q. Charge its memory and all setup.
2. Conditional-image LRF: construct `Solver` and its complete factor base, verify its hash, run one m=3 query on Q, and validate every signed triple against the fixed base and group law. Charge image-table construction, all query work, verification and memory. Record the full relation set and empty outcomes.

The generic method is a conservative reference because its group operations run in Python while the image solver uses native C. For each matched target, compare full process wall and CPU, including failed or censored runs. A 45-second worker cap applies; preserve timeout lower bounds. If the image arm is slower on **every** repetition, report a cold-query lower-bound failure for that cell; do not infer that an unmeasured complete IC solver would be faster. If the image arm is faster, a full relation-collection and linear-algebra campaign remains required. Common target-generation work is performed before arm timing and is the same for both, with its measured cost reported separately.
