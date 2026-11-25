# Frozen three-summand cost pilot (2026-11-25)

## Purpose and scope

Check whether the two-summand conditional-image component gain survives a complete three-summand relation query on the **same targets and factor bases**. This bounded pilot is a falsifier and timing guide for a larger relation-collection campaign; its small panel cannot establish a full index-calculus crossover. Preserve every timeout and empty relation set.

## Fixed inputs

Use the frozen binary curves and factor bases from `manifest.json` for (n,k)=(23,11),(29,13),(29,14). Use the first two random subgroup targets in each existing cell, unchanged. Generate one planted three-point target per cell with `random.Random(20263125+100*n+k)` and three independent draws from the sorted complete factor base, with replacement. The independent Python manifest generator records every point, witness, factor-base hash, source hash and this protocol hash **before timing**. Every arm must use those exact inputs. Repetitions: three per cell, rotating the native LRF/direct order; run the image arm last. Run all arms natively on the same ARM host, single thread, fresh worker processes. Compile the common C binary before timing and record code/compiler/host versions.

## Arms and accounting

Native LRF, direct three-point enumeration, and conditional-image LRF use the repository's `Solver` C kernels. Each worker must construct its own solver and factor-base table, and charge image-table construction in total. Record setup and each query wall/CPU, peak RSS, raw signed relation triples, all failures, and verification wall/CPU. The parent measures common input generation once and reports it separately and added equally to all arms. Worker startup is excluded equally. Charge independent sparse modular elimination over the random-target relation rows as a separate common linear-algebra phase; report coefficient rank, row and column counts, and cost. This phase is diagnostic if rank is insufficient to solve factor-base logs.

Verify every returned triple sums to its target and belongs to the fixed factor base, and verify the planted witness. Completed arms require identical canonical signed triple sets for **every** target. Use a 45-second worker wall cap for this pilot; after a timeout, preserve flushed per-target receipts, elapsed lower bound, and status `timeout_censored`. Do not treat a censored arm as a loss or a completed crossover. Report component wall and CPU as well as common-input/linear-algebra-charged totals. A completed image arm must beat both completed controls on wall and CPU in all three repetitions before advancing that cell to the larger campaign. Repeat close calls on new frozen seeds before attributing a gain.

The next campaign must collect enough independent relations for its declared linear algebra objective and compare full same-target cost, including failed queries and a matched reference method. A pilot crossover alone is not that result.
