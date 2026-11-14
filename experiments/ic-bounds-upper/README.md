# ic-bounds-upper: upper bounds for the two-summand prime chain LFD

Experiments for `docs/proofs/U1-upper-bounds.md` (2026-11-10). Everything is exact arithmetic
modulo p in the chain quotient (group algebra of (Z/2^k)^2, degree = Hamming weight); the engine
`chain_algebra.py` is self-contained (numpy only) and independent of `csd/`.
Run from this directory with `ulimit -v 1500000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

| script | what it does | output | time |
|---|---|---|---|
| `selector_bound_check.py` | executes the proof of Theorem 1 (LFD <= min(k1,k2)+2) with an enforced degree budget and checks m − π(m) ∈ V_max(k1+2, wt m) for every monomial m | `results-selector.jsonl`, `results-selector-summary.json`, `selector.log` | 6 min |
| `lfd_structure.py small` | exact LFD by the filtered-generation criterion (Prop. 4): closures V_D, ideal parts I_D, generation I_{D'+1} = I_{D'} + Σ v I_{D'}; exhaustive planes at k = 1, 2, 3, random at k = 4, 5 | `results-structure-small.jsonl`, `-summary.json` | 12 min |
| `lfd_structure.py sweep --k6 6` | the same on the subgroup two-summand cases of `experiments/ic-curve-sweep` (all with k <= 5, six at k = 6), compared with the recorded exact LFD | `results-structure-sweep.jsonl`, `-summary.json` | 16 min |
| `depth_restriction.py` | Prop. 3: π(V_D(k+1)) ⊆ V_D(k); dimension of the image, of its filtered pieces and of the low-level part, D = 3, 4, 5 | `results-depth-restriction.json`, `depth-restriction.log` | 13 min |
| `lift_of_one.py` | least degree of an element of V_4(6) restricting to 1 at depth 5 | `results-lift-of-one.json` | 1 min |
| `sweep_monotone.py` | LFD nondecreasing in the depth on the certified sweep values | `results-sweep-monotone.json` | seconds |
| `slope_counts.py` | free-extension polynomials Q_D, f_D against \|M_<=D\|; evaluation ranks and the counting inequality of Prop. 5 on Z-heavy instances; thresholds D*(k) | `results-slope-counts.json` | 3 min |

Labels: every number in the outputs is **exact** for the instance computed. A closure is reported
only after a final pass in which every low-degree basis row has been shifted in its final form
(the loop of `chain_algebra.closure` ends exactly then), so the returned space is closed.
In `lfd_structure.py` the value `lfd = 2` means "at most 2" (the bridge to the ordinary ring needs
D >= 2). In `depth_restriction.py` the instance labels `unit-…`/`proper-…` refer to the shallow
depths only; the fields `points_k`, `points_k1` give |Z| at each depth.
