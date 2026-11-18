The lift review found two reproducible boundary defects. This PR fixes both and records the results of all five requested local validation tracks. It stays draft because the full depth-14 certificate, split degree-8 result and depth-18 breakpoint remain open.

## Fixes

- Route the inner row reduction through the bounded exact product. At p=16777213, a 64-by-65 counterexample previously gave 16777148 instead of 16777149. Cap requested product chunks at the safe accumulation length.
- Handle an empty next-degree monomial layer at D=2k. The previous code raised IndexError on k=2, D=4, p=17; both proper and unit boundary cases now pass.
- Neither defect reaches the archived depth-14 cases: their primes are below the affected accumulation bound and D=6<2k.
- Correct the paper README's obsolete Lean audit checklist.

## Completed local evidence

- 32 Python tests pass, including 16 new regression tests. The four-prime arithmetic test and 16 fresh lift runs agree with direct reference and naive closures. Of 18 altered-input controls, 16 are detected; two redundant-generator removals leave the result unchanged.
- A pure-integer derivation DAG checker validates six small closures and a saved certificate. It requires an external target digest and rejects substituted output sets. This is a prototype, not a certificate for the actual depth-14 basis.
- Exact E6 proper-ideal values are 5 at (k,p)=(7,257),(7,769),(8,7681),(9,12289), and 6 at (10,12289). The first and last agree between two closure implementations. An independent involution enumeration checks all five complete zero sets.
- The tiled split engine agrees on 16 control closures. A fresh local ladder checks D3..6 with ranks 1,39,711,5243, then caps at D7 with derived rank 3644. This fresh ladder gives lower bound 7; the older log-supported bound 8 remains separate. D8 is unresolved.
- A closed depth-18 V4 has rank 3821 and excludes 1 on a lemma-backed unit instance, giving LFD at least 5. The predicted value 8 is untested.
- All 439 existing submission artifact consistency checks still pass. The submission PDF is unchanged.

## Candidate proper-ideal proof route

The explicit ideal element h=xy-x-y-8 lies outside the completed lower closures on all five new profiles. A first-derivative commutator passed 61,888 monomial checks. GrowthJetPrototype.lean locally checks two elementary cofactor/commutator implications with only propext. The concrete FEI/growth bridge remains unproved; this is not a new asymptotic theorem and is kept outside the submitted paper.

## Resource and claim boundaries

The owner selected local-only. No paid worker is approved or launched. Initial cell-budget caps, the 240-second reference cap, the three-minute split calibration, and the final local split cap are retained. Full depth-14 derivations and a feasible depth-18 D7/D8 method remain separate obligations.

See README.md, PROTOCOL.md, REVIEW.md and PROPER-IDEAL-ROUTE.md in this directory for source-bound records, chronology, falsifiers and remaining work. Final CI must check this PR's current head, including every Lean module and printed axiom list.
