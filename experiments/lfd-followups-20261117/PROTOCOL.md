# Five LFD follow-ups, 17 November 2026

Tracking: PR #31. This is a follow-up to merged evidence snapshot ce4d785.
The user authorized all five tracks. Paid compute has not been launched.

## Scope and chronology

The certificate prototype, the high-prime arithmetic counterexample, the
initial exceptional-prime probes, the depth-18 degree-4 probe and the
three-minute split retry are exploratory results obtained before this
protocol was written. They are not described as preregistered tests.
The original split prediction (8) and depth-profile law predate this work.
The prospective full split retry below uses that original prediction.

## 1. Derivation certificates

`derivations.py` produces a DAG with normal-form inputs, linear combinations
and variable shifts. `checker.py` uses Python integers and imports no
producer or repository algebra. It binds the parameters and input list,
checks each edge, checks the source degree BEFORE a shift, reconstructs
chain carries independently, and checks the output rows and their digest.
This certifies containment of the provided rows in the normalized chain
closure. It does not certify unprovided rows or raw inputs whose unreduced
degree exceeds the declared presentation's input degree.

Compare the derived independent rows with a direct closure and check
malformed references, inputs, outputs and overdegree carries. The saved
small certificate is replayable without trusting the producer.

The prototype is NOT a depth-14 certificate. The next scaling obligation is
to emit proof operations from the tiled producer, map its orbit basis to
ordinary rows, and bind the checked outputs to the actual large basis.
Both depth-14 inputs require about 42 GiB each just for the archived-style
basis representation; the derivation trace has not been sized. The local
workspace has under 2 GiB free. No full-basis proof run is authorized by the
separate US$100 split-retry question.

## 2. Lift review

Inspect coordinate conversion, regular-top elimination, fall propagation,
fixed-point/saturation exits, degree restrictions and arithmetic bounds.
The high-prime rref counterexample and its correction are recorded in
REVIEW.md. Retain historical run source hashes unchanged; they identify
the older source. Test the corrected implementation against direct small
closures and altered-input controls. An internal review is not an external
review or a Lean proof of the numerical algorithm.

## 3. Prospective full split retry

Fixed instance: m=5, k=3, p=114689, seed=0, frozen sampler from
`scripts/run_ic_split_m.py`. Reconstructed a=52814, b=60265,
target=(52173,107269), normalized system SHA-256
723464b8d858065cf28bbaea733e70fe7c41d6db528d6d4f74dc159c1ee3f627.
The previous log lacks that hash. Therefore first recompute degrees 3..7
and require the recorded ranks 1,39,711,5243,27838, checked closure and
exclusion of 1. This is a fresh run, not checkpoint resumption.

Then compute degree 8 with the tested tiled float32-storage engine. Reading:

- all lower exclusions and 1 derived at degree 8: exact LFD 8;
- degree 8 closed without 1: lower bound 9, refuting the prediction here;
- partial degree 8: retain the verified lower bound 8, upper bound unknown;
- any mismatch in the lower ladder: stop and investigate identity/engine.

Candidate command, only after explicit paid-resource approval:

```
OPENBLAS_NUM_THREADS=32 OMP_NUM_THREADS=32 python3 \
  experiments/lfd-followups-20261117/split_retry.py \
  --large --ladder --seconds 108000 --max-basis-gib 200
```

Request: at least 256 GiB host RAM, at most 30 hours, at most US$100 TOTAL,
including storage. Confirm the live price before reservation; shorten the
time limit if needed to preserve the spend ceiling. One worker and one
instance only. Mirror each completed stage and raw logs during the job,
not only at the end; delete the worker on success, failure or either cap.
No guarantee of completion: the new engine trades memory for more tiled
work. The previous attempt exhausted memory after about 8.4 hours.

Measured local calibration: degree 8 reached rank 1677 in 183.433 seconds
and hit its deadline. It did not complete a closure. Its original runner
source is preserved under `calibration/` with the recorded hash. The
16-case final control suite must pass before the proposed worker runs.

## 4. Proper-ideal exceptional-prime probes

E6: y^2=x^3+6x-7, target (1,0), subgroup H of order 2^k. These are
exceptional admissible primes because the grid has a zero; not every case
is the Fermat-factor diagonal subfamily. Record the complete zero set,
closed spaces, last strict fall, and evaluation-rank/codimension equality
for an exact local-to-global upper bound. Do not infer an all-depth bound.

Cases: (7,257), (7,769), (8,7681), (9,12289), (10,12289).
Preserve initial 30-million-cell caps and the 240-second depth-10 cap.
The larger compressed-reference depth-10 attempt has an 1800-second cap.
A capped partial row space proves no exclusion or upper bound by itself.

## 5. Depth-18 breakpoint

Frozen exploratory instance: p=7340033, curve y^2=x^3+x,
target (3,2205181). Check primality, rationality, Y nonzero, X nonsquare
and 2^19 | p-1. Paper Lemma 'Unit instances' then establishes the unit
ideal; grid enumeration of 2^36 pairs was not performed.

The law predicts LFD 8. A genuine test requires a checked exclusion at
D=7 and derivation of 1 at D=8 (or a contrary certified bound). The degree-4
probe is only an early lower-bound/implementation calibration. Do not
label it confirmation of the depth-18 prediction.

Direct column counts at D=4..8: 66712,443704,2391496,10739176,40999516.
The compressed swap-split worst-case matrix sizes are approximately
2.1 GiB,91.7 GiB,2.6 TiB,52.4 TiB,764.4 TiB, excluding buffers.
These are worst-case storage formulas, not measured usage or runtimes.
A structured method with a defensible resource bound is required before
requesting a paid degree-7/8 run. No such run is included in this task's
split-worker budget question.
