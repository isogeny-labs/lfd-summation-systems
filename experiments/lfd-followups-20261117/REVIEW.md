# Internal review of the independent lift

Reviewed source: pre-fix 8a6a32466677fc11fda4bdf832fd823ff8fef1a7335ebe7dc25cb31a6db62e7a.
Reviewer: Codex in this task; this is not a separate external review.

## Confirmed arithmetic defect and fix

`rref` used a raw matrix product to eliminate the pivots already found
inside its 64-row block. Other multi-term products use `mm`, which limits
the accumulation length to keep every integer sum exactly representable
in float64. The raw product could exceed that limit near p=2^24.

Reproducer p=16777213: a 64 by 65 matrix, identity in its first 63 pivot
columns, p-1 in the last column of those rows, and final row
(p-1 repeated 63 times, 1, p-1). Exact arithmetic gives final free entry
16777149; the pre-fix code gives 16777148. This is a correctness bug,
not a timing difference. It is exercised by the new CI regression.

The fix routes the internal reduction through `mm`. A caller's requested
positive product chunk is capped at the proven safe inner dimension;
nonpositive chunks are refused. Capping preserves the driver's existing
arithmetic test, which intentionally requests chunks beyond that bound.
The four-prime arithmetic self-test passes after the change.

At the archived depth-14 primes, 63*(p-1)^2 is below 2^53-2p, so this
specific defect does not change those runs' arithmetic. Their source hashes
remain historical identifiers and are not rewritten to the corrected hash.

## Other reviewed invariants

- The two orbit parts determine ordinary pivot rows and single rows. The
  self-paired case uses only the plus part. Both-part rows combine using
  the inverse of 2, and the sign changes under swapping are explicit.
- Special tops are exactly those reached from N. A regular top must have
  D+1 pivot divisors; the constructor checks this count. Subtracting the
  corresponding generators cancels that regular coordinate.
- `chunk` reduces candidates against a snapshot, then reduces against any
  rows added since that snapshot before extending K.
- A fall pivot lies in N, so its row is supported there. Shifting each
  changing row once uses the later-row invariant; the final deterministic
  fall-shift pass tests closure when the run is not saturated.
- Saturation is a sound early exit only under the input-containment
  assumption. Random low-row checks do not certify that assumption.
- `include_low` avoids relying on the sampled skipped-generator check;
  the bounded replay exercises both modes and different scheduling seeds.
- Stored float32 residues are exact only below 2^24. Scalar products of
  two residues fit in float64; multi-term products require chunking.

No additional defect was identified in this pass. These checks do not
prove the algorithm, and do not supply a certificate for the depth-14
producer. Input-containment, full input provenance and source-bound output
identities remain explicit requirements for that follow-up.
