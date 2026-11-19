# Complete the three LFD computational obligations

Active full objective: depth-14 input derivation certificate, exact split
five-summand depth-3 degree-8 computation, and the depth-18 breakpoint.
Tracking PR #35. All work is local; no paid resources are authorized.
This is not complete until each full requirement is verified.

## First implementation and evidence

sparse_producer.rs keeps an immutable sparse echelon basis and a queue of
admissible variable shifts. Every inserted pivot logs its source input or
shift, coefficients of prior reducer rows, and normalization. Ordinary
auxiliary exponents and cyclic chain exponents are encoded separately.
The output is a derivation ledger plus the actual produced basis, not
just a dimension or producer status.

replay.py is a separate integer checker. replay.rs is a separate Rust
replayer for larger outputs. Neither imports producer algebra. They
check source degrees before reduction, source/reducer references, every
row identity, normalization and pivot independence. A closedness pass
checks input containment and every admissible shift of the final basis.
A partial basis passing replay proves only membership of those rows.

The fixed split instance over p=114689 (a=52814,b=60265,
target=(52173,107269)) completed D6 with rank 5243 and passed independent
Python replay and closedness. D7 completed with rank 27838; the separate
Rust replay checked all derivations and all admissible shifts, excluding 1.
This supplies a fresh certificate-backed lower bound of 8 for that
instance, beyond the previous log-supported record. It does not prove
membership at D8. The producer used 29.43732 s at D7; replay/closedness
used 50.432899 s. These are local timings, not controlled speedup claims.

The initial gate has ten cases: nine complete and one capped. The unit
control m=5,k=2,D=6 hit its 60-second cap at derived rank 4307; its row
membership replay passed, but its exact closure and predicted unit result
were not established by that attempt. A capped gate is not a passing
complete gate. Other completed gates agree with direct reference closures.
Five malformed-certificate refusal controls pass.

## Live bounded experiments

- D8 on the same fixed split instance: 600-second producer cap, 230000
  rows, 100 million stored terms. This still needs a checked unit witness,
  plus the already checked lower ladder, before calling the exact value 8.
- Depth-10 D5 certificate-scaling probe: hit its 180-second cap at rank
  12287. The separate Rust verifier replayed all retained rows, without a
  closedness claim. This is a scaling probe, not depth 14.
- The runner also stops a child if ledger+basis exceed 2 GiB or free disk
  falls below 2 GiB. External stops retain partial files and make no closure
  claim. These guards are not permission to discard a censored result.

Check actual processes and result.json files before resuming; do not
restart because a tool observation timed out. Initial live handles:
D8 unified session 28633, depth10 session 8269. Successful completion is
not inferred from a stale PID or a log. The results and raw ledger files
are under runs/<tag>/; the small receipts and exact input files are copied
to results/ and inputs/ for version control. Large raw files stay local
and ignored, with their SHA256 values retained in the receipts.

## Provenance

The first gate started before the runner gained the disk watchdog. Its
executed source is preserved as initial-runner.py.txt. Later gate receipts
originally hashed the newer on-disk runner at the end; the receipts now
record that correction explicitly, preserving the previous observed hash.
The producer and Python checker did not change during those runs.

The D7 independent replay was bound to the original input, ledger and
basis hashes before accepting the verifier result. It does not certify
any historical depth-14 W. That requires matching the actual input basis
or a fully checked replacement together with the lift that consumes it.

## Next required work

1. Read the live D8 and depth10 outcomes and preserve every cap/failure.
2. Reduce elimination and replay costs before scaling to depth14 V6.
3. Bind a completed depth14 proof to the actual lift input, including
   its canonical basis identity and source instance. Rank agreement alone
   is insufficient for the requested input certificate.
4. Establish a checked D7 exclusion and D8 witness at depth18; a smaller
   degree closure or theoretical growth result does not decide that cell.

Reproduce refusal controls with:

```
python3 -m pytest -q experiments/lfd-completion-20261118/test_replay.py
```

The committed sample is a complete small certificate. The larger local
raw files are not claimed to be archived remotely. The full objective
remains active.

completion-requirements.json retains the full target identities and proof
requirements. Both recorded depth-14 unit inputs and both originally
planned split instances remain in scope; partial gates do not replace them.
