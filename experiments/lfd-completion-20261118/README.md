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
Ten replay/refusal/custody tests pass across the Python and Rust checkers.

## Bounded experiment outcomes

- D8 on the same fixed split instance: stopped at the free-disk floor
  before its 600-second cap. Last progress record: derived rank 82556,
  13,887,854 stored terms, 442.5 seconds. These are partial-span counts,
  not V8 dimensions. No checked unit witness was obtained.
- Depth-10 D5 certificate-scaling probe: hit its 180-second cap at rank
  12287. The separate Rust verifier replayed all retained rows, without a
  closedness claim. This is a scaling probe, not depth 14.
- The runner also stops a child if ledger+basis exceed 2 GiB or free disk
  falls below 2 GiB. External stops retain partial files and make no closure
  claim. These guards are not permission to discard a censored result.

Check actual processes and result.json files before resuming; do not
restart because a tool observation timed out. Former handles: D8 unified session 28633, depth10 session 8269. Both
are terminal. Successful completion is
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

## Lossless raw archival

After all processes exited, archive_raw.py verified every original raw
SHA256, compressed the ledger and basis, decoded each archive to verify
the same SHA256, then removed only the duplicate uncompressed copy.
Receipts retain both the decoded and compressed identities. Raw .gz files
remain local under runs/. This is lossless storage reduction, not a new
certificate check or experiment completion. Future engines should stream
compressed output to avoid the old uncompressed disk peak.

The runner now refuses to overwrite any existing run directory. It captures
source hashes before launching each child. runner-resource-guard.py.txt
preserves the exact prior runner used for D6/D7/D8 and the depth10 probe;
initial-runner.py.txt remains the source of the first gate suite.

## Heap accumulator and streaming compressed proofs

`heap_producer.rs` retains the same immutable row schedule and proof format,
using a reusable hash accumulator and a lazy leading-term heap. On the
recorded p=114689 split degree-7 input it produced byte-identical ledger
and basis files to `sparse_producer.rs`: 27,838 rows and 106,110 products.
The local comparison measured 27.36 seconds for the baseline and 11.67
seconds for the heap version. These are single local timings, not an
end-to-end solver comparison. `results/heap-D7-comparison.json` binds the
outputs and sources; `heap_accumulator_only.rs.txt` preserves that exact
pre-compression source.

The current producer streams `.gz` output through gzip and waits for both
compressors before reporting success. `run_heap.py INPUT TAG --seconds N`
compiles a source-bound binary into a fresh run directory, records compressed
and decoded SHA256 values separately, and enforces the existing disk floor,
artifact size and external time guards. It changes only the input time budget.
The independent Rust replayer accepts these compressed files directly and
requires successful gzip decoding, including the end-of-stream checksum.

The 21 local tests cover plain and compressed production, exact agreement
with the original producer on both chain and split systems, independent
replay, forged rows, invalid degree shifts, future references, truncated
gzip files and overwrite refusal. Production and full replay of the degree-7
compressed certificate have separate receipts. This engineering work does
not supply the outstanding degree-8 witness, depth-14 input bindings or
depth-18 breakpoint.

The same refusal suite is wired into the Python CI job. Local success is
recorded separately from the CI result on the pushed commit. The second
planned prime, p=163841, now also has a complete degree-7 certificate:
27,838 rows, independently replayed with full closedness and exclusion of 1.
The frozen seed-0 input identity is in split-p163841-instance-identity.json.

## Indexed elimination

`indexed_producer.rs` enumerates the complete degree-bounded monomial domain,
orders it identically to the original producer, and stores sparse basis
rows as `(column_index, coefficient)` pairs. Reductions update a reusable
coefficient array and nonzero bitmap directly. Only input and shift encoding
need a hash lookup. Its degree-7 comparison completed in 4.30 seconds and
matched both previously checked proof files byte for byte; see
`results/indexed-D7-comparison.json`. This timing excludes domain preflight.
It is a separate local measurement, not a matched throughput guarantee.

The engine refuses domains above two million columns before opening proof
outputs. This is an explicit engine resource limit: it does not truncate
the polynomial space and does not resolve the depth-18 requirement. The
original heap engine remains available for sparse exploration beyond that
limit. `run_heap.py --engine indexed` selects the new engine and captures
its source identity. `run_heap-v1.py.txt` preserves the runner already loaded
by the first compressed degree-8 job.

Both engines run through the same byte-agreement and replay tests. There
are now 29 distinct local tests, including an oversized-domain refusal
check. The p=163841 degree-8 computation uses the indexed engine while the
already running p=114689 heap computation retains its original sources.
