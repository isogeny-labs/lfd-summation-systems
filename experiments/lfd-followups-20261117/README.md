# LFD validation follow-ups, 17 November 2026

Draft PR #31 tracks all five requested follow-ups. This directory contains
new experiments and an internal code review; it does not replace the
submission package merged in #24. No paid worker has been launched.

## Results so far

| Track | Completed work | Remaining obligation |
|---|---|---|
| Degree-6 input certificate | Pure-integer DAG checker; six small closures checked against the direct engine; saved replayable certificate and malformed-input controls | Emit and check derivations of the actual depth-14 input basis; prototype is not that certificate |
| Independent lift review | Confirmed high-prime row-reduction defect, corrected it, added CI tests; four-prime arithmetic and 16 fresh lift runs pass | Separate external review / formal algorithm proof; input containment still required |
| Five summands, depth 3 | Tiled engine matches 16 control closures; bounded D8 retry stopped after 183.433 s at rank 1677 | Complete the lower ladder and D8 closure on approved resources |
| Proper-ideal exceptional primes | Exact LFD 5 at (k,p)=(7,257),(7,769),(8,7681),(9,12289) | Depth-10 extended run; general proper-ideal lower bound remains open |
| Depth-18 breakpoint | Unit instance fixed and closed V4 computed: rank 3821, 1 excluded, hence LFD at least 5 | Checked exclusion at D7 and upper witness at D8; prediction 8 is untested |

The depth-7 p=257 profile agrees between the ordinary chain closure and
the compressed reference engine. The other new exact profiles currently
have the compressed reference engine only. Exactness uses the checked
closed space plus evaluation-rank/codimension criterion; these new values
are finite computational evidence, not an all-depth result or Lean replay.

The initial depth-10 reference attempt stopped at 240 seconds with no
completed D6 closure. Its record is retained separately. The extended
attempt has a 1800-second cap. Its output is not an exact value unless
both the full closure check and the local-to-global rank criterion pass.

## Certificate scope

The saved certificate is `certificate-k2-D4.json`. Its checker accepts an
externally supplied generator list and parameters, then checks input,
linear-combination and admissible-shift nodes without importing producer
code. It also verifies every supplied output row against the DAG. It
certifies containment of those rows, not completeness or missing large
basis rows. The direct comparison additionally checks row independence
and equality of rank with the reference closure.

These scripts use the normalized generator as the declared input. They do
not silently admit a higher-degree raw generator at a lower budget.
The inadmissible-carry control specifically rejects a multiplication whose
normal form drops degree but whose source degree already equals D.

## Reproduction

From the repository root:

```
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -m pytest -q tests
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/validate.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/split_retry.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/exceptional_reference.py
python3 scripts/check_publication_artifacts.py
```

`validate.py` and the other experiment scripts rewrite their own output
records and timing fields. Replay in a separate checkout to retain the
saved originals. The archival `calibration/split_retry.py.txt` is the exact
runner source of the initial three-minute retry, before the full-ladder
option was added. Restore it to its original path in a separate checkout
if replaying that calibration. Its hash matches the calibration record.

See [PROTOCOL.md](PROTOCOL.md) for identities, predictions, stopping rules,
chronology and the proposed US$100 split-worker envelope; see
[REVIEW.md](REVIEW.md) for the arithmetic counterexample and review limits.

## Resource boundary

The tiled full-column basis at D8 has a worst-case size of about 191 GiB,
versus about 382 GiB for one float64 full-column basis, before the old
engine's full-size update/copy temporaries. The tiled version also removes
the queued dense row copies. These are storage estimates, not measured
peak usage or a completion guarantee.

The depth-18 direct D8 space has 40,999,516 columns. A compressed split
basis has a worst-case storage estimate near 764 TiB before buffers. A
structured method is needed before allocating a larger run. The proposed
split-worker budget does not cover depth-18 or a full depth-14 certificate.
