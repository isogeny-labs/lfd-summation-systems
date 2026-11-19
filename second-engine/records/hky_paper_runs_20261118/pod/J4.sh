#!/bin/bash
# J4: unit-p786433-a1-s10-k16, D = 6, with the lead-cover skip (research/hky_engine_skip_20261118: PROOF.md,
# gate PASS, README). Binary: pod/BINARY_J4.sha256 = research/hky_engine_skip_20261118/bin/pilot_skip
# (43d5c40d1f, stamped), copied to /workspace/bin/pkm_tower_pilot_skip. J1-J3 keep pod/BINARY.sha256.
# A new job name: a checkpoint of the earlier J4 (no skip, other fingerprint) is never resumed by this one.
# EXPORT=0 (the basis, ~270 GB, does not fit beside two checkpoints on 320 GB); VERIFY=0 (doubles memory).
# Instance file copied to /workspace/instances/.
HERE=$(cd "$(dirname "$0")" && pwd)
export BIN=${BIN:-/workspace/bin/pkm_tower_pilot_skip}
export SHA_FILE=${SHA_FILE:-$HERE/BINARY_J4.sha256}
export ENGINE_FLAGS=${ENGINE_FLAGS:---hky3-skip on --hky3-skip-witness leads}
export EXPORT=${EXPORT:-0}
exec "$HERE/run_job.sh" J4s-unit-p786433-a1-s10-k16-D6 /workspace/instances/unit-p786433-a1-s10-k16.json 6
