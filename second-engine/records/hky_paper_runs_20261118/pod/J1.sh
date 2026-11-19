#!/bin/bash
# J1: unit-p65537-a1-s10-k14, D = 6 (binary: pod/BINARY.sha256 = bin/pkm_tower_pilot_final, 450d5b295b, gated in README section 3).
# Instance file copied to /workspace/instances/.
exec "$(dirname "$0")/run_job.sh" J1-unit-p65537-a1-s10-k14-D6 /workspace/instances/unit-p65537-a1-s10-k14.json 6
