# Relaunch plan: the k = 14 runs lost on 2026-11-11

Status: **launched 2026-11-11** with `relaunch/launch.sh` after the account was funded ($48.91, so the lifetime
caps were shortened to 27,000 s for the GPU pod and 36,000 s for the CPU pod). GPU pod `icb-k14b`
(`mq9aalc4ub1tbl`, H200, $5.29/h, created 06:15 UTC): both self-tests passed through `nopool_run.py`
(06:23 and 06:35 UTC), and the closure of the control instance started at 06:35 UTC. CPU pod `icb-ref11`
(`5uzs5jx9ovkvqf`, 8 vCPU, $0.44/h, created 06:14 UTC): reference engine running. Records are copied to
`relaunch/fetched/` by two detached `mirror.sh` processes. **Outcome (11:22 UTC):** jobs A1, A2 and
the checker on A1's basis are done (control passed: 1,680,848 exactly), job C is done (hashes equal), and both
pods are terminated. Jobs B1 and B2 (p = 163841) were not run: the pod was too slow to finish them before its
cap, so `relaunch/stop_after_cert.sh` stopped it after the control. Results: `RESULTS.md`, section "Relaunch". The text below is the plan as written before the launch:
at that time the status was **not started**. The Runpod balance was -$1.09 when this was written, so no pod can be created until
the account is funded. Nothing below has been run in this form; the command lines are the recorded ones from
`README.md` section 11 and `LIFT.md` section 8 with the instance changed.

## What is missing and what each run decides

| job | input | time on an H200 (from the records) | reading |
|---|---|---|---|
| A1. Closure k = 14, D = 6, basis kept | `instances/proper-p114689-s22-k14.json` (2,370 solutions) | about 1.8 to 2.2 h | expect closed, dim 448,922 |
| A2. Lift D = 6 -> 7 on A1's checkpoint | A1's `--ckpt` directory | about 13 min | **soundness control**: dim V_7 must be at most 1,683,218 - 2,370 = 1,680,848. Above it: the lift has a bug. Exactly on it: LFD = 7 on a random curve over a third prime |
| B1. Closure k = 14, D = 6, basis kept | `instances/unit-p163841-a3-s21-k14.json` | about 1.8 h (6,393 s recorded) | must reproduce hashes e0c4afd7… / 44a89923… |
| B2. Lift D = 6 -> 7 on B1's checkpoint | B1's `--ckpt` directory | about 13 min | 1 in V_7: LFD = 7 exactly on a second unit instance and prime |
| C. Reference engine k = 11, D = 6 (pod CPU) | unit instance over p = 65537 | estimate 5 to 8 h, about 24 GB (scaled from k = 10: 3,129 s, 7.1 GB) | second engine for D = 6 above k = 10 |

A and B in sequence are about 4.5 h of card time, about $24 at $5.29/h; C runs on the pod's CPU at the same
time and sets the total to about 6 to 8 h, $32 to $42. These are estimates from the recorded run times.

## Scripts (added 2026-11-11, `relaunch/`)

`relaunch/launch.sh gpu|ref|both` does everything below: it creates the pod(s), starts `relaunch/mirror.sh`
(the copy loop, detached, which also terminates a pod when its job is done), pushes the files and starts the
jobs. The GPU pod runs `relaunch/chain_k14.py`, which needs no connection to the launching machine: self-tests,
then per instance the closure, the lift (one retry with a smaller work area) and the independent checker, with
progress in `/workspace/out/chain-status.json`. The reference engine runs on its own 8-vCPU pod through
`relaunch/ref_k11.py`, so the card is not billed while it finishes. Records land in `relaunch/fetched/`.

Checked locally before any launch: the chain's three command lines on the k = 10, D = 5 random-curve instance
with the simulated device (closure with checkpoint, lift from the checkpoint with the basis on the host:
60,446, checker: all true), and `ref_k11.py` against `scripts/ic_bounds_d5_verify.py` at k = 7, D = 6 (all
fields and both hashes equal). Not checked anywhere: `relaunch/nopool_run.py` (needs CuPy); the chain runs both
self-tests through it first and stops if either fails. Lifetime caps: 36,000 s for the GPU pod (at most $53)
and 54,000 s for the CPU pod.

## Before creating the pod

1. Start the copy loop first, outside any agent session (`nohup` or launchd): every 10 minutes
   `python3 scripts/ic_bounds_big_pod.py jobs <name>` and `get <name> out/<file> ...` for every file under
   `/workspace/out`, into `experiments/ic-bounds-k14/runs/`. Both earlier losses (2026-11-06 and 2026-11-11)
   were results that existed on a pod and had not been copied off it.
2. Pick `--max-seconds` for the run (36,000 s covers the plan), not the 7-day default.
3. The control token is written to `$POD_STATE_DIR/<name>.token`. Set `POD_STATE_DIR` to a directory that
   survives the session.

## Commands

```
T=scripts/ic_bounds_big_pod.py
python3 $T create --name icb-relaunch --gpu "<H200 type id from `prices`>" --pip "numpy threadpoolctl cupy-cuda12x[ctk]" \
  --volume-gb 0 --disk-gb 200 --max-seconds 36000 --yes-create
python3 $T push icb-relaunch csd/ic_bounds_big.py && python3 $T push icb-relaunch csd/ic_bounds_lift.py
python3 $T push icb-relaunch scripts/ic_bounds_big_run.py && python3 $T push icb-relaunch scripts/ic_bounds_lift_run.py
python3 $T push icb-relaunch csd/ic_bounds_d5_verify.py && python3 $T push icb-relaunch scripts/ic_bounds_d5_verify.py
python3 $T push icb-relaunch experiments/ic-bounds-k14/instances/proper-p114689-s22-k14.json
python3 $T push icb-relaunch experiments/ic-bounds-k14/instances/unit-p163841-a3-s21-k14.json

# self-tests first (both must exit 0)
python3 $T run icb-relaunch st1 -- python3 ic_bounds_big_run.py selftest --device gpu --out /workspace/out/selftest-gpu.json
python3 $T run icb-relaunch st2 -- python3 ic_bounds_lift_run.py selftest --device gpu --out /workspace/out/lift-selftest-gpu.json

# A1 (settings of the recorded second k = 14 run: --dev-gb 105 --work-gb 20 --tile-mb 1024 --batch 1024)
python3 $T run icb-relaunch a1 -- python3 ic_bounds_big_run.py closure --instance /workspace/in/proper-p114689-s22-k14.json --D 6 \
  --device gpu --dev-gb 105 --work-gb 20 --tile-mb 1024 --batch 1024 --ckpt /workspace/ckpt/p114689-k14-D6 --ckpt-minutes 30 --resume \
  --status /workspace/out/proper-p114689-k14-D6.status.json --out /workspace/out/proper-p114689-k14-D6.json
# A2
python3 $T run icb-relaunch a2 -- python3 ic_bounds_lift_run.py lift --ckpt /workspace/ckpt/p114689-k14-D6 \
  --instance /workspace/in/proper-p114689-s22-k14.json --device gpu --basis-on host --dev-gb 130 --work-gb 16 \
  --lift-ckpt /workspace/ckpt/p114689-k14-D6-lift --ckpt-minutes 20 --resume \
  --status /workspace/out/lift-proper-p114689-k14-D6to7.status.json --out /workspace/out/lift-proper-p114689-k14-D6to7.json
# B1, B2: the same two commands with unit-p163841-a3-s21-k14.json and its own checkpoint and output names
# C
python3 $T run icb-relaunch c -- python3 ic_bounds_d5_verify.py closure --k 11 --D 6 --prime 65537 --a 1 --b 0 --seed 10 --full-sweep
```

Not recorded in the repository, so to be settled on the pod: the first lift attempt at k = 14 ran out of card
memory after 6 minutes and the retry "with the memory pool disabled" completed; the switch that disables the
pool is not in the committed driver's options. If the lift fails the same way, lower `--work-gb` and use a new
`--lift-ckpt` directory (`LIFT.md` section 9). The exact invocation of job C should be checked against the
argv of `runs/reference-engine-unit-p65537-k10-D6.json` if that record carries one.

## After the runs

Fetch every record, check the readings in the table, add them to `RESULTS.md` and to the paper's paragraph
"The degree-seven step at depth fourteen" (the sentence that says the k = 14 control has not been completed),
then `terminate` the pod and confirm with `ls`.

## Job B launched on 2026-11-11, 20:26 UTC (owner: "Launch")

Pod `icb-k14d` (`t7xtcu7qokwc72`, H200, US-NC-1, $5.29/h, lifetime cap 12,300 s = 23:51 UTC, fitted to the
balance of $20.22). Chain `relaunch/chain_u163841.py` (job B only: unit instance over p = 163841; closure with the
basis kept, lift, checker); copy loop `mirror.sh` into `relaunch/fetched-d/`, which also terminates the pod when
the chain is finished.

Two findings before the run:
- A first pod (`icb-k14c`, CA-MTL-3, created 20:13 UTC) ran the engine self-test several times slower than the
  recorded runs and was terminated after about nine minutes (about $0.80).
- The cause was not the machine. `relaunch/nopool_run.py` disables the pinned-memory pool as well as the device
  pool; the runs of 2026-11-10 disabled only the device pool. On `icb-k14d` the engine self-test took more than
  230 s with `nopool_run.py` and 100 s with `relaunch/nopool_dev_only.py` (device pool only), which was uploaded
  under the name `nopool_run.py`. This also explains why the relaunch pod of 06:15 UTC was 1.8 times slower in
  the closure and six times slower in the lift. With the device-only wrapper the closure follows the recorded
  trajectory (ranks 26,857 + 27,157 at 270 s, as in the second run of 2026-11-10).
- The launcher has a new option `--datacenters` (allow-list `dataCenterIds`); the pod was created by polling
  US-NC-1 and EU-FR-1 once a minute (fourth try).

Both self-tests passed on the card (20:32 and 20:36 UTC); the closure started at 20:36 UTC. Expected: closure
about 22:26, lift about 22:40, checker about 22:53 UTC. Results will be added to `RESULTS.md`.

### Job B, first attempt failed; second attempt launched 22:05 UTC

The closure on `icb-k14d` exited with code 1 at 21:45 UTC (about 4,100 s, ranks 136,800 + 140,516 at the last
status) and the copy loop terminated the pod, as designed (79 minutes billed). The closure's log was not in the
copy list and is lost; the cause is inferred, not read: `chain_k14.py` asks for a checkpoint every 30 minutes,
the status stood still for more than ten minutes after 2,070 s (a checkpoint of about 120 GB being written), the
exit came at the time of the second one, and the pod had a 150 GB disk (my setting; the earlier relaunch used
250 GB). Two checkpoints near the peak are about 240 GB.

Second attempt: pod `icb-k14e` (`ovbtxrjlxme4fs`, H200, US-NC-1, 200 GB disk, cap 16,200 s), same files, with
`chain_u163841.py` changed to write one checkpoint at the end of the closure phase, and the closure and checker
logs added to the copy list (`relaunch/fetched-e/`). Self-tests passed at 22:07 and 22:11 UTC; closure started
22:11 UTC. The account balance was $112.81 at 22:04 UTC.

### Job B outcome (2026-11-12, 00:27 UTC)

Done on `icb-k14e`: closure 6,634 s (closed, hashes equal to the first run of the instance), lift 677 s (fixed
point, 1 in V_7, dim V_7 = 1,683,218), checker (all true). The copy loop terminated the pod at 00:26:58 UTC;
no pod of this work exists. Records in `runs/` (`rerun-unit-p163841-k14-D6.json`,
`lift-unit-p163841-k14-D6to7.json`, `certcheck-unit-p163841-k14-D6.json`), table in `RESULTS.md`, section "Job B".
