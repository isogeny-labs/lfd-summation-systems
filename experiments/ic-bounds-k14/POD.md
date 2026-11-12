# Runbook: one rented Runpod CPU pod for the k=14 IC-bounds runs

Tool: `scripts/ic_bounds_big_pod.py` (standard library, curl for every call). The
repository is private, so the pod cannot clone it: code and small inputs are pushed
over the pod's proxy URL, guarded by a per-pod token. The pod holds no Runpod, AWS or
GitHub key. Needs `RUNPOD_API_KEY` in the environment only for Runpod calls
(create, ls, status, stop, terminate, prices); remote control needs only the token.

## 1. Create

    T=scripts/ic_bounds_big_pod.py
    python3 $T prices                                  # $/h for cpu5m-32-256, cpu3m-32-256
    python3 $T create --name icb-k14 --dry-run         # inspect the body (token redacted)
    python3 $T create --name icb-k14 --yes-create      # 32 vCPU, 300 GB volume, 40 GB disk, 7-day cap

Options: `--flavors cpu5m,cpu3m` `--vcpu 32` `--volume-gb 300` `--disk-gb 40`
`--max-seconds 604800`. If Runpod rejects the volume on a CPU pod, retry with
`--volume-gb 0 --disk-gb <N>`: `/workspace` then lives on the container disk.
The output gives id, flavor, vCPU, memory, cost per hour and the worst-case cost at the cap.

## 2. Wait until it answers

    python3 $T ping icb-k14        # repeat every ~30 s; boot installs numpy/threadpoolctl first
    python3 $T sys icb-k14         # nproc, memory, /workspace disk, python, numpy, BLAS

On the pod, `/workspace/boot.log` holds the pip output and `/workspace/ctl.log` the
server log (fetch with `get`).

## 3. Push, run, monitor, fetch

    python3 $T push icb-k14 csd/ic_bounds_foo.py            # .py -> /workspace/code/, else /workspace/in/
    python3 $T push icb-k14 data.npy data.npy --dest in
    python3 $T run  icb-k14 k14a -- python3 ic_bounds_foo.py --k 14   # cwd /workspace/code, python3 only
    python3 $T jobs icb-k14                  # state, exit code, elapsed
    python3 $T tail icb-k14 k14a --bytes 8000
    python3 $T signal icb-k14 k14a TERM      # or KILL
    python3 $T get  icb-k14 out/result.json results/result.json [--max-bytes N]   # default 20 MB

`push` and `get` compare sha256 against the server and fail on a mismatch. A job name
that is still running is refused. Jobs survive a restart of the control server; a job
killed with KILL reports exit code -9, one stopped with TERM reports -15. Uploads over
the Runpod proxy may hit an edge size limit (about 100 MB on Cloudflare); split larger
inputs. Results meant to be fetched go under `/workspace/out/` and stay small.

## 4. Cost guard and teardown

- The start command runs `sleep $POD_MAX_SECONDS; kill -TERM 1` (default 604800 s =
  7 days): the pod stops itself at the cap. Worst case = `costPerHr * max_seconds / 3600`,
  printed at creation. Pick `--max-seconds` for the actual run, not the default.
- A stopped pod still bills its volume; only `terminate` removes it.
- When the results are fetched and verified:

      python3 $T terminate <pod-id> --yes
      python3 $T ls                  # must no longer show the pod (or show it gone)

  Then confirm in the Runpod console that no pod named `icb-k14` is left. `stop` and
  `terminate` refuse any pod id not recorded by `create`; other pods in the account
  are never touched.

## 5. What is stored where

| Item | Where |
| --- | --- |
| Control token (never printed, never in the repo) | `scratchpad/pod/<name>.token`, mode 0600 |
| Pod id, name, creation time | `scratchpad/pod/<name>.json` |
| Code, inputs, outputs | pod `/workspace/code`, `/workspace/in`, `/workspace/out` |
| Job logs | pod `/workspace/logs/<job>.log` |
| Job table (pid, command, start, exit code) | pod `/workspace/jobs.json`; exit files in `/workspace/jobs/` |
| Boot and server logs | pod `/workspace/boot.log`, `/workspace/ctl.log` |
| Control server source | embedded in the launcher; `ctl-source` writes it out |

The token is also in the pod environment (`CTL_TOKEN`), so it is visible in the Runpod
console; it is useless once the pod is terminated. Jobs do not inherit it. The state
directory is `$POD_STATE_DIR` (default: the scratchpad path above).

## Smoke test, 2026-11-10 (main session)

A 2-vCPU pod (`cpu3c`, 4 GB, $0.06/h, 30-minute cap, container disk 10 GB, no volume) was
created with `create --yes-create --volume-gb 0`, exercised and terminated after about four
minutes. Verified on the real platform:

- the v1 body is accepted as written with `volumeInGb: 0`; `/workspace` is then the container
  disk and is writable; the 20 KB `CTL_B64` environment variable is accepted;
- plain `pip install numpy threadpoolctl` works in `runpod/base:1.0.2-ubuntu2404`
  (Python 3.12.3, numpy 2.5.4, scipy-openblas 0.3.34, threadpoolctl 3.8.0);
- the control server answered `ping` about 15 s after creation; `sys`, `push` (sha256
  verified), `run`, `tail`, `jobs`, `get` and `terminate` all worked through the proxy URL;
- `os.cpu_count()` inside the pod reports the host's cores (192) while the container is pinned
  to its vCPUs: thread counts must come from `os.sched_getaffinity`.

Not yet verified: the largest container disk a 32-vCPU pod accepts, and uploads over 100 MB
(not needed: code and instance files are small).

## GPU option and smoke test, 2026-11-10 (main session)

`create` now takes `--gpu "<type id>" [--gpu-count N --min-ram-per-gpu G --min-vcpu-per-gpu V]`,
`--image` and `--pip "<packages>"` (installed at boot; default `numpy threadpoolctl`).

A one-GPU test pod (NVIDIA A40, 46 GB, $0.59/h, 25-minute cap) was created with
`--pip "numpy threadpoolctl cupy-cuda12x[ctk]"` on the default image
`runpod/base:1.0.2-ubuntu2404`, exercised and terminated after about five minutes:

- `cupy-cuda12x[ctk]` (CuPy 14.2.0, CUDA runtime 12.9 from pip wheels, driver 595.91) installs
  and works on the base image without a CUDA toolkit in the image; boot to `ping` about 45 s;
- float64 matrix products on the device reduced with `cupy.remainder` agree with Python
  integer arithmetic for p = 65537 (300 x 300, rows and columns of p - 1 included);
- the A40 does 2.5e11 multiply-adds per second in float64 (a consumer-class double-precision
  rate; data-centre cards are far faster, to be measured on the target card before the run).
