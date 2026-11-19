# Launch record, 2026-11-18 06:11 UTC: J1-J4 (paper V_6 closures, depths 14-16)

Owner request: "Do all 4 all at once" after an optimisation round. Every number is `[pcc:unverified]`
(no S3 checkpoint, bench.record or pcc claim reachable from the session container).

Binary `pkm_tower_pilot` 39acdeac... (450d5b295b, merged with engine-opt c51685b6a9, gate PASS),
`run_job.sh`, `BINARY.sha256`; J1/J2 also `hky_lift` 6824a8bd... + `pod_k14.sh` + the recorded lift JSON.
All pushed to `/workspace/in/code/`, SHA-256 verified by each pod. All pods cpu3m (cpu5m sold out),
32 vCPU, cgroup 256 GB (job limit 250 GB), US$1.76/h.

| job | pod | host CPU | instance | flags | disk | lifetime cap (worst case) |
|---|---|---|---|---|---|---|
| J1 | hky-j1 `j5ym8cqln1m237` US-CA-2 | EPYC 7713 | unit-p65537-a1-s10-k14 | VERIFY=1, export, then lift D6->7 | 200 GB | 43,200 s (US$21.12) |
| J2 | hky-j2 `5kmpu4gisqvd24` US-MD-1 | EPYC 7702P | unit-p163841-a3-s21-k14 | VERIFY=1, export, then lift | 200 GB | 43,200 s (US$21.12) |
| J3 | hky-j3 `g9w09ju7t3idas` EU-RO-1 | EPYC 9655 | unit-p65537-a1-s10-k15 | VERIFY=1, export | 320 GB | 216,000 s (US$105.60) |
| J4 | hky-j4 `tonzchna7tb5p8` US-MD-1 | EPYC 7702P | unit-p786433-a1-s10-k16 | EXPORT=0, VERIFY=0 | 320 GB (max) | 864,000 s (US$422.40) |

Predictions (PREDICTIONS.md, before any pod data): J1/J2 12,000 s (7k-24k), 31 GB (x2 with VERIFY);
J3 70,000 s (33k-170k), 75 GB; J4 370,000 s (170k-950k), 170 GB. Expected for J1/J2: profile
0,0,2,154,4451,73587,448922 with 1 not in V_6; lift: 1 in V_7, dim V_7 1,683,218, dim F 50,256.
Note: J4 checkpoints (about 1x RSS, two kept while writing) may exceed the 320 GB disk cap late in
the run; a failed checkpoint is recorded and the closure continues (450d5b295b).

## 2026-11-18 updates

- 06:3x UTC: J4 terminated on owner instruction ("J4 needs to be optimized before we spend 400"); spend under US$0.50.
- 07:16 UTC: hky-j3's container restarted on the host side (RSS was 36 GB of 250 GB); job j3 died with no exit code.
  Found at the 08:41 check; relaunched 08:42 as job `j3r`, resumed from the 30-min checkpoint (rank 526,992,
  sha256 5a657b70...). Lost about 35 min of compute plus 1.5 h idle.
- 08:20 and 09:45 UTC: J1 and J2 finished their closures (last checkpoints at rank 447,644 / 447,570) and
  entered the `--hky3-verify-closed` pass (RSS 38 -> 74 GB).
- 09:46 UTC: **the in-pod lifetime cap does not work.** `pod.py`'s `( sleep N; kill -TERM 1 )` signals PID 1,
  which ignores SIGTERM in these containers. Found while extending the J1/J2 caps: killing the sleep made the
  watchdog fire, and a replacement watchdog misread host uptime and fired too; neither stopped the pod or the
  job. The same holds for every earlier pod launched by this `pod.py`. Spend is bounded only by explicit
  termination; hard-stop reminders are scheduled (J1/J2 at 24 h, J3 at 216,000 s, bl-ladder at 104 h).
- 12:22-13:34 UTC: hky-j3's host hung (control server unreachable from ~12:53), container restarted at 13:34
  and 13:53; job j3r (in verify since ~11:30) died with SIGTERM. Checkpoint (rank 598,501, 40.5 GB,
  sha256 9b620a28...) survived. Resumed 15:4x UTC as job `j3r2` (VERIFY=1; the verify pass restarts).
- 13:37:53 UTC: **J1 done, status ok.** k=14, D=6, unit p=65537: cols 499,178, profile
  [0,0,2,154,4451,73587,448922], 1 not in V_6, rank 448,922; split part ranks 224,688 + 224,234 (cols
  249,824 / 249,354); `--hky3-verify-closed` closed=true, 0 non-zero over 4,120,873 products (18,495 s);
  muladds 1.0729e14; wall 26,723 s; peak RSS 73.9 GB (with verify); basis.bin sha256 06fdd4e4...
  The lift (`hky_lift`, pod_k14.sh) started on the pod at 13:37. Records in `results/J1/`.
- 15:37 UTC: **J4 relaunched with the skip** on owner instruction ("Launch j4 with skip"): pod hky-j4s
  `u2h6guqi1vozbx`, cpu3m, EU-RO-1, EPYC 9655P, 320 GB disk; binary `pilot_skip` 36a047ff... (43d5c40d1f),
  `ENGINE_FLAGS=--hky3-skip on --hky3-skip-witness leads`, EXPORT=0, VERIFY=0, job J4s-unit-p786433-a1-s10-k16-D6.
  Prediction (written before data, PREDICTIONS.md): 230,000 s (100k-950k), 160 GB, ~US$112 central.
  No working in-pod cap: hard-stop reminder set.
- 17:22 UTC: **J1 lift done** (Rust `hky_lift` 6824a8bd..., 32 threads, on the J1 basis 06fdd4e4...): `pod_k14.sh`
  verify exit 0 against the recorded lift `lift-unit-p65537-k14-D6to7.json` — k, D, p, dim V_6 = 448,922,
  1 in V_7 (contains_one true), fixed point, dim V_7 = 1,683,218, profile, dim F = 50,256 and
  dim(V_7 ∩ M_<=6) = 499,178 all **equal**; regular-top count 1,166,824 vs 1,173,002 is informational (a different
  notion: orbit count summed over parts). Input checks: 4 generators in V_6, 84 shifted random combinations stay in V_6.
  madds 5.60e14; peak RSS 234,285 MB (close to the 250 GB job limit); about 3.7 h. The end-of-run fall re-check
  reports 0 rows checked (1 is already in V_7). Records in `results/J1/lift/`.
- 18:15 UTC: J1 pod `j5ym8cqln1m237` terminated (GET 404). Spend about 12.1 h x US$1.76 = about US$21.
- 18:28 UTC: **J4s done** (k=16, D=6, unit p=786433, skip on): 1,149,017 columns, profile
  [0,0,2,162,5079,92323,770594], 1 not in V_6, rank 770,594 (parts 385,594 + 385,000); wall 10,238 s, peak RSS
  51.5 GB; no closedness check (VERIFY=0, memory). Predicted 230,000 s (100k-950k): the run was 22x faster than
  the central prediction and below the band.
- 19:02 UTC: **J3 done** (k=15, D=6, unit p=65537): 768,212 columns, profile [0,0,2,158,4761,82643,600547],
  1 not in V_6, rank 600,547 (parts 300,535 + 300,012); closedness check closed, 0 non-zero over 4,958,581
  products; wall of the last attempt 12,205 s; peak RSS 79.6 GB; basis sha256 5da602a8... (not downloaded).
- **Both profiles equal the frozen recursion** of lfd-summation-systems
  `experiments/ic-bounds-freeext/predictions-2026-11-10.json` ("recursion": p_j(k+1) = p_j(k) + 2 p_{j-1}(k) +
  p_{j-2}(k) for k >= 4D-10), applied from the certified V_6(14) profile: k=15 gives (0,0,2,158,4761,82643,600547)
  and k=16 gives (0,0,2,162,5079,92323,770594), equal entry by entry.
- 19:5x UTC: pods `g9w09ju7t3idas` (J3) and `u2h6guqi1vozbx` (J4s) terminated (GET 404). The J3 basis was not
  kept (a k=15 lift needs about 540 GB, beyond these pods). Spend: J3 about 13.7 h, J4s about 4.3 h.
- 20:02 UTC: **J2 done** (k=14, D=6, unit p=163841): 499,178 columns, profile [0,0,2,154,4451,73587,448922],
  1 not in V_6, rank 448,922 (parts 224,688 + 224,234 in 249,824 / 249,354 columns, as J1 and the paper);
  closedness check closed, 0 non-zero over 4,120,873 products; wall 49,797 s (slow, loaded host);
  peak RSS 74.0 GB; basis sha256 3d1f2996... The Rust lift started on the pod at 20:02.
- 2026-11-19 ~01:50 UTC: **J2 lift done** (Rust `hky_lift`, 32 threads, on the J2 basis 3d1f2996...): verify exit 0
  against `lift-unit-p163841-k14-D6to7.json` — k, D, p, dim V_6 448,922, 1 in V_7, fixed point,
  dim V_7 1,683,218, profile, dim F 50,256, dim(V_7 ∩ M_<=6) 499,178 all **equal**; regular tops 1,166,824 vs
  1,173,002 informational (different notion). Input checks: 4 generators in V_6, 84 shifted combinations stay
  in V_6. madds 5.60e14, peak RSS 234,299 MB. Records in `results/J2/lift/`.
- 2026-11-19 02:41 UTC: the session container restarted and its scratchpad (pod tokens) was wiped; the
  control tokens of hky-j2 and bl-ladder were recovered from the pods' env via the RunPod API. 02:45 UTC: J2 pod
  `5kmpu4gisqvd24` terminated (GET 404). Spend about 20.6 h x US$1.76 = about US$36.
