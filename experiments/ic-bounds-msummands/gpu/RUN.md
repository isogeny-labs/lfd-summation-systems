# m-summand runs on one rented pod (2026-11-12)

Written before the run. Predictions: `../predictions-2026-11-12.json` (sha256 `deb4b304…aad7`, frozen before any closure of
the cases; see `../README.md` for the two smoke closures that preceded it).

## What runs (`chains.json`, wrapper `multi_msum.py`, chain `chain_msum.py`)

One H200 pod. Chains side by side:

| chain | systems | engine | threads |
|---|---|---|---|
| bin | binary r = 16 (G-m4n9k4, 3), r = 18 (G-m3n10k6, 3), r = 21 (G-m3n11k7, 3), r = 20 (G-m4n10k5-koblitz0) | `ic_bounds_bin_big` on the GPU | 2 |
| pm4-a, pm4-b | prime m = 4, k = 4, three instances each (6) | `CompactClosure`, CPU | 8 each |
| pm5 | prime m = 5, k = 3 (6) | `CompactClosure`, CPU | 3 |
| pm3k6 | prime m = 3, k = 6 (2) | `CompactClosure`, CPU | 3 |

The earlier pods of this account had 24 vCPU per H200 (the host reports 192 cores, but the pod's share is 24), so
the threads sum to 24.

Not run (cost): the other r = 20 binary cells, r = 24 (G-m3n12k8, G-m4n12k6), prime m = 3, k = 7. Listed in `chains.json`.

Every instance is rebuilt on the pod from its spec and its hash is compared with the frozen hash in `cost.json`; a
mismatch stops that system. Records leave the pod as `msum-<chain>.tgz`, rewritten every two minutes and copied by the
detached loop `mirror.sh`, which terminates the pod when every chain has ended.

## Change after the freeze

Binary degree caps raised by one (at most r = mk): a local smoke test of the wrapper on G-m3n10k6-koblitz0-ran-ord ended
`capped` at the frozen cap 5 (strict falls at 4 and 5, so LFD >= 5 there, which is above d_in + 1 = 4). The cap is a run
parameter, not a prediction; `cells.json` keeps the frozen value as `max_degree_frozen`. A system still capped at the
raised cap is recorded as capped with its lower bound, never as an exact value.

## Questions the runs decide

- Prime m = 4, k = 4 (d = deg g = 12 < mk = 16): the frozen candidates give 13 (LFD = d + 1) and 17 (bound mk + 1
  attained). This is the discriminating cell; every prime cell run locally had d = mk, where the value is forced.
- Prime m = 5, k = 3 (d = 15 = mk): forced to 16 by Proposition 1 and the kernel bound; a check of the implementation.
- Prime m = 3, k = 6: LFD >= 8 is recorded on one instance; the closure at D = 8 decides 8 or more.
- Binary: the point rule LFD = min(r, d_in + 1) and the bound LFD <= min(r, d_in + 2) at r = 16 to 21.

## Budget

Lifetime cap 6 h (at most $31.74 at $5.29/h).
