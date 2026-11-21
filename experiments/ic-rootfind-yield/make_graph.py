#!/usr/bin/env python3
"""Three-panel figure for the relation-yield study. Reads baseline_summary.jsonl and
frobenius.jsonl in this directory; writes yield_report.svg and yield_report.png."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
base = [json.loads(l) for l in open(os.path.join(here, "baseline_summary.jsonl")) if l.strip()]
frob = [json.loads(l) for l in open(os.path.join(here, "frobenius.jsonl")) if l.strip()]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# Panel 1: cost/relation vs n at the optimal (square) operating point, m=3
opt = {}
for r in base:
    if r["m"] != 3: continue
    if r["delta"] in (0, 1):
        opt.setdefault(r["n"], []).append(r)
ns = sorted(opt)
cost_opt = [min(x["cost_per_relation"] for x in opt[nn]) for nn in ns]
ax[0].semilogy(ns, cost_opt, "o-", color="#1f77b4", label="optimal delta in {0,1}, m=3")
# high-delta and m=4 for contrast
hi = {}
for r in base:
    if r["m"] == 3 and r["delta"] >= 3:
        hi.setdefault(r["n"], []).append(r["cost_per_relation"])
nh = sorted(hi)
ax[0].semilogy(nh, [min(hi[nn]) for nn in nh], "s--", color="#ff7f0e", label="high delta (>=3), m=3")
m4 = [(r["n"], r["cost_per_relation"]) for r in base if r["m"] == 4]
if m4:
    ax[0].semilogy([a for a, _ in m4], [b for _, b in m4], "^", color="#d62728", label="m=4")
ax[0].set_xlabel("field degree n"); ax[0].set_ylabel("cost per relation  (systems+candidates)")
ax[0].set_title("A. Cost/relation: square system is optimal, m=4 worse")
ax[0].legend(fontsize=8); ax[0].grid(True, which="both", alpha=0.3)

# Panel 2: found_per_cand vs n at matched delta (the intrinsic yield decay)
for d, col, mk in ((0, "#2ca02c", "o"), (3, "#9467bd", "s")):
    pts = sorted((r["n"], r["found_per_cand"]) for r in base if r["m"] == 3 and r["delta"] == d)
    if pts:
        ax[1].semilogy([a for a, _ in pts], [b for _, b in pts], mk + "-", color=col, label=f"delta={d}")
ax[1].set_xlabel("field degree n"); ax[1].set_ylabel("found_per_cand")
ax[1].set_title("B. Yield loss is intrinsic: found/cand collapses with n")
ax[1].legend(fontsize=8); ax[1].grid(True, which="both", alpha=0.3)

# Panel 3: Frobenius ideal vs realizable column reduction
fn = sorted(r["n"] for r in frob)
fd = {r["n"]: r for r in frob}
ax[2].plot(fn, [fd[nn]["n_eff_ideal"] for nn in fn], "o-", color="#8c564b", label="ideal |F|/|orbits|")
ax[2].plot(fn, [fd[nn]["real_LA_reduction"] for nn in fn], "s-", color="#e377c2", label="realizable (LRF from V^(0))")
ax[2].axhline(1.0, color="gray", ls=":", lw=1)
ax[2].set_xlabel("field degree n"); ax[2].set_ylabel("column reduction factor")
ax[2].set_title("C. Frobenius no-go: realizable reduction -> 1")
ax[2].legend(fontsize=8); ax[2].grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(os.path.join(here, "yield_report.svg"))
fig.savefig(os.path.join(here, "yield_report.png"), dpi=150)
print("wrote yield_report.svg and yield_report.png")
