#!/usr/bin/env python3
"""Benchmark harness driver (csd/bench).

    python3 scripts/bench.py solvers                                   # what is registered, what is available
    python3 scripts/bench.py run --grid smoke --solvers py_ref,enum,lrf --out experiments/ic-bench/smoke.jsonl
    python3 scripts/bench.py run --cell n=25,k=9,m=3,plant=1 --solvers enum,lrf --aa --pin 0
    python3 scripts/bench.py summarize experiments/ic-bench/smoke.jsonl [--baseline enum]
    python3 scripts/bench.py model --n 131 --m 3 --preset ecc2k130 [--from-records experiments/ic-rootfind/n131.jsonl]

`run` exits non-zero if any row is `mismatch`, `invalid`, `crash` or `error`: a timing is only worth reading from
a row whose answer was checked.  See experiments/ic-bench/PROTOCOL.md for how the sweeps are meant to be used.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from csd.bench import fixtures, grids, model, runner, solvers, stats      # noqa: E402
from csd.bench.spec import Spec                                  # noqa: E402


def fmt_t(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    for unit, scale in (("s", 1), ("ms", 1e-3), ("µs", 1e-6), ("ns", 1e-9)):
        if x >= scale:
            return f"{x / scale:.3g} {unit}"
    return f"{x / 1e-9:.3g} ns"


def cmd_solvers(a):
    for name, s in sorted(solvers.REGISTRY.items()):
        ok, why = s.available()
        print(f"{name:10s} {s.kind:13s} {'available' if ok else 'UNAVAILABLE: ' + why}")


def cmd_run(a):
    names = a.solvers.split(",")
    for n in names:
        if n not in solvers.REGISTRY:
            sys.exit(f"unknown solver {n!r}; see `solvers`")
    if a.fixture:
        specs = []
        for name in a.fixture.split(","):
            specs += fixtures.specs(name)
        if a.filter:
            specs = [s for s in specs if a.filter in s.fixture]
        if a.seeds_max is not None:
            specs = [s for s in specs if s.seed < a.seeds_max]
    elif a.cell:
        specs = [Spec.parse(c) for c in a.cell]
    else:
        specs = grids.GRIDS[a.grid](seeds=a.seeds)
        if a.plant:
            specs = [Spec(**{**s.as_dict(), "plant": True}) for s in specs]
    pin = None if a.pin is None else [int(c) for c in a.pin.split(",")]
    out = Path(a.out or f"experiments/ic-bench/{a.grid or (a.fixture or 'cell').replace(',', '+')}.jsonl")
    bad = runner.run_grid(specs, names, out, reps=a.reps, timeout=a.timeout, threads=a.threads, pin=pin, aa=a.aa,
                          resume=a.resume, order_seed=a.order_seed)
    print(f"wrote {out} and {out.with_suffix('.meta.json')}")
    if bad:
        sys.exit(f"{bad} row(s) failed their checks (mismatch / invalid / crash / error): do not use these timings")


def load(paths):
    rows = []
    for p in paths:
        rows += [json.loads(l) for l in Path(p).read_text().splitlines() if l.strip()]
    return rows


def cmd_summarize(a):
    rows = load(a.records)
    by_status = defaultdict(int)
    for r in rows:
        by_status[(r["solver"], r["status"])] += 1
    print("## Row status\n")
    print("| solver | " + " | ".join(sorted({s for _, s in by_status})) + " |")
    sts = sorted({s for _, s in by_status})
    print("|---|" + "---|" * len(sts))
    for sv in sorted({s for s, _ in by_status}):
        print(f"| {sv} | " + " | ".join(str(by_status.get((sv, s), 0)) for s in sts) + " |")
    failed = [r for r in rows if r["status"] in ("mismatch", "invalid", "crash", "error")]
    if failed:
        print(f"\n**{len(failed)} failed rows; their timings are excluded below.**\n")
    ok = [r for r in rows if r["status"] == "ok" and r.get("arm", "A") == "A"]
    cell = lambda r: (r["n"], r["k"], r["m"], r["basis"], r.get("j"), r["curve"], r.get("plant", False))
    t = defaultdict(lambda: defaultdict(dict))            # cell -> solver -> seed -> t_min
    for r in ok:
        t[cell(r)][r["solver"]][r["seed"]] = r["t_min"]
    names = sorted({s for c in t.values() for s in c})
    print("\n## Per cell: median of t_min over seeds [IQR]; ndec = decompositions summed over seeds\n")
    print("| n | k | m | basis | " + " | ".join(names) + " | ndec |")
    print("|---|---|---|---|" + "---|" * len(names) + "---|")
    ndec = defaultdict(int)
    for r in ok:
        if r.get("n_decs") is not None and r["solver"] == next((s for s in ("py_ref", "enum", "lrf") if s in t[cell(r)]), None):
            ndec[cell(r)] += r["n_decs"]
    for c in sorted(t, key=lambda c: (c[3], str(c[4]), c[2], c[1], c[0])):
        cells = []
        for s in names:
            v = list(t[c].get(s, {}).values())
            if not v:
                cells.append("—")
            else:
                q = stats.summ(v)
                cells.append(f"{fmt_t(q['median'])} [{fmt_t(q['q1'])}–{fmt_t(q['q3'])}]")
        b = f"{c[3]}{'' if c[4] is None else c[4]}{'+planted' if c[6] else ''}"
        print(f"| {c[0]} | {c[1]} | {c[2]} | {b} | " + " | ".join(cells) + f" | {ndec.get(c, '')} |")
    if a.baseline:
        print(f"\n## Speed-up of {a.baseline} over each solver (paired by instance; geometric mean, 95% bootstrap CI)\n")
        print("| n | k | m | solver | slower by | CI | pairs |")
        print("|---|---|---|---|---|---|---|")
        for c in sorted(t, key=lambda c: (c[2], c[1], c[0])):
            for s in names:
                if s == a.baseline or a.baseline not in t[c] or s not in t[c]:
                    continue
                pr = stats.paired_ratio(t[c][s], t[c][a.baseline])
                if pr["n"]:
                    print(f"| {c[0]} | {c[1]} | {c[2]} | {s} | {pr['gmean']:.3g}× | {pr['lo']:.3g}–{pr['hi']:.3g} | {pr['n']} |")
    print("\n## Growth: log2 of median t_min against k, per (solver, m, basis); slope = doublings per step of k\n")
    print("| solver | m | basis | slope | 95% CI | R² | points |")
    print("|---|---|---|---|---|---|---|")
    grp = defaultdict(lambda: defaultdict(list))
    for c in t:
        for s, d in t[c].items():
            grp[(s, c[2], f"{c[3]}{'' if c[4] is None else c[4]}")][c[1]].append(stats.summ(list(d.values()))["median"])
    for key in sorted(grp):
        ks = sorted(grp[key])
        f = stats.fit_log2(ks, [sorted(grp[key][k])[len(grp[key][k]) // 2] for k in ks])
        if "slope" in f:
            print(f"| {key[0]} | {key[1]} | {key[2]} | {f['slope']:.2f} | {f['lo']:.2f}–{f['hi']:.2f} | {f['r2']:.3f} | {f['n']} |")
    refrows = [r for r in ok if r.get("ratio_to_reference")]
    if refrows:
        print("\n## Ratio to the recorded reference (same instance, experiments/ic-rootfind/): median over instances [min–max]\n")
        print("| solver | k | instances | ratio |")
        print("|---|---|---|---|")
        g = defaultdict(list)
        for r in refrows:
            g[(r["solver"], r["k"])].append(r["ratio_to_reference"])
        for key in sorted(g):
            v = g[key]
            print(f"| {key[0]} | {key[1]} | {len(v)} | {sorted(v)[len(v) // 2]:.2f} [{min(v):.2f}–{max(v):.2f}] |")
    arms = defaultdict(dict)
    for r in rows:
        if r["status"] == "ok":
            arms[(r["key"], r["solver"])][r["arm"]] = r["t_min"]
    pairs = [(v["A"], v["B"]) for v in arms.values() if "A" in v and "B" in v]
    if pairs:
        nf = stats.noise_floor(pairs)
        print(f"\n## A/A noise floor (same solver, same instance, two arms): B/A median {nf['median']:.3f}, "
              f"5–95% {nf['p05']:.3f}–{nf['p95']:.3f} over {nf['n']} pairs. Ratios inside this band are not resolved.")
    rss = defaultdict(list)
    for r in ok:
        rss[r["solver"]].append(r["rss_kib"] / 1024)
    print("\n## Peak RSS of the run (MiB, max over instances): " + ", ".join(f"{s} {max(v):.0f}" for s, v in sorted(rss.items())))
    skipped = defaultdict(set)
    for r in rows:
        if r["status"] == "skipped":
            skipped[r["solver"]].add(r.get("reason", ""))
    for s, why in sorted(skipped.items()):
        print(f"\nskipped {s}: " + "; ".join(sorted(why))[:300])


def cmd_fixtures(a):
    if a.what == "list":
        man = json.loads((fixtures.DIR / "MANIFEST.json").read_text())
        desc = {n: d for n, (d, _, _) in fixtures._sets().items()}
        desc["rust-records"] = "Rust-kernel instances up to n = 163 with their recorded cost per target (no expected set)."
        print("| set | rows | complete | with a decomposition | has reference | description |")
        print("|---|---|---|---|---|---|")
        for n in fixtures.available():
            rows = fixtures.load(n)
            print(f"| {n} | {len(rows)} | {sum(1 for r in rows if r['expected'].get('complete'))} | "
                  f"{sum(1 for r in rows if r['expected'].get('count', 0) > 0)} | {sum(1 for r in rows if r['reference'])} | "
                  f"{desc.get(n, '')} |")
    elif a.what == "make":
        for n in (a.set.split(",") if a.set else fixtures.SET_NAMES):
            rows = fixtures.generate(n, log=(lambda *_: None))
            fixtures.write(n, rows)
            print(f"wrote {fixtures.path(n)} ({len(rows)} rows)")
    else:
        bad = fixtures.verify(a.set.split(",") if a.set else None, deep=a.deep)
        drift = fixtures.rng_drift(a.set.split(",") if a.set else None)
        if drift:
            print(f"note: {len(drift)} instance(s) differ from what the current random stream builds (fixtures do not depend on it): {drift[:3]}")
        if bad:
            sys.exit(1)


def cmd_model(a):
    kw = {}
    if a.from_records:
        kw.update(model.from_records(a.from_records))
    for f in ("t_lin_ns", "t_cand_ns", "t_add_ns", "t_la_ns", "slack"):
        if getattr(a, f) is not None:
            kw[f] = getattr(a, f)
    if a.preset == "ecc2k130":
        p = model.ecc2k130(m=a.m, **kw)
    else:
        p = model.Params(n=a.n, m=a.m, frobenius=a.frobenius, cofactor_bits=a.cofactor_bits, **kw)
    ks = range(a.kmin or max(2, p.n // (p.m + 1) - 3), (a.kmax or p.n // 2) + 1)
    print(f"n={p.n} m={p.m}  t_lin={p.t_lin_ns:g} ns  t_cand={p.t_cand_ns:g} ns  t_add={p.t_add_ns:g} ns  "
          f"t_la={p.t_la_ns:g} ns (assumed)  rho: 2^{model.rho_iterations_log2(p):.2f} iterations\n")
    print("| k | log2 mu | log2 relations (s) | log2 algebra (s) | log2 total (s) | total core-years | algebra memory | gap to rho (log2) |")
    print("|---|---|---|---|---|---|---|---|")
    best = model.optimize_k(p)
    for r in model.sweep(p, [k for k in ks if abs(k - best["k"]) <= a.window]):
        mark = " **←**" if r["k"] == best["k"] else ""
        mem = 2.0 ** r["log2_la_mem_bytes"]
        print(f"| {r['k']}{mark} | {r['mu_log2']:.1f} | {r['log2_relations_s']:.1f} | {r['log2_algebra_s']:.1f} | "
              f"{r['log2_total_s']:.1f} | {2.0 ** r['log2_total_core_years']:.3g} | {mem / 2 ** 40:.3g} TiB | {r['log2_gap']:.1f} |")
    print(f"\nrho: {2.0 ** model.rho_stage_log2s(p):.3g} core-seconds = {2.0 ** (model.rho_stage_log2s(p) - model.LOG2_YEAR):.3g} core-years")
    sp = model.required_speedup(p)
    fl = sp["floor"]
    print(f"floor with a FREE relation finder (t_lin = 0): 2^{fl['log2_total_s']:.1f} core-seconds at k = {fl['k']} "
          f"({fl['log2_total_s'] - model.rho_stage_log2s(p):+.1f} vs rho in log2; algebra memory 2^{fl['log2_la_mem_bytes'] - 40:.1f} TiB)")
    if math.isinf(sp["factor_log2"]):
        print("required speed-up of the relation finder to reach rho: none suffices at this m (the floor is above rho)")
    else:
        print(f"required speed-up of the relation finder to reach rho (k re-optimized): 2^{sp['factor_log2']:.1f} "
              f"= {2.0 ** sp['factor_log2']:.3g}x, i.e. t_lin <= {sp.get('t_lin_ns_needed', p.t_lin_ns):.3g} ns, "
              f"break-even at k = {sp.get('k_at_break_even', best['k'])}")
    ob = model.oracle_budget(p)
    if ob["k"] is None:
        print("per-target budget for ANY relation finder: none (the algebra stage alone exceeds rho at every k)")
    else:
        print(f"per-target budget for ANY relation finder (all decompositions of one R, one core), best k = {ob['k']}: "
              f"{fmt_t(ob['scan_s'])}  (mu = 2^{ob['mu_log2']:.1f}, 2^{ob['relations_log2']:.1f} relations, "
              f"{fmt_t(ob['scan_s'] / 2.0 ** ob['mu_log2'])} per relation)"
              + ("" if ob["feasible"] else f"  -- shorter than forming the target ({fmt_t(ob['t_new_s'])}): infeasible"))
    print("Assumptions: |FB| = 2^k, mu = |FB|^m/(m! 2^n), t_la is not measured, kernel limits (e.g. the Rust kernel "
          "needs k >= 35 at n = 131) are not modelled.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("solvers").set_defaults(f=cmd_solvers)
    r = sub.add_parser("run")
    r.add_argument("--grid", choices=sorted(grids.GRIDS), default=None)
    r.add_argument("--fixture", default=None, help="frozen set(s), comma separated: see `fixtures list`")
    r.add_argument("--filter", default=None, help="with --fixture: keep ids containing this text")
    r.add_argument("--seeds-max", type=int, default=None, help="with --fixture: keep seeds below this")
    r.add_argument("--cell", action="append", help="n=25,k=9,m=3[,basis=,j=,curve=,seed=,plant=1] (repeatable)")
    r.add_argument("--solvers", default="py_ref,enum,lrf")
    r.add_argument("--seeds", type=int, default=3)
    r.add_argument("--plant", action="store_true", help="plant a decomposition in every target of the grid")
    r.add_argument("--reps", type=int, default=3)
    r.add_argument("--timeout", type=float, default=120.0)
    r.add_argument("--threads", type=int, default=1)
    r.add_argument("--pin", default=None, help="CPU list, e.g. 0 or 0,1")
    r.add_argument("--aa", action="store_true", help="run each solver twice (arms A and B) to measure the noise floor")
    r.add_argument("--resume", action="store_true")
    r.add_argument("--order-seed", type=int, default=0)
    r.add_argument("--out")
    r.set_defaults(f=cmd_run)
    s = sub.add_parser("summarize")
    s.add_argument("records", nargs="+")
    s.add_argument("--baseline", default=None)
    s.set_defaults(f=cmd_summarize)
    fx = sub.add_parser("fixtures", help="list / make / verify the frozen instance sets")
    fx.add_argument("what", choices=["list", "make", "verify"])
    fx.add_argument("--set", default=None, help="comma separated set names (default: all)")
    fx.add_argument("--deep", action="store_true", help="verify: re-run every solver that verified a row (minutes)")
    fx.set_defaults(f=cmd_fixtures)
    m = sub.add_parser("model")
    m.add_argument("--n", type=int, default=131)
    m.add_argument("--m", type=int, default=3)
    m.add_argument("--preset", choices=["ecc2k130"], default=None)
    m.add_argument("--frobenius", action="store_true")
    m.add_argument("--cofactor-bits", type=int, default=2)
    m.add_argument("--from-records", default=None)
    for f in ("t_lin_ns", "t_cand_ns", "t_add_ns", "t_la_ns", "slack"):
        m.add_argument("--" + f.replace("_", "-"), dest=f, type=float, default=None)
    m.add_argument("--kmin", type=int)
    m.add_argument("--kmax", type=int)
    m.add_argument("--window", type=int, default=4, help="show k within this distance of the optimum")
    m.set_defaults(f=cmd_model)
    a = ap.parse_args()
    if a.cmd == "run" and not a.grid and not a.cell and not a.fixture:
        ap.error("run needs --grid, --cell or --fixture")
    a.f(a)


if __name__ == "__main__":
    main()
