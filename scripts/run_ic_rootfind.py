"""Benchmarks of csd.ic_rootfind (linearized root finding) against its baselines; records in experiments/ic-rootfind/.

    python3 scripts/run_ic_rootfind.py scaling     # m = 3, V = <1, w, ..>, n = 3k - 2: LRF against enumeration
    python3 scripts/run_ic_rootfind.py m4          # m = 4, n = 4k - 3: LRF against meet in the middle and enumeration
    python3 scripts/run_ic_rootfind.py basis       # fixed (n, k): the factor base from <1, w, ..> to random
    python3 scripts/run_ic_rootfind.py regime      # fixed k = 12: n from below 3k to above it
    python3 scripts/run_ic_rootfind.py algebraic   # m = 3: the repository's MutantXL on the Weil descent of S_4
    python3 scripts/run_ic_rootfind.py all

Every row is one seeded instance.  Each method that runs on it must return the same set of decompositions, and
every decomposition is re-added with csd.ic_dreg's group law; a row is written only after both checks.  Times are
single-thread wall clock (the minimum of `reps` runs for the C kernels), excluding the table set-up, which is
reported separately as `t_setup`.
"""
from __future__ import annotations

import argparse
import json
import platform
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from csd import ic_dreg as I          # noqa: E402
from csd.ic_rootfind import Solver    # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "experiments" / "ic-rootfind"


def timed(f, reps):
    best, res = None, None
    for _ in range(reps):
        t = time.perf_counter()
        res = f()
        dt = time.perf_counter() - t
        best = dt if best is None else min(best, dt)
    return res, best


def instance(n, k, seed, basis_kind="standard", j=None):
    rng = random.Random(f"ic-rootfind/{n}/{k}/{seed}/{basis_kind}/{j}")
    F = I.binary_field(n)
    a, b = rng.getrandbits(1), rng.randrange(1, 1 << n)
    E = I.BinaryCurve(F, a, b)
    if basis_kind == "standard":
        basis = I.standard_subspace(F, k)
    elif basis_kind == "random":
        basis = I.random_subspace(F, k, rng)
    else:                                   # <1, .., w^(j-1)> and k - j random vectors
        while True:
            basis = I.standard_subspace(F, j) + [rng.getrandbits(n) for _ in range(k - j)]
            if I._f2_rank(basis) == k:
                break
    R = E.random_point(rng)
    return E, basis, R


def check(E, R, decs):
    for dec in decs:
        acc = None
        for P in dec:
            assert E.on_curve(P)
            acc = E.add(acc, P)
        assert acc == tuple(R), "decomposition does not sum to R"


def one(n, k, m, seed, basis_kind="standard", j=None, brute=True, mitm=False, reps=3):
    E, basis, R = instance(n, k, seed, basis_kind, j)
    t = time.perf_counter()
    S = Solver(E, basis)
    t_setup = time.perf_counter() - t
    (decs, st, rows), t_lrf = timed(lambda: S.lrf(m, R), reps)
    check(E, R, decs)
    row = dict(n=n, k=k, m=m, seed=seed, basis=basis_kind, j=j, a=E.a, b=E.b, xR=R[0], yR=R[1],
               dimV2=S.dimV2, unknowns=k + 1, equations=S.w, nfb=S.nfb, decompositions=len(decs),
               t_setup=round(t_setup, 6), t_lrf=t_lrf, stats=st)
    if brute:
        (bd, _), t_b = timed(lambda: S.brute(m, R), 1 if (m - 1) * k > 18 else reps)
        assert bd == decs, "enumeration and LRF disagree"
        row.update(t_brute=t_b, speedup_vs_brute=t_b / t_lrf)
    if mitm:
        (md, _, tb), t_m = timed(lambda: S.mitm4(R), 1)
        assert md == decs, "meet in the middle and LRF disagree"
        row.update(t_mitm=t_m, mitm_table_bytes=tb, speedup_vs_mitm=t_m / t_lrf)
    return row


def write(name, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w") as f:
        for r in rows:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    print(f"wrote {OUT / name} ({len(rows)} rows)")


def show(r, extra=()):
    keys = ("n", "k", "m", "basis", "j", "dimV2", "nfb", "decompositions", "t_lrf", "t_brute", "t_mitm") + tuple(extra)
    print("  ".join(f"{q}={r[q]:.3g}" if isinstance(r.get(q), float) else f"{q}={r.get(q)}" for q in keys if q in r),
          flush=True)


def scaling(args):
    rows = []
    for k in range(5, args.kmax + 1):
        for seed in range(args.seeds):
            r = one(3 * k - 2, k, 3, seed, brute=k <= args.kbrute)
            rows.append(r)
            show(r)
    write("scaling-m3.jsonl", rows)


def m4(args):
    rows = []
    for k in range(4, 12):
        for seed in range(args.seeds):
            r = one(4 * k - 3, k, 4, seed, brute=k <= 8, mitm=True)
            rows.append(r)
            show(r)
    write("m4.jsonl", rows)


def basis(args):
    rows = []
    n, k = 31, 10
    for j in range(k, -1, -1):
        if j == 1:
            continue                         # <1> + 9 random vectors is as random as j = 0
        for seed in range(args.seeds):
            r = one(n, k, 3, seed, basis_kind="standard" if j == k else "mixed", j=j)
            rows.append(r)
            show(r)
    write("basis-n31-k10.jsonl", rows)


def regime(args):
    """Fixed k = 12, m = 3, V = <1, w, ..>: n from below 3k (free unknowns in the linear system) to above it."""
    rows = []
    k = 12
    for n in range(27, 46, 2):
        for seed in range(args.seeds):
            r = one(n, k, 3, seed)
            rows.append(r)
            show(r)
    write("regime-k12.jsonl", rows)


def algebraic(args):
    """The repository's MutantXL (csd.ic_dreg) on the Weil descent of S_4(x1, x2, x3, xR): time and counted
    operations of the closure at each degree up to the first that gives every linear form vanishing on the
    solution set.  The Boolean system has more solutions than R has decompositions: S_4 = 0 also holds for
    abscissas of points over the quadratic extension and for sign patterns summing to -R."""
    rows = []
    for k in (4, 5):
        for seed in range(args.seeds):
            n = 3 * k - 2
            r = one(n, k, 3, seed)
            E, B, R = instance(n, k, seed)
            t = time.perf_counter()
            eqs = I.binary_descent(E, [B, B, B], R[0])
            t_descent = time.perf_counter() - t
            sols = I.binary_solutions(E, [B, B, B], R[0])
            pts = [c[0] | c[1] << k | c[2] << 2 * k for c in sols]
            v = 3 * k
            target = I._lin_target_f2(pts, v)
            degs = sorted({max(bin(mk).count("1") for mk in eq) for eq in eqs})
            per_degree, d_mxl = [], None
            for D in range(max(degs), v + 1):
                t = time.perf_counter()
                M, cols = I.boolean_macaulay(eqs, v, D)
                rank, piv, ops, rounds, added = I._mxl_closure(M, cols, v, D, 2, I._bool_mulvar)
                dt = time.perf_counter() - t
                lin = sum(1 for c in piv if c >= cols - (v + 1))
                per_degree.append(dict(D=D, cols=cols, rows=int(M.shape[0]), rank=int(rank), t=dt, ops=float(ops)))
                if lin == target:
                    d_mxl = D
                    break
            r.update(boolean_vars=v, boolean_eqs=len(eqs), eq_degrees=degs, boolean_solutions=len(sols),
                     t_descent=t_descent, mxl=per_degree, d_mxl=d_mxl, t_mxl_final=per_degree[-1]["t"],
                     t_mxl_total=sum(d["t"] for d in per_degree), mxl_ops_final=per_degree[-1]["ops"],
                     speedup_vs_mxl_final=per_degree[-1]["t"] / r["t_lrf"])
            rows.append(r)
            show(r, ("boolean_solutions", "d_mxl", "t_mxl_final", "t_mxl_total"))
    write("algebraic-m3.jsonl", rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["scaling", "m4", "basis", "regime", "algebraic", "all"])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=21)
    ap.add_argument("--kbrute", type=int, default=12)
    args = ap.parse_args()
    print(platform.processor() or platform.machine(), platform.python_version())
    for name in (["scaling", "m4", "basis", "regime", "algebraic"] if args.what == "all" else [args.what]):
        globals()[name](args)


if __name__ == "__main__":
    main()
