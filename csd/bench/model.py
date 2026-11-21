"""End-to-end cost model of index calculus on a binary curve against Pollard rho.

The harness *measures* the per-target cost of a relation finder; this module turns measurements into a whole-attack
estimate and says what is assumed.  Nothing here is a measurement of a full attack.

Stages (m summands, factor base V of dimension k, |FB| ~ 2^k points, closed under negation):

  relations  one scan of the factor base for a target R costs |FB|^(m-2) targets of the two-summand problem, each
             t_target = t_lin + (cands - 1) t_cand  (cands = 2^max(0, 3k - 1 - n) for V = <1, w, ..>, dimV2 = 2k-1),
             plus t_new to form the next target R (one group addition by a random walk R <- R + P_j: this is the
             floor that stays when the scan itself is free);
             a scan yields mu = |FB|^m / (m! 2^n) decompositions of a random R (the expectation, not a guarantee)
             and about |FB| (1 + slack) / neg_identify relations are needed.
  algebra    sparse linear algebra on N = |FB| / (neg_identify * fb_reduction) columns, m entries per row, block
             Wiedemann: 2 N matrix-vector products of N m entries each at t_la_ns per entry (an ASSUMED constant: a
             multiply-add modulo the group order; it is not measured here).
  rho        sqrt(pi/2 * ell/sym) iterations, ell = 2^n / cofactor, sym = 2 (negation) or 2 n (Frobenius, Koblitz),
             each costing rho_adds group additions at t_add_ns.

Units: core-seconds; log2 where it is a ratio.  The defaults of t_lin_ns, t_cand_ns and t_add_ns are the n = 131
record of experiments/ic-rootfind/ (1850 ns per target, ~260 ns candidates, 150 ns per addition); `from_records`
re-derives them from a JSONL file.  With them the model reproduces the figures in that README (about one core-year
per scan of 2^44 targets, about 3e13 core-years for the relation stage, rho about 2^61 iterations).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace
from pathlib import Path

LOG2_YEAR = math.log2(365.25 * 24 * 3600)


@dataclass(frozen=True)
class Params:
    n: int
    m: int = 3
    k: int | None = None             # factor-base dimension; None: optimized
    t_lin_ns: float = 1850.0         # per two-summand target, candidates included at cands = 1
    t_cand_ns: float = 260.0         # per extra candidate (ablation at n = 131)
    t_add_ns: float = 150.0          # one group addition with shared inversion
    t_new_ns: float | None = None    # forming a fresh target R; None: one group addition (t_add_ns)
    t_la_ns: float = 10.0            # ASSUMED: one sparse matrix entry update mod the group order
    slack: float = 0.05              # relations beyond |FB|
    neg_identify: int = 2            # P and -P are one column
    fb_reduction: int = 1            # Frobenius orbits of the factor base (Koblitz, invariant V): up to n
    cofactor_bits: int = 2           # log2 of the curve cofactor (ECC2K-130: 4)
    frobenius: bool = False          # rho on Frobenius orbits too
    rho_adds: float = 1.0            # group additions per rho iteration


def ecc2k130(**kw) -> Params:
    """The Certicom ECC2K-130 setting: n = 131, cofactor 4, rho on Frobenius orbits (Bailey et al., ePrint 2009/541)."""
    return Params(**{**dict(n=131, m=3, cofactor_bits=2, frobenius=True), **kw})


def log2fb(k):
    return float(k)


def cands_log2(p: Params, k: int) -> float:
    return max(0.0, 3 * k - 1 - p.n)


def per_target_ns(p: Params, k: int) -> float:
    return (p.t_lin_ns - p.t_cand_ns) + 2.0 ** cands_log2(p, k) * p.t_cand_ns


def mu_log2(p: Params, k: int) -> float:
    return p.m * log2fb(k) - math.log2(math.factorial(p.m)) - p.n


def relation_stage_log2s(p: Params, k: int) -> float:
    """log2 core-seconds of the relation stage."""
    t_new = p.t_add_ns if p.t_new_ns is None else p.t_new_ns
    scan = 2.0 ** ((p.m - 2) * log2fb(k)) * per_target_ns(p, k) * 1e-9 + t_new * 1e-9
    per_rel = math.log2(scan) - mu_log2(p, k)
    need = log2fb(k) + math.log2(1 + p.slack) - math.log2(p.neg_identify)
    return per_rel + need


def la_stage_log2s(p: Params, k: int) -> float:
    N = log2fb(k) - math.log2(p.neg_identify * p.fb_reduction)
    return math.log2(2) + N + (N + math.log2(p.m)) + math.log2(p.t_la_ns * 1e-9)


def la_memory_log2bytes(p: Params, k: int) -> float:
    N = 2.0 ** (log2fb(k) - math.log2(p.neg_identify * p.fb_reduction))
    idx = math.ceil(k / 8)
    vec = math.ceil(p.n / 8)
    return math.log2(N * p.m * (idx + 1) + 3 * N * vec)


def rho_iterations_log2(p: Params) -> float:
    ell = p.n - p.cofactor_bits
    sym = math.log2(2 * (p.n if p.frobenius else 1))
    return 0.5 * (math.log2(math.pi / 2) + ell - sym)


def rho_stage_log2s(p: Params) -> float:
    return rho_iterations_log2(p) + math.log2(p.rho_adds * p.t_add_ns * 1e-9)


def evaluate(p: Params, k: int) -> dict:
    rel, la, rho = relation_stage_log2s(p, k), la_stage_log2s(p, k), rho_stage_log2s(p)
    total = math.log2(2.0 ** rel + 2.0 ** la)
    return dict(n=p.n, m=p.m, k=k, log2_fb=log2fb(k), mu_log2=mu_log2(p, k), cands_log2=cands_log2(p, k),
                log2_relations_s=rel, log2_algebra_s=la, log2_total_s=total, log2_rho_s=rho,
                log2_gap=total - rho, log2_la_mem_bytes=la_memory_log2bytes(p, k),
                log2_total_core_years=total - LOG2_YEAR, log2_rho_core_years=rho - LOG2_YEAR,
                rho_iterations_log2=rho_iterations_log2(p))


def optimize_k(p: Params, kmax: int | None = None) -> dict:
    """Best k for the total cost (k must give mu large enough for the stage to make sense: any k is allowed, mu < 1
    just means several scans per relation)."""
    ks = range(2, (kmax or p.n // 2) + 1)
    return min((evaluate(p, k) for k in ks), key=lambda r: r["log2_total_s"])


def sweep(p: Params, ks) -> list[dict]:
    return [evaluate(p, k) for k in ks]


def free_scan_floor(p: Params) -> dict:
    """The best total when the relation finder is free (t_lin = t_cand = 0): what remains is forming targets,
    t_new / mu per relation, against the algebra stage; the optimum balances them near k = n / (m + 1)."""
    return optimize_k(replace(p, t_lin_ns=1e-12, t_cand_ns=1e-13), kmax=p.n // 2)


def required_speedup(p: Params) -> dict:
    """At the optimal k: by what factor must t_lin and t_cand shrink (all else fixed) for the total to equal rho?
    Bisection on log2 of the factor, re-optimizing k each time.  If even a free relation finder (`floor`) is
    above rho the factor is infinite: the algebra stage and target formation alone exceed rho for this m."""
    rho = rho_stage_log2s(p)
    best = optimize_k(p)
    floor = free_scan_floor(p)
    if best["log2_total_s"] <= rho:
        return dict(factor_log2=0.0, best=best, floor=floor)
    if floor["log2_total_s"] >= rho:
        return dict(factor_log2=math.inf, best=best, floor=floor)
    q = lambda f: optimize_k(replace(p, t_lin_ns=p.t_lin_ns * 2.0 ** -f, t_cand_ns=p.t_cand_ns * 2.0 ** -f))
    lo, hi = 0.0, 80.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if q(mid)["log2_total_s"] > rho:
            lo = mid
        else:
            hi = mid
    r = q(hi)
    return dict(factor_log2=hi, best=best, floor=floor, k_at_break_even=r["k"],
                t_lin_ns_needed=p.t_lin_ns * 2.0 ** -hi)


def oracle_budget(p: Params) -> dict:
    """Solver-agnostic requirement: the longest time (core-seconds) a relation finder may spend on ONE target R, over
    all its decompositions, for the total to equal rho.  For each k: the relations stage may use 2^rho - 2^algebra
    seconds, divided by the relations needed, times mu (decompositions per target, so one relation costs
    budget / mu).  Maximised over k.  Any solver can be compared with this number at its own (m, k): a Groebner basis
    of the Weil descent, a root finder, a SAT call.  `feasible` is False when even forming the target takes longer."""
    rho = 2.0 ** rho_stage_log2s(p)
    t_new = (p.t_add_ns if p.t_new_ns is None else p.t_new_ns) * 1e-9
    best = None
    for k in range(2, p.n // 2 + 1):
        la = 2.0 ** la_stage_log2s(p, k)
        if la >= rho:
            continue
        need = 2.0 ** (log2fb(k) + math.log2(1 + p.slack) - math.log2(p.neg_identify))
        scan = (rho - la) / need * 2.0 ** mu_log2(p, k)
        if best is None or scan > best["scan_s"]:
            best = dict(k=k, scan_s=scan, mu_log2=mu_log2(p, k), relations_log2=math.log2(need),
                        feasible=scan > t_new, t_new_s=t_new)
    return best or dict(k=None, scan_s=0.0, feasible=False, t_new_s=t_new)


def from_records(path, form="rows", threads=1) -> dict:
    """t_lin_ns and t_add_ns from a record in the format of experiments/ic-rootfind/n131.jsonl (median over rows)."""
    lin, add = [], []
    for line in Path(path).read_text().splitlines():
        r = json.loads(line)
        if r.get("form") == form and r.get("threads") == threads:
            lin.append(r["ns_per_target"])
            if "brute_ns_per_add" in r:
                add.append(r["brute_ns_per_add"])
    med = lambda xs: sorted(xs)[len(xs) // 2]
    out = {}
    if lin:
        out["t_lin_ns"] = med(lin)
    if add:
        out["t_add_ns"] = med(add)
    return out
