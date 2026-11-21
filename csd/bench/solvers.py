"""Solver registry.

A solver is anything that, given an Instance, returns the set of decompositions R = P_1 + ... + P_m with every x(P_i)
in V (or, for algebraic solvers, a closing degree).  Each solver separates `prepare` (tables, descent; timed as
t_setup) from `run` (the search; timed `reps` times).  Add a solver by subclassing Solver and decorating it with
@register; the runner needs nothing else.  `available()` reports missing tools with a reason instead of failing.
"""
from __future__ import annotations

import math
import time

from csd import ic_dreg as I

REGISTRY: dict[str, "Solver"] = {}

# Predicted work above which the enumeration-type solvers are skipped instead of run into the timeout.
MAX_ENUM_ADDS = 2 ** 33
MAX_PYREF_BITS = 17
MAX_MXL_VARS = 12             # 1.7 s at 12 Boolean variables, 94 s at 15 (experiments/ic-rootfind/README.md)


def register(cls):
    REGISTRY[cls.name] = cls()
    return cls


def canon(decs):
    """Canonical form of a set of decompositions: sorted tuple of sorted tuples of (x, y)."""
    return tuple(sorted(tuple(sorted((int(x), int(y)) for x, y in d)) for d in decs))


class Solver:
    name = "?"
    kind = "decomposition"        # or "algebraic": returns no decompositions, only a closing degree
    supports_threads = False
    oracle_rank = 99              # lower is more trusted as the reference in cross-checks

    def available(self) -> tuple[bool, str]:
        return True, ""

    def applicable(self, spec) -> tuple[bool, str]:
        return True, ""

    def prepare(self, inst):
        return None

    def run(self, ctx, inst, threads: int = 1):
        """-> (decompositions | None, counters dict)"""
        raise NotImplementedError


def _fb_estimate(spec):
    return 2.0 ** spec.k


CAPS = (1 << 12, 1 << 16, 1 << 20)     # row capacities tried in turn (see _with_cap)


def _with_cap(call):
    """csd.ic_rootfind's methods allocate an output buffer of `cap` rows per call (default 2^20 rows, 48 MB at m = 3).
    That allocation is not part of the algorithm: glibc serves it by mmap on the first call and from the heap with a
    48 MB memset afterwards, 0.07 to 2 ms depending on allocator state, which dominates the 0.1 to 1 ms searches
    of k <= 11 (experiments/ic-bench/RESULTS.md, C1).  Start small and retry on overflow; the cap used is recorded."""
    for cap in CAPS:
        try:
            return call(cap), cap
        except OverflowError:
            continue
    raise OverflowError("more rows than the largest cap")


class _RootfindBase(Solver):
    def available(self):
        try:
            from csd import ic_rootfind
            ic_rootfind._lib()
        except Exception as e:                      # no compiler, no PCLMULQDQ, ...
            return False, f"{type(e).__name__}: {e}"
        return True, ""

    def applicable(self, spec):
        if spec.n > 63:
            return False, "C kernel needs n <= 63 (the Rust kernel in rust/ic_rootfind covers n <= 191)"
        return True, ""

    def prepare(self, inst):
        from csd.ic_rootfind import Solver as RF
        return RF(inst.E, inst.basis)


@register
class LRF(_RootfindBase):
    name = "lrf"
    supports_threads = True
    oracle_rank = 3

    def run(self, ctx, inst, threads=1):
        (decs, st, rows), cap = _with_cap(lambda c: ctx.lrf(inst.spec.m, inst.R, cap=c, threads=threads))
        return decs, dict(st, nfb=ctx.nfb, dimV2=ctx.dimV2, rows_found=rows, cap=cap)


@register
class Enum(_RootfindBase):
    name = "enum"
    supports_threads = True
    oracle_rank = 2

    def applicable(self, spec):
        ok, why = super().applicable(spec)
        if not ok:
            return ok, why
        adds = _fb_estimate(spec) ** (spec.m - 1)
        if adds > MAX_ENUM_ADDS:
            return False, f"predicted ~2^{math.log2(adds):.0f} additions > 2^{math.log2(MAX_ENUM_ADDS):.0f}"
        return True, ""

    def run(self, ctx, inst, threads=1):
        (decs, cnt), cap = _with_cap(lambda c: ctx.brute(inst.spec.m, inst.R, cap=c, threads=threads))
        return decs, dict(nfb=ctx.nfb, rows_found=cnt, cap=cap)


@register
class MITM4(_RootfindBase):
    name = "mitm4"
    oracle_rank = 4

    def applicable(self, spec):
        ok, why = super().applicable(spec)
        if not ok:
            return ok, why
        if spec.m != 4:
            return False, "defined for m = 4 only"
        return True, ""

    def run(self, ctx, inst, threads=1):
        (decs, cnt, table_bytes), cap = _with_cap(lambda c: ctx.mitm4(inst.R, cap=c))
        return decs, dict(nfb=ctx.nfb, rows_found=cnt, table_bytes=table_bytes, cap=cap)


@register
class PyRef(_RootfindBase):
    """Plain enumeration with csd.ic_dreg's group law: slow, independent of the C kernels, the preferred oracle."""
    name = "py_ref"
    oracle_rank = 1

    def applicable(self, spec):
        ok, why = super().applicable(spec)
        if not ok:
            return ok, why
        if spec.k * (spec.m - 1) > MAX_PYREF_BITS:
            return False, f"k(m-1) = {spec.k * (spec.m - 1)} > {MAX_PYREF_BITS}"
        return True, ""

    def run(self, ctx, inst, threads=1):
        return ctx.decompose_reference(inst.spec.m, inst.R), dict(nfb=ctx.nfb)


@register
class MutantXL(Solver):
    """The repository's MutantXL on the Weil descent of S_{m+1} (csd.ic_dreg): the algebraic route.  Returns no
    decompositions (the Boolean system has more solutions than R has decompositions); `closed` says whether the
    closure reached every linear form vanishing on the solution set, and `d_mxl` the degree."""
    name = "mutantxl"
    kind = "algebraic"

    def applicable(self, spec):
        v = spec.m * spec.k
        if v > MAX_MXL_VARS:
            return False, f"{v} Boolean variables > {MAX_MXL_VARS}"
        if spec.basis == "invariant" or spec.n > 32:
            return False, "csd.ic_dreg's descent covers n <= 32"
        return True, ""

    def prepare(self, inst):
        m, k = inst.spec.m, inst.spec.k
        B = [inst.basis] * m
        eqs = I.binary_descent(inst.E, B, inst.R[0])
        sols = I.binary_solutions(inst.E, B, inst.R[0])
        pts = [sum(c[i] << (i * k) for i in range(m)) for c in sols]
        v = m * k
        return dict(eqs=eqs, v=v, target=I._lin_target_f2(pts, v), n_solutions=len(sols),
                    degs=sorted({max(bin(mk).count("1") for mk in eq) for eq in eqs}))

    def run(self, ctx, inst, threads=1):
        eqs, v, target = ctx["eqs"], ctx["v"], ctx["target"]
        per_degree, d_mxl = [], None
        for D in range(max(ctx["degs"]), v + 1):
            t = time.perf_counter()
            M, cols = I.boolean_macaulay(eqs, v, D)
            rank, piv, ops, rounds, added = I._mxl_closure(M, cols, v, D, 2, I._bool_mulvar)
            dt = time.perf_counter() - t
            lin = sum(1 for c in piv if c >= cols - (v + 1))
            per_degree.append(dict(D=D, cols=cols, rows=int(M.shape[0]), rank=int(rank), t=dt, ops=float(ops)))
            if lin == target:
                d_mxl = D
                break
        return None, dict(closed=d_mxl is not None, d_mxl=d_mxl, boolean_vars=v, boolean_solutions=ctx["n_solutions"], mxl_ops_final=per_degree[-1]["ops"],
                          max_matrix_cols=per_degree[-1]["cols"], per_degree=per_degree)


def warmup():
    """Compile and load the C kernel once in the parent, before any child is forked (so no child pays for it)."""
    try:
        from csd import ic_rootfind
        ic_rootfind._lib()
    except Exception:
        pass            # reported per solver by available()
