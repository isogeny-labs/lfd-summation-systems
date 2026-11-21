"""Isolated, cross-checked benchmark runs.

Design (each point answers a way a benchmark number goes wrong):

* **Fork per run.**  The solver runs in a forked child: a hung or crashing kernel costs one row, not the sweep; the
  timeout is enforced by the parent; peak RSS is that of the run alone (/proc/self/clear_refs resets it).
* **Setup is not search.**  `prepare` (tables, descent) is timed once as t_setup; `run` is timed `reps` times
  (t_first is the first, cold one; t_min and t_median use all of them).  GC is off inside the timed region.
* **Order and drift.**  The order of solvers within an instance is a seeded random permutation (position recorded),
  so slow drift of the machine does not favour one solver; `aa=True` runs every solver twice as arms A and B so
  that `stats.noise_floor` can say which differences are below the resolution of the machine.
* **Correctness before speed.**  Every decomposition is re-added to R with csd.ic_dreg's group law; all solvers on
  an instance must return the same set (the trusted reference is the lowest `oracle_rank` that ran).  A row that
  fails either check has status `invalid` / `mismatch` and the driver exits non-zero.
* **Confounders are recorded**: the load average before each run, CPU pinning, threads, and per file a meta
  record (git commit and dirty flag, CPU, compiler, library versions, source hashes of the solvers).
"""
from __future__ import annotations

import gc
import hashlib
import json
import multiprocessing as mp
import os
import platform
import random
import resource
import shutil
import subprocess
import time
import traceback
from pathlib import Path

from . import solvers as S
from .spec import Spec, build

ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------------------------- environment
def _sh(*cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=20, cwd=ROOT).stdout.strip()
    except Exception:
        return ""


def environment() -> dict:
    cpu = ""
    try:
        cpu = next(l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name"))
    except Exception:
        pass
    gov = ""
    try:
        gov = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor").read_text().strip()
    except Exception:
        pass
    import numpy
    srcs = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
            for p in [ROOT / "csd" / "_rootfind.c", ROOT / "csd" / "ic_rootfind.py", ROOT / "csd" / "ic_dreg.py",
                      *sorted(Path(__file__).parent.glob("*.py"))]}
    return dict(commit=_sh("git", "rev-parse", "HEAD"), dirty=bool(_sh("git", "status", "--porcelain")),
                cpu=cpu, nproc=os.cpu_count(), affinity=sorted(os.sched_getaffinity(0)), governor=gov,
                python=platform.python_version(), numpy=numpy.__version__,
                cc=(_sh(shutil.which("cc") or "cc", "--version").splitlines() or [""])[0],
                loadavg_start=os.getloadavg(), utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), sources=srcs)


# ---------------------------------------------------------------------------------------------------- child
def _peak_rss_kib():
    try:
        for line in open("/proc/self/status"):
            if line.startswith("VmHWM:"):
                return int(line.split()[1]), "VmHWM"
    except Exception:
        pass
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, "ru_maxrss"


def _child(conn, solver_name, spec_dict, reps, threads, pin):
    try:
        if pin is not None:
            os.sched_setaffinity(0, set(pin))
        try:
            with open("/proc/self/clear_refs", "w") as f:
                f.write("5")                      # reset the peak-RSS watermark inherited from the parent
        except Exception:
            pass
        solver = S.REGISTRY[solver_name]
        inst = build(Spec(**spec_dict))
        gc.collect()
        gc.disable()
        t0 = time.perf_counter()
        ctx = solver.prepare(inst)
        t_setup = time.perf_counter() - t0
        wall, cpu, decs, counters = [], [], None, {}
        for _ in range(reps):
            t, c = time.perf_counter(), time.process_time()
            decs, counters = solver.run(ctx, inst, threads)
            wall.append(time.perf_counter() - t)
            cpu.append(time.process_time() - c)
        rss, how = _peak_rss_kib()
        conn.send(dict(ok=True, t_setup=t_setup, wall=wall, cpu=cpu, counters=counters, rss_kib=rss, rss_method=how,
                       decs=None if decs is None else S.canon(decs)))
    except BaseException:
        conn.send(dict(ok=False, error=traceback.format_exc()))
    finally:
        conn.close()


def _run_child(solver_name, spec, reps, threads, pin, timeout):
    ctx = mp.get_context("fork")
    parent, child = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_child, args=(child, solver_name, spec.as_dict(), reps, threads, pin))
    p.start()
    child.close()
    try:
        if parent.poll(timeout):
            try:
                msg = parent.recv()
            except EOFError:
                msg = None
            p.join(10)
            if msg is None:
                return dict(status="crash", reason=f"child exited with code {p.exitcode} and no result")
            if not msg["ok"]:
                return dict(status="error", reason=msg["error"].strip().splitlines()[-1], traceback=msg["error"])
            return dict(status="ran", **msg)
        return dict(status="timeout", reason=f"no result in {timeout} s")
    finally:
        if p.is_alive():
            p.kill()
        p.join()
        parent.close()


# ------------------------------------------------------------------------------------------------- checking
def verify(inst, decs) -> str | None:
    """None if every decomposition lies on the curve and sums to R (csd.ic_dreg's group law), else a message."""
    E = inst.E
    for d in decs:
        acc = None
        for P in d:
            if not E.on_curve(P):
                return f"point {P} not on the curve"
            acc = E.add(acc, P)
        if acc is None or tuple(acc) != tuple(inst.R):
            return "decomposition does not sum to R"
    return None


def _digest(decs) -> str:
    return hashlib.sha1(repr(decs).encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------- one cell
def run_instance(spec: Spec, solver_names, *, reps=3, timeout=120.0, threads=1, pin=None, aa=False, order_seed=0):
    """Run the solvers on one instance.  Returns a list of rows (dicts), one per (solver, arm)."""
    inst = build(spec)
    arms = [(s, a) for s in solver_names for a in (("A", "B") if aa else ("A",))]
    random.Random(f"order/{order_seed}/{spec.key()}").shuffle(arms)
    rows, results = [], {}
    for pos, (name, arm) in enumerate(arms):
        solver = S.REGISTRY[name]
        row = dict(spec.as_dict(), key=spec.key(), solver=name, arm=arm, pos=pos, threads=threads, reps=reps,
                   pin=None if pin is None else sorted(pin), load1=os.getloadavg()[0])
        ok, why = solver.available()
        if ok:
            ok, why = solver.applicable(spec)
        if not ok:
            rows.append(dict(row, status="skipped", reason=why))
            continue
        out = _run_child(name, spec, reps, threads if solver.supports_threads else 1, pin, timeout)
        if out["status"] != "ran":
            rows.append(dict(row, **{k: v for k, v in out.items() if k != "traceback"}))
            continue
        wall = out["wall"]
        row.update(status="ok", t_setup=out["t_setup"], times=wall, t_first=wall[0], t_min=min(wall),
                   t_median=sorted(wall)[len(wall) // 2], cpu_min=min(out["cpu"]), rss_kib=out["rss_kib"],
                   rss_method=out["rss_method"], counters=out["counters"], kind=solver.kind)
        if out["decs"] is not None:
            bad = verify(inst, out["decs"])
            row.update(n_decs=len(out["decs"]), decs_sha=_digest(out["decs"]))
            if bad:
                row.update(status="invalid", reason=bad)
            else:
                results.setdefault(name, out["decs"])
                if spec.fixture:
                    from . import fixtures
                    fbad = fixtures.check_result(spec.fixture, out["decs"])
                    row["fixture_check"] = "match" if fbad is None else fbad
                    if fbad:
                        row.update(status="mismatch", reason="fixture: " + fbad)
        if spec.fixture and row["status"] == "ok":
            from . import fixtures
            ref = (fixtures.row(spec.fixture).get("reference") or {})
            rt = ref.get({"lrf": "t_lrf", "enum": "t_brute", "mitm4": "t_mitm"}.get(name, ""))
            if rt:
                row.update(reference_t=rt, ratio_to_reference=row["t_min"] / rt)
        elif solver.kind == "algebraic" and not out["counters"].get("closed", False):
            row.update(status="no_close", reason="closure did not reach the target")
        rows.append(row)
    _cross_check(rows, results)
    return rows


def _cross_check(rows, results):
    """All decomposition sets must be equal; the lowest oracle_rank that ran is the reference."""
    if len(results) < 2:
        for r in rows:
            if r.get("status") == "ok" and r.get("kind") == "decomposition":
                r["checked_against"] = None
        return
    ref_name = min(results, key=lambda s: S.REGISTRY[s].oracle_rank)
    ref = results[ref_name]
    for r in rows:
        if r.get("status") != "ok" or r.get("kind") != "decomposition":
            continue
        r["checked_against"] = ref_name
        if results.get(r["solver"]) != ref:
            r["status"] = "mismatch"
            r["reason"] = f"{r['solver']} and {ref_name} return different sets"


# ------------------------------------------------------------------------------------------------- the grid
def run_grid(specs, solver_names, out: Path, *, reps=3, timeout=120.0, threads=1, pin=None, aa=False, resume=False,
             order_seed=0, log=print):
    """Write one JSONL row per (instance, solver, arm), flushing after each; `resume` skips rows already there."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if resume and out.exists():
        for line in out.read_text().splitlines():
            r = json.loads(line)
            done.add((r["key"], r["solver"], r["arm"]))
    S.warmup()
    meta = dict(environment(), argv_solvers=list(solver_names), reps=reps, timeout=timeout, threads=threads,
                pin=None if pin is None else sorted(pin), aa=aa, order_seed=order_seed, n_specs=len(specs))
    out.with_suffix(".meta.json").write_text(json.dumps(meta, indent=1, sort_keys=True) + "\n")
    bad = 0
    with open(out, "a" if resume else "w") as f:
        for spec in specs:
            arms = [(s, a) for s in solver_names for a in (("A", "B") if aa else ("A",))]
            if resume and all((spec.key(), s, a) in done for s, a in arms):
                continue
            try:
                rows = run_instance(spec, solver_names, reps=reps, timeout=timeout, threads=threads, pin=pin, aa=aa,
                                    order_seed=order_seed)
            except ValueError as e:            # e.g. no invariant subspace of this dimension
                rows = [dict(spec.as_dict(), key=spec.key(), solver="*", arm="A", status="skipped", reason=str(e))]
            for r in rows:
                f.write(json.dumps(r, sort_keys=True) + "\n")
                bad += r["status"] in ("mismatch", "invalid", "crash", "error")
            f.flush()
            log(_line(spec, rows))
    return bad


def _line(spec, rows):
    parts = []
    for r in rows:
        if r["status"] == "ok":
            parts.append(f"{r['solver']}{'' if r['arm'] == 'A' else '/B'}={r['t_min']:.3g}s")
        else:
            parts.append(f"{r['solver']}:{r['status']}")
    return f"{spec.key():44s} " + "  ".join(parts)
