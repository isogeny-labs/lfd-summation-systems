"""Tables of README.md from the records in this directory (median over the seeds of each size)."""
import json
from collections import defaultdict
from pathlib import Path
from statistics import median

D = Path(__file__).parent


def load(name):
    return [json.loads(l) for l in open(D / name)]


def groups(rows, key):
    g = defaultdict(list)
    for r in rows:
        g[key(r)].append(r)
    return sorted(g.items(), key=lambda kv: tuple(-1 if x is None else x for x in kv[0]))


def med(rs, f):
    xs = [r[f] for r in rs if f in r]
    return median(xs) if xs else None


def fmt(x, unit="s"):
    if x is None:
        return "—"
    if unit == "x":
        return f"{x:,.0f}×" if x >= 10 else f"{x:.1f}×"
    return f"{x:.3g} s" if x >= 1e-3 else f"{x * 1e6:.0f} µs"


out = []
rows = load("scaling-m3.jsonl")
out += ["### m = 3, V = <1, w, ..., w^(k-1)>, n = 3k - 2", "",
        "| n | k | \\|FB\\| | decompositions (sum over seeds) | LRF | enumeration | speed-up | LRF per target |",
        "|---|---|---|---|---|---|---|---|"]
for (n, k), rs in groups(rows, lambda r: (r["n"], r["k"])):
    out.append(f"| {n} | {k} | {int(med(rs, 'nfb'))} | {sum(r['decompositions'] for r in rs)} | {fmt(med(rs, 't_lrf'))} | "
               f"{fmt(med(rs, 't_brute'))} | {fmt(med(rs, 'speedup_vs_brute'), 'x')} | "
               f"{med(rs, 't_lrf') / med(rs, 'nfb') * 1e9:.0f} ns |")
rows = load("m4.jsonl")
out += ["", "### m = 4, V = <1, w, ..., w^(k-1)>, n = 4k - 3", "",
        "| n | k | decompositions | LRF | meet in the middle | its table | enumeration | LRF vs MITM | LRF vs enumeration |",
        "|---|---|---|---|---|---|---|---|---|"]
for (n, k), rs in groups(rows, lambda r: (r["n"], r["k"])):
    out.append(f"| {n} | {k} | {sum(r['decompositions'] for r in rs)} | {fmt(med(rs, 't_lrf'))} | {fmt(med(rs, 't_mitm'))} | "
               f"{med(rs, 'mitm_table_bytes') / 2**20:.1f} MiB | {fmt(med(rs, 't_brute'))} | "
               f"{fmt(med(rs, 'speedup_vs_mitm'), 'x')} | {fmt(med(rs, 'speedup_vs_brute'), 'x')} |")
rows = load("basis-n31-k10.jsonl")
out += ["", "### The factor base decides: n = 31, k = 10, m = 3, V = <1, ..., w^(j-1)> + (10 - j) random vectors", "",
        "| j | dim V2 | k + dim V2 - n | candidates per target | LRF | enumeration | speed-up |", "|---|---|---|---|---|---|---|"]
for (j,), rs in reversed(groups(rows, lambda r: (r["j"],))):
    cand = median(r["stats"]["lin_candidates"] / max(1, r["stats"]["targets"]) for r in rs)
    cand = "enumeration fallback" if rs[0]["stats"]["degenerate"] and not rs[0]["stats"]["targets"] else f"{cand:.2f}"
    out.append(f"| {j} | {int(med(rs, 'dimV2'))} | {10 + int(med(rs, 'dimV2')) - 31} | {cand} | {fmt(med(rs, 't_lrf'))} | "
               f"{fmt(med(rs, 't_brute'))} | {fmt(med(rs, 'speedup_vs_brute'), 'x')} |")
rows = load("regime-k12.jsonl")
out += ["", "### Below and above n = 3k: k = 12, m = 3, V = <1, w, ..., w^11>", "",
        "| n | unknowns | equations | targets with a solvable system | candidates per target | decompositions (sum) | LRF | "
        "enumeration | speed-up |", "|---|---|---|---|---|---|---|---|---|"]
for (n,), rs in groups(rows, lambda r: (r["n"],)):
    sv = median(r["stats"]["lin_solvable"] / r["stats"]["targets"] for r in rs)
    cand = median(r["stats"]["lin_candidates"] / r["stats"]["targets"] for r in rs)
    out.append(f"| {n} | {rs[0]['unknowns']} | {rs[0]['equations']} | {sv:.3g} | {cand:.3g} | "
               f"{sum(r['decompositions'] for r in rs)} | {fmt(med(rs, 't_lrf'))} | {fmt(med(rs, 't_brute'))} | "
               f"{fmt(med(rs, 'speedup_vs_brute'), 'x')} |")
rows = load("algebraic-m3.jsonl")
out += ["", "### Against the repository's MutantXL on the Weil descent (m = 3)", "",
        "| n | k | Boolean variables | Boolean solutions | decompositions of R | closing degree | MutantXL at that degree | "
        "MutantXL, all degrees | enumeration | LRF | LRF vs MutantXL |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for (n, k), rs in groups(rows, lambda r: (r["n"], r["k"])):
    out.append(f"| {n} | {k} | {3 * k} | {sorted(r['boolean_solutions'] for r in rs)} | {sorted(r['decompositions'] for r in rs)} | "
               f"{sorted({r['d_mxl'] for r in rs})} | {fmt(med(rs, 't_mxl_final'))} | {fmt(med(rs, 't_mxl_total'))} | "
               f"{fmt(med(rs, 't_brute'))} | {fmt(med(rs, 't_lrf'))} | {fmt(med(rs, 'speedup_vs_mxl_final'), 'x')} |")
(D / "tables.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
