"""One Sage/Singular Boolean Groebner query; stdin and stdout are JSON."""

import json
import resource
import sys
import time

from sage.all import GF, PolynomialRing


def main():
    case = json.load(sys.stdin)
    n, k, delta = (case[key] for key in ("n", "k", "delta"))
    names = [f"g{i}" for i in range(delta)] + [f"z{i}" for i in range(k)]
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    ring = PolynomialRing(GF(2), names=names, order="degrevlex")
    gens = ring.gens()
    g, z = gens[:delta], gens[delta:]
    equations = []
    for bit in range(n):
        eq = ring((case["e2_base"] >> bit) & 1)
        for j in range(delta):
            if case["e2_delta"][j] >> bit & 1:
                eq += g[j]
        for i in range(k):
            if (case["basis_square"][i] ^ case["e1_base_mul"][i]) >> bit & 1:
                eq += z[i]
            for j in range(delta):
                if case["e1_delta_mul"][j][i] >> bit & 1:
                    eq += g[j] * z[i]
        equations.append(eq)
    equations.extend(v * v + v for v in gens)
    construction_wall = time.perf_counter() - start_wall
    construction_cpu = time.process_time() - start_cpu

    gb_start_wall = time.perf_counter()
    gb_start_cpu = time.process_time()
    ideal = ring.ideal(equations)
    basis = ideal.groebner_basis()
    gb_wall = time.perf_counter() - gb_start_wall
    gb_cpu = time.process_time() - gb_start_cpu

    extraction_start_wall = time.perf_counter()
    extraction_start_cpu = time.process_time()
    solutions = ideal.variety()
    pairs = []
    for solution in solutions:
        gmask = sum(int(solution[v]) << j for j, v in enumerate(g))
        zmask = sum(int(solution[v]) << i for i, v in enumerate(z))
        pairs.append([gmask, zmask])
    pairs.sort()
    extraction_wall = time.perf_counter() - extraction_start_wall
    extraction_cpu = time.process_time() - extraction_start_cpu
    json.dump({
        "pairs": pairs,
        "equation_count": len(equations),
        "input_degree_max": max(int(eq.total_degree()) for eq in equations),
        "basis_count": len(basis),
        "max_intermediate_solver_degree": None,
        "construction_wall": construction_wall,
        "construction_cpu": construction_cpu,
        "groebner_wall": gb_wall,
        "groebner_cpu": gb_cpu,
        "extraction_wall": extraction_wall,
        "extraction_cpu": extraction_cpu,
        "query_wall": time.perf_counter() - start_wall,
        "query_cpu": time.process_time() - start_cpu,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }, sys.stdout, sort_keys=True)


if __name__ == "__main__":
    main()
