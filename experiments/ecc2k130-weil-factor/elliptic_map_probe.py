#!/usr/bin/env python3
"""Bounded short elliptic-map probe of the explicit polarization only."""
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import flint


def worker(matrix_path, output_path, gram_mode="exact"):
    here = Path(__file__).parent
    record = json.loads((here / "explicit_ppav.json").read_text())
    matrix_bytes = Path(matrix_path).read_bytes()
    expected = record["polarization"]["principal_matrix_sha256_decimal_rows"]
    assert hashlib.sha256(matrix_bytes).hexdigest() == expected
    principal = [[int(c) for c in row.split(",")]
                 for row in matrix_bytes.decode().splitlines()]
    n, pivot = 260, record["lattice"]["pivot"]
    v, a = record["lattice"]["generator_numerator"], record["lattice"]["degree"]
    omega = flint.fmpz_mat(principal)
    column = [a*principal[i][pivot]-sum(principal[i][j]*v[j]
              for j in range(n) if j != pivot) for i in range(n)]
    for i in range(n):
        omega[i, pivot], omega[pivot, i] = column[i], -column[i]
    omega[pivot, pivot] = 0
    # Symmetric quadratic form q(x)=2*Omega(x,alpha*x).
    gram = flint.fmpz_mat(n, n)
    for i in range(n):
        for j in range(130):
            gram[i, j] = omega[i, j]+2*omega[i, 130+j]
            gram[i, 130+j] = -4*omega[i, j]-omega[i, 130+j]
    assert gram == gram.transpose()
    v263 = [c % 263 for c in v]
    assert v263[pivot] == 1
    mixed = gram*flint.fmpz_mat(n, 1, v263)
    diagonal = sum(v263[i]*int(mixed[i, 0]) for i in range(n))
    assert diagonal % 263**2 == 0
    assert all(int(c) % 263 == 0 for c in mixed.entries())
    for i in range(n):
        gram[i, pivot] = gram[pivot, i] = int(mixed[i, 0])//263
    gram[pivot, pivot] = diagonal//263**2
    started = time.monotonic()
    reduced, transform = gram.lll(transform=True, delta=0.99, eta=0.51,
                                  rep="gram", gram=gram_mode)
    elapsed = time.monotonic()-started
    assert transform*gram*transform.transpose() == reduced
    degrees = [int(reduced[i, i])//2 for i in range(n)]
    assert all(int(reduced[i, i]) % 2 == 0 for i in range(n))
    assert all(d > 0 for d in degrees)
    hit = next((i for i, d in enumerate(degrees) if d in (1, 2, 3, 5, 7, 11)), None)
    result = {
        "schema": "ecc2k130.elliptic-map-probe.v1",
        "principal_matrix_sha256": expected,
        "status": "polarization excluded" if hit is not None else "inconclusive",
        "lll_seconds": elapsed,
        "lll_gram_mode": gram_mode,
        "smallest_observed_basis_degrees": sorted(degrees)[:30],
        "gram_transform_identity_verified": True,
        "scope": "one principal polarization, not the whole isogeny class",
    }
    if hit is not None:
        vector = [int(transform[hit, j]) for j in range(n)]
        exact = sum(vector[i]*int(gram[i, j])*vector[j]
                    for i in range(n) for j in range(n))
        assert exact == 2*degrees[hit]
        result.update({"elliptic_degree": degrees[hit], "vector_in_J_basis": vector,
                       "genus_upper_bound": degrees[hit]**2+1})
    Path(output_path).write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "exact")
    else:
        here = Path(__file__).parent
        approximate = sys.argv[1:] == ["--approx"]
        gram_mode, cap = ("approx", 60) if approximate else ("exact", 120)
        output = here / ("elliptic_map_probe_approx.json" if approximate
                         else "elliptic_map_probe.json")
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "principal.csv"
            subprocess.run([sys.executable, str(here / "construct_ppav.py"),
                            "--matrix-out", str(matrix)], check=True)
            started = time.monotonic()
            try:
                subprocess.run([sys.executable, __file__, "--worker", str(matrix), str(output),
                                gram_mode], check=True, timeout=cap)
            except subprocess.TimeoutExpired:
                output.write_text(json.dumps({
                    "schema": "ecc2k130.elliptic-map-probe.v1", "status": "inconclusive timeout",
                    "worker_cap_seconds": cap, "elapsed_seconds": time.monotonic()-started,
                    "lll_gram_mode": gram_mode,
                    "scope": "one principal polarization, not the whole isogeny class",
                }, indent=2)+"\n")
        print(output.read_text())
