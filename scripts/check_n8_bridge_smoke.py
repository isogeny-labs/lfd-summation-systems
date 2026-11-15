#!/usr/bin/env python3
"""Independent finite smoke tests, NOT a proof of unconditioned LFD growth."""
from itertools import combinations, product
from math import comb


def rank(vectors, p):
    rows = [list(v) for v in vectors]
    if not rows:
        return 0
    count = 0
    for col in range(len(rows[0])):
        pivot = next((j for j in range(count, len(rows)) if rows[j][col] % p), None)
        if pivot is None:
            continue
        rows[count], rows[pivot] = rows[pivot], rows[count]
        inv = pow(rows[count][col] % p, -1, p)
        rows[count] = [x * inv % p for x in rows[count]]
        for j in range(len(rows)):
            if j != count:
                mul = rows[j][col] % p
                rows[j] = [(x - mul * y) % p for x, y in zip(rows[j], rows[count])]
        count += 1
        if count == len(rows):
            break
    return count


def check_splitting():
    p = 7
    vectors = ((1,0,0,0), (0,1,0,0), (1,1,0,0),
               (0,0,1,0), (0,0,0,1), (0,0,1,1))
    T = set(range(3))
    circuits = []
    for size in range(2, len(vectors) + 1):
        for ix in combinations(range(len(vectors)), size):
            if rank([vectors[i] for i in ix], p) == size:
                continue
            if all(rank([vectors[i] for i in ix if i != j], p) == size - 1
                   for j in ix):
                assert set(ix).issubset(T) or set(ix).isdisjoint(T)
                circuits.append(ix)
    assert set(circuits) == {(0,1,2), (3,4,5)}
    count = 0
    for coefficients in product(range(p), repeat=len(vectors)):
        if any(sum(coefficients[i] * vectors[i][j] for i in range(len(vectors))) % p
               for j in range(4)):
            continue
        assert all(sum(coefficients[i] * vectors[i][j] for i in T) % p == 0
                   for j in range(4))
        count += 1
    assert count == p**2
    return count


def check_lambda():
    expected = {1: 25, 2: 36, 3: 48, 4: 59, 5: 71}
    for w, threshold in expected.items():
        upper = 13 * w + 23
        def check(ell):
            C = sum(comb(2 * ell, j) for j in range(min(w, 2 * ell) + 1))
            return 2 ** (ell - 1) < 7776 * C**2
        assert threshold == max(ell for ell in range(1, upper+1) if check(ell))
        assert all(not check(ell) for ell in range(threshold+1, upper+1))
    return expected


def check_unit_mod_two():
    terms = ((1,1,2), (0,1,-2), (1,0,-2), (0,0,-5))
    count = 0
    for N in (1,2,4,8,16):
        for i in range(N):
            for j in range(N):
                col = {}
                for u, v, coeff in terms:
                    idx = ((i + u) % N, (j + v) % N)
                    col[idx] = (col.get(idx, 0) + coeff) % 2
                assert {idx for idx, x in col.items() if x} == {(i,j)}
                count += 1
    return count


if __name__ == "__main__":
    print("PASS circuit splitting:", check_splitting(), "vanishing coefficients")
    print("PASS Lambda counting:", check_lambda())
    print("PASS mod-2 group ring:", check_unit_mod_two(), "columns")
