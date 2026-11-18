"""Bounded derivation producer for chain closures; checker.py is independent.

Every retained row has an input/shift/linear-combination DAG. This is a
reference instrument, not a scalable replacement for the depth-14 engine.
"""
import hashlib
import json
from collections import deque


def cleaned(row, p):
    return {m: int(c) % p for m, c in row.items() if int(c) % p}


def shift(row, var, m, k, p):
    block, bit = divmod(var, k)
    off, mask = block * k, (1 << k) - 1
    return {((e & ~(mask << off)) | (((((e >> off) & mask) + (1 << bit)) & mask) << off)): c
            for e, c in row.items()}


def produce(generators, m, k, p, D, *, max_nodes=200000, max_rank=10000):
    nodes, values, basis, queue = [], [], {}, deque()
    def node(op, row):
        if len(nodes) >= max_nodes:
            raise ValueError('derivation node cap')
        nodes.append(op); values.append(row)
        return len(nodes)-1
    def absorb(row, source):
        terms = [[source, 1]]
        row = dict(row)
        for (_, pivot), (old, ident) in sorted(basis.items(), reverse=True):
            c = row.get(pivot, 0)
            if c:
                terms.append([ident, -c % p])
                for e, v in old.items():
                    row[e] = (row.get(e, 0)-c*v) % p
                row = cleaned(row, p)
        if not row:
            return
        # Degree descending, then mask descending: cancellation cannot raise degree.
        pivot = max(row, key=lambda e: (e.bit_count(), e))
        inv = pow(row[pivot], -1, p)
        row = {e: v*inv % p for e, v in row.items()}
        terms = [[i, c*inv % p] for i, c in terms if c % p]
        ident = node({'op': 'lin', 'terms': terms}, row)
        key = (pivot.bit_count(), pivot)
        basis[key] = (row, ident)
        if len(basis) > max_rank:
            raise ValueError('derivation rank cap')
        if pivot.bit_count() < D:
            queue.append(ident)
    # basis dictionary keys are filtration keys, not raw masks.
    # Convert lookup in elimination explicitly below via a custom mapping loop.
    for i, g in enumerate(generators):
        row = cleaned(g, p)
        if max((e.bit_count() for e in row), default=0) <= D:
            ident = node({'op': 'input', 'index': i}, row)
            absorb(row, ident)
    while queue:
        ident = queue.popleft()
        for var in range(m*k):
            row = shift(values[ident], var, m, k, p)
            source = node({'op': 'shift', 'source': ident, 'variable': var}, row)
            absorb(row, source)
    outputs = [ident for row, ident in basis.values()]
    rows = [sorted(values[i].items()) for i in outputs]
    return dict(format='chain-derivation-v1', m=m, k=k, p=p, D=D,
                generators=[sorted(cleaned(g,p).items()) for g in generators],
                nodes=nodes, outputs=outputs, rows=rows,
                row_sha256=hashlib.sha256(json.dumps(rows,separators=(',', ':')).encode()).hexdigest())
