"""Row echelon forms over GF(2) and F_p with pivot columns and counted operations.

The kernels are in `_elim.c`, compiled once with the system C compiler into
``$CSD_CACHE`` (default ``~/.cache/csd``) and loaded through ctypes. Without a
compiler the numpy fallbacks below are used; they compute the same pivots and
the same operation counts, only slower (tests check that the two agree).

Both functions reduce a copy of the matrix, scanning columns left to right, and
return (rank, pivot columns, ops). `ops` counts row operations as the kernel
performs them (64-bit word XORs over GF(2), multiply-adds over F_p): it is a
counted cost of this elimination, not a timing.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import numpy as np

_SRC = Path(__file__).with_name("_elim.c")


@lru_cache(maxsize=1)
def _lib():
    if os.environ.get("CSD_ELIM_PURE"):
        return None
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        return None
    src = _SRC.read_bytes()
    tag = hashlib.sha256(src).hexdigest()[:16]
    cache = Path(os.environ.get("CSD_CACHE", Path.home() / ".cache" / "csd"))
    cache.mkdir(parents=True, exist_ok=True)
    so = cache / f"elim_{tag}.so"
    if not so.exists():
        tmp = so.with_suffix(f".{os.getpid()}.tmp")
        subprocess.run([cc, "-O3", "-march=native", "-shared", "-fPIC", "-o", str(tmp), str(_SRC)],
                       check=True, capture_output=True)
        os.replace(tmp, so)
    lib = ctypes.CDLL(str(so))
    p = ctypes.c_void_p
    lib.gf2_echelon.argtypes = [p, ctypes.c_long, ctypes.c_long, ctypes.c_long, p, p]
    lib.gf2_echelon.restype = ctypes.c_long
    lib.fp_echelon.argtypes = [p, ctypes.c_long, ctypes.c_long, ctypes.c_uint32, p, p]
    lib.fp_echelon.restype = ctypes.c_long
    return lib


def compiled() -> bool:
    return _lib() is not None


def pack_rows(rows: list[list[int]], cols: int) -> np.ndarray:
    """Bit-pack rows given as lists of column indices (each column at most once)."""
    words = max(1, (cols + 63) // 64)
    M = np.zeros((len(rows), words), dtype=np.uint64)
    for i, r in enumerate(rows):
        if r:
            c = np.asarray(r, dtype=np.int64)
            np.bitwise_xor.at(M[i], c >> 6, np.left_shift(np.uint64(1), (c & 63).astype(np.uint64)))
    return M


def gf2_echelon(M: np.ndarray, cols: int, out: list | None = None) -> tuple[int, list[int], float]:
    """Row echelon form of a bit-packed GF(2) matrix (uint64, shape rows x words).
    If `out` is a list, the echelon rows (rank x words) are appended to it."""
    A = np.ascontiguousarray(M, dtype=np.uint64).copy()
    if out is not None:
        r = _gf2(A, cols)
        out.append(A[:r[0]])
        return r
    return _gf2(A, cols)


def _gf2(A: np.ndarray, cols: int) -> tuple[int, list[int], float]:
    rows, words = A.shape
    if rows == 0:
        return 0, [], 0.0
    lib = _lib()
    if lib is not None:
        piv = np.zeros(min(rows, cols), dtype=np.int64)
        ops = ctypes.c_double(0)
        r = lib.gf2_echelon(A.ctypes.data, rows, words, cols, piv.ctypes.data, ctypes.byref(ops))
        return int(r), piv[:r].tolist(), float(ops.value)
    rank, piv, ops = 0, [], 0.0
    for c in range(cols):
        if rank == rows:
            break
        w, bit = c >> 6, np.uint64(1) << np.uint64(c & 63)
        col = (A[rank:, w] & bit) != 0
        nz = np.flatnonzero(col)
        if nz.size == 0:
            continue
        r = rank + int(nz[0])
        if r != rank:
            A[[rank, r], w:] = A[[r, rank], w:]
        below = rank + 1 + np.flatnonzero((A[rank + 1:, w] & bit) != 0)
        if below.size:
            A[below, w:] ^= A[rank, w:]
            ops += float(below.size * (words - w))
        piv.append(c)
        rank += 1
    return rank, piv, ops


def fp_echelon(M: np.ndarray, p: int, out: list | None = None) -> tuple[int, list[int], float]:
    """Row echelon form over F_p (p < 2^31) of a dense matrix of residues, pivots
    normalised to 1. If `out` is a list, the echelon rows are appended to it."""
    A = np.ascontiguousarray(np.asarray(M, dtype=np.int64) % p, dtype=np.uint32).copy()
    r = _fp(A, p)
    if out is not None:
        out.append(r[3][:r[0]].astype(np.int64))
    return r[:3]


def _fp(A: np.ndarray, p: int):
    rows, cols = A.shape
    if rows == 0:
        return 0, [], 0.0, A
    lib = _lib()
    if lib is not None:
        piv = np.zeros(min(rows, cols), dtype=np.int64)
        ops = ctypes.c_double(0)
        r = lib.fp_echelon(A.ctypes.data, rows, cols, p, piv.ctypes.data, ctypes.byref(ops))
        return int(r), piv[:r].tolist(), float(ops.value), A
    B = A.astype(np.int64)
    rank, piv, ops = 0, [], 0.0
    for c in range(cols):
        if rank == rows:
            break
        nz = np.flatnonzero(B[rank:, c])
        if nz.size == 0:
            continue
        r = rank + int(nz[0])
        if r != rank:
            B[[rank, r], c:] = B[[r, rank], c:]
        B[rank, c:] = B[rank, c:] * pow(int(B[rank, c]), -1, p) % p
        below = rank + 1 + np.flatnonzero(B[rank + 1:, c])
        if below.size:
            f = B[below, c][:, None]
            B[below, c:] = (B[below, c:] - f * B[rank, c:]) % p
            ops += float(below.size * (cols - c))
        piv.append(c)
        rank += 1
    return rank, piv, ops, B
