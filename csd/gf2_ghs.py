"""F_(2^N) on Python integers for ic_dreg / fb_structure: carry-less products, F_2-linear Frobenius maps, trace.

The subset of `csd/binary_ghs.py` (lineage `csd-milestone-1`) that `ic_dreg` and
`fb_structure` use, without its PARI-dependent parts: the modulus is always given.
"""
from __future__ import annotations


def _spread_table():
    t = []
    for v in range(256):
        r = 0
        for i in range(8):
            if v >> i & 1:
                r |= 1 << (2 * i)
        t.append(r)
    return t


_SPREAD = _spread_table()


def clmul(x: int, y: int) -> int:
    """Carry-less product of two bit polynomials."""
    if x.bit_length() < y.bit_length():
        x, y = y, x
    r = 0
    while y:
        low = y & -y
        r ^= x << (low.bit_length() - 1)
        y ^= low
    return r


def poly_mod(x: int, m: int) -> int:
    dm = m.bit_length() - 1
    while x.bit_length() - 1 >= dm:
        x ^= m << (x.bit_length() - 1 - dm)
    return x


class LinearMap:
    """An F_2-linear map on N-bit vectors, applied with byte tables."""

    def __init__(self, images: list[int]):
        self.images = images
        N = len(images)
        self.tables = []
        for k in range(0, N, 8):
            imgs = images[k:k + 8]
            t = [0] * 256
            for v in range(1, 256):
                low = v & -v
                i = low.bit_length() - 1
                t[v] = t[v ^ low] ^ (imgs[i] if i < len(imgs) else 0)
            self.tables.append(t)

    def __call__(self, x: int) -> int:
        r = 0
        for t in self.tables:
            if x & 255:
                r ^= t[x & 255]
            x >>= 8
            if not x:
                break
        return r

    def compose(self, other: "LinearMap") -> "LinearMap":
        """self o other."""
        return LinearMap([self(v) for v in other.images])


class GF2n:
    """F_(2^N) = F_2[w] / (modulus). Elements are ints < 2^N."""

    def __init__(self, N: int, modulus: int):
        self.N = N
        self.modulus = modulus
        assert self.modulus.bit_length() - 1 == N
        self.mask = (1 << N) - 1
        self._frob: dict[int, LinearMap] = {}
        self._frob[1] = LinearMap([self._sq_slow(1 << i) for i in range(N)])
        self._frob[0] = LinearMap([1 << i for i in range(N)])
        self.trace_mask = self._trace_mask()

    def reduce(self, x: int) -> int:
        return poly_mod(x, self.modulus)

    def _sq_slow(self, x: int) -> int:
        r, s = 0, 0
        while x:
            r |= _SPREAD[x & 255] << s
            x >>= 8
            s += 16
        return self.reduce(r)

    def mul(self, x: int, y: int) -> int:
        return self.reduce(clmul(x, y))

    def sq(self, x: int) -> int:
        return self._frob[1](x)

    def frob_map(self, e: int) -> LinearMap:
        """The map x -> x^(2^e), e taken mod N, by binary composition."""
        e %= self.N
        if e not in self._frob:
            if e % 2 == 0:
                h = self.frob_map(e // 2)
                self._frob[e] = h.compose(h)
            else:
                self._frob[e] = self._frob[1].compose(self.frob_map(e - 1))
        return self._frob[e]

    def frob(self, x: int, e: int) -> int:
        return self.frob_map(e)(x)

    def sqrt(self, x: int) -> int:
        return self.frob(x, self.N - 1)

    def inv(self, x: int) -> int:
        """Extended Euclid over F_2[w]."""
        if x == 0:
            raise ZeroDivisionError
        r0, r1, s0, s1 = self.modulus, x, 0, 1
        while r1:
            d = r0.bit_length() - r1.bit_length()
            if d < 0:
                r0, r1, s0, s1 = r1, r0, s1, s0
                continue
            r0 ^= r1 << d
            s0 ^= s1 << d
            if r0.bit_length() < r1.bit_length():
                r0, r1, s0, s1 = r1, r0, s1, s0
        assert r0 == 1
        return self.reduce(s0)

    def _trace_mask(self) -> int:
        """Bit i set iff Tr(w^i) = 1: power sums of the roots of the modulus by Newton's
        identities over F_2."""
        N, m = self.N, self.modulus
        e = [0] * (N + 1)
        for j in range(1, N + 1):
            e[j] = m >> (N - j) & 1
        p = [0] * N
        p[0] = N & 1
        for k in range(1, N):
            s = (k * e[k]) & 1
            for j in range(1, k):
                s ^= e[j] & p[k - j]
            p[k] = s
        return sum(p[i] << i for i in range(N))

    def trace(self, x: int) -> int:
        return bin(x & self.trace_mask).count("1") & 1
