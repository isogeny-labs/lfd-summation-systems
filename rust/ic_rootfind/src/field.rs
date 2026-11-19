//! F_(2^n), n <= 191, three 64-bit words, sparse modulus x^n + sum x^e.

pub const W: usize = 3;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Debug, Hash)]
pub struct Fe(pub [u64; W]);

impl Fe {
    pub const ZERO: Fe = Fe([0; W]);
    pub const ONE: Fe = Fe([1, 0, 0]);
    #[inline]
    pub fn is_zero(&self) -> bool {
        (self.0[0] | self.0[1] | self.0[2]) == 0
    }
    #[inline]
    pub fn bit(&self, i: usize) -> bool {
        (self.0[i >> 6] >> (i & 63)) & 1 == 1
    }
    pub fn unit(i: usize) -> Fe {
        let mut r = Fe::ZERO;
        r.0[i >> 6] = 1u64 << (i & 63);
        r
    }
    /// Index of the highest set bit.
    pub fn top(&self) -> Option<usize> {
        for i in (0..W).rev() {
            if self.0[i] != 0 {
                return Some(64 * i + 63 - self.0[i].leading_zeros() as usize);
            }
        }
        None
    }
    #[inline]
    pub fn and_parity(&self, m: &Fe) -> bool {
        (((self.0[0] & m.0[0]).count_ones() + (self.0[1] & m.0[1]).count_ones() + (self.0[2] & m.0[2]).count_ones()) & 1) == 1
    }
    pub fn hex(&self) -> String {
        let s = format!("{:016x}{:016x}{:016x}", self.0[2], self.0[1], self.0[0]);
        let t = s.trim_start_matches('0');
        if t.is_empty() { "0".into() } else { t.into() }
    }
    pub fn from_hex(s: &str) -> Fe {
        let v = u128::from_str_radix(&format!("{:0>48}", s)[16..], 16).unwrap();
        let hi = u64::from_str_radix(&format!("{:0>48}", s)[..16], 16).unwrap();
        Fe([v as u64, (v >> 64) as u64, hi])
    }
}

impl std::ops::BitXor for Fe {
    type Output = Fe;
    #[inline]
    fn bitxor(self, o: Fe) -> Fe {
        Fe([self.0[0] ^ o.0[0], self.0[1] ^ o.0[1], self.0[2] ^ o.0[2]])
    }
}
impl std::ops::BitXorAssign for Fe {
    #[inline]
    fn bitxor_assign(&mut self, o: Fe) {
        self.0[0] ^= o.0[0];
        self.0[1] ^= o.0[1];
        self.0[2] ^= o.0[2];
    }
}

#[cfg(all(target_arch = "x86_64", target_feature = "pclmulqdq"))]
#[inline(always)]
fn clmul(a: u64, b: u64) -> (u64, u64) {
    use std::arch::x86_64::*;
    unsafe {
        let p = _mm_clmulepi64_si128(_mm_set_epi64x(0, a as i64), _mm_set_epi64x(0, b as i64), 0);
        (_mm_cvtsi128_si64(p) as u64, _mm_extract_epi64(p, 1) as u64)
    }
}

#[cfg(not(all(target_arch = "x86_64", target_feature = "pclmulqdq")))]
#[inline(always)]
fn clmul(a: u64, b: u64) -> (u64, u64) {
    let (mut lo, mut hi) = (0u64, 0u64);
    for i in 0..64 {
        if (b >> i) & 1 == 1 {
            lo ^= a << i;
            if i > 0 {
                hi ^= a >> (64 - i);
            }
        }
    }
    (lo, hi)
}

#[derive(Clone, Debug)]
pub struct Field {
    pub n: usize,
    /// low exponents of the modulus (x^n + sum x^e), at most 4
    terms: [u32; 4],
    nterms: usize,
    pub low: Vec<u32>,
    lowmask: [u64; 3],
}

impl Field {
    /// The lowest-weight irreducible x^n + x^k + 1 with least k, else the least pentanomial
    /// (the rule of csd.ic_dreg.binary_modulus).
    pub fn new(n: usize) -> Field {
        assert!((2..=191).contains(&n));
        for k in 1..n.min(64) {
            let f = Field::with_terms(n, &[k as u32, 0]);
            if f.irreducible() {
                return f;
            }
        }
        for k3 in 3..n.min(64) {
            for k2 in 2..k3 {
                for k1 in 1..k2 {
                    let f = Field::with_terms(n, &[k3 as u32, k2 as u32, k1 as u32, 0]);
                    if f.irreducible() {
                        return f;
                    }
                }
            }
        }
        unreachable!()
    }

    pub fn with_terms(n: usize, low: &[u32]) -> Field {
        assert!(low.iter().all(|&e| e < 64), "modulus low terms must have degree below 64");
        let mut terms = [0u32; 4];
        terms[..low.len()].copy_from_slice(low);
        let mut lowmask = [0u64; 3];
        for &e in low {
            lowmask[(e >> 6) as usize] |= 1u64 << (e & 63);
        }
        Field { n, terms, nterms: low.len(), low: low.to_vec(), lowmask }
    }

    /// Rabin: x^(2^n) = x and gcd(x^(2^(n/q)) - x, f) = 1 for every prime q | n.  The gcd is taken in this
    /// "field" (arithmetic modulo f is a ring for any f), by Euclid on bit vectors.
    fn irreducible(&self) -> bool {
        let x = Fe::unit(1);
        if self.n == 1 {
            return true;
        }
        let mut t = x;
        let mut at = vec![Fe::ZERO; self.n + 1];
        for i in 1..=self.n {
            t = self.sq(t);
            at[i] = t;
        }
        if at[self.n] != x {
            return false;
        }
        let mut m = self.n;
        let mut q = 2;
        let mut primes = vec![];
        while q * q <= m {
            if m % q == 0 {
                primes.push(q);
                while m % q == 0 {
                    m /= q;
                }
            }
            q += 1;
        }
        if m > 1 {
            primes.push(m);
        }
        for q in primes {
            let g = at[self.n / q] ^ x;
            if !self.coprime_to_modulus(g) {
                return false;
            }
        }
        true
    }

    fn coprime_to_modulus(&self, g: Fe) -> bool {
        // polynomials as u256-like [u64; 4]
        let mut a = [0u64; 4];
        a[self.n >> 6] |= 1u64 << (self.n & 63);
        for &e in &self.low {
            a[(e >> 6) as usize] |= 1u64 << (e & 63);
        }
        let mut b = [g.0[0], g.0[1], g.0[2], 0];
        let deg = |p: &[u64; 4]| -> i32 {
            for i in (0..4).rev() {
                if p[i] != 0 {
                    return (64 * i + 63 - p[i].leading_zeros() as usize) as i32;
                }
            }
            -1
        };
        loop {
            let (da, db) = (deg(&a), deg(&b));
            if db < 0 {
                return da == 0;
            }
            if da < db {
                std::mem::swap(&mut a, &mut b);
                continue;
            }
            let s = (da - db) as usize;
            let (ws, bs) = (s >> 6, s & 63);
            for i in (0..4).rev() {
                if i >= ws {
                    let mut v = b[i - ws] << bs;
                    if bs > 0 && i > ws {
                        v |= b[i - ws - 1] >> (64 - bs);
                    }
                    a[i] ^= v;
                }
            }
        }
    }

    /// Fold the part above x^n onto the low part until nothing is left above: H = p >> n, p = (p mod x^n) + H r.
    /// Branch-free inside a pass; H fits three words because the product has degree at most 2n - 2.
    #[inline(always)]
    fn reduce(&self, q: [u64; 6]) -> Fe {
        // constant word offsets keep the product in registers
        match self.n >> 6 {
            0 => self.reduce_w::<0>(q),
            1 => self.reduce_w::<1>(q),
            _ => self.reduce_w::<2>(q),
        }
    }

    #[inline(always)]
    fn reduce_w<const NW: usize>(&self, q: [u64; 6]) -> Fe {
        let (nw, nb) = (NW, (self.n & 63) as u32);
        let mut p = [q[0], q[1], q[2], q[3], q[4], q[5], 0, 0];
        let keep = (1u64 << nb) - 1;
        loop {
            let h0 = (p[nw] >> nb) | ((p[nw + 1] << 1) << (63 - nb));
            let h1 = (p[nw + 1] >> nb) | ((p[nw + 2] << 1) << (63 - nb));
            let h2 = (p[nw + 2] >> nb) | ((p[nw + 3] << 1) << (63 - nb));
            if (h0 | h1 | h2) == 0 {
                break;
            }
            p[nw] &= keep;
            p[nw + 1] = 0;
            p[nw + 2] = 0;
            p[nw + 3] = 0;
            for t in 0..self.nterms {
                let e = self.terms[t];
                p[0] ^= h0 << e;
                p[1] ^= (h1 << e) | ((h0 >> 1) >> (63 - e));
                p[2] ^= (h2 << e) | ((h1 >> 1) >> (63 - e));
                p[3] ^= (h2 >> 1) >> (63 - e);
            }
        }
        Fe([p[0], p[1], p[2]])
    }

    /// a * x (one shift and one conditional fold).
    #[inline(always)]
    pub fn mulx(&self, a: Fe) -> Fe {
        let n = self.n;
        let carry = (a.0[(n - 1) >> 6] >> ((n - 1) & 63)) & 1;
        let mut t = [a.0[0] << 1, (a.0[1] << 1) | (a.0[0] >> 63), (a.0[2] << 1) | (a.0[1] >> 63)];
        if n < 192 {
            t[n >> 6] &= !(1u64 << (n & 63));
        }
        let mk = carry.wrapping_neg();
        t[0] ^= self.lowmask[0] & mk;
        t[1] ^= self.lowmask[1] & mk;
        t[2] ^= self.lowmask[2] & mk;
        Fe(t)
    }

    #[inline(always)]
    pub fn mul(&self, a: Fe, b: Fe) -> Fe {
        let (a, b) = (a.0, b.0);
        let mut p = [0u64; 6];
        if self.n <= 64 {
            let (l, h) = clmul(a[0], b[0]);
            p[0] = l;
            p[1] = h;
        } else if self.n <= 128 {
            let (l0, h0) = clmul(a[0], b[0]);
            let (l2, h2) = clmul(a[1], b[1]);
            let (lm, hm) = clmul(a[0] ^ a[1], b[0] ^ b[1]);
            p[0] = l0;
            p[1] = h0 ^ lm ^ l0 ^ l2;
            p[2] = l2 ^ hm ^ h0 ^ h2;
            p[3] = h2;
        } else {
            for i in 0..3 {
                for j in 0..3 {
                    let (l, h) = clmul(a[i], b[j]);
                    p[i + j] ^= l;
                    p[i + j + 1] ^= h;
                }
            }
        }
        self.reduce(p)
    }

    #[inline(always)]
    pub fn sq(&self, a: Fe) -> Fe {
        let mut p = [0u64; 6];
        for i in 0..3 {
            let (l, h) = clmul(a.0[i], a.0[i]);
            p[2 * i] = l;
            p[2 * i + 1] = h;
        }
        self.reduce(p)
    }

    /// x^(2^n - 2) by Itoh-Tsujii.
    pub fn inv(&self, x: Fe) -> Fe {
        let e = self.n - 1;
        let top = usize::BITS as usize - 1 - e.leading_zeros() as usize;
        let mut have = 1usize;
        let mut t = x;
        for i in (0..top).rev() {
            let mut s = t;
            for _ in 0..have {
                s = self.sq(s);
            }
            t = self.mul(s, t);
            have *= 2;
            if (e >> i) & 1 == 1 {
                t = self.mul(self.sq(t), x);
                have += 1;
            }
        }
        self.sq(t)
    }

    /// out[i] = 1/in[i]; zeros are skipped and stay zero.  One inversion and three multiplications each.
    pub fn batch_inv(&self, input: &[Fe], out: &mut [Fe]) {
        let mut acc = Fe::ONE;
        for i in 0..input.len() {
            out[i] = acc;
            if !input[i].is_zero() {
                acc = self.mul(acc, input[i]);
            }
        }
        acc = self.inv(acc);
        for i in (0..input.len()).rev() {
            if input[i].is_zero() {
                out[i] = Fe::ZERO;
                continue;
            }
            out[i] = self.mul(out[i], acc);
            acc = self.mul(acc, input[i]);
        }
    }

    /// Bit-by-bit product, for the self-test of `mul`.
    pub fn mul_slow(&self, a: Fe, b: Fe) -> Fe {
        let mut r = Fe::ZERO;
        let mut t = a;
        for i in 0..self.n {
            if b.bit(i) {
                r ^= t;
            }
            // t *= x
            let carry = t.bit(self.n - 1);
            t = Fe([t.0[0] << 1, (t.0[1] << 1) | (t.0[0] >> 63), (t.0[2] << 1) | (t.0[1] >> 63)]);
            if carry {
                t.0[self.n >> 6] &= !(1u64 << (self.n & 63));
                for &e in &self.low {
                    t ^= Fe::unit(e as usize);
                }
            }
        }
        r
    }
}
