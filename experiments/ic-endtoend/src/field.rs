//! Prime field F_p (p < 2^62) in Montgomery form, and Z/N helpers.

#[derive(Clone, Copy, Debug)]
pub struct Fp {
    pub p: u64,
    pinv: u64, // -p^{-1} mod 2^64
    r2: u64,   // 2^128 mod p
    pub one: u64,
}

impl Fp {
    pub fn new(p: u64) -> Self {
        assert!(p % 2 == 1 && p < (1u64 << 62));
        let mut inv: u64 = 1;
        for _ in 0..7 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let pinv = inv.wrapping_neg();
        let r = ((1u128 << 64) % p as u128) as u64;
        let r2 = ((r as u128 * r as u128) % p as u128) as u64;
        Fp { p, pinv, r2, one: r }
    }
    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let t = a as u128 * b as u128;
        let m = (t as u64).wrapping_mul(self.pinv);
        let u = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if u >= self.p { u - self.p } else { u }
    }
    #[inline(always)]
    pub fn sqr(&self, a: u64) -> u64 { self.mul(a, a) }
    #[inline(always)]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    #[inline(always)]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    #[inline(always)]
    pub fn neg(&self, a: u64) -> u64 { if a == 0 { 0 } else { self.p - a } }
    pub fn to_m(&self, a: u64) -> u64 { self.mul(a % self.p, self.r2) }
    pub fn from_m(&self, a: u64) -> u64 { self.mul(a, 1) }
    pub fn pow(&self, a: u64, mut e: u64) -> u64 {
        let mut r = self.one;
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 { r = self.mul(r, b); }
            b = self.sqr(b);
            e >>= 1;
        }
        r
    }
    pub fn inv(&self, a: u64) -> u64 { self.pow(a, self.p - 2) }
    pub fn is_square(&self, a: u64) -> bool {
        a == 0 || self.pow(a, (self.p - 1) / 2) == self.one
    }
    /// Tonelli-Shanks; a in Montgomery form, returns a root in Montgomery form.
    pub fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 { return Some(0); }
        if !self.is_square(a) { return None; }
        let mut q = self.p - 1;
        let mut s = 0;
        while q % 2 == 0 { q /= 2; s += 1; }
        let mut zv = 2u64;
        while self.is_square(self.to_m(zv)) { zv += 1; }
        let mut c = self.pow(self.to_m(zv), q);
        let mut x = self.pow(a, (q + 1) / 2);
        let mut t = self.pow(a, q);
        let mut m = s;
        while t != self.one {
            let mut i = 0;
            let mut tt = t;
            while tt != self.one { tt = self.sqr(tt); i += 1; }
            let mut b = c;
            for _ in 0..(m - i - 1) { b = self.sqr(b); }
            x = self.mul(x, b);
            c = self.sqr(b);
            t = self.mul(t, c);
            m = i;
        }
        Some(x)
    }
}

pub fn mulmod(a: u64, b: u64, n: u64) -> u64 { ((a as u128 * b as u128) % n as u128) as u64 }
pub fn addmod(a: u64, b: u64, n: u64) -> u64 { let s = a as u128 + b as u128; (s % n as u128) as u64 }
pub fn submod(a: u64, b: u64, n: u64) -> u64 { addmod(a, n - b % n, n) }
pub fn powmod(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1u64 % n;
    a %= n;
    while e > 0 {
        if e & 1 == 1 { r = mulmod(r, a, n); }
        a = mulmod(a, a, n);
        e >>= 1;
    }
    r
}
pub fn invmod(a: u64, n: u64) -> u64 { powmod(a, n - 2, n) } // n prime

pub fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 { return n == q; }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 { d /= 2; s += 1; }
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 { continue 'w; }
        }
        return false;
    }
    true
}

/// splitmix64
#[derive(Clone)]
pub struct Rng(pub u64);
impl Rng {
    #[inline(always)]
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: u64) -> u64 { ((self.next() as u128 * n as u128) >> 64) as u64 }
}
pub fn derive(seed: u64, tag: u64) -> u64 {
    let mut r = Rng(seed ^ tag.wrapping_mul(0xD1B54A32D192ED03));
    r.next(); r.next()
}
