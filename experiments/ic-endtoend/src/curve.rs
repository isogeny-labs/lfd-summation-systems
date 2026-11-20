//! Short Weierstrass curve over F_p, affine coordinates in Montgomery form.
use crate::field::*;
use crate::hash::HMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Pt { pub x: u64, pub y: u64, pub inf: bool }
pub const O: Pt = Pt { x: 0, y: 0, inf: true };

#[derive(Clone)]
pub struct Curve { pub f: Fp, pub a: u64, pub b: u64, pub n: u64, pub a_raw: u64, pub b_raw: u64 }

impl Curve {
    #[inline(always)]
    pub fn neg(&self, p: Pt) -> Pt { if p.inf { p } else { Pt { x: p.x, y: self.f.neg(p.y), inf: false } } }
    pub fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = &self.f;
        if p.inf { return q; }
        if q.inf { return p; }
        if p.x == q.x {
            if f.add(p.y, q.y) == 0 { return O; }
            return self.dbl(p);
        }
        let l = f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)));
        let x = f.sub(f.sub(f.sqr(l), p.x), q.x);
        let y = f.sub(f.mul(l, f.sub(p.x, x)), p.y);
        Pt { x, y, inf: false }
    }
    pub fn sub(&self, p: Pt, q: Pt) -> Pt { self.add(p, self.neg(q)) }
    pub fn dbl(&self, p: Pt) -> Pt {
        let f = &self.f;
        if p.inf || p.y == 0 { return O; }
        let x2 = f.sqr(p.x);
        let num = f.add(f.add(f.add(x2, x2), x2), self.a);
        let l = f.mul(num, f.inv(f.add(p.y, p.y)));
        let x = f.sub(f.sub(f.sqr(l), p.x), p.x);
        let y = f.sub(f.mul(l, f.sub(p.x, x)), p.y);
        Pt { x, y, inf: false }
    }
    /// scalar multiplication; returns (point, group ops used)
    pub fn mul(&self, p: Pt, k: u64) -> (Pt, u64) {
        let mut r = O;
        let mut ops = 0;
        let bits = 64 - k.leading_zeros();
        for i in (0..bits).rev() {
            if !r.inf { r = self.dbl(r); ops += 1; }
            if (k >> i) & 1 == 1 { r = if r.inf { p } else { ops += 1; self.add(r, p) }; }
        }
        (r, ops)
    }
    pub fn rhs(&self, x: u64) -> u64 {
        let f = &self.f;
        f.add(f.mul(f.add(f.sqr(x), self.a), x), self.b)
    }
    pub fn on_curve(&self, p: Pt) -> bool { p.inf || self.f.sqr(p.y) == self.rhs(p.x) }
    pub fn lift_x(&self, x: u64) -> Option<Pt> { self.f.sqrt(self.rhs(x)).map(|y| Pt { x, y, inf: false }) }
    pub fn random_point(&self, rng: &mut Rng) -> Pt {
        loop {
            let x = self.f.to_m(rng.below(self.f.p));
            if let Some(mut p) = self.lift_x(x) {
                if p.y == 0 { continue; }
                if rng.next() & 1 == 1 { p = self.neg(p); }
                return p;
            }
        }
    }
}

fn isqrt(n: u64) -> u64 { let mut r = (n as f64).sqrt() as u64; while r * r > n { r -= 1; } while (r + 1) * (r + 1) <= n { r += 1; } r }

/// Order of the group if it is prime (BSGS on a random point), else None.
pub fn prime_order(c: &Curve, rng: &mut Rng) -> Option<u64> {
    let p = c.f.p;
    let s = isqrt(p) + 1;
    let lo = p + 1 - 2 * s;
    let hi = p + 1 + 2 * s;
    let pt = c.random_point(rng);
    let w = hi - lo + 1;
    let m = isqrt(w) + 1;
    let mut baby = HMap::with_capacity(m as usize + 1);
    let mut j_p = O;
    let mut list = Vec::with_capacity(m as usize + 1);
    for j in 0..=m {
        if !j_p.inf { baby.insert(j_p.x, j); }
        list.push(j_p);
        j_p = c.add(j_p, pt);
    }
    let step = c.mul(pt, m).0;
    let mut g = c.mul(pt, lo).0;
    let mut base = lo;
    while base <= hi + m {
        let cand: Vec<u64> = if g.inf { vec![base] } else if let Some(j) = baby.get(g.x) {
            let jp = list[j as usize];
            if jp == g { vec![base.wrapping_sub(j)] } else { vec![base + j] }
        } else { vec![] };
        for n in cand {
            if n >= lo && n <= hi && n > 0 {
                if is_prime(n) && c.mul(pt, n).0.inf { return Some(n); }
                return None;
            }
        }
        g = c.add(g, step);
        base += m;
    }
    None
}

pub fn find_prime(b: u32) -> (u64, u32) {
    let kk = (b - 4).min(24);
    let step = 1u64 << kk;
    let mut c = ((1u64 << b) + step - 1) / step;
    loop {
        let p = c * step + 1;
        if p >= (1u64 << b) && is_prime(p) { return (p, kk); }
        c += 1;
    }
}

/// i-th prime-order curve for bit size b (deterministic).
pub fn curves_for(b: u32, count: usize) -> Vec<Curve> {
    let (p, _) = find_prime(b);
    let f = Fp::new(p);
    let mut rng = Rng(20261120u64 ^ b as u64);
    let mut out = Vec::new();
    while out.len() < count {
        let ar = rng.below(p);
        let br = rng.below(p);
        let a = f.to_m(ar);
        let bb = f.to_m(br);
        // nonsingular: 4a^3 + 27b^2 != 0
        let d = f.add(f.mul(f.to_m(4), f.mul(f.sqr(a), a)), f.mul(f.to_m(27), f.sqr(bb)));
        if d == 0 { continue; }
        let mut c = Curve { f, a, b: bb, n: 0, a_raw: ar, b_raw: br };
        if let Some(n) = prime_order(&c, &mut rng) { c.n = n; out.push(c); }
    }
    out
}

/// Elements of t*H_k in Montgomery form (H_k = subgroup of order 2^k).
pub fn coset(f: &Fp, k: u32, t: u64) -> Vec<u64> {
    let mut w = 2u64;
    while f.is_square(f.to_m(w)) { w += 1; }
    let z = f.pow(f.to_m(w), (f.p - 1) >> k); // order exactly 2^k
    let mut v = Vec::with_capacity(1 << k);
    let mut h = t;
    for _ in 0..(1u64 << k) { v.push(h); h = f.mul(h, z); }
    v
}

/// A factor base: class representatives (one per +- pair) and the x-set.
pub struct Base { pub reps: Vec<Pt>, pub xset: HMap, pub label: String }

impl Base {
    pub fn from_xs(c: &Curve, xs: &[u64], label: &str) -> Base {
        let mut reps = Vec::new();
        let mut xset = HMap::with_capacity(xs.len());
        for &x in xs {
            if let Some(p) = c.lift_x(x) {
                if p.y != 0 && xset.insert(x, reps.len() as u64) { reps.push(p); }
            }
        }
        Base { reps, xset, label: label.to_string() }
    }
    pub fn pkm(c: &Curve, k: u32, t: u64) -> Base { Base::from_xs(c, &coset(&c.f, k, t), if t == c.f.one { "subgroup" } else { "coset" }) }
    /// negative control: f random classes
    pub fn random(c: &Curve, f: usize, rng: &mut Rng) -> Base {
        let mut reps = Vec::new();
        let mut xset = HMap::with_capacity(f);
        while reps.len() < f {
            let p = c.random_point(rng);
            if xset.insert(p.x, reps.len() as u64) { reps.push(p); }
        }
        Base { reps, xset, label: "random".into() }
    }
    /// positive control: {+- jG : 1 <= j <= f}
    pub fn ap(c: &Curve, f: usize, g: Pt) -> Base {
        let mut reps = Vec::new();
        let mut xset = HMap::with_capacity(f);
        let mut p = g;
        for _ in 0..f { if xset.insert(p.x, reps.len() as u64) { reps.push(p); } p = c.add(p, g); }
        Base { reps, xset, label: "ap".into() }
    }
    pub fn size(&self) -> usize { 2 * self.reps.len() }
}

/// For a point R and all reps Q_i: x(R - Q_i) and x(R + Q_i) by one batched
/// inversion. Writes into out_minus/out_plus; entries where x_R == x_Q are u64::MAX.
pub fn batch_diff_sum(c: &Curve, r: Pt, reps: &[Pt], scratch: &mut Vec<u64>, om: &mut Vec<u64>, op: &mut Vec<u64>) {
    let f = &c.f;
    let n = reps.len();
    scratch.clear(); om.clear(); op.clear();
    let mut acc = f.one;
    for q in reps {
        let d = f.sub(r.x, q.x);
        scratch.push(acc);
        if d != 0 { acc = f.mul(acc, d); }
    }
    let mut inv = f.inv(acc);
    om.resize(n, 0); op.resize(n, 0);
    for i in (0..n).rev() {
        let q = reps[i];
        let d = f.sub(r.x, q.x);
        if d == 0 { om[i] = u64::MAX; op[i] = u64::MAX; continue; }
        let di = f.mul(inv, scratch[i]);
        inv = f.mul(inv, d);
        let sx = f.add(r.x, q.x);
        let lm = f.mul(f.add(r.y, q.y), di); // R - Q
        let lp = f.mul(f.sub(r.y, q.y), di); // R + Q
        om[i] = f.sub(f.sqr(lm), sx);
        op[i] = f.sub(f.sqr(lp), sx);
    }
}
