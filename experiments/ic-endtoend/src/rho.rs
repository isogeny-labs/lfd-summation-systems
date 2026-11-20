//! Pollard rho: r-adding walk (r = 32), distinguished points, with or without
//! the negation map. Coefficients are recovered by replaying the colliding walks.
use crate::curve::*;
use crate::field::*;
use crate::hash::HMap;
use serde::Serialize;

const R: usize = 32;
const DP_BITS: u32 = 8;
const MAX_WALK: u64 = 20 << DP_BITS;

#[derive(Serialize, Clone, Default)]
pub struct RhoTrial {
    pub b: u32, pub n: u64, pub neg: bool, pub target: u64, pub ok: bool,
    pub walk_adds: u64, pub start_adds: u64, pub replay_adds: u64, pub setup_adds: u64, pub escape_adds: u64,
    pub hash_ops: u64, pub bookkeeping: u64, pub walks: u64, pub abandoned: u64, pub escapes: u64,
}
impl RhoTrial { pub fn group_adds(&self) -> u64 { self.walk_adds + self.start_adds + self.replay_adds + self.setup_adds + self.escape_adds } }

struct Ctx<'a> { c: &'a Curve, m: Vec<Pt>, mc: Vec<(u64, u64)>, neg: bool }

#[inline(always)]
fn h(x: u64) -> usize { (x.wrapping_mul(0x9E3779B97F4A7C15) >> 59) as usize }
#[inline(always)]
fn is_dp(x: u64) -> bool { (x.wrapping_mul(0xA24BAED4963EE407) >> (64 - DP_BITS)) == 0 }

impl<'a> Ctx<'a> {
    #[inline(always)]
    fn canon(&self, p: Pt) -> (Pt, bool) {
        if self.neg && !p.inf { let ny = self.c.f.neg(p.y); if ny < p.y { return (Pt { x: p.x, y: ny, inf: false }, true); } }
        (p, false)
    }
    /// one walk step; returns new point, index used, negation flag, adds spent
    #[inline(always)]
    fn step(&self, p: Pt) -> (Pt, usize, bool, u64) {
        let i = h(p.x);
        if !self.neg { return (self.c.add(p, self.m[i]), i, false, 1); }
        let mut adds = 0;
        for j in 0..R {
            let ii = (i + j) % R;
            let (q, s) = self.canon(self.c.add(p, self.m[ii]));
            adds += 1;
            if q.inf || h(q.x) != ii || j == R - 1 { return (q, ii, s, adds); }
        }
        unreachable!()
    }
    /// Run walk from start (coefficients a,b) for at most `len` steps (or until a DP).
    /// Returns (point, steps taken, a, b) with coefficients only when track.
    fn walk(&self, start: Pt, a0: u64, b0: u64, len: u64, track: bool, t: &mut RhoTrial, replay: bool) -> (Pt, u64, u64, u64, bool) {
        let n = self.c.n;
        let (mut p, s0) = self.canon(start);
        let (mut a, mut b) = if s0 { ((n - a0) % n, (n - b0) % n) } else { (a0, b0) };
        let mut saved = O; let mut saved_at = 0u64;
        let mut steps = 0u64;
        while steps < len {
            if p.inf { return (p, steps, a, b, false); }
            let (q, i, s, ad) = self.step(p);
            if replay { t.replay_adds += ad; } else { t.walk_adds += ad; t.bookkeeping += 2 * ad; }
            p = q; steps += 1;
            if track {
                a = addmod(a, self.mc[i].0, n); b = addmod(b, self.mc[i].1, n);
                if s { a = (n - a) % n; b = (n - b) % n; }
            }
            if self.neg {
                if steps % 64 == 0 { saved = p; saved_at = steps; }
                else if steps - saved_at <= 16 {
                    if !replay { t.bookkeeping += 1; }
                    if p == saved {
                        // fruitless cycle of length L: find the min-x point on it, double it
                        let l = steps - saved_at;
                        let (mut best, mut ba, mut bb) = (p, a, b);
                        let (mut q, mut qa, mut qb) = (p, a, b);
                        for _ in 0..l {
                            let (q2, i2, s2, ad2) = self.step(q);
                            if replay { t.replay_adds += ad2 } else { t.escape_adds += ad2 }
                            q = q2;
                            if track { qa = addmod(qa, self.mc[i2].0, n); qb = addmod(qb, self.mc[i2].1, n); if s2 { qa = (n - qa) % n; qb = (n - qb) % n; } }
                            if q.x < best.x { best = q; ba = qa; bb = qb; }
                        }
                        let (d, sd) = self.canon(self.c.dbl(best));
                        if replay { t.replay_adds += 1 } else { t.escape_adds += 1; t.escapes += 1; }
                        p = d; a = addmod(ba, ba, n); b = addmod(bb, bb, n);
                        if sd { a = (n - a) % n; b = (n - b) % n; }
                        saved = O; saved_at = steps;
                    }
                }
            }
            if is_dp(p.x) { return (p, steps, a, b, true); }
        }
        (p, steps, a, b, false)
    }
}

pub fn rho(c: &Curve, g: Pt, q: Pt, x_true: u64, neg: bool, seed: u64, b_bits: u32) -> RhoTrial {
    let n = c.n;
    let mut t = RhoTrial { b: b_bits, n, neg, target: x_true, ..Default::default() };
    let mut rng = Rng(seed);
    let mut lin = |rng: &mut Rng, t: &mut RhoTrial| -> (Pt, u64, u64) {
        let (a, b) = (1 + rng.below(n - 1), 1 + rng.below(n - 1));
        let (pa, o1) = c.mul(g, a); let (pb, o2) = c.mul(q, b);
        t.setup_adds += o1 + o2 + 1;
        (c.add(pa, pb), a, b)
    };
    let mut m = Vec::new(); let mut mc = Vec::new();
    for _ in 0..R { let (p, a, b) = lin(&mut rng, &mut t); m.push(p); mc.push((a, b)); }
    let ctx = Ctx { c, m, mc, neg };
    let (a0p, a0, b0) = lin(&mut rng, &mut t);
    let (dp_, dc, dd) = lin(&mut rng, &mut t);
    let mut table = HMap::with_capacity(1024);
    let mut start = a0p;
    let mut w: u64 = 0;
    let coefs = |w: u64| (addmod(a0, mulmod(w % n, dc, n), n), addmod(b0, mulmod(w % n, dd, n), n));
    loop {
        if w > 0 { start = c.add(start, dp_); t.start_adds += 1; }
        t.walks += 1;
        let (a, b) = coefs(w);
        let (p, len, _, _, dp) = ctx.walk(start, a, b, MAX_WALK, false, &mut t, false);
        if !dp { t.abandoned += 1; w += 1; continue; }
        t.hash_ops += 1;
        let key = p.x;
        match table.get(key) {
            None => { table.insert(key, (w << 16) | len); t.hash_ops += 1; }
            Some(v) => {
                let (w2, l2) = (v >> 16, v & 0xFFFF);
                // replay both walks with coefficient tracking
                let s1 = c.add(a0p, c.mul(dp_, w % n).0);
                let s2 = c.add(a0p, c.mul(dp_, w2 % n).0);
                t.replay_adds += 2 * (64 - n.leading_zeros() as u64) + 2; // two scalar mults to rebuild starts (upper bound)
                let (c1a, c1b) = coefs(w); let (c2a, c2b) = coefs(w2);
                let (p1, _, a1, b1, _) = ctx.walk(s1, c1a, c1b, len, true, &mut t, true);
                let (p2, _, a2, b2, _) = ctx.walk(s2, c2a, c2b, l2, true, &mut t, true);
                let (num, den) = if p1 == p2 { (submod(a1, a2, n), submod(b2, b1, n)) }
                    else { ((n - addmod(a1, a2, n)) % n, addmod(b1, b2, n)) };
                if den != 0 {
                    let x = mulmod(num, invmod(den, n), n);
                    t.ok = x == x_true && c.mul(g, x).0 == q;
                    return t;
                }
                table.insert(key, (w << 16) | len);
            }
        }
        w += 1;
    }
}
