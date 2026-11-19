//! Monomials and column orders.
//!
//! A monomial is a square-free mask over the 2k chain variables: bit (i-1) k + j is
//! u_{i,j} = X_i^{2^j} (hky.basis/v1, mask_encoding 1).  The exponent of X_i is the k-bit block i.
//! The shift u_{i,j} multiplies by X_i^{2^j}: it adds 2^j to block i modulo 2^k (carry inside the
//! block), so the weight goes up by exactly 1 when bit j of block i is clear and does not go up
//! otherwise.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
pub struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
    fn write_u64(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}
pub type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;

pub fn binom(n: u64, r: u64) -> u64 {
    if r > n {
        return 0;
    }
    let mut c = 1u64;
    for i in 0..r {
        c = c * (n - i) / (i + 1);
    }
    c
}

#[derive(Clone, Copy, Debug)]
pub struct Mon {
    pub k: u32,
}

impl Mon {
    /// shift index s in 0..2k: block i = s / k, bit j = s % k
    #[inline]
    pub fn shift(&self, mask: u64, s: u32) -> u64 {
        let k = self.k;
        let i = s / k;
        let j = s % k;
        let nb = (1u64 << k) - 1;
        let blk = (mask >> (i * k)) & nb;
        let nbk = (blk + (1u64 << j)) & nb;
        (mask & !(nb << (i * k))) | (nbk << (i * k))
    }
    #[inline]
    pub fn sigma(&self, mask: u64) -> u64 {
        let k = self.k;
        let nb = (1u64 << k) - 1;
        (mask >> k) | ((mask & nb) << k)
    }
    pub fn nshifts(&self) -> u32 {
        2 * self.k
    }
    /// all masks of weight w over 2k bits, ascending
    pub fn of_weight(&self, w: u32) -> Vec<u64> {
        let n = 2 * self.k;
        let mut out = Vec::new();
        if w > n {
            return out;
        }
        if w == 0 {
            return vec![0];
        }
        let mut x: u64 = (1u64 << w) - 1;
        let lim = 1u64 << n;
        while x < lim {
            out.push(x);
            // Gosper's hack
            let c = x & x.wrapping_neg();
            let r = x + c;
            x = (((r ^ x) >> 2) / c) | r;
        }
        out
    }
}

/// Columns of M_<=D in a degree-non-increasing order.
pub struct Layout {
    pub k: u32,
    pub d: u32,
    pub mon: Mon,
    pub masks: Vec<u64>,
    pub wt: Vec<u8>,
    pub idx: FxMap<u64, u32>,
}

impl Layout {
    /// hky.basis/v1 column_order 0: degree descending, then mask descending.
    pub fn order0(k: u32, d: u32) -> Layout {
        let mon = Mon { k };
        let mut masks = Vec::new();
        for w in (0..=d).rev() {
            let mut v = mon.of_weight(w);
            v.reverse();
            masks.extend(v);
        }
        Layout::from_masks(k, d, masks).expect("order0")
    }

    /// hky.basis/v1 column_order 1: sigma-orbit grouped (representative r = max(s, sigma s) in
    /// order-0 order, then sigma r).
    pub fn order1(k: u32, d: u32) -> Layout {
        let mon = Mon { k };
        let o0 = Layout::order0(k, d);
        let mut masks = Vec::with_capacity(o0.masks.len());
        for &m in &o0.masks {
            let s = mon.sigma(m);
            if m >= s {
                masks.push(m);
                if s != m {
                    masks.push(s);
                }
            }
        }
        Layout::from_masks(k, d, masks).expect("order1")
    }

    pub fn from_masks(k: u32, d: u32, masks: Vec<u64>) -> Result<Layout, String> {
        let n = 2 * k as u64;
        if n > 64 {
            return Err("2k > 64".into());
        }
        let want: u64 = (0..=d as u64).map(|w| binom(n, w)).sum();
        if masks.len() as u64 != want {
            return Err(format!("{} columns, expected {}", masks.len(), want));
        }
        let mut idx = FxMap::default();
        idx.reserve(masks.len());
        let mut wt = Vec::with_capacity(masks.len());
        let mut prev = u32::MAX;
        for (c, &m) in masks.iter().enumerate() {
            if n < 64 && m >> n != 0 {
                return Err(format!("mask {m:#x} has bits above 2k"));
            }
            let w = m.count_ones();
            if w > d {
                return Err(format!("mask {m:#x} has weight {w} > D"));
            }
            if w > prev {
                return Err(format!("degree increases at column {c}"));
            }
            prev = w;
            wt.push(w as u8);
            if idx.insert(m, c as u32).is_some() {
                return Err(format!("duplicate mask {m:#x}"));
            }
        }
        if masks.last() != Some(&0) {
            return Err("last column is not the monomial 1".into());
        }
        Ok(Layout { k, d, mon: Mon { k }, masks, wt, idx })
    }

    pub fn ncols(&self) -> usize {
        self.masks.len()
    }
    #[inline]
    pub fn col(&self, mask: u64) -> Option<u32> {
        self.idx.get(&mask).copied()
    }
    /// is the order exactly hky.basis/v1 order `o`?
    pub fn is_order(&self, o: u32) -> bool {
        let r = match o {
            0 => Layout::order0(self.k, self.d),
            1 => Layout::order1(self.k, self.d),
            _ => return false,
        };
        r.masks == self.masks
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shift_and_orders() {
        let m = Mon { k: 3 };
        // block 1 = 0b011 + 1 = 0b100
        assert_eq!(m.shift(0b011, 0), 0b100);
        // block 2 bit 1: 0b010<<3 + 0b010<<3 -> 0b100<<3
        assert_eq!(m.shift(0b010 << 3, 4), 0b100 << 3);
        // wrap: 0b111 + 1 = 0
        assert_eq!(m.shift(0b111, 0), 0);
        assert_eq!(m.of_weight(2).len(), 15);
        let l0 = Layout::order0(3, 3);
        let l1 = Layout::order1(3, 3);
        assert_eq!(l0.ncols(), 1 + 6 + 15 + 20);
        assert_eq!(l1.ncols(), l0.ncols());
        assert!(l1.is_order(1) && l0.is_order(0));
    }
}
