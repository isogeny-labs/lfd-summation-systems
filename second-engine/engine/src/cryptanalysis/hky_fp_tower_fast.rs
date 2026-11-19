//! # A faster instrument for the Huang–Kosters–Yeo closure of `hky_fp_tower`.
//!
//! [`hky_fp_tower`](crate::cryptanalysis::hky_fp_tower) is the reference: it
//! computes `V_{F,c}` by the definition and is left untouched. This module
//! computes the same space with the same work-queue algorithm and reports the
//! same [`Closure`] and [`FallDegrees`], through cumulative levels that each
//! change one thing (`research/hky_engine_speed_20261117/README.md` has the
//! measurement of each step and the equality gate they all pass):
//!
//! - level 0: a port of the reference (heap of touched columns, `% p` on
//!   every update, SipHash column index, products materialised when a row
//!   enters). Its `rows_reduced`, `entries` and `muladds` equal the
//!   reference's exactly, which is what shows the port is the same algorithm.
//! - level 1: delayed reduction. The accumulator holds values below `p²`
//!   (`p < 2³¹`, so `acc + f·x < 2⁶³`) and one conditional subtraction of `p²`
//!   replaces the division; a value is reduced modulo `p` once, when its
//!   column is visited.
//! - level 2: no heap and no list of touched columns: the accumulator is
//!   scanned upward between the first and the last column touched (a pivot
//!   row only touches columns after its own), and a row that meets no pivot
//!   added since it was reduced is not reduced twice.
//! - level 3: products `x·g` are formed when they are taken from the queue,
//!   in the parallel section, from the stored row (the queue holds
//!   `(row, variable)` pairs, not polynomials), and the column index uses a
//!   multiplicative hash.
//! - level 4: a stored row with more than a third of its columns non-zero is
//!   kept as a dense `u32` slice from its leading column and subtracted with
//!   a vectorised loop (AVX2 or AVX-512 by run-time detection, the same loop
//!   compiled for the baseline instruction set otherwise).
//!
//! Independently of the level, [`FastOptions::symmetry`] names an involution
//! of the tower variables (for two summands: the swap of the two towers).
//! The engine **proves** it is a symmetry of the system before using it
//! ([`check_symmetry`]: the tower rules are permuted to themselves and every
//! generator is mapped to itself or to its negative, checked coefficient by
//! coefficient) and otherwise runs unsplit and says so in
//! [`FastStats::symmetry`]. With a symmetry `σ` and `p` odd,
//! `V_{F,c} = V⁺ ⊕ V⁻` (the `±1` eigenspaces; `V_{F,c}` is `σ`-stable
//! because it is the least closed space containing a `σ`-stable set). The
//! two parts are computed in orbit coordinates: a column of part `s` is an
//! orbit `{m, σm}` with basis vector `m + s·σm` (`m` alone for a fixed
//! monomial, which has no column in the `−` part). For `g ∈ V^s` and a
//! variable `x` with `σx ≠ x` only `x` of the pair `{x, σx}` is applied and
//! `x·g` is projected to both parts (`2π_q = 1 + q·σ`): `(σx)·g =
//! s·σ(x·g)`, so the two projections span the same space as `x·g` and
//! `(σx)·g`. For a fixed variable `x·g` stays in `V^s`. `dim (V ∩ R_{≤d})` is
//! the sum over the parts (the degree is `σ`-invariant), and `1` is a column
//! of the `+` part. Each part has about half the columns, so an elimination
//! step costs about a quarter.
//!
//! The answer does not depend on the level, the batch size, the order of the
//! queue or the split: `V_{F,c}` is one space. What does depend on them is
//! the work and the censoring by `max_entries` and `budget`, which are limits
//! of an instrument, not of the definition (`max_cols` is checked on the
//! monomial count exactly as in the reference).

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicU64, Ordering as AtOrd};
use std::time::Instant;

use rayon::prelude::*;

use crate::cryptanalysis::hky_fp_tower::{
    count_monomials, Closure, FallDegrees, HkyLimits, HkyRing, NfPoly, MAX_FREE, MAX_TOWER,
};

const FREE_SHIFT: u32 = 56;
const DEG_SHIFT: u32 = 120;
pub(super) const NONE: u32 = u32::MAX;

pub(super) type Key = u128;

fn key(mask: u64, free: &[u8; MAX_FREE]) -> Key {
    let deg = mask.count_ones() + free.iter().map(|&e| u32::from(e)).sum::<u32>();
    let mut k = (u128::from(deg) << DEG_SHIFT) | u128::from(mask);
    for (i, &e) in free.iter().enumerate() {
        k |= u128::from(e) << (FREE_SHIFT + 8 * i as u32);
    }
    k
}

fn key_mask(k: Key) -> u64 {
    (k as u64) & ((1u64 << MAX_TOWER) - 1)
}

fn key_free(k: Key) -> [u8; MAX_FREE] {
    let mut f = [0u8; MAX_FREE];
    for (i, e) in f.iter_mut().enumerate() {
        *e = (k >> (FREE_SHIFT + 8 * i as u32)) as u8;
    }
    f
}

fn key_degree(k: Key) -> u32 {
    (k >> DEG_SHIFT) as u32
}

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    a * b % p
}

fn inv(a: u64, p: u64) -> u64 {
    let (mut r, mut b, mut e) = (1u64, a % p, p - 2);
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, b, p);
        }
        b = mulm(b, b, p);
        e >>= 1;
    }
    r
}

/// [`mul_tower_var`] for the unit tests of `hky_fp_tower_fast2`.
#[cfg(test)]
pub(super) fn mul_tower_var_pub(
    ring: &HkyRing,
    m: u64,
    j: usize,
    coef: u64,
    out: &mut Vec<(u64, u64)>,
) {
    mul_tower_var(ring, m, j, coef, out)
}

/// `y_j · (mask, coef)` in normal form: the reference's `mul_tower_var`.
fn mul_tower_var(ring: &HkyRing, m: u64, j: usize, coef: u64, out: &mut Vec<(u64, u64)>) {
    let p = ring.p;
    let bit = 1u64 << j;
    if m & bit == 0 {
        out.push((m | bit, coef));
        return;
    }
    let r = ring.rules[j];
    let m0 = m & !bit;
    if r.d != 0 {
        out.push((m0, mulm(coef, r.d, p)));
    }
    if r.c != 0 {
        out.push((m, mulm(coef, r.c, p)));
    }
    if let Some(nx) = r.next {
        if r.b != 0 {
            mul_tower_var(ring, m0, nx as usize, mulm(coef, r.b, p), out);
        }
        if r.a != 0 {
            mul_tower_var(ring, m, nx as usize, mulm(coef, r.a, p), out);
        }
    }
}

/// A multiplicative hasher for the 128-bit monomial keys (level 3 on).
#[derive(Default, Clone, Copy)]
pub struct KeyHasher(u64);

impl Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ u64::from(b)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }
    fn write_u128(&mut self, k: u128) {
        let x = (k as u64) ^ ((k >> 64) as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
        let h = x.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        self.0 = h ^ (h >> 29);
    }
}

enum Index {
    Sip(HashMap<Key, u32>),
    Fast(HashMap<Key, u32, BuildHasherDefault<KeyHasher>>),
}

impl Index {
    #[inline]
    fn get(&self, k: &Key) -> u32 {
        match self {
            Index::Sip(m) => m[k],
            Index::Fast(m) => m[k],
        }
    }
}

/// The order the queue of pending rows is served in. The space does not
/// depend on it; the work does, by a factor of several (README, "Schedule").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueueOrder {
    /// The reference's: last in, first out.
    Lifo,
    /// Products of the lowest-degree rows first (last in, first out among
    /// equals). Same space; only the work to reach `1` changes.
    LowDegreeFirst,
    /// Products of the sparsest stored rows first (by the bit length of the
    /// row's number of terms; last in, first out among equals).
    SparseFirst,
}

/// What [`closure_fast`] does differently from the reference.
#[derive(Clone, Debug)]
pub struct FastOptions {
    /// Cumulative optimisation level, 0 (port of the reference) to 4.
    pub level: u8,
    /// A candidate involution of the tower variables (`perm[j]` is the image
    /// of variable `j`), or `None` for no split. Verified before use.
    pub symmetry: Option<Vec<usize>>,
    /// Rows per batch. 0 picks it: the limits' value (the reference's 4096)
    /// under [`QueueOrder::Lifo`], `max(32, 2 · threads)` otherwise.
    pub batch: usize,
    pub order: QueueOrder,
    /// A row is stored dense when `nnz · dense_den > span · dense_num`
    /// (level 4), `span` being the columns from its leading one to the last.
    pub dense_num: usize,
    pub dense_den: usize,
    /// Test hook for the negative controls of the equality gate; 0 in every
    /// measurement. 1: the products by variable 0 are not formed.
    /// 2: the delayed reduction subtracts `p² − 1` instead of `p²`.
    pub fault: u8,
    /// Levels 2 on. A row is reduced only until its leading column is not a
    /// pivot; the reference also clears every later pivot column. Same
    /// leading monomials, same space. **Measured slower** (the rows are
    /// denser: README, "Tried and rejected"); kept off, for the record.
    pub lead_only: bool,
    /// [`QueueOrder::SparseFirst`]: doublings of the number of terms per
    /// priority bucket.
    pub sparse_width: usize,
}

impl Default for FastOptions {
    fn default() -> Self {
        FastOptions {
            level: 4,
            symmetry: None,
            batch: 0,
            order: QueueOrder::SparseFirst,
            dense_num: 1,
            dense_den: 3,
            fault: 0,
            lead_only: false,
            sparse_width: 1,
        }
    }
}

/// Phase timers and structure counts of one closure. Wall-clock seconds of
/// the calling thread except `multiply_cpu_s`, which is summed over threads.
#[derive(Clone, Debug, Default)]
pub struct FastStats {
    /// `off`, `split` or `rejected: <why>`.
    pub symmetry: String,
    /// Columns of each part (one part when unsplit).
    pub part_cols: Vec<u64>,
    /// Pivots of each part when the closure stopped.
    pub part_rank: Vec<u64>,
    pub setup_s: f64,
    /// The parallel section: products (level 3 on) and reduction against
    /// the echelon as it stood when the batch started.
    pub parallel_s: f64,
    /// The sequential section: reduction against the pivots the batch has
    /// added, products (levels 0 to 2) and insertion.
    pub sequential_s: f64,
    pub multiply_cpu_s: f64,
    pub batches: u64,
    pub dense_rows: u64,
    pub sparse_rows: u64,
    /// Stored coefficients (a dense row counts its whole slice).
    pub stored_words: u64,
    pub threads: usize,
    pub batch: usize,
    pub simd: &'static str,
}

/// A stored echelon row, monic, leading column first.
pub(super) enum RowS {
    Sparse(Vec<(u32, u32)>),
    /// `v[i]` is the coefficient of column `lead + i`; `v[0] = 1`.
    Dense {
        lead: u32,
        v: Vec<u32>,
        nnz: u32,
    },
}

impl RowS {
    pub(super) fn nnz(&self) -> usize {
        match self {
            RowS::Sparse(r) => r.len(),
            RowS::Dense { nnz, .. } => *nnz as usize,
        }
    }
    pub(super) fn words(&self) -> usize {
        match self {
            RowS::Sparse(r) => r.len(),
            RowS::Dense { v, .. } => v.len(),
        }
    }
    pub(super) fn for_each(&self, mut f: impl FnMut(u32, u32)) {
        match self {
            RowS::Sparse(r) => {
                for &(c, x) in r {
                    f(c, x);
                }
            }
            RowS::Dense { lead, v, .. } => {
                for (i, &x) in v.iter().enumerate() {
                    if x != 0 {
                        f(lead + i as u32, x);
                    }
                }
            }
        }
    }
}

pub(super) type Row = Vec<(u32, u32)>;

/// One part of the space: its columns (orbit representatives when split),
/// highest degree first.
pub(super) struct Part {
    pub(super) keys: Vec<Key>,
    pub(super) deg: Vec<u8>,
}

struct Echelon {
    pivot: Vec<u32>,
    rows: Vec<RowS>,
    /// Pivots added since the current batch started.
    fresh: u32,
}

/// Per-thread scratch.
struct Scratch {
    acc: Vec<u64>,
    touched: Vec<u32>,
}

impl Scratch {
    fn new(n: usize) -> Self {
        Scratch {
            acc: vec![0u64; n],
            touched: Vec::new(),
        }
    }
}

// ── dense kernels ──────────────────────────────────────────────────────

/// `acc[i] ← acc[i] + f·row[i] mod p²`-lazily: inputs below `p²`, `f < p`,
/// `row[i] < p`, `p < 2³¹`, so the sum is below `2⁶³` and one conditional
/// subtraction of `p²` brings it back below `p²`.
#[inline(always)]
fn axpy_body(acc: &mut [u64], row: &[u32], f: u64, p2: u64) {
    for (a, &x) in acc.iter_mut().zip(row.iter()) {
        let t = (*a + f * u64::from(x)) as i64 - p2 as i64;
        *a = (t + ((t >> 63) & p2 as i64)) as u64;
    }
}

#[inline(never)]
fn axpy_portable(acc: &mut [u64], row: &[u32], f: u64, p2: u64) {
    axpy_body(acc, row, f, p2)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn axpy_avx2(acc: &mut [u64], row: &[u32], f: u64, p2: u64) {
    axpy_body(acc, row, f, p2)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512dq,avx512vl")]
#[inline(never)]
unsafe fn axpy_avx512(acc: &mut [u64], row: &[u32], f: u64, p2: u64) {
    axpy_body(acc, row, f, p2)
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Simd {
    Portable,
    #[cfg_attr(
        not(target_arch = "x86_64"),
        allow(
            dead_code,
            reason = "Preserves the explicit x86 backend vocabulary on other architectures."
        )
    )]
    Avx2,
    #[cfg_attr(
        not(target_arch = "x86_64"),
        allow(
            dead_code,
            reason = "Preserves the explicit x86 backend vocabulary on other architectures."
        )
    )]
    Avx512,
}

pub(super) fn detect_simd() -> Simd {
    #[cfg(target_arch = "x86_64")]
    {
        if std::env::var_os("HKY2_SIMD").is_some_and(|v| v == "portable") {
            return Simd::Portable;
        }
        let force_avx2 = std::env::var_os("HKY2_SIMD").is_some_and(|v| v == "avx2");
        if !force_avx2
            && is_x86_feature_detected!("avx512f")
            && is_x86_feature_detected!("avx512dq")
            && is_x86_feature_detected!("avx512vl")
        {
            return Simd::Avx512;
        }
        if is_x86_feature_detected!("avx2") {
            return Simd::Avx2;
        }
    }
    Simd::Portable
}

#[inline]
pub(super) fn axpy_dense(simd: Simd, acc: &mut [u64], row: &[u32], f: u64, p2: u64) {
    match simd {
        Simd::Portable => axpy_portable(acc, row, f, p2),
        #[cfg(target_arch = "x86_64")]
        // SAFETY: `detect_simd` returned the variant only after the CPU
        // reported the features the function is compiled for.
        Simd::Avx2 => unsafe { axpy_avx2(acc, row, f, p2) },
        #[cfg(target_arch = "x86_64")]
        // SAFETY: as above.
        Simd::Avx512 => unsafe { axpy_avx512(acc, row, f, p2) },
        #[cfg(not(target_arch = "x86_64"))]
        _ => axpy_portable(acc, row, f, p2),
    }
}

/// Reduction of a value below `p²` modulo `p` without a division: `m` is
/// `⌊2⁶⁴ / p⌋`, the quotient estimate is at most 2 too small.
#[inline(always)]
pub(super) fn barrett(v: u64, p: u64, m: u64) -> u64 {
    let q = ((u128::from(v) * u128::from(m)) >> 64) as u64;
    let mut r = v - q * p;
    if r >= p {
        r -= p;
    }
    if r >= p {
        r -= p;
    }
    r
}

struct Ctx {
    p: u64,
    p2: u64,
    barrett_m: u64,
    level: u8,
    simd: Simd,
    dense_num: usize,
    dense_den: usize,
    lead_only: bool,
}

impl Echelon {
    /// The reference's reduction (level 0), or the same with delayed
    /// reduction (level 1). Returns the monic row and the multiply-adds.
    fn reduce_heap(&self, cx: &Ctx, r: &[(u32, u32)], s: &mut Scratch) -> (Row, u64) {
        let lazy = cx.level >= 1;
        let (p, p2) = (cx.p, cx.p2);
        let acc = &mut s.acc;
        let touched = &mut s.touched;
        let mut ops = 0u64;
        let mut heap: BinaryHeap<Reverse<u32>> = BinaryHeap::with_capacity(r.len() * 2);
        for &(col, c) in r {
            if acc[col as usize] == 0 {
                touched.push(col);
                heap.push(Reverse(col));
            }
            acc[col as usize] = (acc[col as usize] + u64::from(c)) % p;
        }
        let mut out: Row = Vec::new();
        while let Some(Reverse(col)) = heap.pop() {
            let v = if lazy {
                acc[col as usize] % p
            } else {
                acc[col as usize]
            };
            if v == 0 {
                acc[col as usize] = 0;
                continue;
            }
            let piv = self.pivot[col as usize];
            if piv == NONE {
                out.push((col, v as u32));
                acc[col as usize] = 0;
                continue;
            }
            let f = p - v;
            let RowS::Sparse(row) = &self.rows[piv as usize] else {
                unreachable!("dense rows exist from level 4 on");
            };
            ops += row.len() as u64;
            if lazy {
                acc[col as usize] = 0;
                for &(c2, x) in &row[1..] {
                    let i = c2 as usize;
                    if acc[i] == 0 {
                        touched.push(c2);
                        heap.push(Reverse(c2));
                    }
                    let t = acc[i] + f * u64::from(x);
                    acc[i] = t.min(t.wrapping_sub(p2));
                }
            } else {
                for &(c2, x) in row {
                    let i = c2 as usize;
                    if acc[i] == 0 && c2 != col {
                        touched.push(c2);
                        heap.push(Reverse(c2));
                    }
                    acc[i] = (acc[i] + f * u64::from(x)) % p;
                }
            }
        }
        for &t in touched.iter() {
            acc[t as usize] = 0;
        }
        touched.clear();
        ops += normalise(&mut out, p);
        (out, ops)
    }

    /// Levels 2 on: no heap and no list of touched columns. The accumulator
    /// is scanned upward from the row's first column to the last column
    /// touched so far (a pivot row only touches columns after its own).
    fn reduce_scan(&self, cx: &Ctx, r: &[(u32, u32)], s: &mut Scratch) -> (Row, u64) {
        let (p, p2, bm) = (cx.p, cx.p2, cx.barrett_m);
        let n = self.pivot.len();
        let acc = &mut s.acc[..n];
        let mut ops = 0u64;
        let (Some(&(first, _)), Some(&(last, _))) = (r.first(), r.last()) else {
            return (Vec::new(), 0);
        };
        for &(col, c) in r {
            acc[col as usize] = u64::from(c);
        }
        let mut out: Row = Vec::new();
        let mut last = last as usize;
        let mut ci = first as usize;
        while ci <= last {
            let raw = acc[ci];
            if raw == 0 {
                ci += 1;
                continue;
            }
            acc[ci] = 0;
            let v = barrett(raw, p, bm);
            if v == 0 {
                ci += 1;
                continue;
            }
            let piv = self.pivot[ci];
            if piv == NONE {
                out.push((ci as u32, v as u32));
                if cx.lead_only {
                    // Lead-only: the row has its leading column; the rest is
                    // gathered as it stands, pivot columns included.
                    for c2 in ci + 1..=last {
                        let raw = acc[c2];
                        if raw != 0 {
                            acc[c2] = 0;
                            let x = barrett(raw, p, bm);
                            if x != 0 {
                                out.push((c2 as u32, x as u32));
                            }
                        }
                    }
                    break;
                }
                ci += 1;
                continue;
            }
            let f = p - v;
            match &self.rows[piv as usize] {
                RowS::Sparse(row) => {
                    ops += row.len() as u64;
                    for &(c2, x) in &row[1..] {
                        let i = c2 as usize;
                        let t = acc[i] + f * u64::from(x);
                        acc[i] = t.min(t.wrapping_sub(p2));
                    }
                    last = last.max(row[row.len() - 1].0 as usize);
                }
                RowS::Dense { v: dv, nnz, .. } => {
                    ops += u64::from(*nnz);
                    axpy_dense(cx.simd, &mut acc[ci + 1..], &dv[1..], f, p2);
                    last = n - 1;
                }
            }
            ci += 1;
        }
        ops += normalise(&mut out, p);
        (out, ops)
    }

    fn reduce(&self, cx: &Ctx, r: &[(u32, u32)], s: &mut Scratch) -> (Row, u64) {
        if cx.level >= 2 {
            self.reduce_scan(cx, r, s)
        } else {
            self.reduce_heap(cx, r, s)
        }
    }

    /// Whether a row reduced before this batch meets a pivot added since.
    fn meets_fresh_pivot(&self, r: &[(u32, u32)], lead_only: bool) -> bool {
        if lead_only {
            return r
                .first()
                .is_some_and(|&(c, _)| self.pivot[c as usize] != NONE);
        }
        self.fresh > 0 && r.iter().any(|&(c, _)| self.pivot[c as usize] != NONE)
    }
}

pub(super) fn normalise(out: &mut Row, p: u64) -> u64 {
    if let Some(&(_, lead)) = out.first() {
        let s = inv(u64::from(lead), p);
        for e in out.iter_mut() {
            e.1 = mulm(u64::from(e.1), s, p) as u32;
        }
        out.len() as u64
    } else {
        0
    }
}

// ── symmetry ───────────────────────────────────────────────────────────

fn sigma_mask(mask: u64, perm: &[usize]) -> u64 {
    let mut out = 0u64;
    let mut bits = mask;
    while bits != 0 {
        let j = bits.trailing_zeros() as usize;
        bits &= bits - 1;
        out |= 1u64 << perm[j];
    }
    out
}

fn sigma_key(k: Key, perm: &[usize]) -> Key {
    let m = key_mask(k);
    (k & !u128::from((1u64 << MAX_TOWER) - 1)) | u128::from(sigma_mask(m, perm))
}

/// Proves that `perm` (an involution of the tower variables; free variables
/// are fixed) is a symmetry of the system: it maps every tower rule to the
/// rule of the image variable, and every generator to itself or to its
/// negative. Returns the sign of each generator, or why the proof fails.
/// Nothing is sampled: every rule and every coefficient is compared.
pub fn check_symmetry(ring: &HkyRing, perm: &[usize], gens: &[NfPoly]) -> Result<Vec<i8>, String> {
    let nt = ring.n_tower();
    let p = ring.p;
    if p.is_multiple_of(2) {
        return Err("p = 2".into());
    }
    if perm.len() != nt {
        return Err(format!(
            "permutation of {} variables, ring has {nt}",
            perm.len()
        ));
    }
    for j in 0..nt {
        if perm[j] >= nt || perm[perm[j]] != j {
            return Err(format!("not an involution at variable {j}"));
        }
    }
    if (0..nt).all(|j| perm[j] == j) {
        return Err("identity".into());
    }
    for j in 0..nt {
        let (r, s) = (ring.rules[j], ring.rules[perm[j]]);
        let next_ok = r.next.map(|x| perm[x as usize]) == s.next.map(|x| x as usize);
        if !(next_ok && r.a == s.a && r.b == s.b && r.c == s.c && r.d == s.d) {
            return Err(format!(
                "tower rule {j} is not the rule of its image {}",
                perm[j]
            ));
        }
    }
    let mut signs = Vec::with_capacity(gens.len());
    for (gi, g) in gens.iter().enumerate() {
        let map: HashMap<(u64, [u8; MAX_FREE]), u64> =
            g.iter().map(|&(m, f, x)| ((m, f), x % p)).collect();
        let mut sign = 0i8;
        for &(m, f, x) in g {
            let x = x % p;
            if x == 0 {
                continue;
            }
            let y = map.get(&(sigma_mask(m, perm), f)).copied().unwrap_or(0);
            let s = if y == x {
                1
            } else if y == p - x {
                -1
            } else {
                return Err(format!("generator {gi} is not mapped to ±itself"));
            };
            if sign == 0 {
                sign = s;
            } else if sign != s {
                return Err(format!("generator {gi} is not mapped to ±itself"));
            }
        }
        signs.push(if sign == 0 { 1 } else { sign });
    }
    Ok(signs)
}

// ── the space: columns, products ───────────────────────────────────────

pub(super) struct Space<'a> {
    ring: &'a HkyRing,
    pub(super) parts: Vec<Part>,
    /// Column of a monomial (unsplit), or of an orbit representative in the
    /// `+` part (split).
    index: Index,
    /// Split only: column in the `−` part of each `+` column, or `NONE`.
    minus_col: Vec<u32>,
    perm: Option<Vec<usize>>,
    nt: usize,
}

fn all_keys(ring: &HkyRing, c: u32) -> Vec<Key> {
    let nt = ring.n_tower();
    let nf = ring.n_free;
    let mut keys = Vec::new();
    let mut masks: Vec<Vec<u64>> = vec![Vec::new(); c as usize + 1];
    fn rec_mask(start: usize, nt: usize, left: u32, cur: u64, out: &mut Vec<Vec<u64>>, depth: u32) {
        out[depth as usize].push(cur);
        if left == 0 {
            return;
        }
        for j in start..nt {
            rec_mask(j + 1, nt, left - 1, cur | 1 << j, out, depth + 1);
        }
    }
    rec_mask(0, nt, c, 0, &mut masks, 0);
    let mut frees: Vec<Vec<[u8; MAX_FREE]>> = vec![Vec::new(); c as usize + 1];
    fn rec_free(
        i: usize,
        nf: usize,
        left: u32,
        cur: &mut [u8; MAX_FREE],
        deg: u32,
        out: &mut Vec<Vec<[u8; MAX_FREE]>>,
    ) {
        if i == nf {
            out[deg as usize].push(*cur);
            return;
        }
        for e in 0..=left {
            cur[i] = e as u8;
            rec_free(i + 1, nf, left - e, cur, deg + e, out);
        }
        cur[i] = 0;
    }
    rec_free(0, nf, c, &mut [0u8; MAX_FREE], 0, &mut frees);
    for (k, ms) in masks.iter().enumerate() {
        for (d, fs) in frees.iter().enumerate() {
            if k + d > c as usize {
                continue;
            }
            for &m in ms {
                for f in fs {
                    keys.push(key(m, f));
                }
            }
        }
    }
    keys.sort_unstable_by_key(|&k| Reverse(k));
    keys
}

impl<'a> Space<'a> {
    pub(super) fn build(
        ring: &'a HkyRing,
        c: u32,
        perm: Option<Vec<usize>>,
        fast_index: bool,
    ) -> Self {
        let all = all_keys(ring, c);
        let nt = ring.n_tower();
        let mk_index = |keys: &[Key]| -> Index {
            if fast_index {
                Index::Fast(
                    keys.iter()
                        .enumerate()
                        .map(|(i, &k)| (k, i as u32))
                        .collect(),
                )
            } else {
                Index::Sip(
                    keys.iter()
                        .enumerate()
                        .map(|(i, &k)| (k, i as u32))
                        .collect(),
                )
            }
        };
        let part_of = |keys: Vec<Key>| -> Part {
            let deg = keys.iter().map(|&k| key_degree(k) as u8).collect();
            Part { keys, deg }
        };
        match perm {
            None => {
                let index = mk_index(&all);
                Space {
                    ring,
                    parts: vec![part_of(all)],
                    index,
                    minus_col: Vec::new(),
                    perm: None,
                    nt,
                }
            }
            Some(perm) => {
                // representatives: the larger key of each orbit
                let plus: Vec<Key> = all
                    .iter()
                    .copied()
                    .filter(|&k| k >= sigma_key(k, &perm))
                    .collect();
                let mut minus = Vec::new();
                let mut minus_col = Vec::with_capacity(plus.len());
                for &k in &plus {
                    if k == sigma_key(k, &perm) {
                        minus_col.push(NONE);
                    } else {
                        minus_col.push(minus.len() as u32);
                        minus.push(k);
                    }
                }
                let index = mk_index(&plus);
                Space {
                    ring,
                    parts: vec![part_of(plus), part_of(minus)],
                    index,
                    minus_col,
                    perm: Some(perm),
                    nt,
                }
            }
        }
    }

    pub(super) fn split(&self) -> bool {
        self.perm.is_some()
    }

    /// The variables whose products are formed (one of each swapped pair).
    pub(super) fn applies(&self, v: usize) -> bool {
        match &self.perm {
            Some(perm) if v < self.nt => perm[v] >= v,
            _ => true,
        }
    }

    /// `x_v · row` for a row of part `s`, as rows of the parts it lands in.
    /// Returns the number of terms formed before merging (the unit the
    /// reference counts) through `terms_out`.
    pub(super) fn multiply(&self, s: usize, row: &RowS, v: usize, out: &mut Vec<(u8, Row)>) -> u64 {
        let ring = self.ring;
        let p = ring.p;
        let nt = self.nt;
        let part = &self.parts[s];
        let mut terms: Vec<(Key, u64)> = Vec::with_capacity(row.nnz() * 2);
        let mut tmp = Vec::new();
        let mut expand = |k: Key, x: u64, terms: &mut Vec<(Key, u64)>| {
            let (m, mut f) = (key_mask(k), key_free(k));
            if v >= nt {
                f[v - nt] += 1;
                terms.push((key(m, &f), x));
            } else {
                tmp.clear();
                mul_tower_var(ring, m, v, x, &mut tmp);
                for &(mm, cc) in &tmp {
                    terms.push((key(mm, &f), cc));
                }
            }
        };
        match &self.perm {
            None => {
                row.for_each(|col, x| expand(part.keys[col as usize], u64::from(x), &mut terms));
                terms.sort_unstable_by_key(|t| Reverse(t.0));
                let mut r: Row = Vec::with_capacity(terms.len());
                let mut i = 0;
                while i < terms.len() {
                    let k = terms[i].0;
                    let mut sum = 0;
                    while i < terms.len() && terms[i].0 == k {
                        sum = (sum + terms[i].1) % p;
                        i += 1;
                    }
                    if sum != 0 {
                        r.push((self.index.get(&k), sum as u32));
                    }
                }
                let n = r.len() as u64;
                out.push((0, r));
                n
            }
            Some(perm) => {
                // lift: column o of part s is m + sgn·σm (m alone if fixed)
                row.for_each(|col, x| {
                    let k = part.keys[col as usize];
                    let x = u64::from(x);
                    expand(k, x, &mut terms);
                    let sk = sigma_key(k, perm);
                    if sk != k {
                        expand(sk, if s == 0 { x } else { p - x }, &mut terms);
                    }
                });
                let fixed_var = v >= nt || perm[v] == v;
                let targets: Vec<usize> = if fixed_var { vec![s] } else { vec![0, 1] };
                let mut total = 0u64;
                for q in targets {
                    // 2π_q: m' = rep → b_rep; m' = σ rep → q·b_rep; fixed → 2·b (q = +) or 0
                    let mut cols: Vec<(u32, u64)> = Vec::with_capacity(terms.len());
                    for &(k, cc) in &terms {
                        let sk = sigma_key(k, perm);
                        let (rep, coef) = if sk == k {
                            if q == 1 {
                                continue;
                            }
                            (k, (2 * cc) % p)
                        } else if k > sk {
                            (k, cc)
                        } else {
                            (sk, if q == 0 { cc } else { (p - cc) % p })
                        };
                        let cp = self.index.get(&rep);
                        let col = if q == 0 {
                            cp
                        } else {
                            self.minus_col[cp as usize]
                        };
                        cols.push((col, coef));
                    }
                    cols.sort_unstable_by_key(|t| t.0);
                    let mut r: Row = Vec::with_capacity(cols.len());
                    let mut i = 0;
                    while i < cols.len() {
                        let c0 = cols[i].0;
                        let mut sum = 0;
                        while i < cols.len() && cols[i].0 == c0 {
                            sum = (sum + cols[i].1) % p;
                            i += 1;
                        }
                        if sum != 0 {
                            r.push((c0, sum as u32));
                        }
                    }
                    total += r.len() as u64;
                    out.push((q as u8, r));
                }
                total
            }
        }
    }

    /// A generator as a row of its part.
    pub(super) fn gen_row(&self, g: &NfPoly, sign: i8) -> (u8, Row) {
        match &self.perm {
            None => {
                let mut r: Row = g
                    .iter()
                    .map(|(m, e, x)| (self.index.get(&key(*m, e)), *x as u32))
                    .collect();
                r.sort_unstable_by_key(|t| t.0);
                (0, r)
            }
            Some(perm) => {
                let q = if sign >= 0 { 0u8 } else { 1u8 };
                let mut r: Row = Vec::new();
                for (m, e, x) in g {
                    let k = key(*m, e);
                    let sk = sigma_key(k, perm);
                    if k < sk || *x == 0 {
                        continue; // the partner term carries the orbit
                    }
                    if k == sk && q == 1 {
                        continue; // zero by antisymmetry (checked)
                    }
                    let cp = self.index.get(&k);
                    let col = if q == 0 {
                        cp
                    } else {
                        self.minus_col[cp as usize]
                    };
                    r.push((col, *x as u32));
                }
                r.sort_unstable_by_key(|t| t.0);
                (q, r)
            }
        }
    }
}

/// The bucket of a row with `nnz` terms in [`QueueOrder::SparseFirst`]:
/// `1 + ⌊log₂ nnz / width⌋` (`width` doublings per bucket).
pub(super) fn sparse_bucket(nnz: usize, width: usize) -> usize {
    let lg = (usize::BITS - 1 - nnz.max(1).leading_zeros()) as usize;
    1 + lg / width.max(1)
}

enum Item {
    Row(u8, Row),
    /// `x_var ·` (row `row` of part `part`), formed when dequeued.
    Prod {
        part: u8,
        row: u32,
        var: u16,
    },
}

/// `V_{F,c}`, as [`closure`](crate::cryptanalysis::hky_fp_tower::closure)
/// computes it, with the instrument of `opts`.
pub fn closure_fast(
    ring: &HkyRing,
    gens: &[NfPoly],
    c: u32,
    limits: &HkyLimits,
    opts: &FastOptions,
) -> (Closure, FastStats) {
    let started = Instant::now();
    let p = ring.p;
    assert!(
        p < 1 << 31 || opts.level == 0,
        "hky_fp_tower_fast: levels 1+ need p < 2^31"
    );
    let n_cols = count_monomials(ring.n_tower(), ring.n_free, c);
    let mut stats = FastStats {
        symmetry: "off".into(),
        threads: rayon::current_num_threads(),
        ..FastStats::default()
    };
    let censored = |why: &str, rows: u64, entries: u64, muladds: u64, stats: FastStats| {
        (
            Closure {
                c,
                cols: n_cols,
                dims: Vec::new(),
                contains_one: false,
                censored: Some(why.to_string()),
                rows_reduced: rows,
                entries,
                muladds,
                ms: started.elapsed().as_secs_f64() * 1e3,
            },
            stats,
        )
    };
    if n_cols > limits.max_cols {
        return censored("max_cols", 0, 0, 0, stats);
    }
    // The symmetry, proved or dropped.
    let live: Vec<&NfPoly> = gens.iter().filter(|g| !g.is_empty()).collect();
    let live_owned: Vec<NfPoly> = live.iter().map(|g| (*g).clone()).collect();
    let (perm, signs) = match &opts.symmetry {
        None => (None, vec![1i8; live.len()]),
        Some(perm) => match check_symmetry(ring, perm, &live_owned) {
            Ok(signs) => {
                stats.symmetry = "split".into();
                (Some(perm.clone()), signs)
            }
            Err(why) => {
                stats.symmetry = format!("rejected: {why}");
                (None, vec![1i8; live.len()])
            }
        },
    };
    let space = Space::build(ring, c, perm, opts.level >= 3);
    let total_cols: usize = space.parts.iter().map(|q| q.keys.len()).sum();
    if !space.split() {
        assert_eq!(
            total_cols as u128, n_cols,
            "hky_fp_tower_fast: column count"
        );
    }
    stats.part_cols = space.parts.iter().map(|q| q.keys.len() as u64).collect();
    let simd = if opts.level >= 4 {
        detect_simd()
    } else {
        Simd::Portable
    };
    stats.simd = match simd {
        Simd::Portable => "portable",
        Simd::Avx2 => "avx2",
        Simd::Avx512 => "avx512",
    };
    let cx = Ctx {
        p,
        p2: if opts.fault == 2 { p * p - 1 } else { p * p },
        barrett_m: (u128::from(u64::MAX) / u128::from(p)) as u64,
        level: opts.level,
        simd,
        dense_num: opts.dense_num,
        dense_den: opts.dense_den,
        lead_only: opts.lead_only && opts.level >= 2,
    };
    let nv = ring.n_vars();
    let applies = |v: usize| space.applies(v) && !(opts.fault == 1 && v == 0);
    let max_n = space.parts.iter().map(|q| q.keys.len()).max().unwrap_or(0);
    // 1 is the last column of part 0.
    let one_col = (space.parts[0].keys.len() - 1) as u32;
    let mut ech: Vec<Echelon> = space
        .parts
        .iter()
        .map(|q| Echelon {
            pivot: vec![NONE; q.keys.len()],
            rows: Vec::new(),
            fresh: 0,
        })
        .collect();
    let lazy_queue = opts.level >= 3;
    // Queue buckets: one (LIFO), or one per degree of the source row.
    let n_buckets = match opts.order {
        QueueOrder::Lifo => 1,
        QueueOrder::LowDegreeFirst => c as usize + 2,
        QueueOrder::SparseFirst => 66,
    };
    let mut queue: Vec<Vec<Item>> = (0..n_buckets).map(|_| Vec::new()).collect();
    for (g, &sign) in live.iter().zip(&signs) {
        if HkyRing::degree(g) <= c {
            let (q, r) = space.gen_row(g, sign);
            if !r.is_empty() {
                queue[0].push(Item::Row(q, r));
            }
        }
    }
    let batch_size = match (opts.batch, opts.order) {
        (0, QueueOrder::Lifo) => limits.batch,
        (0, _) => (2 * rayon::current_num_threads()).max(32),
        (b, _) => b,
    };
    stats.batch = batch_size;
    let mut rows_reduced = 0u64;
    let mut muladds = 0u64;
    let mut entries = 0u64;
    let mut contains_one = false;
    let mul_ns = AtomicU64::new(0);
    stats.setup_s = started.elapsed().as_secs_f64();
    let mut seq = Scratch::new(max_n);
    'outer: while let Some(b) = queue.iter().position(|q| !q.is_empty()) {
        if started.elapsed() > limits.budget {
            stats.multiply_cpu_s = mul_ns.load(AtOrd::Relaxed) as f64 * 1e-9;
            return censored("budget", rows_reduced, entries, muladds, stats);
        }
        let take = queue[b].len().min(batch_size);
        let at = queue[b].len() - take;
        let batch: Vec<Item> = queue[b].drain(at..).collect();
        stats.batches += 1;
        let t0 = Instant::now();
        // Parallel: form the products (level 3 on) and reduce against the
        // echelon as it stands.
        let reduced: Vec<Vec<(u8, Row, u64)>> = batch
            .par_iter()
            .map_init(
                || Scratch::new(max_n),
                |s, item| {
                    let mut cands: Vec<(u8, Row)> = Vec::new();
                    let mut ops = 0u64;
                    match item {
                        Item::Row(q, r) => cands.push((*q, r.clone())),
                        Item::Prod { part, row, var } => {
                            let t = Instant::now();
                            let src = &ech[*part as usize].rows[*row as usize];
                            ops += space.multiply(*part as usize, src, *var as usize, &mut cands);
                            mul_ns.fetch_add(t.elapsed().as_nanos() as u64, AtOrd::Relaxed);
                        }
                    }
                    let mut out = Vec::with_capacity(cands.len());
                    for (i, (q, r)) in cands.into_iter().enumerate() {
                        let (red, o) = ech[q as usize].reduce(&cx, &r, s);
                        out.push((q, red, o + if i == 0 { ops } else { 0 }));
                    }
                    out
                },
            )
            .collect();
        stats.parallel_s += t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        rows_reduced += take as u64;
        for e in ech.iter_mut() {
            e.fresh = 0;
        }
        for group in reduced {
            for (q, r, ops) in group {
                muladds += ops;
                if r.is_empty() {
                    continue;
                }
                let qi = q as usize;
                // Against the pivots this batch has added.
                let r = if cx.level >= 2 && !ech[qi].meets_fresh_pivot(&r, cx.lead_only) {
                    r
                } else {
                    let (r, ops) = ech[qi].reduce(&cx, &r, &mut seq);
                    muladds += ops;
                    r
                };
                let Some(&(lead, _)) = r.first() else {
                    continue;
                };
                let d = space.parts[qi].deg[lead as usize] as u32;
                let row_idx = ech[qi].rows.len() as u32;
                let n_q = ech[qi].pivot.len();
                let span = n_q - lead as usize;
                let stored = if cx.level >= 4 && r.len() * cx.dense_den > span * cx.dense_num {
                    let mut v = vec![0u32; span];
                    for &(col, x) in &r {
                        v[(col - lead) as usize] = x;
                    }
                    stats.dense_rows += 1;
                    RowS::Dense {
                        lead,
                        v,
                        nnz: r.len() as u32,
                    }
                } else {
                    stats.sparse_rows += 1;
                    RowS::Sparse(r)
                };
                if d < c {
                    let bucket = match opts.order {
                        QueueOrder::Lifo => 0,
                        QueueOrder::LowDegreeFirst => d as usize + 1,
                        QueueOrder::SparseFirst => sparse_bucket(stored.nnz(), opts.sparse_width),
                    };
                    if lazy_queue {
                        for v in 0..nv {
                            if applies(v) {
                                queue[bucket].push(Item::Prod {
                                    part: q,
                                    row: row_idx,
                                    var: v as u16,
                                });
                            }
                        }
                    } else {
                        let t = Instant::now();
                        let mut prods = Vec::new();
                        for v in 0..nv {
                            if applies(v) {
                                prods.clear();
                                muladds += space.multiply(qi, &stored, v, &mut prods);
                                for (tq, pr) in prods.drain(..) {
                                    queue[bucket].push(Item::Row(tq, pr));
                                }
                            }
                        }
                        mul_ns.fetch_add(t.elapsed().as_nanos() as u64, AtOrd::Relaxed);
                    }
                }
                entries += stored.nnz() as u64;
                stats.stored_words += stored.words() as u64;
                ech[qi].pivot[lead as usize] = row_idx;
                ech[qi].rows.push(stored);
                ech[qi].fresh += 1;
                if qi == 0 && lead == one_col {
                    contains_one = true;
                    stats.sequential_s += t1.elapsed().as_secs_f64();
                    break 'outer;
                }
                if entries > limits.max_entries {
                    stats.multiply_cpu_s = mul_ns.load(AtOrd::Relaxed) as f64 * 1e-9;
                    return censored("max_entries", rows_reduced, entries, muladds, stats);
                }
            }
        }
        stats.sequential_s += t1.elapsed().as_secs_f64();
    }
    stats.multiply_cpu_s = mul_ns.load(AtOrd::Relaxed) as f64 * 1e-9;
    stats.part_rank = ech.iter().map(|e| e.rows.len() as u64).collect();
    // dims[d] = dim (V ∩ R_{≤d}); with 1 ∈ V, V = R_{≤c}: count monomials.
    let mut dims = vec![0u64; c as usize + 1];
    if contains_one {
        for (qi, part) in space.parts.iter().enumerate() {
            for (col, &dg) in part.deg.iter().enumerate() {
                // an orbit of two monomials has a column in both parts
                let _ = (qi, col);
                dims[dg as usize] += 1;
            }
        }
    } else {
        for (e, part) in ech.iter().zip(&space.parts) {
            for (col, &pv) in e.pivot.iter().enumerate() {
                if pv != NONE {
                    dims[part.deg[col] as usize] += 1;
                }
            }
        }
    }
    for d in 1..dims.len() {
        dims[d] += dims[d - 1];
    }
    if contains_one {
        assert_eq!(
            u128::from(dims[c as usize]),
            n_cols,
            "hky_fp_tower_fast: full space"
        );
    }
    (
        Closure {
            c,
            cols: n_cols,
            dims,
            contains_one,
            censored: None,
            rows_reduced,
            entries,
            muladds,
            ms: started.elapsed().as_secs_f64() * 1e3,
        },
        stats,
    )
}

/// [`fall_degrees_with`](crate::cryptanalysis::hky_fp_tower::fall_degrees_with)
/// on [`closure_fast`]: the same reading of the closures, word for word.
pub fn fall_degrees_fast_with<F: FnMut(&Closure, &FastStats)>(
    ring: &HkyRing,
    gens: &[NfPoly],
    c_min: u32,
    c_max: u32,
    limits: &HkyLimits,
    opts: &FastOptions,
    mut on_closure: F,
) -> FallDegrees {
    let mut closures: Vec<Closure> = Vec::new();
    let mut falls = Vec::new();
    let mut last_fall = None;
    for c in c_min..=c_max {
        let (cl, st) = closure_fast(ring, gens, c, limits, opts);
        let stop = cl.contains_one || cl.censored.is_some();
        if let Some(prev) = closures.last() {
            if cl.censored.is_none()
                && prev.censored.is_none()
                && cl.dims[c as usize - 1] > prev.dim()
            {
                falls.push(c);
            }
        } else if cl.censored.is_none() && c > 0 && cl.dims[c as usize - 1] > 0 && c == c_min {
            let min_deg = gens
                .iter()
                .filter(|g| !g.is_empty())
                .map(HkyRing::degree)
                .min()
                .unwrap_or(0);
            if c_min - 1 < min_deg {
                falls.push(c);
            }
        }
        if cl.contains_one {
            last_fall = Some(c);
        }
        on_closure(&cl, &st);
        closures.push(cl);
        if stop {
            break;
        }
    }
    let first_fall = falls.first().copied();
    FallDegrees {
        closures,
        last_fall,
        first_fall,
        falls,
    }
}

pub fn fall_degrees_fast(
    ring: &HkyRing,
    gens: &[NfPoly],
    c_min: u32,
    c_max: u32,
    limits: &HkyLimits,
    opts: &FastOptions,
) -> FallDegrees {
    fall_degrees_fast_with(ring, gens, c_min, c_max, limits, opts, |_, _| {})
}

/// The swap of the first two of `blocks` equal towers of `t` variables: the
/// candidate symmetry of an `m`-summand system (`blocks = m`).
pub fn swap_first_two_blocks(t: usize, blocks: usize) -> Option<Vec<usize>> {
    if blocks < 2 || t == 0 {
        return None;
    }
    let mut perm: Vec<usize> = (0..t * blocks).collect();
    for j in 0..t {
        perm[j] = t + j;
        perm[t + j] = j;
    }
    Some(perm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::f4_fp_tower::SquareRule;
    use crate::cryptanalysis::hky_fp_tower::fall_degrees;

    const P: u64 = 2013265921;

    fn kummer_rules(t: usize, blocks: usize, top: u64) -> Vec<SquareRule> {
        let mut rules = Vec::new();
        for b in 0..blocks {
            for j in 0..t {
                let last = j + 1 == t;
                rules.push(SquareRule {
                    next: (!last).then_some((b * t + j + 1) as u8),
                    a: 0,
                    b: if last { 0 } else { 1 },
                    c: 0,
                    d: if last { top } else { 0 },
                });
            }
        }
        rules
    }

    /// A general rule `y_j² = a·y_j·y_{j+1} + b·y_{j+1} + c·y_j + d`.
    fn general_rules(t: usize, blocks: usize) -> Vec<SquareRule> {
        let mut rules = Vec::new();
        for b in 0..blocks {
            for j in 0..t {
                let last = j + 1 == t;
                rules.push(SquareRule {
                    next: (!last).then_some((b * t + j + 1) as u8),
                    a: if last { 0 } else { 3 + j as u64 },
                    b: if last { 0 } else { 5 + 2 * j as u64 },
                    c: 7 + j as u64,
                    d: 11 + 3 * j as u64,
                });
            }
        }
        rules
    }

    /// A symmetric quadric in the two top variables and their towers.
    fn symmetric_gen(ring: &HkyRing, t: usize, k: u64) -> NfPoly {
        let n = ring.n_vars();
        let e = |pairs: &[(usize, u32)]| {
            let mut v = vec![0u32; n];
            for &(i, x) in pairs {
                v[i] = x;
            }
            v
        };
        ring.from_raw(&[
            (e(&[(0, 2), (t, 2)]), 1),
            (e(&[(0, 1), (t, 1)]), 17 + k),
            (e(&[(0, 2)]), 5),
            (e(&[(t, 2)]), 5),
            (e(&[(0, 1)]), 9 + k),
            (e(&[(t, 1)]), 9 + k),
            (e(&[]), 23 + 7 * k),
        ])
    }

    fn same(a: &FallDegrees, b: &FallDegrees, what: &str) {
        assert_eq!(a.last_fall, b.last_fall, "{what}: last fall");
        assert_eq!(a.first_fall, b.first_fall, "{what}: first fall");
        assert_eq!(a.falls, b.falls, "{what}: falls");
        assert_eq!(a.closures.len(), b.closures.len(), "{what}: closures");
        for (x, y) in a.closures.iter().zip(&b.closures) {
            assert_eq!(x.dims, y.dims, "{what}: dims at c = {}", x.c);
            assert_eq!(x.contains_one, y.contains_one, "{what}: one at c = {}", x.c);
            assert_eq!(x.censored, y.censored, "{what}: censored at c = {}", x.c);
            assert_eq!(x.cols, y.cols, "{what}: cols at c = {}", x.c);
        }
    }

    fn all_variants(ring: &HkyRing, gens: &[NfPoly], perm: Option<Vec<usize>>, cmax: u32) {
        let limits = HkyLimits::default();
        let reference = fall_degrees(ring, gens, 1, cmax, &limits);
        for level in 0..=4u8 {
            for sym in [false, true] {
                for (batch, order) in [
                    (0, QueueOrder::Lifo),
                    (7, QueueOrder::Lifo),
                    (0, QueueOrder::LowDegreeFirst),
                    (5, QueueOrder::SparseFirst),
                ] {
                    let opts = FastOptions {
                        level,
                        symmetry: if sym { perm.clone() } else { None },
                        batch,
                        order,
                        // dense rows early, so small cases exercise them
                        dense_num: 1,
                        dense_den: 8,
                        fault: 0,
                        lead_only: false,
                        sparse_width: 1,
                    };
                    let fast = fall_degrees_fast(ring, gens, 1, cmax, &limits, &opts);
                    same(
                        &reference,
                        &fast,
                        &format!("level {level} sym {sym} batch {batch} {order:?}"),
                    );
                    if level >= 2 {
                        let opts = FastOptions {
                            lead_only: true,
                            ..opts
                        };
                        let fast = fall_degrees_fast(ring, gens, 1, cmax, &limits, &opts);
                        same(
                            &reference,
                            &fast,
                            &format!("lead-only level {level} sym {sym} batch {batch}"),
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn level_zero_counts_what_the_reference_counts() {
        let ring = HkyRing::new(P, kummer_rules(4, 2, 3), 0);
        let gens = vec![symmetric_gen(&ring, 4, 0)];
        let limits = HkyLimits::default();
        let reference = fall_degrees(&ring, &gens, 1, 5, &limits);
        let opts = FastOptions {
            level: 0,
            order: QueueOrder::Lifo,
            ..FastOptions::default()
        };
        let fast = fall_degrees_fast(&ring, &gens, 1, 5, &limits, &opts);
        same(&reference, &fast, "level 0");
        for (x, y) in reference.closures.iter().zip(&fast.closures) {
            assert_eq!(x.muladds, y.muladds, "muladds at c = {}", x.c);
            assert_eq!(x.rows_reduced, y.rows_reduced);
            assert_eq!(x.entries, y.entries);
        }
    }

    #[test]
    fn every_variant_equals_the_reference_on_kummer_towers() {
        for t in 2..=4 {
            let ring = HkyRing::new(P, kummer_rules(t, 2, 3), 0);
            let gens = vec![symmetric_gen(&ring, t, t as u64)];
            all_variants(&ring, &gens, swap_first_two_blocks(t, 2), 5);
        }
    }

    #[test]
    fn every_variant_equals_the_reference_on_general_rules() {
        for t in 2..=3 {
            let ring = HkyRing::new(P, general_rules(t, 2), 0);
            let gens = vec![symmetric_gen(&ring, t, 1)];
            all_variants(&ring, &gens, swap_first_two_blocks(t, 2), 5);
        }
    }

    #[test]
    fn free_variables_and_three_blocks() {
        // three towers and one free variable; the generators are symmetric
        // under the swap of the first two towers, and one is antisymmetric.
        let t = 2;
        let ring = HkyRing::new(P, kummer_rules(t, 3, 5), 1);
        let n = ring.n_vars();
        let e = |pairs: &[(usize, u32)]| {
            let mut v = vec![0u32; n];
            for &(i, x) in pairs {
                v[i] = x;
            }
            v
        };
        let u = 3 * t;
        let g1 = ring.from_raw(&[
            (e(&[(0, 1), (t, 1)]), 1),
            (e(&[(0, 1)]), 4),
            (e(&[(t, 1)]), 4),
            (e(&[(u, 1)]), P - 1),
            (e(&[]), 6),
        ]);
        let g2 = ring.from_raw(&[
            (e(&[(u, 1), (2 * t, 1)]), 1),
            (e(&[(u, 2)]), 3),
            (e(&[(2 * t, 1)]), 8),
            (e(&[]), 2),
        ]);
        let anti = ring.from_raw(&[(e(&[(0, 1), (u, 1)]), 1), (e(&[(t, 1), (u, 1)]), P - 1)]);
        let perm = swap_first_two_blocks(t, 3);
        assert_eq!(
            check_symmetry(
                &ring,
                perm.as_ref().unwrap(),
                &[g1.clone(), g2.clone(), anti.clone()]
            ),
            Ok(vec![1, 1, -1])
        );
        all_variants(&ring, &[g1.clone(), g2.clone()], perm.clone(), 4);
        all_variants(&ring, &[g1, g2, anti], perm, 4);
    }

    #[test]
    fn a_false_symmetry_is_rejected_and_the_answer_is_unchanged() {
        let t = 3;
        let ring = HkyRing::new(P, kummer_rules(t, 2, 3), 0);
        let mut g = symmetric_gen(&ring, t, 2);
        // break the symmetry in one coefficient
        let i = g
            .iter()
            .position(|&(m, _, _)| m == 1)
            .expect("the term y_0");
        g[i].2 = (g[i].2 + 1) % P;
        let perm = swap_first_two_blocks(t, 2).unwrap();
        assert!(check_symmetry(&ring, &perm, std::slice::from_ref(&g)).is_err());
        let limits = HkyLimits::default();
        let opts = FastOptions {
            symmetry: Some(perm),
            ..FastOptions::default()
        };
        let (cl, st) = closure_fast(&ring, std::slice::from_ref(&g), 4, &limits, &opts);
        assert!(st.symmetry.starts_with("rejected"), "{}", st.symmetry);
        let reference = crate::cryptanalysis::hky_fp_tower::closure(&ring, &[g], 4, &limits);
        assert_eq!(cl.dims, reference.dims);
        // unequal towers are rejected too
        let mut rules = kummer_rules(t, 2, 3);
        rules[t - 1].d = 4;
        let ring2 = HkyRing::new(P, rules, 0);
        let g2 = symmetric_gen(&ring2, t, 2);
        assert!(check_symmetry(&ring2, &swap_first_two_blocks(t, 2).unwrap(), &[g2]).is_err());
    }

    /// A consistent system: a symmetric generator shifted to vanish at a
    /// point of both towers, so `1` is in no closure and every closure is a
    /// proper subspace whose dimensions an arithmetic fault changes.
    fn consistent(t: usize) -> (HkyRing, Vec<NfPoly>) {
        let y0 = 123_456_789u64;
        let mut top = y0;
        for _ in 0..t {
            top = top * top % P;
        }
        let ring = HkyRing::new(P, kummer_rules(t, 2, top), 0);
        let mut point = Vec::new();
        for start in [y0, P - y0] {
            let mut y = start;
            for _ in 0..t {
                point.push(y);
                y = y * y % P;
            }
        }
        let g = symmetric_gen(&ring, t, 3);
        let value = ring.eval(&g, &point);
        let n = ring.n_vars();
        let shifted: Vec<(Vec<u32>, u64)> = g
            .iter()
            .map(|&(m, _, x)| ((0..n).map(|j| (m >> j & 1) as u32).collect(), x))
            .chain(std::iter::once((vec![0u32; n], P - value)))
            .collect();
        let h = ring.from_raw(&shifted);
        assert_eq!(ring.eval(&h, &point), 0);
        (ring, vec![h])
    }

    #[test]
    fn a_consistent_system_never_contains_one() {
        let t = 4;
        let (ring, gens) = consistent(t);
        let reference = fall_degrees(&ring, &gens, 1, 5, &HkyLimits::default());
        assert_eq!(reference.last_fall, None);
        assert!(reference
            .closures
            .iter()
            .all(|c| u128::from(c.dim()) < c.cols));
        all_variants(&ring, &gens, swap_first_two_blocks(t, 2), 5);
    }

    #[test]
    fn the_comparison_detects_a_faulty_engine() {
        // Negative controls of the equality gate: an engine that skips the
        // products by one variable, and one whose delayed reduction is off
        // by one, must not agree with the reference. The system is
        // consistent: on a closure that contains 1 the dimensions are those
        // of the whole space whatever the arithmetic did, so only closures
        // without 1 give the comparison its power.
        let t = 4;
        let (ring, gens) = consistent(t);
        let limits = HkyLimits::default();
        let a = fall_degrees(&ring, &gens, 1, 5, &limits);
        for fault in [1u8, 2] {
            for sym in [false, true] {
                for level in [1u8, 2, 4] {
                    let opts = FastOptions {
                        fault,
                        level,
                        symmetry: if sym {
                            swap_first_two_blocks(t, 2)
                        } else {
                            None
                        },
                        dense_num: 1,
                        dense_den: 8,
                        ..FastOptions::default()
                    };
                    let b = fall_degrees_fast(&ring, &gens, 1, 5, &limits, &opts);
                    let differ = a.last_fall != b.last_fall
                        || a.closures.len() != b.closures.len()
                        || a.closures
                            .iter()
                            .zip(&b.closures)
                            .any(|(x, y)| x.dims != y.dims);
                    assert!(
                        differ,
                        "fault {fault} sym {sym} level {level} must not pass"
                    );
                }
            }
        }
    }

    #[test]
    fn dense_kernel_matches_the_definition() {
        let p2 = P * P;
        let row: Vec<u32> = (0..1000u64)
            .map(|i| ((i * 2654435761) % P) as u32)
            .collect();
        let f = P - 12345;
        for simd in [Simd::Portable, detect_simd()] {
            let mut acc: Vec<u64> = (0..1000u64).map(|i| (i * i * 40503 + 7) % p2).collect();
            let want: Vec<u64> = acc
                .iter()
                .zip(&row)
                .map(|(&a, &x)| {
                    ((u128::from(a) + u128::from(f) * u128::from(x)) % u128::from(P)) as u64
                })
                .collect();
            axpy_dense(simd, &mut acc, &row, f, p2);
            assert!(acc.iter().all(|&a| a < p2));
            let got: Vec<u64> = acc.iter().map(|&a| a % P).collect();
            assert_eq!(got, want, "{simd:?}");
        }
        let m = (u128::from(u64::MAX) / u128::from(P)) as u64;
        for v in [0, 1, P - 1, P, P + 1, p2 - 1, p2 / 2, 123456789012345678] {
            assert_eq!(barrett(v, P, m), v % P);
        }
    }
}
