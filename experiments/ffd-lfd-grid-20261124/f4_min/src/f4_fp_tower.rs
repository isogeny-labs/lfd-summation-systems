//! # F4 over a tower quotient ring: the sparse engine for PKM systems.
//!
//! The decomposition systems of Petit–Kosters–Messeng tower factor bases
//! (`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`) live in
//! the quotient `R = F_p[y, u] / (towers)`. Every tower variable `y_j`
//! satisfies one quadratic relation
//!
//! ```text
//!     y_j² = a·y_j·y_{j+1} + b·y_{j+1} + c·y_j + d
//! ```
//!
//! with `y_{j+1}` the next level of the same tower. The relation for each
//! family:
//! - Kummer: `y_j² = y_{j+1}`;
//! - Dickson: `y_j² = y_{j+1} + 2`;
//! - a 2-isogeny step: `y_j² = y_j·y_{j+1} − ξ·y_{j+1} + ξ·y_j − τ`.
//!
//! On the top level `a = b = 0`. Under grevlex with every level larger than
//! the next, these relations are a Gröbner basis with leading monomials
//! `y_j²` (Lemma 1 of the note). Every element of `R` therefore has a
//! unique normal form that is square-free in the tower variables. The
//! chain presentation's free unknowns `u_k`, the abscissae of partial sums,
//! stay ordinary polynomial variables.
//!
//! This engine is Faugère's F4 on `R`, built like
//! [`pq_f4_f2`](crate::cryptanalysis::pq_f4_f2), which is F4 on the boolean
//! ring. Its parts:
//! - **Monomials** are a 64-bit tower mask plus up to four free exponents,
//!   packed into one `u128` whose integer order is grevlex with variable 0
//!   largest, the order [`f4_fp`](crate::cryptanalysis::f4_fp) uses.
//! - **Products are rewritten by the tower as they are formed**, so the
//!   tower equations never enter a matrix.
//! - **A tower pair `(g, y)`** exists for every tower variable `y` of
//!   `LM(g)`. It is the S-polynomial of `g` against `y² − …`, at degree
//!   `deg LM(g) + 1`, and enters the matrix as the single row `NF(y·g)`.
//!   Without these pairs the basis is not closed in `R`.
//! - **Critical pairs** go through Gebauer–Möller (Becker–Weispfenning
//!   `UPDATE`). Tower pairs are never pruned and never used to prune, as the
//!   boolean engine treats its field pairs.
//! - **The elimination follows Faugère–Lachartre.** Symbolic preprocessing
//!   gives reducer rows with distinct leading columns, and the shortest
//!   S-row at each lcm column joins them. The columns no reducer leads are
//!   then exactly the monomials without a divisor. The reducer block is
//!   reduced to dense rows over those columns (`B' = A⁻¹B`, the rightmost
//!   lead first), every S-row is reduced by it in one pass (`D − C·B'`), and
//!   the residues are echelonized in chunks, each reduced in parallel by the
//!   pivots found before it. Every pivot is a new leading monomial.
//! - **A full echelon ends the step's elimination.** Once the residues'
//!   echelon has a pivot on every column without a divisor, every later
//!   S-row reduces to zero, and a pivot is never rewritten once found: the
//!   remaining S-rows are skipped (`full_rank_exit`, on by default). The
//!   basis, the trace and the solving degree are the same with it off; only
//!   the multiply-adds and the count of residual rows fall. It matters most
//!   in a refutation's last step, which in the rounds' systems ends with a
//!   full echelon (note §15).
//! - **Rows are formed on demand.** A row is kept as the product that makes
//!   it and formed again, a block at a time, when it is eliminated, so a
//!   step's memory is `B'` rather than the matrix. The scan keeps the
//!   products it forms for the elimination to read back, up to `2²³` terms a
//!   step (64 MiB of stored term payload; the vectors' spare capacity comes
//!   on top, see `PRODUCT_CACHE_TERMS`); the rest of a larger step is
//!   formed again.
//! - **The basis keeps the echelon's rows.** A new element is the pivot row
//!   the elimination leaves, dense over the step's columns without a divisor
//!   from its lead on, and the elements a step adds share that column list:
//!   four bytes an entry, zeros included, where a sparse term with its
//!   `u128` monomial costs twenty. The new elements' tails are mostly
//!   non-zero, and late in a large run the basis is most of the memory.
//! - **Arithmetic** below `2³¹` is lazy: an accumulator is reduced only when
//!   it crosses `2⁶³`, and the loop is dispatched at run time to AVX2 or
//!   AVX-512. Above `2³¹` every step is reduced (Barrett).
//!
//! What it measures:
//! - `solving_degree_max`: the highest degree of a step at which the basis
//!   gained an element. This is the solving degree of `f4_fp` and of the
//!   framework (`docs/ic/FRAMEWORK.md` §3).
//! - The unit: `F_p` multiply-adds in the eliminations.
//!
//! The input polynomials must already be tower normal forms. [`TowerRing::from_raw`]
//! reduces arbitrary ones. The engine therefore measures the note's
//! `reduced` presentation (§2.2): a summation polynomial of raw degree 4
//! enters at its normal form's degree, which is 2 on a Kummer tower.
//!
//! Limits: at most 64 tower variables, 4 free variables with exponents
//! below 256, and `p < 2³²`.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use rayon::prelude::*;

// ── Arithmetic in F_p, p < 2^32 ──────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub(crate) struct Fp {
    pub(crate) p: u64,
    m: u64,
}

impl Fp {
    pub(crate) fn new(p: u64) -> Self {
        assert!(
            (3..1 << 32).contains(&p),
            "f4_fp_tower: the prime must be odd and below 2^32"
        );
        Fp {
            p,
            m: (u128::from(u64::MAX) + 1).div_euclid(u128::from(p)) as u64,
        }
    }
    /// `x mod p` for any `x < 2⁶⁴` (Barrett; the quotient estimate is low
    /// by at most one).
    #[inline(always)]
    pub(crate) fn reduce(self, x: u64) -> u64 {
        let q = ((u128::from(x) * u128::from(self.m)) >> 64) as u64;
        let r = x - q * self.p;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }
    #[inline(always)]
    pub(crate) fn mul(self, a: u64, b: u64) -> u64 {
        self.reduce(a * b)
    }
    #[inline(always)]
    pub(crate) fn add(self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline(always)]
    pub(crate) fn neg(self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    pub(crate) fn pow(self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        a %= self.p;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    pub(crate) fn inv(self, a: u64) -> u64 {
        assert!(!a.is_multiple_of(self.p), "f4_fp_tower: inverse of zero");
        self.pow(a, self.p - 2)
    }

    /// Montgomery's trick batch inversion over `F_p`: computes `vals[i]^{-1} mod p`
    /// for every entry using a single modular inversion and `3*(k - 1)` multiplications.
    /// Returns `None` immediately if any entry is `0 mod p` (Schwartz–Zippel pivot check).
    pub(crate) fn batch_inv(self, vals: &[u64]) -> Option<Vec<u64>> {
        if vals.is_empty() {
            return Some(Vec::new());
        }
        let mut prefix = Vec::with_capacity(vals.len());
        let mut acc = 1u64;
        for &v in vals {
            let vr = v % self.p;
            if vr == 0 {
                return None;
            }
            prefix.push(acc);
            acc = self.mul(acc, vr);
        }
        let mut inv_acc = self.inv(acc);
        let mut out = vec![0u64; vals.len()];
        for i in (0..vals.len()).rev() {
            let vr = vals[i] % self.p;
            out[i] = self.mul(inv_acc, prefix[i]);
            inv_acc = self.mul(inv_acc, vr);
        }
        Some(out)
    }
}

// ── Monomials ────────────────────────────────────────────────────────

/// A monomial of `R`, stored as its grevlex sort key: the larger key is
/// the larger monomial.
///
/// Layout, least significant first:
/// - bits `0..64`: the complement of the tower mask. Bit `j` of the mask is
///   `y_j`, variable `j`.
/// - bits `64..96`: for free variable `f` (global index `n_tower + f`), the
///   byte `255 − e_f` at bit `64 + 8f`.
/// - bits `96..104`: the total degree.
///
/// Grevlex compares degrees first, then the exponent of the *last*
/// variable, smaller winning, then the one before it, and so on. The last
/// free variable is the most significant byte of the free field, and the
/// free field sits above the tower field, so integer comparison is
/// grevlex. Within the tower field the highest differing bit decides,
/// which is the last differing variable, and the complement makes a 0
/// there, the smaller exponent, the larger key.
pub type Mono = u128;

const FREE_SHIFT: u32 = 64;
const DEG_SHIFT: u32 = 96;
/// Free variables a monomial can carry.
pub const MAX_FREE: usize = 4;

/// Positive integer weights `1..=8` on the tower variables: the grading of
/// a weighted run ([`TowerF4Options::weights`]).
///
/// The weight of variable `j` is `1 + Σ_b 2^b·bit_j(planes[b])`, so the
/// default (all planes zero) is the unit grading, under which every key is
/// the one [`mono`] makes. A weighted degree with positive weights, broken
/// by the same reverse-lexicographic comparison of the key's low half, is a
/// monomial order. Its use here: a radix-`2^k` Kummer tower `u^{2^k} = u'`
/// is the radix-2 tower `y_0 = u, y_{i+1} = y_i²` with weights
/// `1, 2, …, 2^{k−1}` on each block of `k` levels. The ring, its rules and its
/// standard monomials are the same; only the grading, and so the order in
/// which F4 meets the monomials, changes. [`TowerWeights::check`] says when
/// the rules still have leading monomials `y_j²`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct TowerWeights {
    planes: [u64; 3],
}

impl TowerWeights {
    /// The unit grading.
    pub const UNIT: TowerWeights = TowerWeights { planes: [0; 3] };

    /// Weights given per tower variable, each in `1..=8`; variables past the
    /// slice's end get weight 1.
    pub fn from_slice(w: &[u32]) -> Self {
        assert!(w.len() <= 64, "f4_fp_tower: at most 64 tower weights");
        let mut planes = [0u64; 3];
        for (j, &x) in w.iter().enumerate() {
            assert!(
                (1..=8).contains(&x),
                "f4_fp_tower: weight {x} of y_{j} is not in 1..=8"
            );
            let e = x - 1;
            for (b, plane) in planes.iter_mut().enumerate() {
                if (e >> b) & 1 == 1 {
                    *plane |= 1u64 << j;
                }
            }
        }
        TowerWeights { planes }
    }

    /// Radix-`2^k` blocks on `n` tower variables laid out as `blocks` towers
    /// of `t` consecutive levels each (`n = blocks·t`): level `i` of every
    /// tower gets weight `2^(i mod k)`. `k = 1` is the unit grading.
    pub fn radix_blocks(blocks: usize, t: usize, k: u32) -> Self {
        assert!(
            (1..=4).contains(&k),
            "f4_fp_tower: radix exponent k must be 1..=4"
        );
        let w: Vec<u32> = (0..blocks)
            .flat_map(|_| (0..t).map(move |i| 1u32 << (i as u32 % k)))
            .collect();
        Self::from_slice(&w)
    }

    pub fn is_unit(self) -> bool {
        self.planes == [0; 3]
    }

    /// The weight of tower variable `j`.
    #[inline]
    pub fn weight(self, j: u32) -> u32 {
        1 + self
            .planes
            .iter()
            .enumerate()
            .map(|(b, p)| (((p >> j) & 1) as u32) << b)
            .sum::<u32>()
    }

    /// The weighted degree of a tower mask.
    #[inline]
    pub fn mask_degree(self, mask: u64) -> u32 {
        mask.count_ones()
            + (mask & self.planes[0]).count_ones()
            + 2 * (mask & self.planes[1]).count_ones()
            + 4 * (mask & self.planes[2]).count_ones()
    }

    /// Panics unless every rule of `ring` keeps `y_j²` as its leading
    /// monomial in this grading, so that the tower stays a Gröbner basis
    /// with pairwise coprime leads (Lemma 1 of the note):
    /// - a `b·y_next` term needs `w(next) ≤ 2·w(j)` (at equality the reverse
    ///   lexicographic tie-break puts `y_next`, the later variable, below);
    /// - an `a·y_j·y_next` term needs `w(next) ≤ w(j)`, by the same argument;
    /// - `c·y_j` and `d` are always below.
    pub fn check(self, ring: &TowerRing) {
        assert!(
            self.planes
                .iter()
                .all(|p| ring.n_tower() == 64 || p >> ring.n_tower() == 0),
            "f4_fp_tower: weights set on variables past the tower"
        );
        for (j, r) in ring.rules.iter().enumerate() {
            let wj = self.weight(j as u32);
            if let Some(nx) = r.next {
                let wn = self.weight(u32::from(nx));
                if r.b != 0 {
                    assert!(
                        wn <= 2 * wj,
                        "f4_fp_tower: weight of y_{nx} ({wn}) above twice that of y_{j} ({wj}): y_{j}² is not the leading monomial"
                    );
                }
                if r.a != 0 {
                    assert!(
                        wn <= wj,
                        "f4_fp_tower: rule {j} has a y_{j}·y_{nx} term and w(y_{nx}) = {wn} > w(y_{j}) = {wj}"
                    );
                }
            }
        }
    }
}

/// The key of degree `deg`, tower mask `mask` and free exponents `free`.
#[inline]
fn key_with_degree(deg: u32, mask: u64, free: [u8; MAX_FREE]) -> Mono {
    let mut key = (u128::from(deg) << DEG_SHIFT) | u128::from(!mask);
    for (k, &e) in free.iter().enumerate() {
        key |= u128::from(255 - e) << (FREE_SHIFT + 8 * k as u32);
    }
    key
}

/// The monomial with tower mask `mask` and free exponents `free`, in the
/// unit grading.
#[inline]
pub fn mono(mask: u64, free: [u8; MAX_FREE]) -> Mono {
    mono_w(mask, free, TowerWeights::UNIT)
}

/// The monomial with tower mask `mask` and free exponents `free`, in the
/// grading `w`.
#[inline]
pub fn mono_w(mask: u64, free: [u8; MAX_FREE], w: TowerWeights) -> Mono {
    let deg = w.mask_degree(mask) + free.iter().map(|&e| u32::from(e)).sum::<u32>();
    key_with_degree(deg, mask, free)
}

/// The tower variables of `m`, as a mask.
#[inline]
pub fn mask_of(m: Mono) -> u64 {
    !(m as u64)
}

/// The free exponents of `m`.
#[inline]
pub fn free_of(m: Mono) -> [u8; MAX_FREE] {
    let mut f = [0u8; MAX_FREE];
    for (k, e) in f.iter_mut().enumerate() {
        *e = 255 - ((m >> (FREE_SHIFT + 8 * k as u32)) & 0xff) as u8;
    }
    f
}

/// Total degree of `m`.
#[inline]
pub fn degree_of(m: Mono) -> u32 {
    (m >> DEG_SHIFT) as u32
}

/// The monomial `1`.
pub const ONE: Mono = (0xffff_ffffu128 << FREE_SHIFT) | (u64::MAX as u128);

/// `a | b` in the polynomial ring.
#[inline]
pub(crate) fn divides(a: Mono, b: Mono) -> bool {
    if mask_of(a) & !mask_of(b) != 0 || degree_of(a) > degree_of(b) {
        return false;
    }
    let (fa, fb) = (free_of(a), free_of(b));
    fa.iter().zip(&fb).all(|(x, y)| x <= y)
}

/// `a | b` from the masks and free exponents alone, whatever grading
/// either key is in.
#[inline]
pub(crate) fn divides_ungraded(a: Mono, b: Mono) -> bool {
    if mask_of(a) & !mask_of(b) != 0 {
        return false;
    }
    let (fa, fb) = (free_of(a), free_of(b));
    fa.iter().zip(&fb).all(|(x, y)| x <= y)
}

#[inline]
pub(crate) fn lcm(a: Mono, b: Mono) -> Mono {
    lcm_w(a, b, TowerWeights::UNIT)
}

/// The lcm of `a` and `b`, keyed in the grading `w` (both must be).
#[inline]
pub(crate) fn lcm_w(a: Mono, b: Mono, w: TowerWeights) -> Mono {
    let (fa, fb) = (free_of(a), free_of(b));
    let mut f = [0u8; MAX_FREE];
    for k in 0..MAX_FREE {
        f[k] = fa[k].max(fb[k]);
    }
    mono_w(mask_of(a) | mask_of(b), f, w)
}

/// `a / b`, for `b | a`. Its degree is `deg a − deg b` in any grading, so
/// the quotient is keyed in the grading of its arguments.
#[inline]
pub(crate) fn quotient(a: Mono, b: Mono) -> Mono {
    let (fa, fb) = (free_of(a), free_of(b));
    let mut f = [0u8; MAX_FREE];
    for k in 0..MAX_FREE {
        f[k] = fa[k] - fb[k];
    }
    key_with_degree(degree_of(a) - degree_of(b), mask_of(a) & !mask_of(b), f)
}

#[inline]
pub(crate) fn coprime(a: Mono, b: Mono) -> bool {
    if mask_of(a) & mask_of(b) != 0 {
        return false;
    }
    let (fa, fb) = (free_of(a), free_of(b));
    fa.iter().zip(&fb).all(|(x, y)| *x == 0 || *y == 0)
}

#[inline]
pub(crate) fn add_free(a: [u8; MAX_FREE], b: [u8; MAX_FREE]) -> [u8; MAX_FREE] {
    let mut f = [0u8; MAX_FREE];
    for k in 0..MAX_FREE {
        f[k] = a[k]
            .checked_add(b[k])
            .expect("f4_fp_tower: a free exponent passed 255");
    }
    f
}

/// splitmix64 over both halves of a monomial key: the low half is a
/// complemented mask whose low bits barely vary between the monomials of
/// one system, so a multiplicative hash would crowd them into few buckets
/// (see `fx_hash::MaskHasher`).
#[derive(Default, Clone, Copy)]
pub(crate) struct MonoHasher(u64);

impl Hasher for MonoHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut w = [0u8; 8];
            w[..chunk.len()].copy_from_slice(chunk);
            self.write_u64(u64::from_le_bytes(w));
        }
    }
    #[inline]
    fn write_u64(&mut self, value: u64) {
        let mut z = (self.0 ^ value).wrapping_add(0x9e37_79b9_7f4a_7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        self.0 = z ^ (z >> 31);
    }
    #[inline]
    fn write_u128(&mut self, value: u128) {
        self.write_u64(value as u64);
        self.write_u64((value >> 64) as u64);
    }
    #[inline]
    fn write_usize(&mut self, value: usize) {
        self.write_u64(value as u64);
    }
}

pub(crate) type MonoMap<V> = HashMap<Mono, V, BuildHasherDefault<MonoHasher>>;

// ── The ring ─────────────────────────────────────────────────────────

/// How one tower variable's square rewrites:
/// `y_j² = a·y_j·y_next + b·y_next + c·y_j + d`, all coefficients in
/// `[0, p)`. On a tower's top level `next` is `None` and `a = b = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SquareRule {
    pub next: Option<u8>,
    pub a: u64,
    pub b: u64,
    pub c: u64,
    pub d: u64,
}

/// Combine repeated masks and drop zero coefficients, in place.
fn merge_masks(fp: Fp, terms: &mut Vec<(u64, u64)>) {
    if terms.len() < 2 {
        return;
    }
    terms.sort_unstable_by_key(|t| t.0);
    let mut out = 0usize;
    let mut i = 0usize;
    while i < terms.len() {
        let m = terms[i].0;
        let mut c = 0u64;
        while i < terms.len() && terms[i].0 == m {
            c = fp.add(c, terms[i].1);
            i += 1;
        }
        if c != 0 {
            terms[out] = (m, c);
            out += 1;
        }
    }
    terms.truncate(out);
}

/// The two working lists of [`TowerRing::mask_product_into`].
#[derive(Default)]
struct MaskScratch {
    cur: Vec<(u64, u64)>,
    next: Vec<(u64, u64)>,
}

/// `F_p[y_0, …, y_{n−1}, u_0, …] / (y_j² − rule_j)`. The tower variables
/// come first; the free variables follow, at global indices
/// `n_tower .. n_tower + n_free`.
#[derive(Clone, Debug)]
pub struct TowerRing {
    pub p: u64,
    /// One rule per tower variable.
    pub rules: Vec<SquareRule>,
    pub n_free: usize,
}

impl TowerRing {
    pub fn new(p: u64, rules: Vec<SquareRule>, n_free: usize) -> Self {
        let ring = Self { p, rules, n_free };
        ring.check();
        ring
    }

    pub fn n_tower(&self) -> usize {
        self.rules.len()
    }

    pub fn n_vars(&self) -> usize {
        self.rules.len() + self.n_free
    }

    /// Panics unless the rules describe a tower this engine handles: every
    /// `next` is a later tower variable, so that each rewrite lowers the
    /// monomial in grevlex, and a top rule has `a = b = 0`.
    pub fn check(&self) {
        assert!(
            self.rules.len() <= 64,
            "f4_fp_tower: at most 64 tower variables"
        );
        assert!(
            self.n_free <= MAX_FREE,
            "f4_fp_tower: at most {MAX_FREE} free variables"
        );
        for (j, r) in self.rules.iter().enumerate() {
            for &x in &[r.a, r.b, r.c, r.d] {
                assert!(x < self.p, "f4_fp_tower: rule {j} has a coefficient ≥ p");
            }
            match r.next {
                Some(nx) => assert!(
                    (nx as usize) > j && (nx as usize) < self.rules.len(),
                    "f4_fp_tower: rule {j} must point at a later tower variable"
                ),
                None => assert!(
                    r.a == 0 && r.b == 0,
                    "f4_fp_tower: top rule {j} cannot use a next level"
                ),
            }
        }
    }

    /// `y_j · (mask, coef)` in normal form, pushed onto `out` as tower
    /// masks with coefficients (repeats allowed).
    fn mul_var_into(&self, fp: Fp, m: u64, j: usize, coef: u64, out: &mut Vec<(u64, u64)>) {
        let bit = 1u64 << j;
        if m & bit == 0 {
            out.push((m | bit, coef));
            return;
        }
        let r = self.rules[j];
        let m0 = m & !bit;
        if r.d != 0 {
            out.push((m0, fp.mul(coef, r.d)));
        }
        if r.c != 0 {
            out.push((m, fp.mul(coef, r.c)));
        }
        if let Some(nx) = r.next {
            let nx = nx as usize;
            if r.b != 0 {
                self.mul_var_into(fp, m0, nx, fp.mul(coef, r.b), out);
            }
            if r.a != 0 {
                self.mul_var_into(fp, m, nx, fp.mul(coef, r.a), out);
            }
        }
    }

    /// `coef · y^ma · y^mb` in normal form (tower parts only), pushed onto
    /// `out` with repeats allowed. `scratch` holds the working lists, so a
    /// caller that forms many products allocates them once.
    fn mask_product_into(
        &self,
        fp: Fp,
        ma: u64,
        mb: u64,
        coef: u64,
        out: &mut Vec<(u64, u64)>,
        scratch: &mut MaskScratch,
    ) {
        let common = ma & mb;
        if common == 0 {
            out.push((ma | mb, coef));
            return;
        }
        // y^ma · y^mb = y^(ma ∪ mb) · Π_{j ∈ ma ∩ mb} y_j.
        let MaskScratch { cur, next } = scratch;
        cur.clear();
        cur.push((ma | mb, coef));
        let mut bits = common;
        while bits != 0 {
            let j = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            next.clear();
            for &(m, c) in cur.iter() {
                self.mul_var_into(fp, m, j, c, next);
            }
            // Without merging, a chain of rewrites repeats monomials and
            // the list grows exponentially in the number of common bits.
            merge_masks(fp, next);
            std::mem::swap(cur, next);
        }
        out.extend_from_slice(cur);
    }

    /// `u · g` in normal form (unit grading).
    pub fn mul_mono(&self, g: &RPoly, u: Mono) -> RPoly {
        let fp = Fp::new(self.p);
        self.mul_mono_fp(fp, g, u, TowerWeights::UNIT)
    }

    /// `u · g` in normal form, with `g`, `u` and the result keyed in `w`.
    pub fn mul_mono_w(&self, g: &RPoly, u: Mono, w: TowerWeights) -> RPoly {
        let fp = Fp::new(self.p);
        self.mul_mono_fp(fp, g, u, w)
    }

    fn mul_mono_fp(&self, fp: Fp, g: &RPoly, u: Mono, w: TowerWeights) -> RPoly {
        let terms = g.monos.iter().copied().zip(g.coefs.iter().copied());
        self.mul_terms_fp(fp, terms, g.len(), u, w)
    }

    fn mul_elem(&self, fp: Fp, g: &Elem, u: Mono, w: TowerWeights) -> RPoly {
        self.mul_terms_fp(fp, g.terms(), g.len(), u, w)
    }

    /// `u · g` in normal form, for `g` given by its terms, descending and
    /// with non-zero coefficients (`len` of them). `g`, `u` and the result
    /// are keyed in the grading `w`.
    pub(crate) fn mul_terms_fp(
        &self,
        fp: Fp,
        g: impl Iterator<Item = (Mono, u32)>,
        len: usize,
        u: Mono,
        w: TowerWeights,
    ) -> RPoly {
        let (mu, fu) = (mask_of(u), free_of(u));
        // When `u` has no free variable, the key of a plain product is the
        // key of the term plus a constant: the degree rises by that of `u`,
        // the complemented mask loses `mu` (the masks are disjoint, so those
        // bits are set in it and nothing borrows) and the free bytes stay.
        let plain_shift =
            (fu == [0; MAX_FREE]).then(|| (u128::from(degree_of(u)) << DEG_SHIFT) - u128::from(mu));
        // A term that shares no tower variable with `u` is a plain product,
        // injective and order-preserving: those terms stay distinct and
        // descending, and go straight into the result. Only the terms that
        // rewrite (`carried`) come out in no order, and only they are sorted.
        let mut monos: Vec<Mono> = Vec::with_capacity(len);
        let mut coefs: Vec<u32> = Vec::with_capacity(len);
        let mut carried: Vec<(Mono, u64)> = Vec::new();
        let mut tmp: Vec<(u64, u64)> = Vec::new();
        let mut scratch = MaskScratch::default();
        for (m, c) in g {
            let mm = mask_of(m);
            if mm & mu == 0 {
                monos.push(match plain_shift {
                    Some(shift) => m + shift,
                    None => mono_w(mm | mu, add_free(free_of(m), fu), w),
                });
                coefs.push(c);
            } else {
                let f = add_free(free_of(m), fu);
                tmp.clear();
                self.mask_product_into(fp, mm, mu, u64::from(c), &mut tmp, &mut scratch);
                carried.extend(tmp.iter().map(|&(x, cc)| (mono_w(x, f, w), cc)));
            }
        }
        if carried.is_empty() {
            return RPoly { monos, coefs };
        }
        carried.sort_unstable_by_key(|t| std::cmp::Reverse(t.0));
        // Merge the two descending lists, adding the coefficients of a
        // repeated monomial and dropping a zero sum: what `from_terms` does
        // to the concatenation. Both lists hold coefficients below `p`.
        let mut out = RPoly {
            monos: Vec::with_capacity(monos.len() + carried.len()),
            coefs: Vec::with_capacity(monos.len() + carried.len()),
        };
        let (mut i, mut j) = (0, 0);
        loop {
            let m = match (monos.get(i), carried.get(j)) {
                (Some(&a), Some(&(b, _))) => a.max(b),
                (Some(&a), None) => a,
                (None, Some(&(b, _))) => b,
                (None, None) => break,
            };
            let mut c = 0u64;
            if i < monos.len() && monos[i] == m {
                c = u64::from(coefs[i]);
                i += 1;
            }
            while j < carried.len() && carried[j].0 == m {
                c = fp.add(c, carried[j].1);
                j += 1;
            }
            if c != 0 {
                out.monos.push(m);
                out.coefs.push(c as u32);
            }
        }
        out
    }

    /// The normal form of a polynomial given by arbitrary exponent vectors
    /// (length [`n_vars`](Self::n_vars)) and coefficients.
    pub fn from_raw(&self, terms: &[(Vec<u32>, u64)]) -> RPoly {
        let fp = Fp::new(self.p);
        let n_t = self.n_tower();
        let mut out: Vec<(Mono, u64)> = Vec::new();
        let mut cur: Vec<(u64, u64)> = Vec::new();
        let mut next: Vec<(u64, u64)> = Vec::new();
        for (e, c) in terms {
            assert_eq!(e.len(), self.n_vars(), "f4_fp_tower: exponent length");
            let c = c % self.p;
            if c == 0 {
                continue;
            }
            let mut f = [0u8; MAX_FREE];
            for k in 0..self.n_free {
                f[k] = u8::try_from(e[n_t + k]).expect("f4_fp_tower: free exponent above 255");
            }
            cur.clear();
            cur.push((0, c));
            for (j, &ej) in e.iter().take(n_t).enumerate() {
                for _ in 0..ej {
                    next.clear();
                    for &(m, cc) in &cur {
                        self.mul_var_into(fp, m, j, cc, &mut next);
                    }
                    merge_masks(fp, &mut next);
                    std::mem::swap(&mut cur, &mut next);
                }
            }
            for &(m, cc) in &cur {
                out.push((mono(m, f), cc));
            }
        }
        RPoly::from_terms(fp, out)
    }

    /// Parent-to-child basis specialization in the same [`TowerRing`]: evaluates tower
    /// variable `var` at `val mod p` in every polynomial of `basis`, merges terms in
    /// grading `w`, and prepends the linear generator `y_{var} - val`.
    /// Returns `(child_inputs, fp_mults)`.
    pub fn specialize_basis_at_var(
        &self,
        basis: &[RPoly],
        var: usize,
        val: u64,
        w: TowerWeights,
    ) -> (Vec<RPoly>, u64) {
        assert!(var < self.n_tower(), "f4_fp_tower: var out of range");
        let fp = Fp::new(self.p);
        let val = val % self.p;
        let bit = 1u64 << var;
        let mut mults = 0u64;
        let mut out = Vec::with_capacity(basis.len() + 1);

        // Linear branch generator `y_{var} - val`.
        let mut lin_terms = vec![(mono_w(bit, [0; MAX_FREE], w), 1u64)];
        if val != 0 {
            lin_terms.push((ONE, fp.neg(val)));
        }
        out.push(RPoly::from_terms(fp, lin_terms));

        for g in basis {
            let mut terms = Vec::with_capacity(g.len());
            for (&m, &c) in g.monos.iter().zip(&g.coefs) {
                let mask = mask_of(m);
                let free = free_of(m);
                if (mask & bit) != 0 {
                    if val == 0 {
                        continue;
                    }
                    let cc = if val == 1 {
                        u64::from(c)
                    } else {
                        mults += 1;
                        fp.mul(u64::from(c), val)
                    };
                    terms.push((mono_w(mask & !bit, free, w), cc));
                } else {
                    terms.push((mono_w(mask, free, w), u64::from(c)));
                }
            }
            let sp = RPoly::from_terms(fp, terms);
            if sp.is_unit() {
                return (vec![sp], mults);
            }
            if !sp.is_empty() {
                out.push(sp);
            }
        }
        (out, mults)
    }

    /// Surjective homomorphism specializing a top-level tower variable `var`
    /// (`rules[var].next == None`) to `val mod p` and dropping `var` from the ring
    /// (`n_tower - 1` tower variables). Any lower rule pointing to `var` becomes a
    /// top-level rule with `c' = c + a·val, d' = d + b·val (mod p)`.
    /// Returns `(child_ring, child_basis, child_weights, fp_mults)`.
    pub fn specialize_drop_top_var(
        &self,
        basis: &[RPoly],
        var: usize,
        val: u64,
        w: TowerWeights,
    ) -> (TowerRing, Vec<RPoly>, TowerWeights, u64) {
        assert!(var < self.n_tower(), "f4_fp_tower: var out of range");
        assert!(
            self.rules[var].next.is_none(),
            "f4_fp_tower: specialize_drop_top_var requires a top-level variable"
        );
        let fp = Fp::new(self.p);
        let val = val % self.p;
        let mut mults = 0u64;

        let mut new_rules = Vec::with_capacity(self.n_tower() - 1);
        for (j, r) in self.rules.iter().enumerate() {
            if j == var {
                continue;
            }
            match r.next {
                Some(nxt) if usize::from(nxt) == var => {
                    let av = if r.a % self.p == 0 || val == 0 {
                        0
                    } else {
                        mults += 1;
                        fp.mul(r.a % self.p, val)
                    };
                    let bv = if r.b % self.p == 0 || val == 0 {
                        0
                    } else {
                        mults += 1;
                        fp.mul(r.b % self.p, val)
                    };
                    new_rules.push(SquareRule {
                        next: None,
                        a: 0,
                        b: 0,
                        c: fp.add(r.c % self.p, av),
                        d: fp.add(r.d % self.p, bv),
                    });
                }
                Some(nxt) if usize::from(nxt) > var => {
                    new_rules.push(SquareRule {
                        next: Some(nxt - 1),
                        ..*r
                    });
                }
                _ => new_rules.push(*r),
            }
        }
        let child_ring = TowerRing::new(self.p, new_rules, self.n_free);
        let child_w = if w.is_unit() {
            TowerWeights::UNIT
        } else {
            let ws: Vec<u32> = (0..self.n_tower())
                .filter(|&j| j != var)
                .map(|j| w.weight(j as u32))
                .collect();
            TowerWeights::from_slice(&ws)
        };

        let low_mask = if var == 0 { 0u64 } else { (1u64 << var) - 1 };
        let bit = 1u64 << var;
        let high_mask = if var == 63 {
            0u64
        } else {
            !((1u64 << (var + 1)) - 1)
        };

        let mut child_basis = Vec::with_capacity(basis.len());
        for g in basis {
            let mut terms = Vec::with_capacity(g.len());
            for (&m, &c) in g.monos.iter().zip(&g.coefs) {
                let mask = mask_of(m);
                let free = free_of(m);
                let new_mask = (mask & low_mask) | ((mask & high_mask) >> 1);
                if (mask & bit) != 0 {
                    if val == 0 {
                        continue;
                    }
                    let cc = if val == 1 {
                        u64::from(c)
                    } else {
                        mults += 1;
                        fp.mul(u64::from(c), val)
                    };
                    terms.push((mono_w(new_mask, free, child_w), cc));
                } else {
                    terms.push((mono_w(new_mask, free, child_w), u64::from(c)));
                }
            }
            let sp = RPoly::from_terms(fp, terms);
            if sp.is_unit() {
                return (child_ring, vec![sp], child_w, mults);
            }
            if !sp.is_empty() {
                child_basis.push(sp);
            }
        }
        (child_ring, child_basis, child_w, mults)
    }

    /// `g` at a point: one value per variable, tower variables first.
    pub fn eval(&self, g: &RPoly, x: &[u64]) -> u64 {
        let fp = Fp::new(self.p);
        let n_t = self.n_tower();
        let mut acc = 0u64;
        for (&m, &c) in g.monos.iter().zip(&g.coefs) {
            let mut t = u64::from(c);
            let mut bits = mask_of(m);
            while bits != 0 {
                let j = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                t = fp.mul(t, x[j] % self.p);
            }
            for (k, &e) in free_of(m).iter().enumerate().take(self.n_free) {
                if e > 0 {
                    t = fp.mul(t, fp.pow(x[n_t + k], u64::from(e)));
                }
            }
            acc = fp.add(acc, t);
        }
        acc
    }
}

// ── Polynomials ──────────────────────────────────────────────────────

/// An element of `R` in normal form: monomials strictly descending,
/// coefficients non-zero and below `p`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RPoly {
    pub monos: Vec<Mono>,
    pub coefs: Vec<u32>,
}

impl RPoly {
    /// Sort, merge repeats and drop zeros.
    pub(crate) fn from_terms(fp: Fp, mut terms: Vec<(Mono, u64)>) -> RPoly {
        terms.sort_unstable_by_key(|t| std::cmp::Reverse(t.0));
        let mut monos = Vec::with_capacity(terms.len());
        let mut coefs: Vec<u32> = Vec::with_capacity(terms.len());
        let mut i = 0;
        while i < terms.len() {
            let m = terms[i].0;
            let mut c = 0u64;
            while i < terms.len() && terms[i].0 == m {
                c = fp.add(c, fp.reduce(terms[i].1));
                i += 1;
            }
            if c != 0 {
                monos.push(m);
                coefs.push(c as u32);
            }
        }
        RPoly { monos, coefs }
    }

    pub fn len(&self) -> usize {
        self.monos.len()
    }

    pub fn is_empty(&self) -> bool {
        self.monos.is_empty()
    }

    /// Leading monomial.
    pub fn lm(&self) -> Option<Mono> {
        self.monos.first().copied()
    }

    pub fn degree(&self) -> u32 {
        self.monos.iter().map(|&m| degree_of(m)).max().unwrap_or(0)
    }

    /// Whether this is a non-zero constant.
    pub fn is_unit(&self) -> bool {
        self.monos.len() == 1 && self.monos[0] == ONE
    }
}

/// A basis element, kept as the elimination leaves it: the coefficients of
/// the monomials `cols[start..]`, zeros included, from the leading monomial
/// on. The elements one step adds share that step's list of columns without
/// a divisor, so an entry costs four bytes where a sparse term costs twenty,
/// and a new element is the echelon's own row, not a copy of it.
#[derive(Clone, Debug)]
pub(crate) struct Elem {
    /// Monomials, descending.
    pub(crate) cols: Arc<[Mono]>,
    pub(crate) start: u32,
    /// `vals[i]` is the coefficient of `cols[start + i]`; `vals[0]` is not 0.
    pub(crate) vals: Vec<u32>,
    /// Non-zero entries of `vals`.
    pub(crate) len: u32,
}

impl Elem {
    pub(crate) fn new(cols: Arc<[Mono]>, start: u32, vals: Vec<u32>) -> Self {
        debug_assert!(vals.first().is_some_and(|&v| v != 0));
        debug_assert!(start as usize + vals.len() <= cols.len());
        let len = vals.iter().filter(|&&v| v != 0).count() as u32;
        Elem {
            cols,
            start,
            vals,
            len,
        }
    }

    pub(crate) fn from_rpoly(p: RPoly) -> Self {
        let len = p.len() as u32;
        Elem {
            cols: p.monos.into(),
            start: 0,
            vals: p.coefs,
            len,
        }
    }

    /// Leading monomial.
    pub(crate) fn lm(&self) -> Mono {
        self.cols[self.start as usize]
    }

    /// Number of terms.
    pub(crate) fn len(&self) -> usize {
        self.len as usize
    }

    /// The terms, descending.
    pub(crate) fn terms(&self) -> impl Iterator<Item = (Mono, u32)> + '_ {
        self.cols[self.start as usize..]
            .iter()
            .zip(&self.vals)
            .filter(|&(_, &v)| v != 0)
            .map(|(&m, &v)| (m, v))
    }

    pub(crate) fn to_rpoly(&self) -> RPoly {
        let (monos, coefs) = self.terms().unzip();
        RPoly { monos, coefs }
    }
}

// ── Options and report ───────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct TowerF4Options {
    /// Pairs above this degree are not processed; the run then stops and
    /// reports how many were left.
    pub max_degree: u32,
    /// Cooperative stop: the run returns `timed_out` once past it.
    pub deadline: Option<Instant>,
    /// Stop once the leading monomials of the basis leave at most this many
    /// standard monomials (with pairs pending: the basis is then not
    /// certified). The count bounds the number of solutions.
    pub stop_staircase: Option<usize>,
    /// Stop, reporting `oversize`, before building a matrix with more than
    /// this many non-zero entries, or a reduced reducer block with more
    /// than this many dense entries.
    pub max_nnz: Option<u64>,
    /// When set, the bound on the reduced reducer block alone, in place of
    /// `max_nnz`. The block is held in memory, 4 bytes an entry; a
    /// matrix's rows are formed when they are reduced, so its nonzeros
    /// never are.
    pub max_dense: Option<u64>,
    /// Called with each step's number (from 1) and trace as soon as the step
    /// ends, so that a run that dies still leaves its trace.
    pub on_step: Option<fn(usize, &StepTrace)>,
    /// Called once if a stop cuts a step off, at the moment it does, with
    /// the step's number, whether the stop is `oversize` (else the
    /// deadline), and the step as far as it got
    /// ([`TowerF4Report::stopped_step`]). It comes before the report's basis
    /// is built, which for a large stopped run can be the largest
    /// allocation left.
    pub on_stop: Option<fn(usize, bool, &StepTrace)>,
    /// The same for the signature engine's steps.
    pub on_sig_step: Option<fn(usize, &super::sig_fp_tower::SigStep)>,
    /// Stop reducing a step's S-rows once the echelon of their residues
    /// has a pivot on every column the residues live on: every later row
    /// then reduces to zero, and no pivot changes. On by default; off
    /// reproduces the elimination row for row, multiply-adds included.
    pub full_rank_exit: bool,
    /// Learn an [`F4Trace`]: on a run that ends in a refutation, record for
    /// each step the S-rows that gave a new pivot (and the promoted rows),
    /// and the leading monomials of the new elements. Costs nothing in
    /// multiply-adds; the result is [`TowerF4Report::learned`].
    pub record_trace: bool,
    /// Replay a learned trace on another specialisation of the same system
    /// shape: each step reduces only the recorded rows (and only the
    /// reducers they reach), and stops with a failure as soon as the new
    /// elements' leads differ from the trace's. Only a refutation (`1` in
    /// the ideal) is a result of a replay; it is sound whatever the trace
    /// was learned on, because every row computed lies in the ideal. A
    /// replay that does not refute has not decided anything.
    pub replay: Option<Arc<F4Trace>>,
    /// The grading of the run: weights on the tower variables
    /// ([`TowerWeights`]). `None` is the unit grading, and leaves every key,
    /// every step and every count as before. With weights the inputs are
    /// re-keyed at entry, the steps and `solving_degree_max` are weighted
    /// degrees, and the basis comes back keyed in the grading
    /// ([`degree_of`] on it is the weighted degree; [`mask_of`] and
    /// [`free_of`] are unchanged).
    pub weights: Option<TowerWeights>,
}

/// What an F4 run that ended in a refutation learned about its system
/// shape: per step, which rows to reduce, plus a compiled minimal-1-witness
/// linear-combination DAG for zero-overhead Schwartz–Zippel replay.
/// See [`TowerF4Options::replay`].
#[derive(Clone, Debug, Default)]
pub struct F4Trace {
    pub steps: Vec<TraceStep>,
    pub(crate) compiled: Option<Arc<CompiledRefutation>>,
}

/// One step of an [`F4Trace`].
#[derive(Clone, Debug, Default)]
pub struct TraceStep {
    pub degree: u32,
    /// `(multiplier, element)` of the rows kept: S-rows that gave a pivot,
    /// and the S-rows promoted to reducers.
    pub kept: std::collections::HashSet<(Mono, u32)>,
    /// Leading monomials of the step's new elements, ascending.
    pub fresh_lms: Vec<Mono>,
}

impl F4Trace {
    /// Rows kept over all steps (a measure of the trace's size).
    pub fn kept_rows(&self) -> usize {
        self.steps.iter().map(|s| s.kept.len()).sum()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ScatterPlan {
    elem_id: u32,
    expected_lead: u32,
    unit_moves: Vec<(u32, u32)>,
    scaled_moves: Vec<(u32, u32, u32)>,
}

impl ScatterPlan {
    fn eval_row(&self, fp: Fp, elem_vals: &[u32], n_cols: usize, muladds: &mut u64) -> Row {
        if self.scaled_moves.is_empty() {
            // Fast path when tower masks are disjoint: unit_moves are strictly ordered if
            // monomials stay distinct, or we accumulate into a small sparse/dense buffer.
            let mut pairs: Vec<(u32, u64)> = Vec::with_capacity(self.unit_moves.len());
            for &(src, dst) in &self.unit_moves {
                let v = elem_vals[src as usize];
                if v != 0 {
                    pairs.push((dst, u64::from(v)));
                }
            }
            pairs.sort_unstable_by_key(|&(c, _)| c);
            let mut cols = Vec::with_capacity(pairs.len());
            let mut vals = Vec::with_capacity(pairs.len());
            let mut i = 0;
            while i < pairs.len() {
                let c = pairs[i].0;
                let mut acc = pairs[i].1;
                i += 1;
                while i < pairs.len() && pairs[i].0 == c {
                    acc = fp.add(acc, pairs[i].1);
                    i += 1;
                }
                if acc != 0 {
                    cols.push(c);
                    vals.push(acc as u32);
                }
            }
            let _ = n_cols;
            return Row { cols, vals };
        }
        let mut pairs: Vec<(u32, u64)> =
            Vec::with_capacity(self.unit_moves.len() + self.scaled_moves.len());
        for &(src, dst) in &self.unit_moves {
            let v = elem_vals[src as usize];
            if v != 0 {
                pairs.push((dst, u64::from(v)));
            }
        }
        for &(src, dst, factor) in &self.scaled_moves {
            let v = elem_vals[src as usize];
            if v != 0 {
                *muladds += 1;
                pairs.push((dst, fp.mul(u64::from(v), u64::from(factor))));
            }
        }
        pairs.sort_unstable_by_key(|&(c, _)| c);
        let mut cols = Vec::with_capacity(pairs.len());
        let mut vals = Vec::with_capacity(pairs.len());
        let mut i = 0;
        while i < pairs.len() {
            let c = pairs[i].0;
            let mut acc = pairs[i].1;
            i += 1;
            while i < pairs.len() && pairs[i].0 == c {
                acc = fp.add(acc, pairs[i].1);
                i += 1;
            }
            if acc != 0 {
                cols.push(c);
                vals.push(acc as u32);
            }
        }
        Row { cols, vals }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CompiledSRow {
    plan: ScatterPlan,
    expected_q_lead: u32,
    out_elem: Option<(u32, Vec<bool>)>,
}

#[derive(Clone, Debug)]
pub(crate) struct CompiledStep {
    degree: u32,
    n_cols: usize,
    nq: usize,
    q_of: Vec<u32>,
    q_start: Vec<u32>,
    pivot_of: Vec<u32>,
    piv_leads: Vec<u32>,
    piv_order: Vec<u32>,
    pivs: Vec<ScatterPlan>,
    s_rows: Vec<CompiledSRow>,
    is_final: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct InitElemSpec {
    expected_lm: Mono,
    needed: bool,
    monos: Vec<Mono>,
}

#[derive(Clone, Debug)]
pub(crate) struct CompiledRefutation {
    n_elems: usize,
    init_elems: Vec<InitElemSpec>,
    steps: Vec<CompiledStep>,
}

pub(crate) struct UsefulRecord {
    pub(crate) mult: Mono,
    pub(crate) g: u32,
    pub(crate) q_lead_l: u32,
    pub(crate) used_q_leads: Vec<u32>,
    pub(crate) elem_id: Option<u32>,
}

pub(crate) struct RecordedStep {
    pub(crate) degree: u32,
    pub(crate) cols: Vec<Mono>,
    pub(crate) col_of: MonoMap<u32>,
    pub(crate) pivot_of: Vec<u32>,
    pub(crate) q_col: Vec<u32>,
    pub(crate) n_reducers: usize,
    pub(crate) piv: Vec<(Mono, u32, Mono)>,
    pub(crate) s_useful: Vec<UsefulRecord>,
}

fn build_scatter_plan(
    ring: &TowerRing,
    fp: Fp,
    w: TowerWeights,
    elem_id: u32,
    supp: &[Mono],
    mult: Mono,
    expected_lead_mono: Mono,
    c_col_of: &MonoMap<u32>,
) -> ScatterPlan {
    let mut unit_moves = Vec::with_capacity(supp.len());
    let mut scaled_moves = Vec::new();
    for (src_idx, &m_src) in supp.iter().enumerate() {
        let single = RPoly {
            monos: vec![m_src],
            coefs: vec![1],
        };
        let prod = ring.mul_mono_fp(fp, &single, mult, w);
        for (&m_dst, &coef) in prod.monos.iter().zip(&prod.coefs) {
            let dst_col = c_col_of[&m_dst];
            if coef == 1 {
                unit_moves.push((src_idx as u32, dst_col));
            } else if coef != 0 {
                scaled_moves.push((src_idx as u32, dst_col, coef));
            }
        }
    }
    ScatterPlan {
        elem_id,
        expected_lead: c_col_of[&expected_lead_mono],
        unit_moves,
        scaled_moves,
    }
}

pub(crate) fn compile_refutation(
    ring: &TowerRing,
    fp: Fp,
    w: TowerWeights,
    init_lms: &[Mono],
    init_monos: &[Vec<Mono>],
    elem_monos: &[Vec<Mono>],
    rec_steps: &[RecordedStep],
) -> Option<(CompiledRefutation, Vec<TraceStep>)> {
    if rec_steps.is_empty() {
        return None;
    }
    let n_steps = rec_steps.len();
    let total_elems = elem_monos.len();
    let mut need_elem = vec![false; total_elems];
    let mut step_kept_u: Vec<Vec<usize>> = vec![Vec::new(); n_steps];

    // Phase 1: Backward reachability from the single `ONE` row at `n_steps - 1`.
    for t in (0..n_steps).rev() {
        let rs = &rec_steps[t];
        let nu = rs.s_useful.len();
        let mut u_of_ql: HashMap<u32, usize> = HashMap::with_capacity(nu);
        for (u, ur) in rs.s_useful.iter().enumerate() {
            u_of_ql.insert(ur.q_lead_l, u);
        }
        let mut need_u = vec![false; nu];
        if t == n_steps - 1 {
            let u_one = (0..nu).find(|&u| {
                let ql = rs.s_useful[u].q_lead_l as usize;
                rs.cols[rs.q_col[ql] as usize] == ONE
            })?;
            need_u[u_one] = true;
        } else {
            for (u, ur) in rs.s_useful.iter().enumerate() {
                if ur.elem_id.is_some_and(|id| need_elem[id as usize]) {
                    need_u[u] = true;
                }
            }
        }
        for u in (0..nu).rev() {
            if !need_u[u] {
                continue;
            }
            for &ql in &rs.s_useful[u].used_q_leads {
                if let Some(&up) = u_of_ql.get(&ql) {
                    need_u[up] = true;
                }
            }
        }
        let kept_u: Vec<usize> = (0..nu).filter(|&u| need_u[u]).collect();
        if kept_u.is_empty() {
            continue;
        }
        // Mark basis dependencies of kept S-rows and transitively reachable reducers.
        let mut need_piv = vec![false; rs.piv.len()];
        let mut work: Vec<usize> = Vec::new();
        let mark_poly = |mult: Mono,
                         g: u32,
                         skip_lead: Option<Mono>,
                         need_piv: &mut Vec<bool>,
                         work: &mut Vec<usize>| {
            let dummy = RPoly {
                monos: elem_monos[g as usize].clone(),
                coefs: vec![1; elem_monos[g as usize].len()],
            };
            let prod = ring.mul_mono_fp(fp, &dummy, mult, w);
            for &m in &prod.monos {
                if skip_lead == Some(m) {
                    continue;
                }
                if let Some(&c) = rs.col_of.get(&m) {
                    let pj = rs.pivot_of[c as usize];
                    if pj != NONE && !need_piv[pj as usize] {
                        need_piv[pj as usize] = true;
                        work.push(pj as usize);
                    }
                }
            }
        };
        for &u in &kept_u {
            let ur = &rs.s_useful[u];
            need_elem[ur.g as usize] = true;
            mark_poly(ur.mult, ur.g, None, &mut need_piv, &mut work);
        }
        while let Some(pj) = work.pop() {
            let (pmult, pg, plead) = rs.piv[pj];
            need_elem[pg as usize] = true;
            mark_poly(pmult, pg, Some(plead), &mut need_piv, &mut work);
        }
        step_kept_u[t] = kept_u;
    }

    // Phase 2: Forward compilation over tight reachable supports.
    let mut cur_supp: Vec<Vec<Mono>> = vec![Vec::new(); total_elems];
    let mut init_elems = Vec::with_capacity(init_lms.len());
    for (g, &lm) in init_lms.iter().enumerate() {
        if need_elem[g] {
            cur_supp[g] = init_monos[g].clone();
        }
        init_elems.push(InitElemSpec {
            expected_lm: lm,
            needed: need_elem[g],
            monos: init_monos[g].clone(),
        });
    }

    let mut compiled_steps = Vec::new();
    let mut minimal_trace_steps = Vec::new();
    for t in 0..n_steps {
        let kept_u = &step_kept_u[t];
        if kept_u.is_empty() {
            continue;
        }
        let rs = &rec_steps[t];
        let is_final = t == n_steps - 1;
        let mut need_piv = vec![false; rs.piv.len()];
        let mut work: Vec<usize> = Vec::new();
        let mut touched: MonoMap<()> = MonoMap::default();
        if is_final {
            touched.insert(ONE, ());
        }

        let scan_into = |mult: Mono,
                         g: u32,
                         skip_lead: Option<Mono>,
                         need_piv: &mut Vec<bool>,
                         work: &mut Vec<usize>,
                         touched: &mut MonoMap<()>| {
            for &m_src in &cur_supp[g as usize] {
                let single = RPoly {
                    monos: vec![m_src],
                    coefs: vec![1],
                };
                let prod = ring.mul_mono_fp(fp, &single, mult, w);
                for &m in &prod.monos {
                    touched.insert(m, ());
                    if skip_lead == Some(m) {
                        continue;
                    }
                    if let Some(&c) = rs.col_of.get(&m) {
                        let pj = rs.pivot_of[c as usize];
                        if pj != NONE && !need_piv[pj as usize] {
                            need_piv[pj as usize] = true;
                            work.push(pj as usize);
                        }
                    }
                }
            }
        };

        for &u in kept_u {
            let ur = &rs.s_useful[u];
            let orig_lm = rs.cols[rs.q_col[ur.q_lead_l as usize] as usize];
            touched.insert(orig_lm, ());
            scan_into(ur.mult, ur.g, None, &mut need_piv, &mut work, &mut touched);
        }
        while let Some(pj) = work.pop() {
            let (pmult, pg, plead) = rs.piv[pj];
            touched.insert(plead, ());
            scan_into(
                pmult,
                pg,
                Some(plead),
                &mut need_piv,
                &mut work,
                &mut touched,
            );
        }

        let mut c_cols: Vec<Mono> = touched.keys().copied().collect();
        c_cols.sort_unstable_by(|a, b| b.cmp(a));
        let c_col_of: MonoMap<u32> = c_cols
            .iter()
            .enumerate()
            .map(|(i, &m)| (m, i as u32))
            .collect();
        let n_cols = c_cols.len();

        let kept_piv_idx: Vec<usize> = (0..rs.piv.len()).filter(|&j| need_piv[j]).collect();
        let mut pivot_of = vec![NONE; n_cols];
        let mut piv_leads = Vec::with_capacity(kept_piv_idx.len());
        let mut pivs = Vec::with_capacity(kept_piv_idx.len());
        for (new_k, &orig_j) in kept_piv_idx.iter().enumerate() {
            let (pmult, pg, plead) = rs.piv[orig_j];
            let plan = build_scatter_plan(
                ring,
                fp,
                w,
                pg,
                &cur_supp[pg as usize],
                pmult,
                plead,
                &c_col_of,
            );
            pivot_of[plan.expected_lead as usize] = new_k as u32;
            piv_leads.push(plan.expected_lead);
            pivs.push(plan);
        }
        let mut piv_order: Vec<u32> = (0..pivs.len() as u32).collect();
        piv_order.sort_unstable_by_key(|&k| std::cmp::Reverse(piv_leads[k as usize]));

        let mut q_of = vec![NONE; n_cols];
        let mut q_start = vec![0u32; n_cols];
        let mut c_q_monos: Vec<Mono> = Vec::new();
        for c in 0..n_cols {
            q_start[c] = c_q_monos.len() as u32;
            if pivot_of[c] == NONE {
                q_of[c] = c_q_monos.len() as u32;
                c_q_monos.push(c_cols[c]);
            }
        }
        let nq = c_q_monos.len();
        if is_final && c_q_monos.last() != Some(&ONE) {
            return None;
        }

        let mut s_rows = Vec::with_capacity(kept_u.len());
        let mut kept_set: std::collections::HashSet<(Mono, u32)> =
            std::collections::HashSet::with_capacity(kept_u.len() + kept_piv_idx.len());
        let mut fresh_lms = Vec::with_capacity(kept_u.len());
        for &orig_j in &kept_piv_idx {
            if orig_j >= rs.n_reducers {
                kept_set.insert((rs.piv[orig_j].0, rs.piv[orig_j].1));
            }
        }
        for &u in kept_u {
            let ur = &rs.s_useful[u];
            kept_set.insert((ur.mult, ur.g));
            let orig_lm = rs.cols[rs.q_col[ur.q_lead_l as usize] as usize];
            fresh_lms.push(orig_lm);
            let lead_col = *c_col_of.get(&orig_lm)?;
            let ql = q_of[lead_col as usize];
            if ql == NONE {
                return None;
            }
            // Note: `ur.mult * LM(g)` can differ from `orig_lm` (because S-rows cancel up to `orig_lm`),
            // so `expected_lead` of the unreduced product is meaningless; `expected_q_lead` is what matters!
            let plan = build_scatter_plan(
                ring,
                fp,
                w,
                ur.g,
                &cur_supp[ur.g as usize],
                ur.mult,
                orig_lm,
                &c_col_of,
            );
            let out_elem = if !is_final && ur.elem_id.is_some_and(|id| need_elem[id as usize]) {
                let id = ur.elem_id.unwrap();
                let target_supp = &elem_monos[id as usize];
                let q_tail = &c_q_monos[ql as usize..];
                if target_supp.first() != q_tail.first() {
                    return None;
                }
                let mut keep_mask = vec![false; q_tail.len()];
                let mut ts_idx = 0usize;
                for (q_rel, &qm) in q_tail.iter().enumerate() {
                    if ts_idx < target_supp.len() && qm == target_supp[ts_idx] {
                        keep_mask[q_rel] = true;
                        ts_idx += 1;
                    }
                }
                if ts_idx != target_supp.len() {
                    return None;
                }
                cur_supp[id as usize] = target_supp.clone();
                Some((id, keep_mask))
            } else {
                None
            };
            s_rows.push(CompiledSRow {
                plan,
                expected_q_lead: ql,
                out_elem,
            });
        }
        fresh_lms.sort_unstable();
        minimal_trace_steps.push(TraceStep {
            degree: rs.degree,
            kept: kept_set,
            fresh_lms,
        });
        compiled_steps.push(CompiledStep {
            degree: rs.degree,
            n_cols,
            nq,
            q_of,
            q_start,
            pivot_of,
            piv_leads,
            piv_order,
            pivs,
            s_rows,
            is_final,
        });
    }

    if compiled_steps.last().is_none_or(|st| !st.is_final) {
        return None;
    }
    Some((
        CompiledRefutation {
            n_elems: total_elems,
            init_elems,
            steps: compiled_steps,
        },
        minimal_trace_steps,
    ))
}

pub(crate) fn replay_compiled(
    start: &[RPoly],
    fp: Fp,
    opts: &TowerF4Options,
    dag: &CompiledRefutation,
    rep: &mut TowerF4Report,
) {
    if start.len() != dag.init_elems.len() {
        rep.replay_failure = Some("initial basis count differs");
        return;
    }
    let mut polys: Vec<Vec<u32>> = vec![Vec::new(); dag.n_elems];
    for (g, spec) in dag.init_elems.iter().enumerate() {
        if !spec.needed {
            continue;
        }
        if start[g].lm() != Some(spec.expected_lm) {
            rep.replay_failure = Some("initial lead differs");
            return;
        }
        let mut dense = vec![0u32; spec.monos.len()];
        for (&m, &c) in start[g].monos.iter().zip(&start[g].coefs) {
            // `spec.monos` is strictly descending (`b.cmp(&a)`).
            match spec.monos.binary_search_by(|probe| m.cmp(probe)) {
                Ok(idx) => dense[idx] = c,
                Err(_) => {
                    rep.replay_failure = Some("input support outside calibration trace");
                    return;
                }
            }
        }
        polys[g] = dense;
    }

    let mut max_cols_so_far = 0usize;
    let mut max_deg = 0u32;
    for (step_idx, st) in dag.steps.iter().enumerate() {
        if opts.deadline.is_some_and(|dl| Instant::now() >= dl) {
            rep.timed_out = true;
            rep.replay_failure = Some("deadline");
            return;
        }
        rep.steps = step_idx + 1;
        rep.degree_reached = rep.degree_reached.max(st.degree);
        max_deg = max_deg.max(st.degree);
        max_cols_so_far = max_cols_so_far.max(st.n_cols);
        rep.max_cols = rep.max_cols.max(st.n_cols);
        rep.max_rows = rep.max_rows.max(st.pivs.len() + st.s_rows.len());

        let kern = Kernel::with_budget(fp, 2 * st.n_cols + 4);

        // 1. Form kept reducer/promoted rows and batch-normalize any promoted rows whose lead != 1.
        let n_piv = st.pivs.len();
        let mut raw_pivs: Vec<Row> = Vec::with_capacity(n_piv);
        let mut to_inv_idx: Vec<usize> = Vec::new();
        let mut to_inv_vals: Vec<u64> = Vec::new();
        for (k, plan) in st.pivs.iter().enumerate() {
            let r = plan.eval_row(
                fp,
                &polys[plan.elem_id as usize],
                st.n_cols,
                &mut rep.muladds,
            );
            if r.cols.first().copied() != Some(plan.expected_lead) {
                rep.replay_failure = Some("reducer pivot vanished or shifted");
                return;
            }
            if r.vals[0] != 1 {
                to_inv_idx.push(k);
                to_inv_vals.push(u64::from(r.vals[0]));
            }
            raw_pivs.push(r);
        }
        if !to_inv_vals.is_empty() {
            let Some(invs) = fp.batch_inv(&to_inv_vals) else {
                rep.replay_failure = Some("promoted pivot zero");
                return;
            };
            rep.muladds += 3 * (to_inv_vals.len().saturating_sub(1) as u64);
            for (&k, &inv) in to_inv_idx.iter().zip(&invs) {
                raw_pivs[k].vals[0] = 1;
                for v in &mut raw_pivs[k].vals[1..] {
                    rep.muladds += 1;
                    *v = fp.mul(u64::from(*v), inv) as u32;
                }
            }
        }

        // Compute B' = A^{-1} B in pre-sorted rightmost-lead-first order.
        let mut bred: Vec<Vec<u32>> = vec![Vec::new(); n_piv];
        for &k_u32 in &st.piv_order {
            let k = k_u32 as usize;
            let s0 = st.q_start[st.piv_leads[k] as usize] as usize;
            let mut acc = vec![0u64; st.nq - s0];
            let row = &raw_pivs[k];
            for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                let q = st.q_of[c as usize];
                if q != NONE {
                    acc[q as usize - s0] = u64::from(v);
                }
            }
            for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                let cu = c as usize;
                if st.q_of[cu] != NONE {
                    continue;
                }
                let j = st.pivot_of[cu] as usize;
                let off = st.q_start[cu] as usize - s0;
                rep.muladds += kern.sub(&mut acc[off..], u64::from(v), &bred[j]);
            }
            bred[k] = acc.iter().map(|&x| fp.reduce(x) as u32).collect();
        }

        // 2. Reduce kept S-rows by B' and run minimal intra-step echelon.
        let mut ech: Vec<(u32, Vec<u32>)> = Vec::new();

        for (si, sr) in st.s_rows.iter().enumerate() {
            let row = sr.plan.eval_row(
                fp,
                &polys[sr.plan.elem_id as usize],
                st.n_cols,
                &mut rep.muladds,
            );
            let mut acc = vec![0u64; st.nq];
            for (&c, &v) in row.cols.iter().zip(&row.vals) {
                let q = st.q_of[c as usize];
                if q != NONE {
                    acc[q as usize] = u64::from(v);
                }
            }
            for (&c, &v) in row.cols.iter().zip(&row.vals) {
                let cu = c as usize;
                if st.q_of[cu] == NONE {
                    let off = st.q_start[cu] as usize;
                    let j = st.pivot_of[cu] as usize;
                    rep.muladds += kern.sub(&mut acc[off..], u64::from(v), &bred[j]);
                }
            }
            for (ql, erow) in &ech {
                let l = *ql as usize;
                let a = fp.reduce(acc[l]);
                if a != 0 {
                    rep.muladds += kern.sub(&mut acc[l..], a, erow);
                }
            }

            let exp_l = sr.expected_q_lead as usize;
            // Verify all columns strictly to the left of `exp_l` cancelled to 0,
            // and `acc[exp_l]` is non-zero!
            if acc[..exp_l].iter().any(|&x| fp.reduce(x) != 0) {
                rep.replay_failure = Some("new leads differ from the trace");
                rep.replay_mismatch = Some((rep.steps, si, st.s_rows.len()));
                return;
            }
            let lead_val = fp.reduce(acc[exp_l]);
            if lead_val == 0 {
                rep.replay_failure = Some("witness pivot vanished");
                rep.replay_mismatch = Some((rep.steps, si, st.s_rows.len()));
                return;
            }

            if st.is_final && si + 1 == st.s_rows.len() {
                // Exact algebraic soundness check: `exp_l == st.nq - 1` corresponds to `ONE`,
                // all columns `0 .. st.nq - 1` are verified 0, and `lead_val != 0`, so
                // `lead_val · 1 ∈ I(x_R)` with zero division needed!
                if exp_l + 1 != st.nq {
                    rep.replay_failure = Some("final row not at unit column");
                    return;
                }
                rep.inconsistent = true;
                rep.basis = vec![RPoly {
                    monos: vec![ONE],
                    coefs: vec![1],
                }];
                rep.solving_degree_max = max_deg;
                rep.max_cols_to_solution = max_cols_so_far;
                rep.steps_to_solution = rep.steps;
                rep.last_productive_degree = st.degree;
                return;
            }

            let inv = fp.inv(lead_val);
            let mut full_monic = Vec::with_capacity(st.nq - exp_l);
            full_monic.push(1u32);
            for &x in &acc[exp_l + 1..] {
                let rx = fp.reduce(x);
                if rx == 0 {
                    full_monic.push(0);
                } else {
                    rep.muladds += 1;
                    full_monic.push(fp.mul(rx, inv) as u32);
                }
            }
            if let Some((eid, ref keep_mask)) = sr.out_elem {
                let mut compact = Vec::with_capacity(keep_mask.len());
                for (q_rel, &val) in full_monic.iter().enumerate() {
                    if keep_mask[q_rel] {
                        compact.push(val);
                    } else if val != 0 {
                        rep.replay_failure = Some("element support widened");
                        return;
                    }
                }
                polys[eid as usize] = compact;
            }
            if si + 1 < st.s_rows.len() {
                let pos = ech.partition_point(|&(ql, _)| ql < sr.expected_q_lead);
                ech.insert(pos, (sr.expected_q_lead, full_monic));
            }
        }
    }
    rep.replay_failure = Some("trace exhausted");
}

impl TowerF4Options {
    pub fn new(max_degree: u32) -> Self {
        TowerF4Options {
            max_degree,
            deadline: None,
            stop_staircase: None,
            max_nnz: None,
            max_dense: None,
            on_step: None,
            on_stop: None,
            on_sig_step: None,
            full_rank_exit: true,
            record_trace: false,
            replay: None,
            weights: None,
        }
    }
    /// Run in the grading `w` (see [`TowerF4Options::weights`]).
    pub fn with_weights(mut self, w: TowerWeights) -> Self {
        self.weights = (!w.is_unit()).then_some(w);
        self
    }
    pub fn on_step(mut self, f: fn(usize, &StepTrace)) -> Self {
        self.on_step = Some(f);
        self
    }
    pub fn with_budget(mut self, budget: Duration) -> Self {
        self.deadline = Some(Instant::now() + budget);
        self
    }
    pub fn stopping_below(mut self, standard_monomials: usize) -> Self {
        self.stop_staircase = Some(standard_monomials);
        self
    }
    pub fn with_max_nnz(mut self, nnz: u64) -> Self {
        self.max_nnz = Some(nnz);
        self
    }
    pub fn with_max_dense(mut self, entries: u64) -> Self {
        self.max_dense = Some(entries);
        self
    }
}

/// One F4 step, for traces.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct StepTrace {
    pub degree: u32,
    pub critical_pairs: usize,
    pub tower_pairs: usize,
    pub s_rows: usize,
    pub reducer_rows: usize,
    /// S-rows moved into the reducers, one per lcm column.
    pub promoted_rows: usize,
    pub cols: usize,
    pub nnz: u64,
    /// S-rows the reducers alone do not reduce to zero, and the columns
    /// their residues live on (those without a divisor).
    pub residual_rows: usize,
    pub residual_cols: usize,
    /// How many S-rows had been reduced when the echelon of the residues
    /// reached a pivot on every residual column, if it did; and the S-rows
    /// never reduced because of it (`full_rank_exit`). Skipped rows are
    /// not counted in `residual_rows`.
    pub full_rank_after: Option<usize>,
    pub skipped_rows: usize,
    /// Entries of the reduced reducer block, dense over those columns.
    pub dense_entries: u64,
    pub fresh: usize,
    pub fresh_min_degree: u32,
    pub muladds: u64,
    pub basis_active: usize,
    /// Elements kept after the step, active or not (an inactive one stays
    /// while a pending pair uses it), and the entries they hold.
    pub basis_kept: usize,
    pub basis_entries: u64,
    pub pairs_left: usize,
    pub ms: f64,
    /// Where `ms` went: symbolic preprocessing, the reduced reducer block
    /// `B'`, the S-rows and their echelon (rows formed on demand count
    /// where they are used), and the basis update.
    pub ms_rows: f64,
    pub ms_reduce: f64,
    pub ms_echelon: f64,
    pub ms_update: f64,
}

#[derive(Clone, Debug, Default)]
pub struct TowerF4Report {
    /// The active basis when the run ended (`[1]` when inconsistent). It
    /// is a Gröbner basis only if the run ended with no pairs left.
    pub basis: Vec<RPoly>,
    pub inconsistent: bool,
    /// Highest step degree processed.
    pub degree_reached: u32,
    /// Highest step degree at which the basis gained an element, or at
    /// which `1` appeared: the solving degree.
    pub solving_degree_max: u32,
    /// Degree of the last such step.
    pub last_productive_degree: u32,
    /// Widest matrix, and steps, up to and including the last productive
    /// step.
    pub max_cols_to_solution: usize,
    pub steps_to_solution: usize,
    pub steps: usize,
    pub max_rows: usize,
    pub max_cols: usize,
    pub max_nnz: u64,
    pub max_residual_rows: usize,
    /// S-rows never reduced because their step's echelon was already full
    /// (`full_rank_exit`).
    pub rows_skipped: u64,
    pub max_dense_entries: u64,
    /// The unit: `F_p` multiply-adds in the eliminations.
    pub muladds: u64,
    pub critical_pairs_reduced: u64,
    pub tower_pairs_reduced: u64,
    pub pairs_product_skipped: u64,
    pub pairs_chain_skipped: u64,
    pub reducer_rows: u64,
    /// Pairs left above `max_degree` when the run stopped there.
    pub pairs_above_bound: usize,
    /// Set when [`TowerF4Options::stop_staircase`] ended the run.
    pub staircase_at_stop: Option<usize>,
    pub timed_out: bool,
    pub oversize: bool,
    pub ms: f64,
    pub trace: Vec<StepTrace>,
    /// The step a stop (`timed_out` or `oversize`) cut off, as far as it
    /// got. It is not in `trace`, and `on_step` never sees it. A stop in
    /// symbolic preprocessing leaves only the degree, the pairs and the
    /// nonzeros scanned; a stop in the elimination also leaves the matrix,
    /// its columns without a divisor (`residual_cols`) and `B'`'s planned
    /// entries (`dense_entries`). The basis fields are those before the
    /// step, and `pairs_left` counts the step's pairs again.
    pub stopped_step: Option<StepTrace>,
    /// With [`TowerF4Options::record_trace`] and a refutation: the trace.
    pub learned: Option<F4Trace>,
    /// With [`TowerF4Options::replay`]: why the replay stopped without a
    /// refutation (a replay that refuted leaves `None` and `inconsistent`).
    pub replay_failure: Option<&'static str>,
    /// For a lead mismatch: (step, new elements found, new elements in the trace).
    pub replay_mismatch: Option<(usize, usize, usize)>,
}

// ── Sparse rows and the two eliminations ─────────────────────────────

/// A sparse row: ascending column indices (column 0 is the largest
/// monomial) and non-zero coefficients.
#[derive(Clone, Debug, Default)]
pub(crate) struct Row {
    pub(crate) cols: Vec<u32>,
    pub(crate) vals: Vec<u32>,
}

pub(crate) const NONE: u32 = u32::MAX;

/// A pivot row prepared for [`Acc::sub_tail`]: the row, lead first and
/// columns ascending, and the words of the accumulator's bitset that its tail
/// has columns in, with the bits to set in each.
///
/// Setting the bits of a long row one entry at a time makes the entries that
/// share a word wait, each for the store of the one before it. The words a
/// row touches do not change, so the list is made once with the pivot and
/// read every time the pivot is subtracted.
pub(crate) struct Pivot {
    row: Row,
    /// `(word, bits)`, in the order of the tail's columns.
    runs: Vec<(u32, u64)>,
}

impl Pivot {
    pub(crate) fn new(row: Row) -> Self {
        let mut runs: Vec<(u32, u64)> = Vec::new();
        for &c in &row.cols[1..] {
            let (w, b) = (c / 64, 1u64 << (c % 64));
            match runs.last_mut() {
                Some(r) if r.0 == w => r.1 |= b,
                _ => runs.push((w, b)),
            }
        }
        Pivot { row, runs }
    }

    pub(crate) fn into_row(self) -> Row {
        self.row
    }
}

/// A dense accumulator with a bitset of the columns it may be non-zero
/// in, scanned in ascending order while rows are subtracted. `hi` is the
/// last bitset word anything was written to, so a scan stops there instead
/// of at the end of the matrix.
///
/// While the matrix is narrow enough that no entry can reach `2⁶⁴`, an entry
/// is held unreduced (`lazy`): a subtraction adds its product and leaves the
/// reduction to the read, which goes through [`take`](Acc::take). A column
/// left holding a multiple of `p` reads as zero there, as it did when every
/// subtraction reduced.
pub(crate) struct Acc {
    pub(crate) vals: Vec<u64>,
    pub(crate) bits: Vec<u64>,
    pub(crate) hi: usize,
    lazy: bool,
}

impl Acc {
    /// An accumulator for rows of `n_cols` columns over `fp`.
    ///
    /// A column starts below `p` and gains at most one product `(p − 1)²`
    /// from each pivot row subtracted, and a scan meets each column's pivot
    /// once, at a lead to the left of it: fewer than `n_cols` products. Their
    /// sum stays below `2⁶⁴` when `n_cols · (p − 1)² + p` does, and any wider
    /// matrix reduces at every subtraction.
    pub(crate) fn new(fp: Fp, n_cols: usize) -> Self {
        let worst = n_cols as u128 * u128::from(fp.p - 1).pow(2) + u128::from(fp.p);
        Acc {
            vals: vec![0; n_cols],
            bits: vec![0; n_cols.div_ceil(64).max(1)],
            hi: 0,
            lazy: worst < 1 << 64,
        }
    }

    #[inline]
    pub(crate) fn load(&mut self, row: &Row) {
        for (&c, &v) in row.cols.iter().zip(&row.vals) {
            self.vals[c as usize] = u64::from(v);
            self.bits[c as usize / 64] |= 1u64 << (c % 64);
        }
        self.hi = row.cols.last().map_or(0, |&c| c as usize / 64);
    }

    /// Column `c`, reduced, leaving it zero.
    #[inline(always)]
    pub(crate) fn take(&mut self, fp: Fp, c: usize) -> u64 {
        let v = std::mem::take(&mut self.vals[c]);
        if self.lazy {
            fp.reduce(v)
        } else {
            v
        }
    }

    /// Subtract `f · pivot` (without its lead, which the caller zeroes), for
    /// `f < p`. The columns written are to the right of the lead, so those
    /// in the lead's own word, which the caller holds in `word` with the
    /// accumulator's copy cleared, are collected from the accumulator once
    /// the row is done.
    #[inline]
    pub(crate) fn sub_tail(
        &mut self,
        fp: Fp,
        f: u64,
        piv: &Pivot,
        cur_word: usize,
        word: &mut u64,
    ) -> u64 {
        let nf = fp.neg(f);
        let (cols, vals) = (&piv.row.cols[1..], &piv.row.vals[1..]);
        if self.lazy {
            for (&c, &v) in cols.iter().zip(vals) {
                self.vals[c as usize] += nf * u64::from(v);
            }
        } else {
            for (&c, &v) in cols.iter().zip(vals) {
                let x = &mut self.vals[c as usize];
                *x = fp.reduce(*x + nf * u64::from(v));
            }
        }
        for &(w, b) in &piv.runs {
            self.bits[w as usize] |= b;
        }
        *word |= std::mem::take(&mut self.bits[cur_word]);
        if let Some(&c) = piv.row.cols.last() {
            self.hi = self.hi.max(c as usize / 64);
        }
        cols.len() as u64
    }

    /// Move every non-zero entry from word `w` on into `out`, scaled by
    /// `scale`, leaving the accumulator clean.
    pub(crate) fn drain_from(&mut self, fp: Fp, w: usize, scale: u64, out: &mut Row) {
        let mut w = w;
        while w <= self.hi {
            let mut word = self.bits[w];
            self.bits[w] = 0;
            while word != 0 {
                let b = word.trailing_zeros() as usize;
                word &= word - 1;
                let c = w * 64 + b;
                let v = self.take(fp, c);
                if v != 0 {
                    out.cols.push(c as u32);
                    out.vals.push(fp.mul(v, scale) as u32);
                }
            }
            w += 1;
        }
    }

    /// Leave the accumulator clean from word `w` on, discarding what it
    /// holds there.
    pub(crate) fn clear_from(&mut self, w: usize) {
        let mut w = w;
        while w <= self.hi {
            let mut word = std::mem::take(&mut self.bits[w]);
            while word != 0 {
                let b = word.trailing_zeros() as usize;
                word &= word - 1;
                self.vals[w * 64 + b] = 0;
            }
            w += 1;
        }
    }
}

/// Echelon form of `rows`: each row is reduced, at its lead, by the
/// pivots so far and kept, monic, if anything survives. Rows are
/// semi-reduced: once the lead is a new pivot column the rest is kept as it
/// stands. Returns the pivots, or `None` past the deadline.
pub(crate) fn echelon(
    fp: Fp,
    mut rows: Vec<Row>,
    n_cols: usize,
    muladds: &mut u64,
    deadline: Option<Instant>,
) -> Option<Vec<Row>> {
    rows.retain(|r| !r.cols.is_empty());
    rows.sort_by_key(|r| (r.cols[0], r.cols.len()));
    let mut pivot_of = vec![NONE; n_cols];
    let mut pivots: Vec<Pivot> = Vec::new();
    let mut acc = Acc::new(fp, n_cols);
    for (k, row) in rows.into_iter().enumerate() {
        if k % 256 == 0 && deadline.is_some_and(|d| Instant::now() >= d) {
            return None;
        }
        // A lead nobody holds: the row is a pivot as it stands.
        if pivot_of[row.cols[0] as usize] == NONE {
            let inv = fp.inv(u64::from(row.vals[0]));
            let vals = row
                .vals
                .iter()
                .map(|&v| fp.mul(u64::from(v), inv) as u32)
                .collect();
            pivot_of[row.cols[0] as usize] = pivots.len() as u32;
            pivots.push(Pivot::new(Row {
                cols: row.cols,
                vals,
            }));
            continue;
        }
        acc.load(&row);
        let mut w = row.cols[0] as usize / 64;
        let mut found: Option<Row> = None;
        'scan: while w <= acc.hi {
            let mut word = acc.bits[w];
            acc.bits[w] = 0;
            while word != 0 {
                let b = word.trailing_zeros() as usize;
                word &= word - 1;
                let c = w * 64 + b;
                let v = acc.take(fp, c);
                if v == 0 {
                    continue;
                }
                let p = pivot_of[c];
                if p != NONE {
                    *muladds += acc.sub_tail(fp, v, &pivots[p as usize], w, &mut word);
                    continue;
                }
                // A new lead: keep the rest of the row as it stands.
                let inv = fp.inv(v);
                let mut out = Row {
                    cols: vec![c as u32],
                    vals: vec![1],
                };
                acc.bits[w] = word;
                acc.drain_from(fp, w, inv, &mut out);
                found = Some(out);
                break 'scan;
            }
            w += 1;
        }
        if let Some(out) = found {
            pivot_of[out.cols[0] as usize] = pivots.len() as u32;
            pivots.push(Pivot::new(out));
        }
    }
    Some(pivots.into_iter().map(Pivot::into_row).collect())
}

/// `acc += a · b` with lazy reduction: every entry stays below `2⁶³`
/// because `a · b < 2⁶²` and `m`, a multiple of `p` at least `2⁶³ − p`, is
/// subtracted whenever an entry reaches `2⁶³`.
#[inline(always)]
fn lazy_axpy(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    let n = acc.len().min(b.len());
    for (x, &y) in acc[..n].iter_mut().zip(&b[..n]) {
        let s = *x + u64::from(a) * u64::from(y);
        *x = s - ((s >> 63).wrapping_neg() & m);
    }
}

/// `acc += a · b` with no reduction at all: the caller has bounded the
/// number of products an entry can take (see [`Kernel::with_budget`]) so
/// the sum stays below `2⁶⁴`, and reads every entry through
/// [`Fp::reduce`].
#[inline(always)]
fn exact_axpy(acc: &mut [u64], a: u32, b: &[u32]) {
    let n = acc.len().min(b.len());
    for (x, &y) in acc[..n].iter_mut().zip(&b[..n]) {
        *x += u64::from(a) * u64::from(y);
    }
}

/// [`exact_axpy`] compiled for AVX2.
///
/// # Safety
/// The CPU must support AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn exact_axpy_avx2(acc: &mut [u64], a: u32, b: &[u32]) {
    exact_axpy(acc, a, b)
}

/// [`exact_axpy`] compiled for AVX-512.
///
/// # Safety
/// The CPU must support AVX-512F.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn exact_axpy_avx512(acc: &mut [u64], a: u32, b: &[u32]) {
    exact_axpy(acc, a, b)
}

/// [`lazy_axpy`] compiled for AVX2.
///
/// # Safety
/// The CPU must support AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn lazy_axpy_avx2(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    lazy_axpy(acc, a, b, m)
}

/// [`lazy_axpy`] compiled for AVX-512.
///
/// # Safety
/// The CPU must support AVX-512F.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn lazy_axpy_avx512(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    lazy_axpy(acc, a, b, m)
}

/// The widest vector unit the lazy kernel can use on this CPU.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
enum Lanes {
    Scalar,
    Avx2,
    Avx512,
}

impl Lanes {
    fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512f") {
                return Lanes::Avx512;
            }
            if std::is_x86_feature_detected!("avx2") {
                return Lanes::Avx2;
            }
        }
        Lanes::Scalar
    }
}

/// `F_p` accumulation for the dense kernels. Below `2³¹` a product fits in
/// 62 bits, so an accumulator is kept below `2⁶³` by subtracting a multiple
/// of `p` only when it crosses `2⁶³`, and read through [`Fp::reduce`]
/// (lazy reduction); that loop vectorizes. Above, every step is reduced.
#[derive(Clone, Copy)]
struct Kernel {
    fp: Fp,
    /// The largest multiple of `p` at most `2⁶³`, when `p < 2³¹`.
    lazy: Option<u64>,
    lanes: Lanes,
    /// No entry can take enough products to reach `2⁶⁴`: nothing is
    /// subtracted at all and every read goes through [`Fp::reduce`].
    exact: bool,
}

impl Kernel {
    /// A kernel for accumulators that take at most `max_products` products
    /// each between two full reductions, all read through [`Fp::reduce`]:
    /// when `max_products · (p − 1)² + p < 2⁶⁴` the entries are only added
    /// to, with no conditional subtraction of a multiple of `p`.
    fn with_budget(fp: Fp, max_products: usize) -> Self {
        let mut k = Self::with_lanes(fp, Lanes::detect());
        let worst = max_products as u128 * u128::from(fp.p - 1).pow(2) + u128::from(fp.p);
        k.exact = worst < 1 << 64;
        k
    }

    fn with_lanes(fp: Fp, lanes: Lanes) -> Self {
        Kernel {
            fp,
            lazy: (fp.p < 1 << 31).then(|| ((1u64 << 63) / fp.p) * fp.p),
            lanes,
            exact: false,
        }
    }

    /// `acc += a · b` entrywise, for `a < p`. Every entry of `acc` must be
    /// below `2⁶³` (lazy) or below `p` (otherwise), and stays so.
    #[inline]
    fn axpy(self, acc: &mut [u64], a: u32, b: &[u32]) {
        if self.exact {
            match self.lanes {
                #[cfg(target_arch = "x86_64")]
                // SAFETY: `Lanes::detect` found AVX-512F on this CPU.
                Lanes::Avx512 => unsafe { exact_axpy_avx512(acc, a, b) },
                #[cfg(target_arch = "x86_64")]
                // SAFETY: `Lanes::detect` found AVX2 on this CPU.
                Lanes::Avx2 => unsafe { exact_axpy_avx2(acc, a, b) },
                _ => exact_axpy(acc, a, b),
            }
            return;
        }
        match (self.lazy, self.lanes) {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: `Lanes::detect` found AVX-512F on this CPU.
            (Some(m), Lanes::Avx512) => unsafe { lazy_axpy_avx512(acc, a, b, m) },
            #[cfg(target_arch = "x86_64")]
            // SAFETY: `Lanes::detect` found AVX2 on this CPU.
            (Some(m), Lanes::Avx2) => unsafe { lazy_axpy_avx2(acc, a, b, m) },
            (Some(m), _) => lazy_axpy(acc, a, b, m),
            (None, _) => {
                let fp = self.fp;
                for (x, &y) in acc.iter_mut().zip(b) {
                    *x = fp.reduce(*x + u64::from(a) * u64::from(y));
                }
            }
        }
    }

    /// `acc −= v · b`, the step every elimination takes.
    #[inline]
    fn sub(self, acc: &mut [u64], v: u64, b: &[u32]) -> u64 {
        self.axpy(acc, self.fp.neg(self.fp.reduce(v)) as u32, b);
        b.len() as u64
    }
}

/// Why an elimination stopped early.
pub(crate) enum Stop {
    Deadline,
    Oversize,
}

/// Counters of one elimination, for the trace.
#[derive(Default)]
pub(crate) struct ElimStats {
    /// S-rows that the reducers alone do not reduce to zero, among those
    /// reduced.
    pub(crate) residual_rows: usize,
    /// S-rows reduced when the echelon had a pivot on every `Q` column, if
    /// it did, and the S-rows then skipped.
    pub(crate) full_rank_after: Option<usize>,
    pub(crate) skipped_rows: usize,
    /// Columns no reducer leads: those without a divisor.
    pub(crate) q_cols: usize,
    /// Entries of the reduced reducer block `B'`.
    pub(crate) dense_entries: u64,
    pub(crate) muladds: u64,
    /// Time spent on `B'`, and on the S-rows and their echelon.
    pub(crate) ms_reducers: f64,
    pub(crate) ms_residues: f64,
    /// Input: form `B'` only for the reducers the S-rows reach, directly or
    /// through other reducers (a replay; the plain run forms all of it).
    pub(crate) prune_reducers: bool,
    /// Input: record intra-step echelon dependencies for minimal-1-witness trace compilation.
    pub(crate) record_deps: bool,
    /// Output: the S-rows (indices) that gave a new pivot, in order.
    pub(crate) useful: Vec<u32>,
    /// Output (when `record_deps`): `(q_lead_l, used_earlier_q_leads)` per entry in `useful`.
    pub(crate) useful_deps: Vec<(u32, Vec<u32>)>,
}

/// Forms row `k` of a set of rows on demand.
pub(crate) type FormRow<'a> = &'a (dyn Fn(usize) -> Row + Sync);

/// One step's matrix, eliminated in the manner of Faugère–Lachartre.
///
/// The pivot rows, formed by `piv_row`, have distinct, monic leads
/// `piv_leads` (`pivot_of[col]` indexes them): the `A|B` rows. The `n_s`
/// S-rows, formed by `s_row`, are the `C|D` rows. Rows are formed a block
/// at a time, when they are reduced. The columns no pivot row leads, `Q`,
/// are the ones without a divisor, and are handled densely:
/// 1. every pivot row is reduced to its lead plus a dense row over the `Q`
///    columns to its right, the rightmost lead first (`B' = A⁻¹B`), in
///    blocks whose rows are independent but for the few that point inside
///    the block;
/// 2. every S-row is reduced by those rows in one pass, leaving a dense row
///    over `Q` (`D − C·B'`), with no fill-in on the pivot columns;
/// 3. the residues are echelonized in chunks: each chunk is reduced in
///    parallel by the pivots found so far, then finished one row at a time.
///
/// Once the echelon has a pivot on every `Q` column, every later S-row
/// reduces to zero and no pivot changes (a pivot is never rewritten after
/// it is found). With `full_rank_exit` the rest of the S-rows are then
/// skipped: the rows of the current chunk not yet finished, and every
/// later chunk. The echelon is the same either way.
///
/// Returns the pivots of the echelon, monic and semi-reduced, as the dense
/// rows over `Q` that step 3 leaves. `B'` larger than `max_dense` entries
/// stops it as oversize.
#[allow(clippy::too_many_arguments)]
pub(crate) fn eliminate(
    fp: Fp,
    n_cols: usize,
    piv_leads: &[u32],
    piv_row: FormRow,
    pivot_of: &[u32],
    n_s: usize,
    s_row: FormRow,
    deadline: Option<Instant>,
    max_dense: Option<u64>,
    full_rank_exit: bool,
    st: &mut ElimStats,
) -> Result<Echelon, Stop> {
    // An entry of an accumulator takes at most one product per column of
    // its row at a pivot column (at most `n_cols`), and later one per
    // earlier lead (at most `nq ≤ n_cols`, the old and the chunk's new
    // leads together), the rest being reduced in between: `2 · n_cols`
    // products bound every entry between two reductions.
    let kern = Kernel::with_budget(fp, 2 * n_cols + 4);
    // Q, ascending; `q_start[c]` is the first Q index at or after `c`.
    let mut q_of = vec![NONE; n_cols];
    let mut q_start = vec![0u32; n_cols];
    let mut q_col: Vec<u32> = Vec::new();
    for c in 0..n_cols {
        q_start[c] = q_col.len() as u32;
        if pivot_of[c] == NONE {
            q_of[c] = q_col.len() as u32;
            q_col.push(c as u32);
        }
    }
    let nq = q_col.len();
    st.q_cols = nq;
    let tail_len = |c: u32| nq - q_start[c as usize] as usize;
    st.dense_entries = piv_leads.iter().map(|&c| tail_len(c) as u64).sum();
    if max_dense.is_some_and(|cap| st.dense_entries > cap) {
        return Err(Stop::Oversize);
    }
    let expired = AtomicBool::new(false);
    let past = || {
        if expired.load(AtomicOrdering::Relaxed) {
            return true;
        }
        if deadline.is_some_and(|d| Instant::now() >= d) {
            expired.store(true, AtomicOrdering::Relaxed);
            return true;
        }
        false
    };

    // 1. B' = A⁻¹B, the rightmost lead first. A pivot row's tail points
    //    only at pivot rows further right, so at rows already reduced, or
    //    at rows earlier in its own block, which wait for the sequential
    //    pass over the block.
    let started = Instant::now();
    let n_piv = piv_leads.len();
    let mut order: Vec<u32> = (0..n_piv as u32).collect();
    order.sort_unstable_by_key(|&k| std::cmp::Reverse(piv_leads[k as usize]));
    let mut pos = vec![0u32; n_piv];
    for (i, &k) in order.iter().enumerate() {
        pos[k as usize] = i as u32;
    }
    let mut bred: Vec<Vec<u32>> = vec![Vec::new(); n_piv];
    // A replay forms `B'` only for the reducers its S-rows reach.
    let needed: Option<Vec<bool>> = st.prune_reducers.then(|| {
        let mut need = vec![false; n_piv];
        let mut work: Vec<u32> = Vec::new();
        let mark = |row: &Row, from: usize, need: &mut Vec<bool>, work: &mut Vec<u32>| {
            for &c in &row.cols[from..] {
                if q_of[c as usize] == NONE {
                    let j = pivot_of[c as usize];
                    if !need[j as usize] {
                        need[j as usize] = true;
                        work.push(j);
                    }
                }
            }
        };
        for k in 0..n_s {
            mark(&s_row(k), 0, &mut need, &mut work);
        }
        while let Some(j) = work.pop() {
            mark(&piv_row(j as usize), 1, &mut need, &mut work);
        }
        need
    });
    const BLOCK: usize = 64;
    for (b, block) in order.chunks(BLOCK).enumerate() {
        let first = (b * BLOCK) as u32;
        type Partial = (Vec<u64>, Vec<(u32, u32)>, u64);
        let partial: Vec<Partial> = block
            .par_iter()
            .map(|&k| {
                let s0 = q_start[piv_leads[k as usize] as usize] as usize;
                let mut acc = vec![0u64; nq - s0];
                let mut later = Vec::new();
                let mut work = 0u64;
                if past() || needed.as_ref().is_some_and(|n| !n[k as usize]) {
                    return (Vec::new(), later, work);
                }
                let row = piv_row(k as usize);
                debug_assert_eq!(row.cols[0], piv_leads[k as usize]);
                debug_assert_eq!(row.vals[0], 1, "a pivot row must be monic");
                // The row's own Q entries first: the subtractions below
                // write to Q columns anywhere to their right.
                for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                    if q_of[c as usize] != NONE {
                        acc[q_of[c as usize] as usize - s0] = u64::from(v);
                    }
                }
                for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                    let cu = c as usize;
                    if q_of[cu] != NONE {
                        continue;
                    }
                    let j = pivot_of[cu];
                    if pos[j as usize] < first {
                        let off = q_start[cu] as usize - s0;
                        work += kern.sub(&mut acc[off..], u64::from(v), &bred[j as usize]);
                    } else {
                        later.push((j, v));
                    }
                }
                (acc, later, work)
            })
            .collect();
        if expired.load(AtomicOrdering::Relaxed) {
            return Err(Stop::Deadline);
        }
        for (&k, (mut acc, later, work)) in block.iter().zip(partial) {
            st.muladds += work;
            let s0 = q_start[piv_leads[k as usize] as usize] as usize;
            for (j, v) in later {
                let off = q_start[piv_leads[j as usize] as usize] as usize - s0;
                st.muladds += kern.sub(&mut acc[off..], u64::from(v), &bred[j as usize]);
            }
            bred[k as usize] = acc.iter().map(|&x| fp.reduce(x) as u32).collect();
        }
    }
    let residues_started = Instant::now();
    st.ms_reducers = (residues_started - started).as_secs_f64() * 1e3;

    // 2 and 3. The S-rows, in chunks. `ech[l]` is the pivot with lead `l`
    // (a Q index), stored from `l` on, with 1 at `l`; empty when there is
    // none. `leads` is ascending.
    let mut ech: Vec<Vec<u32>> = vec![Vec::new(); nq];
    let mut leads: Vec<u32> = Vec::new();
    // A task reduces a block of rows together, one reducer at a time, so
    // that each dense row it reads is read once per block rather than once
    // per row. The block's accumulators are sized to stay in cache.
    let per_task = ((1usize << 18) / (8 * nq.max(1))).clamp(4, 64);
    let chunk = per_task * 2 * rayon::current_num_threads().max(1);
    type Reduced = (Vec<(Option<Vec<u64>>, bool, Vec<u32>)>, u64);
    let all: Vec<usize> = (0..n_s).collect();
    // S-rows reduced so far. With no residual column the echelon is full
    // before the first row.
    let mut done = 0usize;
    let rec_deps = st.record_deps;
    if nq == 0 {
        st.full_rank_after = Some(0);
    }
    for rows in all.chunks(chunk) {
        if full_rank_exit && leads.len() == nq {
            // Added to the rows the previous chunk skipped after its echelon
            // filled.
            st.skipped_rows += n_s - done;
            break;
        }
        let reduced: Vec<Reduced> = rows
            .par_chunks(per_task)
            .map(|ks| {
                if past() {
                    return (Vec::new(), 0);
                }
                let blk: Vec<Row> = ks.iter().map(|&k| s_row(k)).collect();
                let mut acc = vec![0u64; blk.len() * nq];
                let mut work = 0u64;
                // (column, row, coefficient) of the pivot-column entries,
                // grouped by column: every row that meets a reducer takes
                // it while it is in cache.
                let mut ents: Vec<(u32, u32, u32)> = Vec::new();
                for (r, row) in blk.iter().enumerate() {
                    for (&c, &v) in row.cols.iter().zip(&row.vals) {
                        let q = q_of[c as usize];
                        if q != NONE {
                            acc[r * nq + q as usize] = u64::from(v);
                        } else {
                            ents.push((c, r as u32, v));
                        }
                    }
                }
                // (column, row) is unique, so the key orders the entries
                // exactly as the whole tuple did
                ents.sort_unstable_by_key(|&(c, r, _)| u64::from(c) << 32 | u64::from(r));
                for &(c, r, v) in &ents {
                    let (off, end) = (
                        r as usize * nq + q_start[c as usize] as usize,
                        (r as usize + 1) * nq,
                    );
                    work += kern.sub(
                        &mut acc[off..end],
                        u64::from(v),
                        &bred[pivot_of[c as usize] as usize],
                    );
                }
                let residual: Vec<bool> = (0..blk.len())
                    .map(|r| acc[r * nq..(r + 1) * nq].iter().any(|&x| fp.reduce(x) != 0))
                    .collect();
                let mut live: Vec<usize> = (0..blk.len()).filter(|&r| residual[r]).collect();
                let mut used_old: Vec<Vec<u32>> = if rec_deps {
                    vec![Vec::new(); blk.len()]
                } else {
                    Vec::new()
                };
                for &l in &leads {
                    let lu = l as usize;
                    for &r in &live {
                        let a = fp.reduce(acc[r * nq + lu]);
                        if a != 0 {
                            if rec_deps {
                                used_old[r].push(l);
                            }
                            work += kern.sub(&mut acc[r * nq + lu..(r + 1) * nq], a, &ech[lu]);
                        }
                    }
                }
                let mut out: Vec<(Option<Vec<u64>>, bool, Vec<u32>)> =
                    residual.iter().map(|&x| (None, x, Vec::new())).collect();
                live.retain(|&r| {
                    let row = &mut acc[r * nq..(r + 1) * nq];
                    for x in row.iter_mut() {
                        *x = fp.reduce(*x);
                    }
                    row.iter().any(|&x| x != 0)
                });
                for r in live {
                    out[r].0 = Some(acc[r * nq..(r + 1) * nq].to_vec());
                    if rec_deps {
                        out[r].2 = std::mem::take(&mut used_old[r]);
                    }
                }
                (out, work)
            })
            .collect();
        if expired.load(AtomicOrdering::Relaxed) {
            return Err(Stop::Deadline);
        }
        // What this chunk adds, reduced one row at a time by the pivots it
        // added before it (ascending); every residue already vanishes on
        // the older leads.
        let mut new_leads: Vec<u32> = Vec::new();
        let mut residues: Vec<(Option<Vec<u64>>, bool, Vec<u32>)> = Vec::with_capacity(rows.len());
        for (out, work) in reduced {
            st.muladds += work;
            residues.extend(out);
        }
        for (i, (acc, residual, mut used)) in residues.into_iter().enumerate() {
            if full_rank_exit && leads.len() + new_leads.len() == nq {
                // The echelon filled up earlier in this chunk: the rest of
                // it reduces to zero.
                st.skipped_rows += rows.len() - i;
                break;
            }
            st.residual_rows += usize::from(residual);
            let Some(mut acc) = acc else { continue };
            for &l in &new_leads {
                let lu = l as usize;
                let a = fp.reduce(acc[lu]);
                if a != 0 {
                    if rec_deps {
                        used.push(l);
                    }
                    st.muladds += kern.sub(&mut acc[lu..], a, &ech[lu]);
                }
            }
            let Some(l) = acc.iter().position(|&x| fp.reduce(x) != 0) else {
                continue;
            };
            let inv = fp.inv(fp.reduce(acc[l]));
            ech[l] = acc[l..]
                .iter()
                .map(|&x| fp.mul(fp.reduce(x), inv) as u32)
                .collect();
            let at = new_leads.partition_point(|&m| m < l as u32);
            new_leads.insert(at, l as u32);
            st.useful.push(rows[i] as u32);
            if rec_deps {
                st.useful_deps.push((l as u32, used));
            }
            if leads.len() + new_leads.len() == nq {
                st.full_rank_after = Some(done + i + 1);
            }
        }
        done += rows.len();
        leads.extend(new_leads);
        leads.sort_unstable();
    }
    st.ms_residues = residues_started.elapsed().as_secs_f64() * 1e3;

    Ok(Echelon {
        q_col,
        leads,
        rows: ech,
    })
}

/// The pivots an elimination finds, as it leaves them. `rows[l]`, for each
/// `l` of `leads` (ascending), is the pivot whose lead is the `Q` column `l`,
/// which is column `q_col[l]` of the step: monic, dense over the `Q` columns
/// from `l` on. The other `rows` are empty.
pub(crate) struct Echelon {
    pub(crate) q_col: Vec<u32>,
    pub(crate) leads: Vec<u32>,
    pub(crate) rows: Vec<Vec<u32>>,
}

impl Echelon {
    /// The pivots as sparse rows over the step's columns, leads ascending.
    #[cfg(test)]
    fn sparse_rows(&self) -> Vec<Row> {
        self.leads
            .iter()
            .map(|&l| {
                let mut row = Row::default();
                for (i, &v) in self.rows[l as usize].iter().enumerate() {
                    if v != 0 {
                        row.cols.push(self.q_col[l as usize + i]);
                        row.vals.push(v);
                    }
                }
                row
            })
            .collect()
    }
}

// ── Pairs and the basis ──────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
enum PairKind {
    /// The S-polynomial of basis elements `i < j`.
    Critical(u32, u32),
    /// `y_v · g_i` for a tower variable `v` of `LM(g_i)`.
    Tower(u32, u8),
}

#[derive(Clone, Copy, Debug)]
struct Pair {
    kind: PairKind,
    /// The lcm of the leading monomials (critical pairs).
    lcm: Mono,
    deg: u32,
}

struct State {
    polys: Vec<Elem>,
    lm: Vec<Mono>,
    active: Vec<bool>,
    /// Active leading monomial → element.
    lm_index: MonoMap<u32>,
    pairs: Vec<Pair>,
    /// The grading every key of this run is in.
    w: TowerWeights,
}

/// Every monomial `L'` with `lo | L' | L`, `L'` ≠ `L`, is passed to `f`
/// until it returns `true`; returns whether one did.
pub(crate) fn any_between(lo: Mono, l: Mono, f: impl FnMut(Mono) -> bool) -> bool {
    any_between_w(lo, l, TowerWeights::UNIT, f)
}

/// [`any_between`] with every monomial keyed in the grading `w`.
pub(crate) fn any_between_w(
    lo: Mono,
    l: Mono,
    w: TowerWeights,
    mut f: impl FnMut(Mono) -> bool,
) -> bool {
    let (lo_mask, l_mask) = (mask_of(lo), mask_of(l));
    let extra = l_mask & !lo_mask;
    let (lo_f, l_f) = (free_of(lo), free_of(l));
    // enumerate submasks s of `extra` and free vectors in [lo_f, l_f]
    let mut s = extra;
    loop {
        let mut fv = lo_f;
        loop {
            let cand = mono_w(lo_mask | s, fv, w);
            if cand != l && f(cand) {
                return true;
            }
            // next free vector (odometer)
            let mut k = 0;
            loop {
                if k == MAX_FREE {
                    break;
                }
                if fv[k] < l_f[k] {
                    fv[k] += 1;
                    break;
                }
                fv[k] = lo_f[k];
                k += 1;
            }
            if k == MAX_FREE {
                break;
            }
        }
        if s == 0 {
            break;
        }
        s = (s - 1) & extra;
    }
    false
}

impl State {
    /// Becker–Weispfenning `UPDATE` for the critical pairs of a new
    /// element, plus its tower pairs.
    fn insert(&mut self, h_poly: Elem, rep: &mut TowerF4Report) {
        let h = self.polys.len();
        let lh = h_poly.lm();
        self.polys.push(h_poly);
        self.lm.push(lh);
        self.active.push(false);

        let dh = degree_of(lh);
        let mut bits = mask_of(lh);
        while bits != 0 {
            let v = bits.trailing_zeros();
            bits &= bits - 1;
            self.pairs.push(Pair {
                kind: PairKind::Tower(h as u32, v as u8),
                lcm: lh,
                deg: dh + self.w.weight(v),
            });
        }

        // New pairs (g, h). Keep one per lcm that no other new lcm divides
        // strictly (criterion M), none for an lcm a coprime pair also has
        // (criterion F with the product criterion), and drop the coprime
        // ones: what `UPDATE`'s first two loops leave.
        // lcm → (the first element with it, whether any pair with it is
        // coprime). A group with a coprime pair is dropped whole, so its
        // representative never matters.
        let mut groups: MonoMap<(u32, bool)> = MonoMap::default();
        for g in 0..h {
            if !self.active[g] {
                continue;
            }
            let l = lcm_w(lh, self.lm[g], self.w);
            let e = groups.entry(l).or_insert((g as u32, false));
            if coprime(lh, self.lm[g]) {
                e.1 = true;
            }
        }
        let mut kept: Vec<(u32, Mono)> = Vec::new();
        for (&l, &(g, has_coprime)) in &groups {
            if any_between_w(lh, l, self.w, |c| groups.contains_key(&c)) {
                rep.pairs_chain_skipped += 1;
                continue;
            }
            if has_coprime {
                rep.pairs_product_skipped += 1;
                continue;
            }
            kept.push((g, l));
        }
        // Old critical pairs whose lcm `LM(h)` divides, strictly on both
        // sides (criterion B).
        let lm = &self.lm;
        let w = self.w;
        let mut dropped = 0u64;
        self.pairs.retain(|p| match p.kind {
            PairKind::Tower(..) => true,
            PairKind::Critical(i, j) => {
                let keep = !divides(lh, p.lcm)
                    || lcm_w(lm[i as usize], lh, w) == p.lcm
                    || lcm_w(lm[j as usize], lh, w) == p.lcm;
                if !keep {
                    dropped += 1;
                }
                keep
            }
        });
        rep.pairs_chain_skipped += dropped;
        kept.sort_unstable_by_key(|&(g, _)| g);
        for (g, l) in kept {
            self.pairs.push(Pair {
                kind: PairKind::Critical(g, h as u32),
                lcm: l,
                deg: degree_of(l),
            });
        }
        for g in 0..h {
            if self.active[g] && divides(lh, self.lm[g]) {
                self.active[g] = false;
                if self.lm_index.get(&self.lm[g]) == Some(&(g as u32)) {
                    self.lm_index.remove(&self.lm[g]);
                }
            }
        }
        self.active[h] = true;
        self.lm_index.insert(lh, h as u32);
    }

    /// An active element whose leading monomial divides `m`, the shortest.
    fn reducer_for(&self, m: Mono) -> Option<u32> {
        let mask = mask_of(m);
        let f = free_of(m);
        let mut best: Option<u32> = None;
        let mut s = mask;
        loop {
            let mut fv = [0u8; MAX_FREE];
            loop {
                if let Some(&g) = self.lm_index.get(&mono_w(s, fv, self.w)) {
                    if best
                        .is_none_or(|b| self.polys[g as usize].len() < self.polys[b as usize].len())
                    {
                        best = Some(g);
                    }
                }
                let mut k = 0;
                loop {
                    if k == MAX_FREE {
                        break;
                    }
                    if fv[k] < f[k] {
                        fv[k] += 1;
                        break;
                    }
                    fv[k] = 0;
                    k += 1;
                }
                if k == MAX_FREE {
                    break;
                }
            }
            if s == 0 {
                break;
            }
            s = (s - 1) & mask;
        }
        best
    }

    fn active_lms(&self) -> Vec<Mono> {
        (0..self.polys.len())
            .filter(|&g| self.active[g])
            .map(|g| self.lm[g])
            .collect()
    }
}

/// The number of standard monomials of `R / (lms)`, or `None` if it
/// exceeds `bound` or is infinite (a free variable without a pure power
/// among the leading monomials).
pub(crate) fn staircase_at_most(
    lms: &[Mono],
    n_tower: usize,
    n_free: usize,
    bound: usize,
) -> Option<usize> {
    let mut cap = [0u8; MAX_FREE];
    for (k, c) in cap.iter_mut().enumerate().take(n_free) {
        *c = lms
            .iter()
            .filter(|&&m| {
                mask_of(m) == 0
                    && free_of(m)
                        .iter()
                        .enumerate()
                        .all(|(i, &e)| (i == k) == (e > 0))
            })
            .map(|&m| free_of(m)[k])
            .min()?;
    }
    fn walk(
        k: usize,
        mask: u64,
        free: &mut [u8; MAX_FREE],
        ctx: (&[Mono], usize, usize, &[u8; MAX_FREE]),
        count: &mut usize,
        bound: usize,
    ) -> bool {
        let (lms, n_tower, n_free, cap) = ctx;
        // Divisibility on masks and exponents alone: the leads may be keyed
        // in a weighted grading, which `divides`'s degree test would mix
        // with this unit-keyed `cur`.
        let cur = mono(mask, *free);
        if lms.iter().any(|&l| divides_ungraded(l, cur)) {
            return true;
        }
        if k == n_tower + n_free {
            *count += 1;
            return *count <= bound;
        }
        if k < n_tower {
            walk(k + 1, mask, free, ctx, count, bound)
                && walk(k + 1, mask | 1 << k, free, ctx, count, bound)
        } else {
            let f = k - n_tower;
            for e in 0..cap[f] {
                free[f] = e;
                if !walk(k + 1, mask, free, ctx, count, bound) {
                    free[f] = 0;
                    return false;
                }
            }
            free[f] = 0;
            true
        }
    }
    let mut count = 0usize;
    let mut free = [0u8; MAX_FREE];
    walk(
        0,
        0,
        &mut free,
        (lms, n_tower, n_free, &cap),
        &mut count,
        bound,
    )
    .then_some(count)
}

/// A row of a step's matrix, kept as the product that forms it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RowDesc {
    pub(crate) mult: Mono,
    pub(crate) g: u32,
    /// The product's leading monomial and number of terms.
    pub(crate) lead: Mono,
    pub(crate) len: u32,
    /// Where the product itself is held in the step's [`ProductCache`], or
    /// [`NONE`] when it is not and must be formed again.
    pub(crate) at: u32,
}

/// The most terms a step keeps of the products its symbolic preprocessing
/// forms. A term costs eight bytes (a `u32` column number and a `u32`
/// coefficient), so the stored payload is at most 64 MiB however large the
/// step. This bounds `len`, not allocated memory: the two vectors grow by
/// doubling, so their capacity (and, while one reallocates, the old and new
/// buffers together) can exceed that, by up to about twice in the worst
/// case. That is still small beside the working set of the steps whose
/// products do not fit in memory together, which is why a row is formed on
/// demand at all. A product that would pass the limit is not kept and is
/// formed again when its row is eliminated.
pub(crate) const PRODUCT_CACHE_TERMS: usize = 1 << 23;

/// The products a step's scan has formed, one after another, so that the
/// elimination reads each row back instead of multiplying again. A row is
/// eliminated once, and the products are deterministic: a cached row is the
/// row the product would give.
///
/// A term is kept as the number its monomial was given when the scan first
/// met it (the value of `examined`), not as the monomial: the columns are
/// not known until the scan is over, and a number is turned into its column
/// by an array read where a monomial needs a hash lookup.
#[derive(Default)]
pub(crate) struct ProductCache {
    pub(crate) ids: Vec<u32>,
    pub(crate) coefs: Vec<u32>,
    /// Terms this cache may hold.
    limit: usize,
}

impl ProductCache {
    pub(crate) fn new(limit: usize) -> Self {
        ProductCache {
            // A position is a `u32` below `NONE`.
            limit: limit.min(NONE as usize - 1),
            ..Default::default()
        }
    }

    /// Whether a product of `terms` terms fits.
    pub(crate) fn has_room(&self, terms: usize) -> bool {
        self.ids.len() + terms <= self.limit
    }

    /// Keeps a product, which must fit: its position.
    pub(crate) fn keep(&mut self, ids: &[u32], coefs: &[u32]) -> u32 {
        debug_assert!(self.has_room(ids.len()) && ids.len() == coefs.len());
        let at = self.ids.len() as u32;
        self.ids.extend_from_slice(ids);
        self.coefs.extend_from_slice(coefs);
        at
    }
}

/// The number `examined` holds for `m`: a monomial is numbered when it is
/// first met, in order of arrival, and then joins `frontier`.
#[inline]
pub(crate) fn number_of(examined: &mut MonoMap<u32>, frontier: &mut Vec<Mono>, m: Mono) -> u32 {
    let next = examined.len() as u32;
    *examined.entry(m).or_insert_with(|| {
        frontier.push(m);
        next
    })
}

/// The columns of a matrix over the monomials of `examined`, numbered as
/// [`number_of`] does: the monomials, descending (column `i` is `cols[i]`),
/// and the column of each number.
pub(crate) fn order_columns(examined: &MonoMap<u32>) -> (Vec<Mono>, Vec<u32>) {
    let mut numbered: Vec<(Mono, u32)> = examined.iter().map(|(&m, &n)| (m, n)).collect();
    numbered.par_sort_unstable_by(|a, b| b.0.cmp(&a.0));
    let mut col_of_number = vec![0u32; numbered.len()];
    for (i, &(_, n)) in numbered.iter().enumerate() {
        col_of_number[n as usize] = i as u32;
    }
    (
        numbered.into_iter().map(|(m, _)| m).collect(),
        col_of_number,
    )
}

/// Forms the products `mult · polys[g]` of `todo` a batch at a time and
/// records each non-empty one in `out`. Every monomial of a product, past
/// its first `skip`, that `examined` has not seen joins `frontier` and is
/// numbered there, in order of arrival. A product is kept in `cache` while it
/// has room and dropped otherwise. Returns the number of terms formed.
///
/// The first `skip` monomials are ones `examined` has seen.
#[allow(clippy::too_many_arguments)]
fn scan_products(
    ring: &TowerRing,
    fp: Fp,
    w: TowerWeights,
    polys: &[Elem],
    todo: &[(Mono, u32)],
    skip: usize,
    out: &mut Vec<RowDesc>,
    examined: &mut MonoMap<u32>,
    frontier: &mut Vec<Mono>,
    cache: &mut ProductCache,
) -> u64 {
    let mut terms = 0u64;
    let mut ids: Vec<u32> = Vec::new();
    for batch in todo.chunks(1024) {
        let prods: Vec<RPoly> = batch
            .par_iter()
            .map(|&(mult, g)| ring.mul_elem(fp, &polys[g as usize], mult, w))
            .collect();
        for (&(mult, g), p) in batch.iter().zip(&prods) {
            let Some(lead) = p.lm() else { continue };
            terms += p.len() as u64;
            let keep = cache.has_room(p.len());
            ids.clear();
            if keep {
                ids.extend(p.monos.iter().take(skip).map(|m| examined[m]));
            }
            for &m in p.monos.iter().skip(skip) {
                let id = number_of(examined, frontier, m);
                if keep {
                    ids.push(id);
                }
            }
            out.push(RowDesc {
                mult,
                g,
                lead,
                len: p.len() as u32,
                at: if keep {
                    cache.keep(&ids, &p.coefs)
                } else {
                    NONE
                },
            });
        }
    }
    terms
}

// ── The algorithm ────────────────────────────────────────────────────

/// **F4 on `R`** for the ideal generated by `input` (tower normal forms),
/// with the normal selection strategy. See the module notes.
pub fn f4_tower(input: &[RPoly], ring: &TowerRing, opts: &TowerF4Options) -> TowerF4Report {
    f4_tower_cached(input, ring, opts, PRODUCT_CACHE_TERMS)
}

/// [`f4_tower`] with a step keeping at most `cache_terms` terms of its
/// products (see [`ProductCache`]). The result does not depend on it.
fn f4_tower_cached(
    input: &[RPoly],
    ring: &TowerRing,
    opts: &TowerF4Options,
    cache_terms: usize,
) -> TowerF4Report {
    ring.check();
    let w = opts.weights.unwrap_or(TowerWeights::UNIT);
    w.check(ring);
    let started = Instant::now();
    let fp = Fp::new(ring.p);
    let mut rep = TowerF4Report::default();
    let finish = |mut rep: TowerF4Report| {
        rep.ms = started.elapsed().as_secs_f64() * 1e3;
        rep
    };
    let unit = || {
        vec![RPoly {
            monos: vec![ONE],
            coefs: vec![1],
        }]
    };

    let mut inputs: Vec<RPoly> = input.iter().filter(|p| !p.is_empty()).cloned().collect();
    if !w.is_unit() {
        // Re-key every term in the grading. The map is injective on
        // (mask, free), so the terms stay distinct; only their order moves.
        for p in &mut inputs {
            let mut t: Vec<(Mono, u32)> = p
                .monos
                .iter()
                .zip(&p.coefs)
                .map(|(&m, &c)| (mono_w(mask_of(m), free_of(m), w), c))
                .collect();
            t.sort_unstable_by_key(|x| std::cmp::Reverse(x.0));
            p.monos = t.iter().map(|x| x.0).collect();
            p.coefs = t.iter().map(|x| x.1).collect();
        }
    }
    if inputs.iter().any(|p| p.is_unit()) {
        rep.inconsistent = true;
        rep.learned = opts.record_trace.then(F4Trace::default);
        rep.basis = unit();
        return finish(rep);
    }

    // The initial echelon: distinct leading monomials to start from.
    let (start, _) = {
        let mut all: Vec<Mono> = inputs
            .iter()
            .flat_map(|p| p.monos.iter().copied())
            .collect();
        all.sort_unstable_by(|a, b| b.cmp(a));
        all.dedup();
        let col_of: MonoMap<u32> = all
            .iter()
            .enumerate()
            .map(|(i, &m)| (m, i as u32))
            .collect();
        let rows: Vec<Row> = inputs
            .iter()
            .map(|p| Row {
                cols: p.monos.iter().map(|m| col_of[m]).collect(),
                vals: p.coefs.clone(),
            })
            .collect();
        match echelon(fp, rows, all.len(), &mut rep.muladds, opts.deadline) {
            Some(piv) => (
                piv.into_iter()
                    .map(|r| RPoly {
                        monos: r.cols.iter().map(|&c| all[c as usize]).collect(),
                        coefs: r.vals,
                    })
                    .collect::<Vec<_>>(),
                (),
            ),
            None => {
                rep.timed_out = true;
                rep.basis = inputs;
                return finish(rep);
            }
        }
    };
    if start.iter().any(|p| p.is_unit()) {
        rep.inconsistent = true;
        rep.learned = opts.record_trace.then(F4Trace::default);
        rep.basis = unit();
        return finish(rep);
    }
    let mut start = start;
    start.sort_by_key(|p| std::cmp::Reverse(p.lm().unwrap()));

    if let Some(trace) = &opts.replay {
        if let Some(dag) = &trace.compiled {
            replay_compiled(&start, fp, opts, dag, &mut rep);
            return finish(rep);
        }
    }

    let mut s = State {
        polys: Vec::new(),
        lm: Vec::new(),
        active: Vec::new(),
        lm_index: MonoMap::default(),
        pairs: Vec::new(),
        w,
    };
    let mut init_lms: Vec<Mono> = Vec::new();
    let mut init_monos: Vec<Vec<Mono>> = Vec::new();
    let mut elem_monos: Vec<Vec<Mono>> = Vec::new();
    if opts.record_trace {
        for p in &start {
            init_lms.push(p.lm().unwrap());
            init_monos.push(p.monos.clone());
            elem_monos.push(p.monos.clone());
        }
    }
    for p in start {
        s.insert(Elem::from_rpoly(p), &mut rep);
    }
    // (solving_degree_max, max_cols_to_solution, steps_to_solution)
    let mut learned = (0u32, 0usize, 0usize);
    let mut trace_steps: Vec<TraceStep> = Vec::new();
    let mut rec_steps: Vec<RecordedStep> = Vec::new();
    let mut max_cols_so_far = 0usize;

    let active_basis = |s: &State| -> Vec<RPoly> {
        (0..s.polys.len())
            .filter(|&g| s.active[g])
            .map(|g| s.polys[g].to_rpoly())
            .collect()
    };
    let kept = |s: &State, tr: &mut StepTrace| {
        tr.basis_active = s.active.iter().filter(|&&a| a).count();
        tr.basis_kept = s.polys.len();
        tr.basis_entries = s.polys.iter().map(|p| p.vals.len() as u64).sum();
        tr.pairs_left = s.pairs.len();
    };
    let record = |rep: &mut TowerF4Report, learned: (u32, usize, usize)| {
        rep.solving_degree_max = learned.0;
        rep.max_cols_to_solution = learned.1;
        rep.steps_to_solution = learned.2;
    };

    while !s.pairs.is_empty() {
        if opts.deadline.is_some_and(|d| Instant::now() >= d) {
            rep.timed_out = true;
            break;
        }
        let d = s.pairs.iter().map(|p| p.deg).min().unwrap();
        if d > opts.max_degree {
            rep.pairs_above_bound = s.pairs.len();
            break;
        }
        let step_started = Instant::now();
        let (selected, rest): (Vec<Pair>, Vec<Pair>) = s.pairs.drain(..).partition(|p| p.deg == d);
        s.pairs = rest;
        rep.steps += 1;
        rep.degree_reached = rep.degree_reached.max(d);
        let mut tr = StepTrace {
            degree: d,
            ..Default::default()
        };

        // The rows the step reduces: both halves of every critical pair and
        // every tower product, each (multiplier, element) once.
        let mut seen: HashMap<(Mono, u32), (), BuildHasherDefault<MonoHasher>> = HashMap::default();
        let mut jobs: Vec<(Mono, u32)> = Vec::new();
        let mut lcm_cols: MonoMap<()> = MonoMap::default();
        for p in &selected {
            match p.kind {
                PairKind::Critical(i, j) => {
                    tr.critical_pairs += 1;
                    lcm_cols.insert(p.lcm, ());
                    for g in [i, j] {
                        let mult = quotient(p.lcm, s.lm[g as usize]);
                        if seen.insert((mult, g), ()).is_none() {
                            jobs.push((mult, g));
                        }
                    }
                }
                PairKind::Tower(g, v) => {
                    tr.tower_pairs += 1;
                    // A tower multiplier lies inside LM(g), a half's outside it.
                    let mult = mono_w(1u64 << v, [0; MAX_FREE], s.w);
                    if seen.insert((mult, g), ()).is_none() {
                        jobs.push((mult, g));
                    }
                }
            }
        }
        // Replay: keep only the rows the trace kept at this step. The lcm
        // columns stay as they are: a column of `lcm_cols` is never searched
        // for a reducer, and the row promoted there (kept) is what leads it;
        // dropping a column whose promoted row is a tower product would send
        // it to the reducer search and change the matrix.
        let rstep: Option<&TraceStep> = match &opts.replay {
            None => None,
            Some(trace) => match trace.steps.get(rep.steps - 1) {
                None => {
                    rep.replay_failure = Some("trace exhausted");
                    break;
                }
                Some(ts) if ts.degree != d => {
                    rep.replay_failure = Some("step degree differs");
                    break;
                }
                Some(ts) => Some(ts),
            },
        };
        if let Some(ts) = rstep {
            jobs.retain(|j| ts.kept.contains(j));
        }
        rep.critical_pairs_reduced += tr.critical_pairs as u64;
        rep.tower_pairs_reduced += tr.tower_pairs as u64;

        // Symbolic preprocessing, one frontier at a time. A row is kept as
        // the product that makes it and formed again when it is eliminated:
        // the products of a large step do not fit in memory together, their
        // monomials do. The first `cache_terms` terms of products are kept
        // as well, and the elimination reads those back.
        // Every monomial the matrix will have a column for, each with the
        // number it was first met under: the lcm columns first.
        let mut examined: MonoMap<u32> = lcm_cols
            .keys()
            .enumerate()
            .map(|(i, &m)| (m, i as u32))
            .collect();
        let mut no_divisor: MonoMap<()> = MonoMap::default();
        let mut frontier: Vec<Mono> = Vec::new();
        let mut s_desc: Vec<RowDesc> = Vec::with_capacity(jobs.len());
        let mut cache = ProductCache::new(cache_terms);
        let mut nnz = scan_products(
            ring,
            fp,
            s.w,
            &s.polys,
            &jobs,
            0,
            &mut s_desc,
            &mut examined,
            &mut frontier,
            &mut cache,
        );
        let mut r_desc: Vec<RowDesc> = Vec::new();
        let mut aborted = false;
        while !frontier.is_empty() {
            if opts.deadline.is_some_and(|dl| Instant::now() >= dl) {
                aborted = true;
                break;
            }
            let found: Vec<(Mono, Option<u32>)> = frontier
                .par_iter()
                .map(|&m| (m, s.reducer_for(m)))
                .collect();
            frontier.clear();
            let mut todo: Vec<(Mono, u32)> = Vec::new();
            for &(m, g) in &found {
                match g {
                    Some(g) => todo.push((quotient(m, s.lm[g as usize]), g)),
                    None => {
                        no_divisor.insert(m, ());
                    }
                }
            }
            let before = r_desc.len();
            nnz += scan_products(
                ring,
                fp,
                s.w,
                &s.polys,
                &todo,
                1,
                &mut r_desc,
                &mut examined,
                &mut frontier,
                &mut cache,
            );
            debug_assert!(
                r_desc[before..]
                    .iter()
                    .zip(found.iter().filter(|f| f.1.is_some()))
                    .all(|(d, f)| d.lead == f.0),
                "a reducer must lead with its monomial"
            );
            if opts.max_nnz.is_some_and(|cap| nnz > cap) {
                rep.oversize = true;
                aborted = true;
                break;
            }
        }
        if aborted {
            if !rep.oversize {
                rep.timed_out = true;
            }
            s.pairs.extend(selected);
            tr.nnz = nnz;
            kept(&s, &mut tr);
            tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
            if let Some(f) = opts.on_stop {
                f(rep.steps, rep.oversize, &tr);
            }
            rep.stopped_step = Some(tr);
            break;
        }
        rep.reducer_rows += r_desc.len() as u64;

        // Columns, descending.
        let (cols, col_of_id) = order_columns(&examined);
        // The map of the scan becomes the column map: each monomial's number
        // is replaced by its column. The numbers kept in the cache are turned
        // into columns by `col_of_id`.
        let mut col_of = examined;
        for c in col_of.values_mut() {
            *c = col_of_id[*c as usize];
        }
        let n_cols = cols.len();
        let mut pivot_of = vec![NONE; n_cols];
        let mut piv: Vec<RowDesc> = r_desc;
        for (k, d) in piv.iter().enumerate() {
            pivot_of[col_of[&d.lead] as usize] = k as u32;
        }
        let n_reducers = piv.len();
        // The shortest S-row at each lcm column joins the reducers, made
        // monic, and the other rows there are reduced by it: the A|B and
        // C|D split. The residues then live on the columns without a
        // divisor alone, where every pivot is a new leading monomial. (The
        // sort is stable: ties keep the order of the jobs.)
        s_desc.sort_by_key(|d| (std::cmp::Reverse(d.lead), d.len));
        let mut s_list: Vec<RowDesc> = Vec::with_capacity(s_desc.len());
        for d in s_desc {
            let c = col_of[&d.lead] as usize;
            if pivot_of[c] == NONE && lcm_cols.contains_key(&d.lead) {
                pivot_of[c] = piv.len() as u32;
                piv.push(d);
            } else {
                s_list.push(d);
            }
        }
        tr.promoted_rows = piv.len() - n_reducers;
        tr.s_rows = s_list.len();
        tr.reducer_rows = n_reducers;
        tr.cols = n_cols;
        tr.nnz = nnz;
        rep.max_rows = rep.max_rows.max(s_list.len() + piv.len());
        rep.max_cols = rep.max_cols.max(n_cols);
        rep.max_nnz = rep.max_nnz.max(nnz);
        max_cols_so_far = max_cols_so_far.max(n_cols);

        tr.ms_rows = step_started.elapsed().as_secs_f64() * 1e3;

        // The elimination; the columns without a pivot row are exactly the
        // ones without a divisor, so every pivot it returns is new.
        debug_assert!(
            (0..n_cols).all(|c| (pivot_of[c] == NONE) == no_divisor.contains_key(&cols[c])),
            "the columns without a reducer must be those without a divisor"
        );
        let polys = &s.polys;
        let form = |d: &RowDesc| -> Row {
            if d.at != NONE {
                let at = d.at as usize..d.at as usize + d.len as usize;
                return Row {
                    cols: cache.ids[at.clone()]
                        .iter()
                        .map(|&id| col_of_id[id as usize])
                        .collect(),
                    vals: cache.coefs[at].to_vec(),
                };
            }
            let p = ring.mul_elem(fp, &polys[d.g as usize], d.mult, s.w);
            debug_assert_eq!(
                p.len(),
                d.len as usize,
                "a row must be formed as it was scanned"
            );
            Row {
                cols: p.monos.iter().map(|m| col_of[m]).collect(),
                vals: p.coefs,
            }
        };
        let piv_row = |k: usize| -> Row {
            let mut r = form(&piv[k]);
            if r.vals[0] != 1 {
                // Only a promoted tower row can lead with another
                // coefficient: products of monic elements stay monic.
                debug_assert!(k >= n_reducers);
                let inv = fp.inv(u64::from(r.vals[0]));
                for v in &mut r.vals {
                    *v = fp.mul(u64::from(*v), inv) as u32;
                }
            }
            r
        };
        let s_row = |k: usize| form(&s_list[k]);
        let piv_leads: Vec<u32> = piv.iter().map(|d| col_of[&d.lead]).collect();
        let mut st = ElimStats {
            prune_reducers: rstep.is_some(),
            record_deps: opts.record_trace,
            ..Default::default()
        };
        let result = eliminate(
            fp,
            n_cols,
            &piv_leads,
            &piv_row,
            &pivot_of,
            s_list.len(),
            &s_row,
            opts.deadline,
            opts.max_dense.or(opts.max_nnz),
            opts.full_rank_exit,
            &mut st,
        );
        // Nothing reads the cached products past the elimination.
        drop(cache);
        drop(col_of_id);
        tr.muladds = st.muladds;
        rep.muladds += st.muladds;
        tr.residual_rows = st.residual_rows;
        tr.full_rank_after = st.full_rank_after;
        tr.skipped_rows = st.skipped_rows;
        rep.rows_skipped += st.skipped_rows as u64;
        tr.residual_cols = st.q_cols;
        tr.dense_entries = st.dense_entries;
        tr.ms_reduce = st.ms_reducers;
        tr.ms_echelon = st.ms_residues;
        rep.max_residual_rows = rep.max_residual_rows.max(st.residual_rows);
        rep.max_dense_entries = rep.max_dense_entries.max(st.dense_entries);
        let ech = match result {
            Ok(ech) => ech,
            Err(stop) => {
                match stop {
                    Stop::Deadline => rep.timed_out = true,
                    Stop::Oversize => rep.oversize = true,
                }
                s.pairs.extend(selected);
                kept(&s, &mut tr);
                tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
                if let Some(f) = opts.on_stop {
                    f(rep.steps, rep.oversize, &tr);
                }
                rep.stopped_step = Some(tr);
                break;
            }
        };
        let update_started = Instant::now();

        // The new elements are the echelon's rows, over the step's columns
        // without a divisor.
        let Echelon {
            q_col,
            leads,
            rows: mut ech_rows,
        } = ech;
        let q_monos: Arc<[Mono]> = q_col.iter().map(|&c| cols[c as usize]).collect();
        let mut fresh: Vec<Elem> = leads
            .iter()
            .map(|&l| {
                Elem::new(
                    q_monos.clone(),
                    l,
                    std::mem::take(&mut ech_rows[l as usize]),
                )
            })
            .collect();
        drop(ech_rows);
        tr.fresh = fresh.len();
        tr.fresh_min_degree = fresh.iter().map(|p| degree_of(p.lm())).min().unwrap_or(0);
        if !fresh.is_empty() {
            learned = (learned.0.max(d), max_cols_so_far, rep.steps);
            rep.last_productive_degree = d;
        }
        let refuted = fresh.iter().any(|p| p.lm() == ONE);
        let mut fresh_lms: Vec<Mono> = Vec::new();
        if opts.record_trace || rstep.is_some() {
            fresh_lms = fresh.iter().map(|p| p.lm()).collect();
            fresh_lms.sort_unstable();
        }
        if opts.record_trace {
            let mut kept_rows_set: std::collections::HashSet<(Mono, u32)> = st
                .useful
                .iter()
                .map(|&k| (s_list[k as usize].mult, s_list[k as usize].g))
                .collect();
            kept_rows_set.extend(piv[n_reducers..].iter().map(|d| (d.mult, d.g)));
            trace_steps.push(TraceStep {
                degree: d,
                kept: kept_rows_set,
                fresh_lms: fresh_lms.clone(),
            });

            let base_elem_id = s.polys.len() as u32;
            let s_useful: Vec<UsefulRecord> = st
                .useful
                .iter()
                .zip(st.useful_deps)
                .map(|(&s_idx, (ql, used_ql))| {
                    let elem_id = if !refuted {
                        let rank = leads.binary_search(&ql).expect("useful lead in leads");
                        Some(base_elem_id + rank as u32)
                    } else {
                        None
                    };
                    UsefulRecord {
                        mult: s_list[s_idx as usize].mult,
                        g: s_list[s_idx as usize].g,
                        q_lead_l: ql,
                        used_q_leads: used_ql,
                        elem_id,
                    }
                })
                .collect();
            if !refuted {
                for f in &fresh {
                    elem_monos.push(f.to_rpoly().monos);
                }
            }
            rec_steps.push(RecordedStep {
                degree: d,
                cols,
                col_of,
                pivot_of,
                q_col,
                n_reducers,
                piv: piv.iter().map(|d| (d.mult, d.g, d.lead)).collect(),
                s_useful,
            });
        }
        if let Some(ts) = rstep {
            // A refutation stands on its own; anything else must follow the
            // trace, or the replay stops.
            if !refuted && fresh_lms != ts.fresh_lms {
                rep.replay_failure = Some("new leads differ from the trace");
                rep.replay_mismatch = Some((rep.steps, fresh_lms.len(), ts.fresh_lms.len()));
                break;
            }
        }
        if refuted {
            if opts.record_trace {
                // Preserve the complete step schedule for independent file
                // replay. The upstream optimized trace prunes zero steps and
                // requires a private compiled DAG that cannot be serialized.
                rep.learned = Some(F4Trace {
                    steps: trace_steps,
                    compiled: None,
                });
            }
            rep.inconsistent = true;
            rep.basis = unit();
            record(&mut rep, learned);
            kept(&s, &mut tr);
            tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
            if let Some(f) = opts.on_step {
                f(rep.steps, &tr);
            }
            rep.trace.push(tr);
            return finish(rep);
        }
        // Largest leading monomial first, so that a later, smaller one
        // retires the earlier ones it divides (`UPDATE`'s G_new).
        fresh.sort_by_key(|p| std::cmp::Reverse(p.lm()));
        for p in fresh {
            debug_assert_eq!(p.vals[0], 1, "an echelon row is monic");
            s.insert(p, &mut rep);
        }
        kept(&s, &mut tr);
        tr.ms_update = update_started.elapsed().as_secs_f64() * 1e3;
        tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
        if let Some(f) = opts.on_step {
            f(rep.steps, &tr);
        }
        rep.trace.push(tr);

        if let Some(bound) = opts.stop_staircase {
            if !s.pairs.is_empty() {
                if let Some(count) =
                    staircase_at_most(&s.active_lms(), ring.n_tower(), ring.n_free, bound)
                {
                    rep.staircase_at_stop = Some(count);
                    break;
                }
            }
        }
    }
    record(&mut rep, learned);
    rep.basis = active_basis(&s);
    finish(rep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    fn fp(p: u64) -> Fp {
        Fp::new(p)
    }

    /// grevlex on exponent vectors, variable 0 largest.
    fn grevlex_cmp(a: &[u32], b: &[u32]) -> std::cmp::Ordering {
        let (da, db): (u32, u32) = (a.iter().sum(), b.iter().sum());
        if da != db {
            return da.cmp(&db);
        }
        for i in (0..a.len()).rev() {
            if a[i] != b[i] {
                return b[i].cmp(&a[i]);
            }
        }
        std::cmp::Ordering::Equal
    }

    fn key_of(e: &[u32], n_tower: usize) -> Mono {
        let mut mask = 0u64;
        for (j, &x) in e.iter().take(n_tower).enumerate() {
            assert!(x <= 1);
            if x == 1 {
                mask |= 1 << j;
            }
        }
        let mut f = [0u8; MAX_FREE];
        for (k, &x) in e.iter().skip(n_tower).enumerate() {
            f[k] = x as u8;
        }
        mono(mask, f)
    }

    #[test]
    fn keys_order_like_grevlex() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let (n_tower, n_free) = (9usize, 2usize);
        for _ in 0..20_000 {
            let mut e1 = vec![0u32; n_tower + n_free];
            let mut e2 = vec![0u32; n_tower + n_free];
            for e in [&mut e1, &mut e2] {
                for (i, x) in e.iter_mut().enumerate() {
                    *x = if i < n_tower {
                        rng.gen_range(0..2)
                    } else {
                        rng.gen_range(0..4)
                    };
                }
            }
            let (k1, k2) = (key_of(&e1, n_tower), key_of(&e2, n_tower));
            assert_eq!(k1.cmp(&k2), grevlex_cmp(&e1, &e2), "{e1:?} {e2:?}");
            assert_eq!(degree_of(k1), e1.iter().sum::<u32>());
            // lcm, divisibility and quotient agree with exponent vectors
            let l: Vec<u32> = e1.iter().zip(&e2).map(|(a, b)| *a.max(b)).collect();
            assert_eq!(lcm(k1, k2), key_of(&l, n_tower));
            let div = e1.iter().zip(&e2).all(|(a, b)| a <= b);
            assert_eq!(divides(k1, k2), div);
            if div {
                let q: Vec<u32> = e2.iter().zip(&e1).map(|(b, a)| b - a).collect();
                assert_eq!(quotient(k2, k1), key_of(&q, n_tower));
            }
        }
        assert_eq!(mono(0, [0; MAX_FREE]), ONE);
    }

    /// Every kernel this CPU can run against exact arithmetic, over enough
    /// steps that the lazy accumulators cross `2⁶³` many times.
    #[test]
    fn the_kernels_agree_with_exact_arithmetic() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        #[cfg_attr(not(target_arch = "x86_64"), allow(unused_mut))]
        let mut lanes = vec![Lanes::Scalar];
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                lanes.push(Lanes::Avx2);
            }
            if std::is_x86_feature_detected!("avx512f") {
                lanes.push(Lanes::Avx512);
            }
        }
        for &p in &[
            3u64,
            786433,
            2013265921,
            (1 << 31) - 1,
            3221225473,
            4294967291,
        ] {
            let f = fp(p);
            for _ in 0..20 {
                let n = rng.gen_range(0..200usize);
                let init: Vec<u64> = (0..n).map(|_| rng.gen_range(0..p)).collect();
                let steps: Vec<(u32, Vec<u32>)> = (0..rng.gen_range(1..40))
                    .map(|_| {
                        let a = rng.gen_range(0..p) as u32;
                        (a, (0..n).map(|_| rng.gen_range(0..p) as u32).collect())
                    })
                    .collect();
                let mut exact: Vec<u128> = init.iter().map(|&x| u128::from(x)).collect();
                for (a, b) in &steps {
                    for (x, &y) in exact.iter_mut().zip(b) {
                        *x = (*x + u128::from(*a) * u128::from(y)) % u128::from(p);
                    }
                }
                for &l in &lanes {
                    let kern = Kernel::with_lanes(f, l);
                    let mut acc = init.clone();
                    for (a, b) in &steps {
                        kern.axpy(&mut acc, *a, b);
                    }
                    let got: Vec<u128> = acc.iter().map(|&x| u128::from(f.reduce(x))).collect();
                    assert_eq!(got, exact, "p = {p}, {l:?}");
                }
                // the unreduced kernel, where the budget allows it
                let kern = Kernel::with_budget(f, steps.len() + 1);
                assert_eq!(kern.exact, p == 3 || p == 786433, "p = {p}");
                if kern.exact {
                    let mut acc = init.clone();
                    for (a, b) in &steps {
                        kern.axpy(&mut acc, *a, b);
                    }
                    let got: Vec<u128> = acc.iter().map(|&x| u128::from(f.reduce(x))).collect();
                    assert_eq!(got, exact, "p = {p}, exact");
                }
            }
        }
    }

    /// The step's elimination against the plain sequential echelon, on
    /// random matrices that cross its block and chunk boundaries, with both
    /// kernels (below and above `2³¹`).
    #[test]
    fn the_elimination_matches_a_plain_echelon() {
        use rand::seq::SliceRandom;
        use rand::{Rng, SeedableRng};
        use std::collections::BTreeSet;
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        for &p in &[786433u64, 2013265921, 3221225473] {
            let f = fp(p);
            for _ in 0..6 {
                let n_cols = rng.gen_range(150..400usize);
                let density = rng.gen_range(0.02..0.15);
                let random_row = |lead: usize, rng: &mut rand::rngs::StdRng| {
                    let mut row = Row {
                        cols: vec![lead as u32],
                        vals: vec![rng.gen_range(1..p) as u32],
                    };
                    for c in lead + 1..n_cols {
                        if rng.gen_bool(density) {
                            row.cols.push(c as u32);
                            row.vals.push(rng.gen_range(1..p) as u32);
                        }
                    }
                    row
                };
                // Pivot rows: distinct, monic leads, in a shuffled order.
                let mut leads: Vec<usize> = (0..n_cols).filter(|_| rng.gen_bool(0.6)).collect();
                leads.shuffle(&mut rng);
                let mut pivot_of = vec![NONE; n_cols];
                let mut piv_rows = Vec::new();
                for &c in &leads {
                    let mut row = random_row(c, &mut rng);
                    row.vals[0] = 1;
                    pivot_of[c] = piv_rows.len() as u32;
                    piv_rows.push(row);
                }
                // Up to several chunks of S-rows (a chunk holds at most
                // 128 rows per thread).
                let s_rows: Vec<Row> = (0..rng.gen_range(50..2000))
                    .map(|_| {
                        let lead = rng.gen_range(0..n_cols);
                        random_row(lead, &mut rng)
                    })
                    .collect();
                let mut st = ElimStats::default();
                let piv_leads: Vec<u32> = piv_rows.iter().map(|r| r.cols[0]).collect();
                let Ok(ech) = eliminate(
                    f,
                    n_cols,
                    &piv_leads,
                    &|k| piv_rows[k].clone(),
                    &pivot_of,
                    s_rows.len(),
                    &|k| s_rows[k].clone(),
                    None,
                    None,
                    false,
                    &mut st,
                ) else {
                    panic!("no deadline and no cap were set");
                };
                let got = ech.sparse_rows();
                let mut all = piv_rows.clone();
                all.extend(s_rows.iter().cloned());
                let mut w = 0;
                let reference = echelon(f, all.clone(), n_cols, &mut w, None).unwrap();
                let ref_leads: BTreeSet<u32> = reference.iter().map(|r| r.cols[0]).collect();
                let piv_leads: BTreeSet<u32> = leads.iter().map(|&c| c as u32).collect();
                let new_leads: BTreeSet<u32> = got.iter().map(|r| r.cols[0]).collect();
                assert_eq!(new_leads.len(), got.len(), "p = {p}: repeated leads");
                assert!(new_leads.is_disjoint(&piv_leads), "p = {p}");
                assert_eq!(
                    ref_leads,
                    piv_leads.union(&new_leads).copied().collect(),
                    "p = {p}: a different row space"
                );
                // The rows lie in the row space: adding them keeps the rank.
                all.extend(got.iter().cloned());
                assert_eq!(
                    echelon(f, all, n_cols, &mut w, None).unwrap().len(),
                    reference.len(),
                    "p = {p}: a row outside the row space"
                );
                for r in &got {
                    assert_eq!(r.vals[0], 1);
                    assert!(r.cols.windows(2).all(|w| w[0] < w[1]));
                    assert!(r.cols.iter().all(|&c| pivot_of[c as usize] == NONE));
                    assert!(r.vals.iter().all(|&v| v != 0 && u64::from(v) < p));
                }
            }
        }
    }

    /// The full-rank exit leaves the echelon as it was: the same pivots,
    /// row for row, and no more multiply-adds, on random matrices whose
    /// S-rows outnumber the columns without a pivot row (so that most
    /// fill the echelon, some in the middle of a chunk), and on matrices
    /// where they do not.
    #[test]
    fn the_full_rank_exit_changes_no_pivot() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(23);
        let (mut exits, mut skipped_total) = (0, 0usize);
        for &p in &[786433u64, 2013265921, 3221225473] {
            let f = fp(p);
            for trial in 0..10 {
                let n_cols = rng.gen_range(60..300usize);
                let density = rng.gen_range(0.05..0.4);
                let random_row = |lead: usize, rng: &mut rand::rngs::StdRng| {
                    let mut row = Row {
                        cols: vec![lead as u32],
                        vals: vec![rng.gen_range(1..p) as u32],
                    };
                    for c in lead + 1..n_cols {
                        if rng.gen_bool(density) {
                            row.cols.push(c as u32);
                            row.vals.push(rng.gen_range(1..p) as u32);
                        }
                    }
                    row
                };
                let mut pivot_of = vec![NONE; n_cols];
                let mut piv_rows = Vec::new();
                for c in 0..n_cols {
                    if rng.gen_bool(0.5) {
                        let mut row = random_row(c, &mut rng);
                        row.vals[0] = 1;
                        pivot_of[c] = piv_rows.len() as u32;
                        piv_rows.push(row);
                    }
                }
                let nq = pivot_of.iter().filter(|&&x| x == NONE).count();
                // Even trials: several times as many S-rows as Q columns.
                // Odd trials: fewer, so the echelon cannot fill.
                let n_s = if trial % 2 == 0 {
                    rng.gen_range(2 * nq + 1..6 * nq + 40)
                } else {
                    rng.gen_range(1..nq.max(2))
                };
                let s_rows: Vec<Row> = (0..n_s)
                    .map(|_| {
                        let lead = rng.gen_range(0..n_cols);
                        random_row(lead, &mut rng)
                    })
                    .collect();
                let piv_leads: Vec<u32> = piv_rows.iter().map(|r| r.cols[0]).collect();
                let run = |exit: bool| {
                    let mut st = ElimStats::default();
                    let Ok(ech) = eliminate(
                        f,
                        n_cols,
                        &piv_leads,
                        &|k| piv_rows[k].clone(),
                        &pivot_of,
                        s_rows.len(),
                        &|k| s_rows[k].clone(),
                        None,
                        None,
                        exit,
                        &mut st,
                    ) else {
                        panic!("no deadline and no cap were set");
                    };
                    (ech, st)
                };
                let (full, a) = run(false);
                let (cut, b) = run(true);
                assert_eq!(full.leads, cut.leads, "p = {p}, trial {trial}");
                assert_eq!(full.rows, cut.rows, "p = {p}, trial {trial}");
                assert_eq!(full.q_col, cut.q_col);
                assert_eq!(
                    a.full_rank_after, b.full_rank_after,
                    "p = {p}, trial {trial}"
                );
                assert_eq!(a.skipped_rows, 0);
                assert!(b.muladds <= a.muladds, "p = {p}, trial {trial}");
                assert!(b.residual_rows <= a.residual_rows);
                match a.full_rank_after {
                    None => {
                        assert!(full.leads.len() < nq);
                        assert_eq!(b.skipped_rows, 0);
                        assert_eq!(b.muladds, a.muladds);
                        assert_eq!(b.residual_rows, a.residual_rows);
                    }
                    Some(k) => {
                        assert_eq!(full.leads.len(), nq, "p = {p}, trial {trial}");
                        assert!(k <= n_s);
                        // Every row after the one that filled the echelon
                        // is skipped: the rest of its chunk, and every chunk
                        // after it.
                        assert_eq!(b.skipped_rows, n_s - k, "p = {p}, trial {trial}");
                        if b.skipped_rows > 0 {
                            exits += 1;
                            skipped_total += b.skipped_rows;
                        }
                    }
                }
            }
        }
        assert!(exits >= 5, "the exit fired in only {exits} trials");
        eprintln!("the exit fired in {exits} trials and skipped {skipped_total} S-rows");
    }

    /// The whole engine with the exit on and off: the same basis, the same
    /// report and the same trace, but for the multiply-adds, the residual
    /// rows and the skipped rows, on Kummer systems like the pilot's and
    /// on towers with every rule coefficient non-zero.
    #[test]
    fn the_full_rank_exit_changes_no_basis() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(2027);
        let (mut exits, mut saved) = (0, 0u64);
        for trial in 0..40 {
            let general = trial % 2 == 1;
            let (ring, g) = if general {
                let ring = random_ring(&mut rng, 1_000_003, 2, 2 + trial % 3, trial % 2, true);
                let g = biquadratic(&mut rng, &ring, 2 + trial % 3, None);
                (ring, g)
            } else {
                let p = 97u64;
                let t = 1 + trial % 4;
                let (r0, _) = kummer_block(p, t, 0);
                let (r1, _) = kummer_block(p, t, t);
                let ring = TowerRing {
                    p,
                    rules: [r0, r1].concat(),
                    n_free: 0,
                };
                let g = biquadratic(&mut rng, &ring, t, None);
                (ring, g)
            };
            if g.is_empty() {
                continue;
            }
            let on = TowerF4Options::new(64);
            let mut off = TowerF4Options::new(64);
            off.full_rank_exit = false;
            let a = f4_tower(std::slice::from_ref(&g), &ring, &off);
            let b = f4_tower(std::slice::from_ref(&g), &ring, &on);
            assert_eq!(a.basis, b.basis, "trial {trial}");
            assert_eq!(a.inconsistent, b.inconsistent);
            assert_eq!(a.solving_degree_max, b.solving_degree_max);
            assert_eq!(a.last_productive_degree, b.last_productive_degree);
            assert_eq!(a.degree_reached, b.degree_reached);
            assert_eq!(
                (
                    a.steps,
                    a.max_rows,
                    a.max_cols,
                    a.max_nnz,
                    a.max_dense_entries
                ),
                (
                    b.steps,
                    b.max_rows,
                    b.max_cols,
                    b.max_nnz,
                    b.max_dense_entries
                )
            );
            assert_eq!(
                (
                    a.critical_pairs_reduced,
                    a.tower_pairs_reduced,
                    a.reducer_rows
                ),
                (
                    b.critical_pairs_reduced,
                    b.tower_pairs_reduced,
                    b.reducer_rows
                )
            );
            assert_eq!(a.rows_skipped, 0);
            assert!(b.muladds <= a.muladds, "trial {trial}");
            assert_eq!(a.trace.len(), b.trace.len());
            for (x, y) in a.trace.iter().zip(&b.trace) {
                assert_eq!(
                    (
                        x.degree,
                        x.s_rows,
                        x.reducer_rows,
                        x.promoted_rows,
                        x.cols,
                        x.nnz
                    ),
                    (
                        y.degree,
                        y.s_rows,
                        y.reducer_rows,
                        y.promoted_rows,
                        y.cols,
                        y.nnz
                    )
                );
                assert_eq!(
                    (
                        x.residual_cols,
                        x.dense_entries,
                        x.fresh,
                        x.fresh_min_degree
                    ),
                    (
                        y.residual_cols,
                        y.dense_entries,
                        y.fresh,
                        y.fresh_min_degree
                    )
                );
                assert_eq!(
                    (x.basis_active, x.basis_kept, x.basis_entries, x.pairs_left),
                    (y.basis_active, y.basis_kept, y.basis_entries, y.pairs_left)
                );
                assert_eq!(x.full_rank_after, y.full_rank_after);
                assert_eq!(x.skipped_rows, 0);
                // Every S-row after the one that filled the echelon.
                assert_eq!(
                    y.skipped_rows,
                    y.full_rank_after.map_or(0, |k| y.s_rows - k)
                );
                if y.skipped_rows == 0 {
                    assert_eq!(x.muladds, y.muladds);
                    assert_eq!(x.residual_rows, y.residual_rows);
                } else {
                    assert!(y.muladds <= x.muladds);
                    assert!(y.residual_rows <= x.residual_rows);
                }
            }
            if b.rows_skipped > 0 {
                exits += 1;
                saved += a.muladds - b.muladds;
            }
        }
        assert!(exits > 0, "the exit never fired");
        eprintln!("the exit fired on {exits} systems and saved {saved} multiply-adds");
    }

    /// Generic reduction by `y_j² − rule_j` on exponent vectors, as an
    /// independent reference for the rewriting.
    fn reference_nf(ring: &TowerRing, terms: &[(Vec<u32>, u64)]) -> BTreeMap<Vec<u32>, u64> {
        let p = ring.p;
        let f = fp(p);
        let mut work: BTreeMap<Vec<u32>, u64> = BTreeMap::new();
        for (e, c) in terms {
            let v = work.entry(e.clone()).or_insert(0);
            *v = f.add(*v, c % p);
        }
        loop {
            let hit = work
                .iter()
                .find(|(e, c)| **c != 0 && e.iter().take(ring.n_tower()).any(|&x| x >= 2))
                .map(|(e, c)| (e.clone(), *c));
            let Some((e, c)) = hit else { break };
            work.remove(&e);
            let j = e.iter().position(|&x| x >= 2).unwrap();
            let r = ring.rules[j];
            let mut base = e.clone();
            base[j] -= 2;
            let mut push = |add: &[(usize, u32)], coef: u64| {
                if coef == 0 {
                    return;
                }
                let mut t = base.clone();
                for &(v, k) in add {
                    t[v] += k;
                }
                let v = work.entry(t).or_insert(0);
                *v = f.add(*v, f.mul(c, coef));
            };
            push(&[], r.d);
            push(&[(j, 1)], r.c);
            if let Some(nx) = r.next {
                push(&[(nx as usize, 1)], r.b);
                push(&[(j, 1), (nx as usize, 1)], r.a);
            }
        }
        work.retain(|_, c| *c != 0);
        work
    }

    fn random_ring(
        rng: &mut impl rand::Rng,
        p: u64,
        blocks: usize,
        t: usize,
        n_free: usize,
        general: bool,
    ) -> TowerRing {
        let mut rules = Vec::new();
        for b in 0..blocks {
            for j in 0..t {
                let top = j + 1 == t;
                let r = if general {
                    SquareRule {
                        next: (!top).then_some((b * t + j + 1) as u8),
                        a: if top { 0 } else { rng.gen_range(0..p) },
                        b: if top { 0 } else { rng.gen_range(0..p) },
                        c: rng.gen_range(0..p),
                        d: rng.gen_range(0..p),
                    }
                } else {
                    SquareRule {
                        next: (!top).then_some((b * t + j + 1) as u8),
                        a: 0,
                        b: if top { 0 } else { 1 },
                        c: 0,
                        d: if top { rng.gen_range(1..p) } else { 0 },
                    }
                };
                rules.push(r);
            }
        }
        TowerRing { p, rules, n_free }
    }

    #[test]
    fn products_match_a_generic_reduction() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let p = 1_000_003u64;
        for trial in 0..60 {
            let general = trial % 2 == 0;
            let ring = random_ring(&mut rng, p, 2, 4, 1, general);
            let n = ring.n_vars();
            // Sparse raw terms: the reference rewrites whole exponent
            // vectors, and a dense vector like (2, 2, 2, 2) on one tower
            // swells through millions of intermediate monomials before it
            // settles (the engine never sees them: it stays square-free
            // after every variable).
            let rand_raw =
                |rng: &mut rand::rngs::StdRng, terms: usize, max_e: u32| -> Vec<(Vec<u32>, u64)> {
                    (0..terms)
                        .map(|_| {
                            let mut e = vec![0u32; n];
                            let mut budget = 4u32;
                            for x in e.iter_mut() {
                                if budget > 0 && rng.gen_range(0..10) < 3 {
                                    *x = rng.gen_range(1..=max_e.min(budget));
                                    budget -= *x;
                                }
                            }
                            (e, rng.gen_range(1..p))
                        })
                        .collect()
                };
            // from_raw against the reference
            let raw = rand_raw(&mut rng, 6, 3);
            let got = ring.from_raw(&raw);
            let want = reference_nf(&ring, &raw);
            let got_map: BTreeMap<Vec<u32>, u64> = got
                .monos
                .iter()
                .zip(&got.coefs)
                .map(|(&m, &c)| {
                    let mut e = vec![0u32; n];
                    let mask = mask_of(m);
                    for (j, x) in e.iter_mut().enumerate().take(ring.n_tower()) {
                        *x = ((mask >> j) & 1) as u32;
                    }
                    for k in 0..ring.n_free {
                        e[ring.n_tower() + k] = u32::from(free_of(m)[k]);
                    }
                    (e, u64::from(c))
                })
                .collect();
            assert_eq!(got_map, want, "trial {trial}");
            // monomial × polynomial against from_raw of the raw product
            let g = ring.from_raw(&rand_raw(&mut rng, 5, 1));
            let u_raw: Vec<u32> = (0..n)
                .map(|i| {
                    if i < ring.n_tower() {
                        rng.gen_range(0..2)
                    } else {
                        rng.gen_range(0..3)
                    }
                })
                .collect();
            let u = key_of(&u_raw, ring.n_tower());
            let prod = ring.mul_mono(&g, u);
            let raw_prod: Vec<(Vec<u32>, u64)> = g
                .monos
                .iter()
                .zip(&g.coefs)
                .map(|(&m, &c)| {
                    let mut e = vec![0u32; n];
                    let mask = mask_of(m);
                    for (j, x) in e.iter_mut().enumerate().take(ring.n_tower()) {
                        *x = ((mask >> j) & 1) as u32 + u_raw[j];
                    }
                    for k in 0..ring.n_free {
                        e[ring.n_tower() + k] =
                            u32::from(free_of(m)[k]) + u_raw[ring.n_tower() + k];
                    }
                    (e, u64::from(c))
                })
                .collect();
            assert_eq!(prod, ring.from_raw(&raw_prod), "trial {trial}");
        }
    }

    /// A basis element kept as an echelon row, zeros and all, over a column
    /// list it shares with others, has its polynomial's terms and products.
    #[test]
    fn a_stored_element_is_its_polynomial() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let p = 1_000_003u64;
        for trial in 0..40 {
            let ring = random_ring(&mut rng, p, 2, 4, 1, trial % 2 == 0);
            let (n, nt) = (ring.n_vars(), ring.n_tower());
            let exponents = |rng: &mut rand::rngs::StdRng| -> Vec<u32> {
                (0..n)
                    .map(|i| rng.gen_range(0..if i < nt { 2 } else { 3 }))
                    .collect()
            };
            let raw: Vec<(Vec<u32>, u64)> = (0..6)
                .map(|_| (exponents(&mut rng), rng.gen_range(1..p)))
                .collect();
            let g = ring.from_raw(&raw);
            if g.is_empty() {
                continue;
            }
            // Its monomials among others, some before its lead, some between
            // its terms and some after the last.
            let mut cols: Vec<Mono> = g.monos.clone();
            cols.extend((0..20).map(|_| key_of(&exponents(&mut rng), nt)));
            cols.sort_unstable_by(|a, b| b.cmp(a));
            cols.dedup();
            let start = cols.iter().position(|&m| m == g.monos[0]).unwrap();
            let vals: Vec<u32> = cols[start..]
                .iter()
                .map(|m| {
                    g.monos
                        .iter()
                        .position(|x| x == m)
                        .map_or(0, |i| g.coefs[i])
                })
                .collect();
            let e = Elem::new(cols.into(), start as u32, vals);
            assert_eq!(e.lm(), g.monos[0]);
            assert_eq!(e.len(), g.len());
            assert_eq!(e.to_rpoly(), g, "trial {trial}");
            let u = key_of(&exponents(&mut rng), nt);
            assert_eq!(
                ring.mul_elem(fp(p), &e, u, TowerWeights::UNIT),
                ring.mul_mono(&g, u),
                "trial {trial}"
            );
        }
    }

    // ── Small towers with points, for end-to-end checks ──

    /// A Kummer tower of `t` levels over `p` (needs `2^t | p − 1`): the
    /// rules, and every point `(y_0 … y_{t−1})` with `y_0` in `μ_{2^t}`.
    fn kummer_block(p: u64, t: usize, base: usize) -> (Vec<SquareRule>, Vec<Vec<u64>>) {
        let f = fp(p);
        let zeta = (2..p)
            .map(|h| f.pow(h, (p - 1) >> t))
            .find(|&z| f.pow(z, 1 << (t - 1)) != 1)
            .unwrap();
        let rules: Vec<SquareRule> = (0..t)
            .map(|j| SquareRule {
                next: (j + 1 < t).then_some((base + j + 1) as u8),
                a: 0,
                b: if j + 1 < t { 1 } else { 0 },
                c: 0,
                d: if j + 1 < t { 0 } else { 1 },
            })
            .collect();
        let points = (0..1u64 << t)
            .map(|k| {
                let mut y = vec![f.pow(zeta, k)];
                for j in 1..t {
                    y.push(f.mul(y[j - 1], y[j - 1]));
                }
                y
            })
            .collect();
        (rules, points)
    }

    /// A random polynomial of bidegree at most (2, 2) in `x_1 = y_{0,0}`,
    /// `x_2 = y_{1,0}`, the shape of `S₃`, in tower normal form.
    fn biquadratic(
        rng: &mut impl rand::Rng,
        ring: &TowerRing,
        t: usize,
        zero_at: Option<&[u64]>,
    ) -> RPoly {
        let n = ring.n_vars();
        let mut raw = Vec::new();
        for a in 0..3u32 {
            for b in 0..3u32 {
                let mut e = vec![0u32; n];
                e[0] = a;
                e[t] = b;
                raw.push((e, rng.gen_range(1..ring.p)));
            }
        }
        let mut g = ring.from_raw(&raw);
        if let Some(pt) = zero_at {
            let v = ring.eval(&g, pt);
            if v != 0 {
                g = ring.from_raw(&[raw, vec![(vec![0; n], ring.p - v)]].concat());
            }
        }
        g
    }

    #[test]
    fn refutes_exactly_the_empty_systems_and_counts_the_rest() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(2026);
        let p = 97u64; // 2^5 | 96
        for t in 1..=4usize {
            let (r0, pts) = kummer_block(p, t, 0);
            let (r1, _) = kummer_block(p, t, t);
            let ring = TowerRing {
                p,
                rules: [r0, r1].concat(),
                n_free: 0,
            };
            for trial in 0..12 {
                let plant: Option<Vec<u64>> = (trial % 3 == 0).then(|| {
                    let a = &pts[rng.gen_range(0..pts.len())];
                    let b = &pts[rng.gen_range(0..pts.len())];
                    [a.clone(), b.clone()].concat()
                });
                let g = biquadratic(&mut rng, &ring, t, plant.as_deref());
                // brute force over V × V
                let mut zeros = Vec::new();
                for a in &pts {
                    for b in &pts {
                        let x = [a.clone(), b.clone()].concat();
                        if ring.eval(&g, &x) == 0 {
                            zeros.push(x);
                        }
                    }
                }
                let rep = f4_tower(&[g], &ring, &TowerF4Options::new(64));
                assert!(!rep.timed_out && rep.pairs_above_bound == 0);
                assert_eq!(rep.inconsistent, zeros.is_empty(), "t {t} trial {trial}");
                if !rep.inconsistent {
                    for z in &zeros {
                        for b in &rep.basis {
                            assert_eq!(ring.eval(b, z), 0, "a solution was lost");
                        }
                    }
                    // A zero-dimensional radical ideal: its staircase counts
                    // its points.
                    let lms: Vec<Mono> = rep.basis.iter().map(|b| b.lm().unwrap()).collect();
                    assert_eq!(
                        staircase_at_most(&lms, ring.n_tower(), 0, 1 << 12),
                        Some(zeros.len())
                    );
                }
            }
        }
    }

    // ── Weighted gradings (radix-2^k presentations of the same tower) ──

    /// The key `mono` made before gradings existed, written out again.
    fn unit_key_reference(mask: u64, free: [u8; MAX_FREE]) -> Mono {
        let deg = mask.count_ones() + free.iter().map(|&e| u32::from(e)).sum::<u32>();
        let mut key = (u128::from(deg) << DEG_SHIFT) | u128::from(!mask);
        for (k, &e) in free.iter().enumerate() {
            key |= u128::from(255 - e) << (FREE_SHIFT + 8 * k as u32);
        }
        key
    }

    #[test]
    fn unit_weights_leave_every_key_unchanged() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        assert!(TowerWeights::radix_blocks(3, 5, 1).is_unit());
        assert!(TowerWeights::from_slice(&[1; 64]).is_unit());
        for _ in 0..20_000 {
            let (ma, mb): (u64, u64) = (rng.gen(), rng.gen());
            let fa = [
                rng.gen_range(0..9u8),
                rng.gen_range(0..9),
                0,
                rng.gen_range(0..3),
            ];
            let fb = [rng.gen_range(0..9u8), 0, rng.gen_range(0..4), 0];
            let (a, b) = (mono(ma, fa), mono(mb, fb));
            assert_eq!(a, unit_key_reference(ma, fa));
            assert_eq!(mono_w(ma, fa, TowerWeights::UNIT), a);
            let l = lcm(a, b);
            assert_eq!(lcm_w(a, b, TowerWeights::UNIT), l);
            // The quotient now subtracts degrees; it must equal the old
            // recomputation exactly.
            let (qa, qb) = (free_of(l), free_of(b));
            let mut qf = [0u8; MAX_FREE];
            for k in 0..MAX_FREE {
                qf[k] = qa[k] - qb[k];
            }
            assert_eq!(quotient(l, b), unit_key_reference(mask_of(l) & !mb, qf));
        }
    }

    #[test]
    fn a_weighted_key_order_is_a_monomial_order() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        for k in 2..=3u32 {
            let w = TowerWeights::radix_blocks(4, 6, k);
            assert_eq!(w.weight(0), 1);
            assert_eq!(w.weight(1), 2);
            assert_eq!(w.weight(6), 1);
            for _ in 0..20_000 {
                let n = 24u32;
                let full = (1u64 << n) - 1;
                let (ma, mb, mu) = (
                    rng.gen::<u64>() & full,
                    rng.gen::<u64>() & full,
                    rng.gen::<u64>() & full,
                );
                let (ma, mb) = (ma & !mu, mb & !mu);
                let z = [0u8; MAX_FREE];
                let (a, b, u) = (mono_w(ma, z, w), mono_w(mb, z, w), mono_w(mu, z, w));
                let (au, bu) = (mono_w(ma | mu, z, w), mono_w(mb | mu, z, w));
                // Multiplicative, and the plain-product shortcut is exact.
                assert_eq!(a < b, au < bu);
                let shift = (u128::from(degree_of(u)) << DEG_SHIFT) - u128::from(mu);
                assert_eq!(a + shift, au);
                // Divisibility implies order, 1 is least, and the degree is
                // the weighted one.
                assert!(au >= a && a >= ONE);
                assert_eq!(
                    degree_of(a),
                    (0..n)
                        .filter(|j| ma >> j & 1 == 1)
                        .map(|j| w.weight(j))
                        .sum::<u32>()
                );
                assert!(divides_ungraded(a, au) && divides(a, au));
                assert_eq!(quotient(au, a), u);
                let _ = (b, bu);
            }
        }
    }

    #[test]
    fn weighted_runs_refute_exactly_the_empty_systems_and_count_the_rest() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(2027);
        let p = 97u64; // 2^5 | 96
        for (t, k) in [(2usize, 2u32), (4, 2), (3, 3), (4, 3), (5, 2)] {
            let (r0, pts) = kummer_block(p, t, 0);
            let (r1, _) = kummer_block(p, t, t);
            let ring = TowerRing {
                p,
                rules: [r0, r1].concat(),
                n_free: 0,
            };
            let w = TowerWeights::radix_blocks(2, t, k);
            for trial in 0..10 {
                let plant: Option<Vec<u64>> = (trial % 2 == 0).then(|| {
                    let a = &pts[rng.gen_range(0..pts.len())];
                    let b = &pts[rng.gen_range(0..pts.len())];
                    [a.clone(), b.clone()].concat()
                });
                let g = biquadratic(&mut rng, &ring, t, plant.as_deref());
                let mut zeros = Vec::new();
                for a in &pts {
                    for b in &pts {
                        let x = [a.clone(), b.clone()].concat();
                        if ring.eval(&g, &x) == 0 {
                            zeros.push(x);
                        }
                    }
                }
                let unit = f4_tower(std::slice::from_ref(&g), &ring, &TowerF4Options::new(64));
                let rep = f4_tower(&[g], &ring, &TowerF4Options::new(64).with_weights(w));
                assert!(!rep.timed_out && rep.pairs_above_bound == 0);
                assert_eq!(
                    rep.inconsistent,
                    zeros.is_empty(),
                    "t {t} k {k} trial {trial}"
                );
                assert_eq!(rep.inconsistent, unit.inconsistent);
                if !rep.inconsistent {
                    for z in &zeros {
                        for b in &rep.basis {
                            assert_eq!(ring.eval(b, z), 0, "a solution was lost (t {t} k {k})");
                        }
                    }
                    // Every ideal of F_p^{V×V} is radical: the staircase of
                    // the weighted basis counts the points too.
                    let lms: Vec<Mono> = rep.basis.iter().map(|b| b.lm().unwrap()).collect();
                    assert_eq!(
                        staircase_at_most(&lms, ring.n_tower(), 0, 1 << 12),
                        Some(zeros.len()),
                        "t {t} k {k} trial {trial}"
                    );
                    // The basis is keyed in the grading.
                    for b in &rep.basis {
                        for &m in &b.monos {
                            assert_eq!(m, mono_w(mask_of(m), free_of(m), w));
                        }
                    }
                }
            }
        }
    }

    #[test]
    #[should_panic(expected = "is not the leading monomial")]
    fn a_grading_that_unseats_the_square_is_refused() {
        let (rules, _) = kummer_block(97, 2, 0);
        let ring = TowerRing {
            p: 97,
            rules,
            n_free: 0,
        };
        // w(y_1) = 3 > 2·w(y_0): y_1 would lead y_0² − y_1.
        TowerWeights::from_slice(&[1, 3]).check(&ring);
    }

    #[test]
    fn a_free_variable_is_eliminated_like_a_chain_link() {
        // x1 in a Kummer set V, u free: u² − x1 = 0 and u − c = 0 for a c
        // whose square is or is not in V.
        let p = 97u64;
        let t = 3usize;
        let (rules, pts) = kummer_block(p, t, 0);
        let ring = TowerRing {
            p,
            rules,
            n_free: 1,
        };
        let n = ring.n_vars();
        let f = fp(p);
        for c in 1..p {
            let sq = f.mul(c, c);
            let in_v = pts.iter().any(|y| y[0] == sq);
            let mut e_u2 = vec![0u32; n];
            e_u2[t] = 2;
            let mut e_x = vec![0u32; n];
            e_x[0] = 1;
            let mut e_u = vec![0u32; n];
            e_u[t] = 1;
            let g1 = ring.from_raw(&[(e_u2, 1), (e_x, p - 1)]);
            let g2 = ring.from_raw(&[(e_u, 1), (vec![0; n], p - c)]);
            let rep = f4_tower(&[g1, g2], &ring, &TowerF4Options::new(64));
            assert_eq!(rep.inconsistent, !in_v, "c = {c}");
        }
    }

    #[test]
    fn the_stop_keeps_every_solution() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        let p = 97u64;
        let t = 4usize;
        let (r0, pts) = kummer_block(p, t, 0);
        let (r1, _) = kummer_block(p, t, t);
        let ring = TowerRing {
            p,
            rules: [r0, r1].concat(),
            n_free: 0,
        };
        for _ in 0..10 {
            let a = pts[rng.gen_range(0..pts.len())].clone();
            let b = pts[rng.gen_range(0..pts.len())].clone();
            let x = [a, b].concat();
            let g = biquadratic(&mut rng, &ring, t, Some(&x));
            let full = f4_tower(std::slice::from_ref(&g), &ring, &TowerF4Options::new(64));
            let stopped = f4_tower(&[g], &ring, &TowerF4Options::new(64).stopping_below(4));
            assert!(!stopped.inconsistent);
            for b in &stopped.basis {
                assert_eq!(ring.eval(b, &x), 0);
            }
            assert!(stopped.solving_degree_max <= full.solving_degree_max);
            if let Some(k) = stopped.staircase_at_stop {
                assert!((1..=4).contains(&k));
            }
        }
    }

    #[test]
    fn the_degree_bound_stops_and_reports() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        let p = 97u64;
        let t = 4usize;
        let (r0, _) = kummer_block(p, t, 0);
        let (r1, _) = kummer_block(p, t, t);
        let ring = TowerRing {
            p,
            rules: [r0, r1].concat(),
            n_free: 0,
        };
        let g = biquadratic(&mut rng, &ring, t, None);
        let rep = f4_tower(&[g], &ring, &TowerF4Options::new(2));
        assert!(rep.degree_reached <= 2);
        assert!(rep.inconsistent || rep.pairs_above_bound > 0);
    }

    /// A run the size cap stops reports the step it cut off: the step the
    /// uncapped run took there, as far as the stop let it get, after the
    /// same earlier steps. A run that ends reports none.
    #[test]
    fn the_size_cap_reports_the_stopped_step() {
        use rand::SeedableRng;
        use std::sync::atomic::AtomicUsize;
        // `on_stop` is a plain function, so it reports through a static that
        // only this test touches.
        static STOPS: AtomicUsize = AtomicUsize::new(0);
        fn on_stop(k: usize, oversize: bool, st: &StepTrace) {
            assert!(k > 0 && oversize && st.fresh == 0);
            STOPS.fetch_add(1, AtomicOrdering::Relaxed);
        }
        let mut rng = rand::rngs::StdRng::seed_from_u64(31);
        let (mut in_elimination, mut in_preprocessing, mut dense_only) = (0, 0, 0);
        for trial in 0..24 {
            let p = 97u64;
            let t = 2 + trial % 3;
            let (r0, _) = kummer_block(p, t, 0);
            let (r1, _) = kummer_block(p, t, t);
            let ring = TowerRing {
                p,
                rules: [r0, r1].concat(),
                n_free: 0,
            };
            let g = biquadratic(&mut rng, &ring, t, None);
            let full = f4_tower(std::slice::from_ref(&g), &ring, &TowerF4Options::new(64));
            assert!(full.stopped_step.is_none() && !full.oversize);
            let caps: Vec<u64> = full
                .trace
                .iter()
                .flat_map(|s| [s.nnz, s.dense_entries])
                .filter(|&x| x > 0)
                .map(|x| x - 1)
                .collect();
            for cap in caps {
                // `max_dense` alone bounds `B'` and nothing else: the run
                // stops at the first step whose `B'` passes it.
                let mut dense = TowerF4Options::new(64).with_max_dense(cap);
                dense.on_stop = Some(on_stop);
                let rep = f4_tower(std::slice::from_ref(&g), &ring, &dense);
                let first = full.trace.iter().position(|s| s.dense_entries > cap);
                match (&rep.stopped_step, first) {
                    (None, None) => assert_eq!(rep.basis, full.basis),
                    (Some(st), Some(k)) => {
                        dense_only += 1;
                        assert!(rep.oversize && rep.trace.len() == k);
                        let want = &full.trace[k];
                        assert_eq!(
                            (st.cols, st.nnz, st.residual_cols, st.dense_entries),
                            (want.cols, want.nnz, want.residual_cols, want.dense_entries)
                        );
                    }
                    _ => panic!("trial {trial} cap {cap}: max_dense stopped elsewhere"),
                }
                let mut opts = TowerF4Options::new(64).with_max_nnz(cap);
                opts.on_stop = Some(on_stop);
                let rep = f4_tower(std::slice::from_ref(&g), &ring, &opts);
                let Some(st) = &rep.stopped_step else {
                    assert!(!rep.oversize);
                    assert_eq!(rep.basis, full.basis, "trial {trial} cap {cap}");
                    continue;
                };
                assert!(rep.oversize && !rep.inconsistent, "trial {trial} cap {cap}");
                let k = rep.trace.len();
                assert_eq!(rep.steps, k + 1);
                let want = &full.trace[k];
                assert_eq!(st.degree, want.degree);
                for (x, y) in rep.trace.iter().zip(&full.trace) {
                    assert_eq!((x.cols, x.nnz, x.fresh), (y.cols, y.nnz, y.fresh));
                }
                // The basis before the step, and its pairs given back.
                if let Some(b) = k.checked_sub(1).map(|j| &full.trace[j]) {
                    assert_eq!(
                        (
                            st.basis_active,
                            st.basis_kept,
                            st.basis_entries,
                            st.pairs_left
                        ),
                        (b.basis_active, b.basis_kept, b.basis_entries, b.pairs_left)
                    );
                }
                assert!(st.pairs_left >= want.critical_pairs + want.tower_pairs);
                assert_eq!((st.fresh, st.residual_rows, st.muladds), (0, 0, 0));
                if st.cols > 0 {
                    // Stopped by `B'`: the whole matrix was built.
                    in_elimination += 1;
                    assert_eq!(
                        (
                            st.s_rows,
                            st.reducer_rows,
                            st.promoted_rows,
                            st.cols,
                            st.nnz
                        ),
                        (
                            want.s_rows,
                            want.reducer_rows,
                            want.promoted_rows,
                            want.cols,
                            want.nnz
                        )
                    );
                    assert_eq!(
                        (st.residual_cols, st.dense_entries),
                        (want.residual_cols, want.dense_entries)
                    );
                    assert!(st.dense_entries > cap);
                    assert_eq!(
                        rep.max_dense_entries,
                        full.trace[..=k]
                            .iter()
                            .map(|s| s.dense_entries)
                            .max()
                            .unwrap()
                    );
                } else {
                    // Stopped while the rows were scanned.
                    in_preprocessing += 1;
                    assert!(st.nnz > cap && st.nnz <= want.nnz);
                    assert_eq!(st.dense_entries, 0);
                }
            }
        }
        assert!(
            in_elimination > 0 && in_preprocessing > 0 && dense_only > 0,
            "{in_elimination} {in_preprocessing} {dense_only}"
        );
        assert_eq!(
            STOPS.load(AtomicOrdering::Relaxed),
            in_elimination + in_preprocessing + dense_only
        );
    }

    /// `Fp::reduce` is the remainder for every `x < 2⁶⁴`, at the edges of
    /// each multiple of `p` as well as at random.
    #[test]
    fn reduce_is_the_remainder() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(61);
        for &p in &[
            3u64,
            5,
            97,
            786433,
            (1 << 31) - 1,
            (1 << 31) + 11,
            4294967291,
            4294967295,
        ] {
            let f = fp(p);
            let mut xs = vec![0, 1, p - 1, p, p + 1, u64::MAX, u64::MAX - 1, u64::MAX - p];
            for _ in 0..2000 {
                let k = rng.gen_range(0..=u64::MAX / p);
                xs.extend([k * p, (k * p).saturating_sub(1), rng.gen::<u64>()]);
            }
            for x in xs {
                assert_eq!(f.reduce(x), x % p, "p = {p}, x = {x}");
            }
        }
    }

    /// The product of a polynomial and a monomial, formed by sorting every
    /// term's product together, as it was before the terms that carry were
    /// told from the ones that do not.
    fn product_by_full_sort(ring: &TowerRing, f: Fp, g: &RPoly, u: Mono) -> (RPoly, usize) {
        let (mu, fu) = (mask_of(u), free_of(u));
        let mut terms = Vec::new();
        let mut tmp = Vec::new();
        for (&m, &c) in g.monos.iter().zip(&g.coefs) {
            let free = add_free(free_of(m), fu);
            tmp.clear();
            ring.mask_product_into(
                f,
                mask_of(m),
                mu,
                u64::from(c),
                &mut tmp,
                &mut MaskScratch::default(),
            );
            terms.extend(tmp.iter().map(|&(x, cc)| (mono(x, free), cc)));
        }
        let before = terms.len();
        (RPoly::from_terms(f, terms), before)
    }

    /// Merging the plain terms with the sorted carried ones is `from_terms`
    /// of all of them: on general and Kummer towers, over small and large
    /// primes, with repeated monomials and cancelling sums.
    #[test]
    fn a_product_merges_plain_and_carried_terms_like_a_full_sort() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(31);
        let (mut mixed, mut merged, mut cancelled) = (0, 0, 0);
        // Products with plain terms, by whether `u` has a free variable: the
        // key sum or `mono` forms their keys.
        let (mut by_sum, mut by_mono) = (0, 0);
        for &p in &[3u64, 7, 97, 1_000_003, 4294967291] {
            let f = fp(p);
            for trial in 0..120 {
                let ring = random_ring(&mut rng, p, 2, 4, 2, trial % 2 == 0);
                let (n, nt) = (ring.n_vars(), ring.n_tower());
                let exps = |rng: &mut rand::rngs::StdRng, density: f64| -> Vec<u32> {
                    (0..n)
                        .map(|i| {
                            if !rng.gen_bool(density) {
                                0
                            } else if i < nt {
                                1
                            } else {
                                rng.gen_range(1..3)
                            }
                        })
                        .collect()
                };
                let raw: Vec<(Vec<u32>, u64)> = (0..rng.gen_range(1..16))
                    .map(|_| (exps(&mut rng, 0.4), rng.gen_range(1..p)))
                    .collect();
                let g = ring.from_raw(&raw);
                let density = rng.gen_range(0.05..0.8);
                let u = key_of(&exps(&mut rng, density), nt);
                let (want, before) = product_by_full_sort(&ring, f, &g, u);
                assert_eq!(ring.mul_mono(&g, u), want, "p = {p}, trial {trial}");
                let carries = g.monos.iter().any(|&m| mask_of(m) & mask_of(u) != 0);
                let plain = g.monos.iter().any(|&m| mask_of(m) & mask_of(u) == 0);
                mixed += usize::from(carries && plain);
                merged += usize::from(want.len() < before);
                cancelled += usize::from(want.len() + 2 <= before);
                by_sum += usize::from(plain && free_of(u) == [0; MAX_FREE]);
                by_mono += usize::from(plain && free_of(u) != [0; MAX_FREE]);
            }
        }
        assert!(
            mixed > 100 && merged > 100 && cancelled > 50,
            "the cases that merge were not reached: {mixed} {merged} {cancelled}"
        );
        assert!(
            by_sum > 100 && by_mono > 100,
            "both ways to key a plain product are to be reached: {by_sum} {by_mono}"
        );
    }

    /// The same on the widest towers, 63 and 64 tower variables as one chain,
    /// as blocks of three or four and as one block each: terms and multipliers
    /// with `y_63`, the top bit of the mask, in the key sum and in the merge.
    #[test]
    fn a_product_on_a_64_variable_tower_merges_like_a_full_sort() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(43);
        let (mut mixed, mut top_bit, mut by_sum, mut by_mono) = (0, 0, 0, 0);
        for &p in &[3u64, 97, 1_000_003, 4294967291] {
            let f = fp(p);
            for &(blocks, t) in &[(1usize, 64usize), (1, 63), (16, 4), (21, 3), (64, 1)] {
                for n_free in [0, 1, 4] {
                    for _ in 0..12 {
                        let ring = random_ring(&mut rng, p, blocks, t, n_free, true);
                        ring.check();
                        let (n, nt) = (ring.n_vars(), ring.n_tower());
                        // Terms spread over the whole mask; a few common bits
                        // with the multiplier, so that the rewrites stay small.
                        let raw: Vec<(Vec<u32>, u64)> = (0..rng.gen_range(1..20))
                            .map(|_| {
                                let density = rng.gen_range(0.02..0.25);
                                let mut e: Vec<u32> = (0..n)
                                    .map(|i| {
                                        if i < nt {
                                            u32::from(rng.gen_bool(density))
                                        } else {
                                            rng.gen_range(0..4)
                                        }
                                    })
                                    .collect();
                                if nt == 64 && rng.gen_bool(0.3) {
                                    e[63] = 1;
                                }
                                (e, rng.gen_range(1..p))
                            })
                            .collect();
                        let g = ring.from_raw(&raw);
                        let mut mask = 0u64;
                        for _ in 0..rng.gen_range(0..6) {
                            mask |= 1 << rng.gen_range(0..nt);
                        }
                        if nt == 64 && rng.gen_bool(0.4) {
                            mask |= 1 << 63;
                        }
                        let mut free = [0u8; MAX_FREE];
                        if rng.gen_bool(0.5) {
                            for e in free.iter_mut().take(n_free) {
                                *e = rng.gen_range(0..3);
                            }
                        }
                        let u = mono(mask, free);
                        let (want, _) = product_by_full_sort(&ring, f, &g, u);
                        assert_eq!(ring.mul_mono(&g, u), want, "p = {p}, {blocks} x {t}");
                        let carries = g.monos.iter().any(|&m| mask_of(m) & mask != 0);
                        let plain = g.monos.iter().any(|&m| mask_of(m) & mask == 0);
                        mixed += usize::from(carries && plain);
                        top_bit += usize::from(mask >> 63 == 1 && nt == 64);
                        by_sum += usize::from(plain && free == [0; MAX_FREE]);
                        by_mono += usize::from(plain && free != [0; MAX_FREE]);
                    }
                }
            }
        }
        assert!(
            mixed > 50 && top_bit > 20 && by_sum > 50 && by_mono > 20,
            "{mixed} {top_bit} {by_sum} {by_mono}"
        );
    }

    /// The key of a plain product, one without a tower variable in common,
    /// is the term's key plus the constant `mul_terms_fp` adds when the
    /// multiplier has no free variable: at the edges of the mask (bit 63,
    /// the empty mask) and of the free bytes.
    #[test]
    fn a_plain_product_key_is_the_term_key_plus_a_constant() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(17);
        for _ in 0..20_000 {
            let (mm, mu) = match rng.gen_range(0..4) {
                0 => (rng.gen::<u64>(), 0),
                1 => (0, rng.gen::<u64>()),
                2 => (1 << 63, rng.gen::<u64>() >> 1),
                _ => {
                    let a = rng.gen::<u64>();
                    (a, !a & rng.gen::<u64>())
                }
            };
            let free: [u8; MAX_FREE] = std::array::from_fn(|_| rng.gen_range(0..=255));
            let (m, u) = (mono(mm, free), mono(mu, [0; MAX_FREE]));
            let shift = (u128::from(degree_of(u)) << DEG_SHIFT) - u128::from(mask_of(u));
            assert_eq!(mask_of(m) & mask_of(u), 0);
            assert_eq!(m + shift, mono(mm | mu, free), "mask {mm:#x} by {mu:#x}");
        }
    }

    /// The counters and basis of a run, and its trace without the timings.
    fn the_run(r: &TowerF4Report) -> (Vec<RPoly>, Vec<u64>, serde_json::Value) {
        let mut trace = serde_json::to_value(&r.trace).unwrap();
        for step in trace.as_array_mut().unwrap() {
            step.as_object_mut()
                .unwrap()
                .retain(|k, _| !k.starts_with("ms"));
        }
        let counters = vec![
            u64::from(r.inconsistent),
            u64::from(r.degree_reached),
            u64::from(r.solving_degree_max),
            u64::from(r.last_productive_degree),
            r.max_cols_to_solution as u64,
            r.steps_to_solution as u64,
            r.steps as u64,
            r.max_rows as u64,
            r.max_cols as u64,
            r.max_nnz,
            r.max_residual_rows as u64,
            r.max_dense_entries,
            r.muladds,
            r.critical_pairs_reduced,
            r.tower_pairs_reduced,
            r.pairs_product_skipped,
            r.pairs_chain_skipped,
            r.reducer_rows,
            r.pairs_above_bound as u64,
            r.staircase_at_stop.map_or(u64::MAX, |c| c as u64),
            u64::from(r.timed_out),
            u64::from(r.oversize),
        ];
        (r.basis.clone(), counters, trace)
    }

    /// A step keeps the products it formed while they fit the cache and
    /// forms the rest again: with no room, with room for some rows only and
    /// with room for all, the run is the same, to the last counter.
    #[test]
    fn the_run_does_not_depend_on_the_product_cache() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(73);
        let limits = [0, 1, 25, 400, 6000, 100_000, usize::MAX];
        let mut big_runs = 0;
        for trial in 0..14 {
            let (ring, input, max_degree) = if trial % 2 == 0 {
                // Two Kummer towers and a biquadratic in their first levels.
                let p = 97u64;
                let t = 3 + trial / 8;
                let (r0, pts) = kummer_block(p, t, 0);
                let (r1, _) = kummer_block(p, t, t);
                let ring = TowerRing {
                    p,
                    rules: [r0, r1].concat(),
                    n_free: 0,
                };
                let x = [pts[rng.gen_range(0..pts.len())].clone(), pts[0].clone()].concat();
                let plant = (trial % 4 == 0).then_some(&x[..]);
                let g = biquadratic(&mut rng, &ring, t, plant);
                (ring, vec![g], 64)
            } else {
                // A general tower, a free variable and several generators.
                let p = 1_000_003u64;
                let ring = random_ring(&mut rng, p, 2, 2 + trial % 3, trial % 3, true);
                let n = ring.n_vars();
                let input: Vec<RPoly> = (0..1 + trial % 3)
                    .map(|_| {
                        let raw: Vec<(Vec<u32>, u64)> = (0..5)
                            .map(|_| {
                                let mut left = 3u32;
                                let e: Vec<u32> = (0..n)
                                    .map(|_| {
                                        if left > 0 && rng.gen_range(0..3) == 0 {
                                            left -= 1;
                                            1
                                        } else {
                                            0
                                        }
                                    })
                                    .collect();
                                (e, rng.gen_range(1..p))
                            })
                            .collect();
                        ring.from_raw(&raw)
                    })
                    .collect();
                (ring, input, 40)
            };
            let opts = TowerF4Options::new(max_degree);
            let want = the_run(&f4_tower_cached(&input, &ring, &opts, usize::MAX));
            for &limit in &limits {
                let got = the_run(&f4_tower_cached(&input, &ring, &opts, limit));
                assert_eq!(got, want, "trial {trial}, cache of {limit} terms");
            }
            let full = f4_tower_cached(&input, &ring, &opts, usize::MAX);
            big_runs += usize::from(full.max_nnz > 400);
        }
        assert!(
            big_runs >= 4,
            "too few runs have steps larger than the small caches: {big_runs}"
        );
    }

    /// Monomials are numbered in order of arrival, once each, and
    /// `order_columns` turns the numbers into the descending order of the
    /// monomials.
    #[test]
    fn monomials_are_numbered_by_arrival_and_ordered_into_columns() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        let mut examined: MonoMap<u32> = MonoMap::default();
        let mut frontier = Vec::new();
        let mut seen: Vec<Mono> = Vec::new();
        for _ in 0..3000 {
            let m = mono(rng.gen_range(0..64u64), [rng.gen_range(0..3), 0, 0, 0]);
            let n = number_of(&mut examined, &mut frontier, m);
            match seen.iter().position(|&x| x == m) {
                Some(k) => assert_eq!(n as usize, k, "a repeat keeps its number"),
                None => {
                    assert_eq!(n as usize, seen.len(), "a new one takes the next");
                    assert_eq!(frontier.last(), Some(&m));
                    seen.push(m);
                }
            }
        }
        assert_eq!(frontier, seen, "the frontier is the new monomials in order");
        assert!(seen.len() > 100);
        let (cols, col_of_number) = order_columns(&examined);
        let mut want = seen.clone();
        want.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(cols, want);
        for (k, m) in seen.iter().enumerate() {
            assert_eq!(cols[col_of_number[k] as usize], *m);
        }
    }

    #[test]
    fn the_product_cache_keeps_what_fits_and_nothing_else() {
        // A product of `k` terms, numbered from `from`, with coefficient `from`.
        let keep = |c: &mut ProductCache, k: usize, from: u32| -> u32 {
            if !c.has_room(k) {
                return NONE;
            }
            let ids: Vec<u32> = (from..from + k as u32).collect();
            c.keep(&ids, &vec![from; k])
        };
        let mut c = ProductCache::new(10);
        assert_eq!(keep(&mut c, 4, 100), 0);
        assert_eq!(keep(&mut c, 7, 200), NONE, "one past the room left");
        assert_eq!(keep(&mut c, 6, 300), 4, "an exact fit");
        assert_eq!(keep(&mut c, 1, 400), NONE, "full");
        assert_eq!(keep(&mut c, 0, 500), 10, "an empty product always fits");
        assert_eq!(c.ids.len(), 10);
        assert_eq!(c.coefs.len(), 10);
        assert_eq!(c.ids[4..10], [300, 301, 302, 303, 304, 305]);
        assert_eq!(c.coefs[..4], [100; 4]);
        assert_eq!(c.coefs[4..], [300; 6]);
        // No limit still leaves `NONE` unused as a position.
        assert_eq!(ProductCache::new(usize::MAX).limit, NONE as usize - 1);
        assert!(!ProductCache::new(0).has_room(1));
    }

    /// Every sparse row ascending, with entries below `p`, over `n_cols`
    /// columns: `lead` first, then each later column with probability
    /// `density`.
    fn random_sparse_row(
        rng: &mut rand::rngs::StdRng,
        p: u64,
        n_cols: usize,
        lead: usize,
        density: f64,
    ) -> Row {
        use rand::Rng;
        let mut row = Row {
            cols: vec![lead as u32],
            vals: vec![rng.gen_range(1..p) as u32],
        };
        for c in lead + 1..n_cols {
            if rng.gen_bool(density) {
                row.cols.push(c as u32);
                row.vals.push(rng.gen_range(1..p) as u32);
            }
        }
        row
    }

    /// [`echelon`] by dense, exact arithmetic: the same rows in the same
    /// order, each reduced at its lead and kept as it stands past a new one.
    /// Returns the pivots and the multiply-adds.
    fn dense_echelon(p: u64, mut rows: Vec<Row>, n_cols: usize) -> (Vec<Row>, u64) {
        let f = fp(p);
        rows.retain(|r| !r.cols.is_empty());
        rows.sort_by_key(|r| (r.cols[0], r.cols.len()));
        let mut pivot_of = vec![NONE; n_cols];
        let mut pivots: Vec<Row> = Vec::new();
        let mut work = 0u64;
        for row in rows {
            let mut v = vec![0u64; n_cols];
            for (&c, &x) in row.cols.iter().zip(&row.vals) {
                v[c as usize] = u64::from(x);
            }
            for c in row.cols[0] as usize..n_cols {
                if v[c] == 0 {
                    continue;
                }
                let q = pivot_of[c];
                if q != NONE {
                    let a = v[c];
                    v[c] = 0;
                    let piv = &pivots[q as usize];
                    for (&cc, &x) in piv.cols[1..].iter().zip(&piv.vals[1..]) {
                        let sub = a * u64::from(x) % p;
                        v[cc as usize] = (v[cc as usize] + p - sub) % p;
                    }
                    work += piv.cols.len() as u64 - 1;
                    continue;
                }
                let inv = f.inv(v[c]);
                let mut out = Row {
                    cols: vec![c as u32],
                    vals: vec![1],
                };
                for (cc, &x) in v.iter().enumerate().skip(c + 1) {
                    if x != 0 {
                        out.cols.push(cc as u32);
                        out.vals
                            .push((u128::from(x) * u128::from(inv) % u128::from(p)) as u32);
                    }
                }
                pivot_of[c] = pivots.len() as u32;
                pivots.push(out);
                break;
            }
        }
        (pivots, work)
    }

    /// A pivot lists exactly the bitset words its tail has columns in, with
    /// the bits of those columns, whatever the order and the density.
    #[test]
    fn a_pivot_lists_the_words_its_tail_touches() {
        use rand::seq::SliceRandom;
        use rand::{Rng, SeedableRng};
        use std::collections::BTreeMap;
        let mut rng = rand::rngs::StdRng::seed_from_u64(83);
        for _ in 0..300 {
            let n = rng.gen_range(1..700usize);
            let density = rng.gen_range(0.01..1.0);
            let mut row = random_sparse_row(&mut rng, 97, n, 0, density);
            if rng.gen_bool(0.3) {
                // The list does not rely on the columns being ascending.
                row.cols[1..].shuffle(&mut rng);
            }
            let mut want: BTreeMap<u32, u64> = BTreeMap::new();
            for &c in &row.cols[1..] {
                *want.entry(c / 64).or_insert(0) |= 1u64 << (c % 64);
            }
            let piv = Pivot::new(row);
            let mut got: BTreeMap<u32, u64> = BTreeMap::new();
            for &(w, b) in &piv.runs {
                *got.entry(w).or_insert(0) |= b;
            }
            assert_eq!(got, want);
            assert!(piv.runs.iter().all(|&(_, b)| b != 0));
        }
    }

    /// The columns at which the accumulator turns lazy: `n_cols · (p − 1)² +
    /// p < 2⁶⁴`, at the last width that holds and the first that does not.
    #[test]
    fn the_lazy_accumulator_guard_is_exact() {
        for &p in &[
            3u64, 97, 786433, 469762049, 998244353, 2013265921, 4294967291,
        ] {
            let f = fp(p);
            let sq = u128::from(p - 1).pow(2);
            let n_max = ((1u128 << 64) - 1 - u128::from(p)) / sq;
            if n_max > 4096 {
                assert!(Acc::new(f, 4096).lazy, "p = {p}");
                continue;
            }
            let n = n_max as usize;
            assert!(n as u128 * sq + u128::from(p) < 1 << 64);
            assert!((n as u128 + 1) * sq + u128::from(p) >= 1 << 64);
            assert!(Acc::new(f, n).lazy, "p = {p}, {n} columns");
            assert!(!Acc::new(f, n + 1).lazy, "p = {p}, {} columns", n + 1);
        }
    }

    /// The largest value a column can reach: every pivot row full of `p − 1`,
    /// every multiplier `p − 1`, at the widest lazy matrix. The sum is held
    /// exactly, and reads as its remainder.
    #[test]
    fn a_lazy_column_at_the_guard_holds_the_exact_sum() {
        for &p in &[786433u64, 469762049, 998244353, 2013265921, 4294967291] {
            let f = fp(p);
            let sq = u128::from(p - 1).pow(2);
            let n = (((1u128 << 64) - 1 - u128::from(p)) / sq).min(700) as usize;
            let mut acc = Acc::new(f, n);
            assert!(acc.lazy, "p = {p}");
            let start = Row {
                cols: (0..n as u32).collect(),
                vals: vec![(p - 1) as u32; n],
            };
            acc.load(&start);
            let mut want = vec![u128::from(p - 1); n];
            for lead in 0..n {
                let row = Row {
                    cols: (lead as u32..n as u32).collect(),
                    vals: vec![(p - 1) as u32; n - lead],
                };
                let mut word = 0u64;
                acc.sub_tail(f, 1, &Pivot::new(row), lead / 64, &mut word);
                for w in &mut want[lead + 1..] {
                    *w += sq;
                }
            }
            let top = *want.iter().max().unwrap();
            assert!(top < 1 << 64, "p = {p}: the guard let a sum pass 2^64");
            for (c, &w) in want.iter().enumerate() {
                assert_eq!(u128::from(acc.vals[c]), w, "p = {p}, column {c}");
                assert_eq!(
                    u128::from(acc.take(f, c)),
                    w % u128::from(p),
                    "p = {p}, column {c}"
                );
            }
        }
    }

    /// A column that holds a multiple of `p` reads as zero: not pushed out
    /// by a drain, and not a pivot candidate for a scan.
    #[test]
    fn a_multiple_of_p_in_a_lazy_column_reads_as_zero() {
        let p = 786433u64;
        let f = fp(p);
        let mut acc = Acc::new(f, 3);
        assert!(acc.lazy);
        // Column 1 starts at 5; subtracting 7 · row[1] with v chosen so that
        // 5 + (p − 7)·v ≡ 0 leaves the raw value 5 + (p − 7)·v, a multiple of
        // p that is not zero.
        let v = f.mul(p - 5, f.inv(p - 7));
        let lead = Row {
            cols: vec![0, 1, 2],
            vals: vec![1, v as u32, 3],
        };
        acc.load(&Row {
            cols: vec![0, 1],
            vals: vec![7, 5],
        });
        let mut word = 0;
        acc.take(f, 0);
        acc.bits[0] = 0;
        acc.sub_tail(f, 7, &Pivot::new(lead), 0, &mut word);
        assert_ne!(acc.vals[1], 0, "the raw value is a multiple of p, not 0");
        assert_eq!(acc.vals[1] % p, 0);
        assert_eq!(word, 0b110);
        acc.bits[0] = word;
        let mut out = Row::default();
        acc.drain_from(f, 0, 1, &mut out);
        assert_eq!(out.cols, vec![2], "column 1 reads as zero");
        assert_eq!(out.vals, vec![f.mul(3, p - 7) as u32]);
        assert!(acc.vals.iter().all(|&x| x == 0) && acc.bits.iter().all(|&b| b == 0));
    }

    /// [`Acc::clear_from`] leaves the accumulator as a drain does.
    #[test]
    fn clearing_the_accumulator_leaves_it_clean() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(81);
        for &p in &[97u64, 4294967291] {
            let f = fp(p);
            for _ in 0..30 {
                let n = rng.gen_range(1..400usize);
                let mut acc = Acc::new(f, n);
                let row = random_sparse_row(&mut rng, p, n, 0, 0.3);
                acc.load(&row);
                for _ in 0..4 {
                    let lead = rng.gen_range(0..n);
                    let pivot = random_sparse_row(&mut rng, p, n, lead, 0.3);
                    let mut word = 0;
                    acc.sub_tail(
                        f,
                        rng.gen_range(1..p),
                        &Pivot::new(pivot),
                        lead / 64,
                        &mut word,
                    );
                    acc.bits[lead / 64] |= word;
                }
                acc.clear_from(0);
                assert!(acc.vals.iter().all(|&x| x == 0), "p = {p}");
                assert!(acc.bits.iter().all(|&b| b == 0), "p = {p}");
            }
        }
    }

    /// The accumulator, lazy or not, echelonizes as exact arithmetic does:
    /// the same pivots, entry for entry, and the same count of multiply-adds,
    /// at widths on both sides of the guard.
    #[test]
    fn an_accumulator_echelon_matches_dense_exact_arithmetic() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(97);
        let (mut lazy_runs, mut eager_runs, mut cancelled) = (0, 0, 0u64);
        for (p, widths) in [
            (3u64, vec![150usize, 330]),
            (5, vec![200]),
            (97, vec![64, 301]),
            (786433, vec![190, 400]),
            (469762049, vec![83, 84, 130]),
            (998244353, vec![18, 19, 70]),
            (2013265921, vec![4, 5, 90]),
            (4294967291, vec![1, 2, 60]),
        ] {
            let f = fp(p);
            for &n_cols in &widths {
                for round in 0..5 {
                    let density = [0.9, 0.5, 0.2, 0.05, 1.0][round];
                    let n_rows = rng.gen_range(n_cols / 2 + 1..2 * n_cols + 2);
                    let mut rows: Vec<Row> = Vec::new();
                    for _ in 0..n_rows {
                        let lead = rng.gen_range(0..n_cols);
                        rows.push(random_sparse_row(&mut rng, p, n_cols, lead, density));
                    }
                    // Rows that are sums of earlier ones: cancellations.
                    for _ in 0..n_rows / 4 {
                        let (a, b) = (rng.gen_range(0..rows.len()), rng.gen_range(0..rows.len()));
                        let mut v = vec![0u64; n_cols];
                        for r in [&rows[a], &rows[b]] {
                            for (&c, &x) in r.cols.iter().zip(&r.vals) {
                                v[c as usize] = (v[c as usize] + u64::from(x)) % p;
                            }
                        }
                        let sum = Row {
                            cols: (0..n_cols as u32).filter(|&c| v[c as usize] != 0).collect(),
                            vals: v.iter().filter(|&&x| x != 0).map(|&x| x as u32).collect(),
                        };
                        if !sum.cols.is_empty() {
                            rows.push(sum);
                        }
                    }
                    let (want, want_work) = dense_echelon(p, rows.clone(), n_cols);
                    let mut work = 0u64;
                    let got = echelon(f, rows, n_cols, &mut work, None).unwrap();
                    assert_eq!(work, want_work, "p = {p}, {n_cols} columns, round {round}");
                    assert_eq!(
                        got.len(),
                        want.len(),
                        "p = {p}, {n_cols} columns, round {round}"
                    );
                    for (g, w) in got.iter().zip(&want) {
                        assert_eq!(g.cols, w.cols, "p = {p}, {n_cols} columns, round {round}");
                        assert_eq!(g.vals, w.vals, "p = {p}, {n_cols} columns, round {round}");
                    }
                    cancelled += (n_rows as u64).saturating_sub(want.len() as u64);
                    if Acc::new(f, n_cols).lazy {
                        lazy_runs += 1;
                    } else {
                        eager_runs += 1;
                    }
                }
            }
        }
        assert!(lazy_runs > 20 && eager_runs > 20 && cancelled > 100);
    }

    #[test]
    fn compiled_refutation_and_basis_specialization_sound_and_fast() {
        // Kummer tower of height t = 3 for m = 2 summands over p = 65537 (2^16 + 1).
        let p = 65_537u64;
        let t = 3usize;
        let mut rules = Vec::new();
        for _ in 0..2 {
            let base = rules.len();
            for j in 0..(t - 1) {
                rules.push(SquareRule {
                    next: Some((base + j + 1) as u8),
                    a: 0,
                    b: 1,
                    c: 0,
                    d: 0,
                });
            }
            rules.push(SquareRule {
                next: None,
                a: 0,
                b: 0,
                c: 0,
                d: 1,
            });
        }
        let ring = TowerRing::new(p, rules, 0);
        // S_3(3·y_0, 3·y_3, xr) on y^2 = x^3 + 7x + 11.
        let (a, b, g) = (7u64, 11u64, 3u64);
        let build_sys = |xr: u64| -> Vec<RPoly> {
            // (g y0 - g y3)^2 xr^2 - 2[(g y0 + g y3)(g^2 y0 y3 + a) + 2b] xr + (g^2 y0 y3 - a)^2 - 4b(g y0 + g y3)
            let n = ring.n_vars();
            let var = |i: usize, c: u64| {
                let mut e = vec![0u32; n];
                e[i] = 1;
                vec![(e, c)]
            };
            let x1 = ring.from_raw(&var(0, g));
            let x2 = ring.from_raw(&var(t, g));
            let mut terms = Vec::new();
            // Sample on exact symbolic evaluation through from_raw:
            let g2 = g * g % p;
            let xr2 = xr * xr % p;
            let z = vec![0u32; n];
            let mut e0 = z.clone();
            e0[0] = 1;
            let mut e1 = z.clone();
            e1[t] = 1;
            let mut e01 = z.clone();
            e01[0] = 1;
            e01[t] = 1;
            let mut e0_2 = z.clone();
            e0_2[0] = 2;
            let mut e1_2 = z.clone();
            e1_2[t] = 2;
            let mut e0_2_1 = z.clone();
            e0_2_1[0] = 2;
            e0_2_1[t] = 1;
            let mut e0_1_2 = z.clone();
            e0_1_2[0] = 1;
            e0_1_2[t] = 2;
            let mut e01_2 = z.clone();
            e01_2[0] = 2;
            e01_2[t] = 2;

            // (x1 - x2)^2 xr^2 = g^2 xr^2 (y0^2 - 2 y0 y3 + y3^2)
            let c_sq = g2 * xr2 % p;
            terms.push((e0_2, c_sq));
            terms.push((e1_2, c_sq));
            terms.push((e01.clone(), (p - 2 * c_sq % p) % p));

            // -2[(g(y0+y3))(g^2 y0 y3 + a) + 2b] xr
            let neg2xr = (p - 2 * xr % p) % p;
            let cg3 = neg2xr * (g * g2 % p) % p;
            terms.push((e0_2_1, cg3));
            terms.push((e0_1_2, cg3));
            let cga = neg2xr * (g * a % p) % p;
            terms.push((e0.clone(), cga));
            terms.push((e1.clone(), cga));
            terms.push((z.clone(), neg2xr * (2 * b % p) % p));

            // (g^2 y0 y3 - a)^2 - 4b g(y0 + y3)
            terms.push((e01_2, g2 * g2 % p));
            terms.push((e01, (p - 2 * g2 % p * a % p) % p));
            terms.push((z, a * a % p));
            let neg4bg = (p - 4 * b % p * g % p) % p;
            terms.push((e0, neg4bg));
            terms.push((e1, neg4bg));
            let _ = (x1, x2);
            vec![ring.from_raw(&terms)]
        };

        // Calibration on xr = 12345 (generic UNSAT).
        let mut learn_opts = TowerF4Options::new(12);
        learn_opts.record_trace = true;
        let cal = f4_tower(&build_sys(12_345), &ring, &learn_opts);
        assert!(cal.inconsistent);
        let trace = Arc::new(cal.learned.expect("learned trace"));
        assert!(trace.compiled.is_some(), "compiled minimal-1-witness DAG");

        // Replay across 32 random targets: every replay refutation MUST agree with plain F4,
        // and replay muladds must be strictly smaller than plain F4 muladds.
        let mut replay_opts = TowerF4Options::new(12);
        replay_opts.replay = Some(trace);
        let mut refuted = 0usize;
        for xr in 100..132u64 {
            let sys = build_sys(xr);
            let plain = f4_tower(&sys, &ring, &TowerF4Options::new(12));
            let fast = f4_tower(&sys, &ring, &replay_opts);
            if fast.inconsistent {
                assert!(plain.inconsistent, "soundness violated at xr = {xr}");
                assert!(
                    fast.muladds <= plain.muladds,
                    "fast muladds {} > plain {}",
                    fast.muladds,
                    plain.muladds
                );
                refuted += 1;
            }
        }
        assert!(refuted >= 28, "Schwartz-Zippel completeness: {refuted}/32");

        // Test parent-to-child specialization helpers (`specialize_basis_at_var` & `specialize_drop_top_var`).
        let sys_sat_probe = build_sys(500);
        let mut low_opts = TowerF4Options::new(3);
        low_opts.max_degree = 3;
        let parent_rep = f4_tower(&sys_sat_probe, &ring, &low_opts);
        if !parent_rep.inconsistent {
            let top_var = t - 1;
            for root in [1u64, p - 1] {
                let (child_in, _) = ring.specialize_basis_at_var(
                    &parent_rep.basis,
                    top_var,
                    root,
                    TowerWeights::UNIT,
                );
                let r_keep = f4_tower(&child_in, &ring, &TowerF4Options::new(12));
                let (c_ring, c_basis, _, _) = ring.specialize_drop_top_var(
                    &parent_rep.basis,
                    top_var,
                    root,
                    TowerWeights::UNIT,
                );
                let r_drop = f4_tower(&c_basis, &c_ring, &TowerF4Options::new(12));
                assert_eq!(r_keep.inconsistent, r_drop.inconsistent);
            }
        }
    }
}
