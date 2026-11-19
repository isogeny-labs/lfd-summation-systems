//! # The Huang–Kosters–Yeo closure over a PKM tower quotient ring.
//!
//! `research/prime_fall_degree_20261027/PREREGISTRATION.md` asks for the
//! first and last fall degree of prime-field Semaev decomposition systems on
//! Petit–Kosters–Messeng tower factor bases. This module computes them
//! exactly, by the definition, with no Gröbner basis.
//!
//! The ring is the one [`f4_fp_tower`](crate::cryptanalysis::f4_fp_tower)
//! works in, `R = F_p[y, u] / (y_j² − a·y_j·y_{j+1} − b·y_{j+1} − c·y_j − d)`.
//! The tower relations are a Gröbner basis with leading monomials `y_j²`
//! (note Lemma 1), so every element has a normal form square-free in the
//! tower variables, and rewriting a product into it never raises the
//! degree. The chain's free unknowns `u_k` stay ordinary variables.
//!
//! Definitions (`research/notes/index-calculus/RESEARCH_FALL_DEGREE_BOUNDS.md`
//! §0, the binary program's convention, with the tower in place of the
//! field equations):
//! - `V_{F,c}` is the smallest `F_p`-subspace of `R_{≤c}` that contains
//!   `F ∩ R_{≤c}` and, with every `g` of degree `< c`, every `x·g` for a
//!   variable `x`. Products are formed only when `deg g + 1 ≤ c`, then
//!   rewritten. Rewriting is degree-non-increasing and stays inside the
//!   closure of `F ∪ towers`, so this is HKY's `V_{F ∪ T, c}` read in `R`
//!   (the boolean argument of the ledger §0, word for word).
//! - The **last fall degree** of an inconsistent system is the least `c`
//!   with `1 ∈ V_{F,c}`.
//! - The **operational first fall degree** is the least `c` with
//!   `V_{F,c} ∩ R_{≤c−1} ≠ V_{F,c−1}`.
//!
//! The closure is a work queue (the binary engine's,
//! `examples/hky_last_fall_degree.rs`): every row that enters the echelon
//! with degree `< c` has its multiples by every variable queued. Columns are
//! the monomials of degree `≤ c` sorted by degree, highest first, so a
//! reduced row's first column is its leading monomial and its degree, and
//! `V_{F,c} ∩ R_{≤d}` is spanned by the echelon rows that lead at degree
//! `≤ d`. Rows are fully reduced against the pivots before them.
//!
//! Limits: at most 56 tower variables, 8 free variables with exponents below
//! 256, `p < 2³²`.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cryptanalysis::f4_fp_tower::SquareRule;

/// Tower variables a monomial can carry.
pub const MAX_TOWER: usize = 56;
/// Free variables a monomial can carry.
pub const MAX_FREE: usize = 8;

const FREE_SHIFT: u32 = 56;
const DEG_SHIFT: u32 = 120;

type Key = u128;

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

/// The ring: tower rules (one per tower variable, as in `f4_fp_tower`) and
/// the number of free variables, which follow the tower variables.
#[derive(Clone, Debug)]
pub struct HkyRing {
    pub p: u64,
    pub rules: Vec<SquareRule>,
    pub n_free: usize,
}

/// A polynomial of `R` in normal form, as `(monomial, coefficient)` with
/// the monomial as tower mask plus free exponents.
pub type NfPoly = Vec<(u64, [u8; MAX_FREE], u64)>;

impl HkyRing {
    pub fn new(p: u64, rules: Vec<SquareRule>, n_free: usize) -> Self {
        assert!(p < 1 << 32, "hky_fp_tower: p must be below 2^32");
        assert!(
            rules.len() <= MAX_TOWER,
            "hky_fp_tower: at most {MAX_TOWER} tower variables"
        );
        assert!(
            n_free <= MAX_FREE,
            "hky_fp_tower: at most {MAX_FREE} free variables"
        );
        for (j, r) in rules.iter().enumerate() {
            if let Some(nx) = r.next {
                assert!((nx as usize) > j && (nx as usize) < rules.len());
            } else {
                assert!(r.a == 0 && r.b == 0);
            }
        }
        HkyRing { p, rules, n_free }
    }

    pub fn n_tower(&self) -> usize {
        self.rules.len()
    }

    pub fn n_vars(&self) -> usize {
        self.rules.len() + self.n_free
    }

    /// `y_j · (mask, coef)` in normal form (`f4_fp_tower::mul_var_into`).
    fn mul_tower_var(&self, m: u64, j: usize, coef: u64, out: &mut Vec<(u64, u64)>) {
        let p = self.p;
        let bit = 1u64 << j;
        if m & bit == 0 {
            out.push((m | bit, coef));
            return;
        }
        let r = self.rules[j];
        let m0 = m & !bit;
        if r.d != 0 {
            out.push((m0, mulm(coef, r.d, p)));
        }
        if r.c != 0 {
            out.push((m, mulm(coef, r.c, p)));
        }
        if let Some(nx) = r.next {
            if r.b != 0 {
                self.mul_tower_var(m0, nx as usize, mulm(coef, r.b, p), out);
            }
            if r.a != 0 {
                self.mul_tower_var(m, nx as usize, mulm(coef, r.a, p), out);
            }
        }
    }

    fn merge(&self, terms: &mut Vec<(u64, u64)>) {
        terms.sort_unstable_by_key(|t| t.0);
        let mut out = 0;
        let mut i = 0;
        while i < terms.len() {
            let m = terms[i].0;
            let mut c = 0;
            while i < terms.len() && terms[i].0 == m {
                c = (c + terms[i].1) % self.p;
                i += 1;
            }
            if c != 0 {
                terms[out] = (m, c);
                out += 1;
            }
        }
        terms.truncate(out);
    }

    /// The normal form of a polynomial given by exponent vectors of length
    /// [`n_vars`](Self::n_vars) (tower variables first).
    pub fn from_raw(&self, terms: &[(Vec<u32>, u64)]) -> NfPoly {
        let nt = self.n_tower();
        let mut acc: HashMap<(u64, [u8; MAX_FREE]), u64> = HashMap::new();
        let (mut cur, mut next) = (Vec::new(), Vec::new());
        for (e, c) in terms {
            assert_eq!(e.len(), self.n_vars(), "hky_fp_tower: exponent length");
            let c = c % self.p;
            if c == 0 {
                continue;
            }
            let mut f = [0u8; MAX_FREE];
            for k in 0..self.n_free {
                f[k] = u8::try_from(e[nt + k]).expect("hky_fp_tower: free exponent above 255");
            }
            cur.clear();
            cur.push((0u64, c));
            for (j, &ej) in e.iter().take(nt).enumerate() {
                for _ in 0..ej {
                    next.clear();
                    for &(m, cc) in &cur {
                        self.mul_tower_var(m, j, cc, &mut next);
                    }
                    self.merge(&mut next);
                    std::mem::swap(&mut cur, &mut next);
                }
            }
            for &(m, cc) in &cur {
                let v = acc.entry((m, f)).or_insert(0);
                *v = (*v + cc) % self.p;
            }
        }
        let mut out: NfPoly = acc
            .into_iter()
            .filter(|&(_, c)| c != 0)
            .map(|((m, f), c)| (m, f, c))
            .collect();
        out.sort_unstable_by_key(|&(m, f, _)| Reverse(key(m, &f)));
        out
    }

    /// The degree of a normal form.
    pub fn degree(f: &NfPoly) -> u32 {
        f.iter()
            .map(|(m, e, _)| key_degree(key(*m, e)))
            .max()
            .unwrap_or(0)
    }

    /// Evaluate a normal form at a point (one value per variable).
    pub fn eval(&self, f: &NfPoly, x: &[u64]) -> u64 {
        let p = self.p;
        let nt = self.n_tower();
        let mut acc = 0;
        for (m, e, c) in f {
            let mut t = *c;
            let mut bits = *m;
            while bits != 0 {
                let j = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                t = mulm(t, x[j] % p, p);
            }
            for k in 0..self.n_free {
                for _ in 0..e[k] {
                    t = mulm(t, x[nt + k] % p, p);
                }
            }
            acc = (acc + t) % p;
        }
        acc
    }
}

/// The columns of `R_{≤c}`: every normal-form monomial of degree `≤ c`,
/// highest degree first.
struct Columns {
    keys: Vec<Key>,
    index: HashMap<Key, u32>,
    /// `degree_start[d]` is the first column of degree `d` (columns of
    /// degree `d` are `degree_start[d] .. degree_start[d − 1]`).
    deg: Vec<u8>,
}

/// The number of normal-form monomials of degree `≤ c`.
pub fn count_monomials(n_tower: usize, n_free: usize, c: u32) -> u128 {
    let binom = |n: u128, k: u128| -> u128 {
        if k > n {
            return 0;
        }
        let mut r = 1u128;
        for i in 0..k {
            r = r * (n - i) / (i + 1);
        }
        r
    };
    let mut total = 0u128;
    for k in 0..=c.min(n_tower as u32) {
        // free monomials in n_free variables of degree ≤ c − k
        let rest = u128::from(c - k);
        let free = if n_free == 0 {
            1
        } else {
            binom(n_free as u128 + rest, rest)
        };
        total += binom(n_tower as u128, u128::from(k)) * free;
    }
    total
}

impl Columns {
    fn build(ring: &HkyRing, c: u32) -> Columns {
        let nt = ring.n_tower();
        let nf = ring.n_free;
        let mut keys = Vec::new();
        // tower masks by popcount
        let mut masks: Vec<Vec<u64>> = vec![Vec::new(); c as usize + 1];
        fn rec_mask(
            start: usize,
            nt: usize,
            left: u32,
            cur: u64,
            out: &mut Vec<Vec<u64>>,
            depth: u32,
        ) {
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
        let index = keys
            .iter()
            .enumerate()
            .map(|(i, &k)| (k, i as u32))
            .collect();
        let deg = keys.iter().map(|&k| key_degree(k) as u8).collect();
        Columns { keys, index, deg }
    }
}

/// A sparse row over [`Columns`]: increasing column (decreasing monomial),
/// non-zero coefficients.
type Row = Vec<(u32, u32)>;

/// Limits on one closure.
#[derive(Clone, Debug)]
pub struct HkyLimits {
    pub max_cols: u128,
    /// Total stored row entries in the echelon.
    pub max_entries: u64,
    pub budget: Duration,
    pub batch: usize,
}

impl Default for HkyLimits {
    fn default() -> Self {
        HkyLimits {
            max_cols: 2_000_000,
            max_entries: 600_000_000,
            budget: Duration::from_secs(3600),
            batch: 4096,
        }
    }
}

/// One closure `V_{F,c}`.
#[derive(Clone, Debug)]
pub struct Closure {
    pub c: u32,
    pub cols: u128,
    /// `dims[d] = dim (V_{F,c} ∩ R_{≤d})`, `d = 0..=c`; empty if censored.
    pub dims: Vec<u64>,
    pub contains_one: bool,
    /// Stopped by a limit before the fixpoint (and before `1`).
    pub censored: Option<String>,
    pub rows_reduced: u64,
    pub entries: u64,
    /// `F_p` multiply-adds: one per entry of a pivot row subtracted during
    /// reduction, plus one per term of a variable product before rewriting
    /// merges it (the unit `f4_fp_tower` reports as `muladds`).
    pub muladds: u64,
    pub ms: f64,
}

impl Closure {
    pub fn dim(&self) -> u64 {
        self.dims.last().copied().unwrap_or(0)
    }
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

struct Echelon {
    p: u64,
    /// `pivot[col]` = index into `rows` of the row leading at `col`.
    pivot: Vec<u32>,
    rows: Vec<Row>,
    entries: u64,
}

const NONE: u32 = u32::MAX;

impl Echelon {
    /// Reduce `r` by every pivot (full reduction) and make it monic; also
    /// the multiply-adds it took.
    fn reduce(&self, r: &Row, acc: &mut [u64], touched: &mut Vec<u32>) -> (Row, u64) {
        let mut ops = 0u64;
        let p = self.p;
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
            let v = acc[col as usize];
            if v == 0 {
                continue;
            }
            let piv = self.pivot[col as usize];
            if piv == NONE {
                out.push((col, v as u32));
                acc[col as usize] = 0;
                continue;
            }
            // acc -= v · row (rows are monic)
            let f = p - v;
            ops += self.rows[piv as usize].len() as u64;
            for &(c2, x) in &self.rows[piv as usize] {
                let i = c2 as usize;
                if acc[i] == 0 && c2 != col {
                    touched.push(c2);
                    heap.push(Reverse(c2));
                }
                acc[i] = (acc[i] + f * u64::from(x)) % p;
            }
            debug_assert_eq!(acc[col as usize], 0);
        }
        for &t in touched.iter() {
            acc[t as usize] = 0;
        }
        touched.clear();
        if let Some(&(_, lead)) = out.first() {
            let s = inv(u64::from(lead), p);
            for e in out.iter_mut() {
                e.1 = mulm(u64::from(e.1), s, p) as u32;
            }
            ops += out.len() as u64;
        }
        (out, ops)
    }
}

/// `V_{F,c}` for the normal-form generators `gens`.
pub fn closure(ring: &HkyRing, gens: &[NfPoly], c: u32, limits: &HkyLimits) -> Closure {
    let started = Instant::now();
    let p = ring.p;
    let n_cols = count_monomials(ring.n_tower(), ring.n_free, c);
    let censored = |why: &str, rows: u64, entries: u64, muladds: u64| Closure {
        c,
        cols: n_cols,
        dims: Vec::new(),
        contains_one: false,
        censored: Some(why.to_string()),
        rows_reduced: rows,
        entries,
        muladds,
        ms: started.elapsed().as_secs_f64() * 1e3,
    };
    if n_cols > limits.max_cols {
        return censored("max_cols", 0, 0, 0);
    }
    let cols = Columns::build(ring, c);
    assert_eq!(
        cols.keys.len() as u128,
        n_cols,
        "hky_fp_tower: column count"
    );
    let one_col = (cols.keys.len() - 1) as u32;
    let to_row = |f: &[(u64, [u8; MAX_FREE], u64)]| -> Row {
        let mut r: Row = f
            .iter()
            .map(|(m, e, x)| (cols.index[&key(*m, e)], *x as u32))
            .collect();
        r.sort_unstable_by_key(|t| t.0);
        r
    };
    let nv = ring.n_vars();
    let nt = ring.n_tower();
    // x_v · row, in normal form, for deg row < c.
    let multiply = |r: &Row, v: usize| -> Row {
        let mut terms: Vec<(Key, u64)> = Vec::with_capacity(r.len() * 2);
        let mut tmp = Vec::new();
        for &(col, x) in r {
            let k = cols.keys[col as usize];
            let (m, mut f) = (key_mask(k), key_free(k));
            if v >= nt {
                f[v - nt] += 1;
                terms.push((key(m, &f), u64::from(x)));
            } else {
                tmp.clear();
                ring.mul_tower_var(m, v, u64::from(x), &mut tmp);
                for &(mm, cc) in &tmp {
                    terms.push((key(mm, &f), cc));
                }
            }
        }
        terms.sort_unstable_by_key(|t| Reverse(t.0));
        let mut out: Row = Vec::with_capacity(terms.len());
        let mut i = 0;
        while i < terms.len() {
            let k = terms[i].0;
            let mut s = 0;
            while i < terms.len() && terms[i].0 == k {
                s = (s + terms[i].1) % p;
                i += 1;
            }
            if s != 0 {
                out.push((cols.index[&k], s as u32));
            }
        }
        out
    };
    let mut ech = Echelon {
        p,
        pivot: vec![NONE; cols.keys.len()],
        rows: Vec::new(),
        entries: 0,
    };
    let mut queue: Vec<Row> = gens
        .iter()
        .filter(|g| !g.is_empty() && HkyRing::degree(g) <= c)
        .map(|g| to_row(g))
        .collect();
    let mut rows_reduced = 0u64;
    let mut muladds = 0u64;
    let mut contains_one = false;
    let n = cols.keys.len();
    'outer: while !queue.is_empty() {
        if started.elapsed() > limits.budget {
            return censored("budget", rows_reduced, ech.entries, muladds);
        }
        let take = queue.len().min(limits.batch);
        let batch: Vec<Row> = queue.drain(queue.len() - take..).collect();
        // Reduce the batch in parallel against the echelon as it stands.
        let reduced: Vec<(Row, u64)> = batch
            .par_iter()
            .map_init(
                || (vec![0u64; n], Vec::new()),
                |(acc, touched), r| ech.reduce(r, acc, touched),
            )
            .collect();
        rows_reduced += take as u64;
        let mut acc = vec![0u64; n];
        let mut touched = Vec::new();
        for (r, ops) in reduced {
            muladds += ops;
            if r.is_empty() {
                continue;
            }
            // Against the pivots this batch has added.
            let (r, ops) = ech.reduce(&r, &mut acc, &mut touched);
            muladds += ops;
            let Some(&(lead, _)) = r.first() else {
                continue;
            };
            let d = cols.deg[lead as usize] as u32;
            if d < c {
                for v in 0..nv {
                    let prod = multiply(&r, v);
                    muladds += prod.len() as u64;
                    queue.push(prod);
                }
            }
            ech.entries += r.len() as u64;
            ech.pivot[lead as usize] = ech.rows.len() as u32;
            ech.rows.push(r);
            if lead == one_col {
                contains_one = true;
                break 'outer;
            }
            if ech.entries > limits.max_entries {
                return censored("max_entries", rows_reduced, ech.entries, muladds);
            }
        }
    }
    let mut dims = vec![0u64; c as usize + 1];
    for (col, &pv) in ech.pivot.iter().enumerate() {
        if pv != NONE {
            dims[cols.deg[col] as usize] += 1;
        }
    }
    for d in 1..dims.len() {
        dims[d] += dims[d - 1];
    }
    if contains_one {
        // 1 ∈ V gives V = R_{≤c}: every monomial is a multiple of 1 by
        // variables, one at a time, each product of degree ≤ c.
        let mut full = vec![0u64; c as usize + 1];
        for &dg in &cols.deg {
            full[dg as usize] += 1;
        }
        for d in 1..full.len() {
            full[d] += full[d - 1];
        }
        dims = full;
    }
    Closure {
        c,
        cols: n_cols,
        dims,
        contains_one,
        censored: None,
        rows_reduced,
        entries: ech.entries,
        muladds,
        ms: started.elapsed().as_secs_f64() * 1e3,
    }
}

/// The closures `V_{F,c}` for `c = c_min ..= c_max`, stopping after the
/// first that contains `1` or is censored, with the two fall degrees read
/// off them.
#[derive(Clone, Debug)]
pub struct FallDegrees {
    pub closures: Vec<Closure>,
    /// Least `c` with `1 ∈ V_{F,c}`, if reached.
    pub last_fall: Option<u32>,
    /// Least `c` with `V_{F,c} ∩ R_{≤c−1} ≠ V_{F,c−1}`, among the closures
    /// computed (each fall degree compares two consecutive closures).
    pub first_fall: Option<u32>,
    /// Every degree `c` at which a fall occurred.
    pub falls: Vec<u32>,
}

pub fn fall_degrees(
    ring: &HkyRing,
    gens: &[NfPoly],
    c_min: u32,
    c_max: u32,
    limits: &HkyLimits,
) -> FallDegrees {
    fall_degrees_with(ring, gens, c_min, c_max, limits, |_| {})
}

/// [`fall_degrees`], calling `on_closure` with each closure `V_{F,c}` as soon
/// as it is computed, before the next degree starts. A caller that records
/// these keeps every decided degree if the process is stopped during a later
/// closure; a later run can resume at the first undecided degree with
/// `c_min` set to it (each closure is computed independently).
pub fn fall_degrees_with<F: FnMut(&Closure)>(
    ring: &HkyRing,
    gens: &[NfPoly],
    c_min: u32,
    c_max: u32,
    limits: &HkyLimits,
    mut on_closure: F,
) -> FallDegrees {
    let mut closures: Vec<Closure> = Vec::new();
    let mut falls = Vec::new();
    let mut last_fall = None;
    for c in c_min..=c_max {
        let cl = closure(ring, gens, c, limits);
        let stop = cl.contains_one || cl.censored.is_some();
        if let Some(prev) = closures.last() {
            if cl.censored.is_none()
                && prev.censored.is_none()
                && cl.dims[c as usize - 1] > prev.dim()
            {
                falls.push(c);
            }
        } else if cl.censored.is_none() && c > 0 && cl.dims[c as usize - 1] > 0 && c == c_min {
            // No V_{F,c−1} computed: a fall at c_min is recorded only when
            // c_min − 1 is below every generator's degree, so V_{F,c−1} = 0.
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
        on_closure(&cl);
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

#[cfg(test)]
mod tests {
    use super::*;

    // x² − 1 = 0, x² − 4 = 0 and x − 1 = 0 over F_101 have no common root:
    // the closure reaches 1 only once the product x·(x − 1) is admitted.
    fn inconsistent() -> (HkyRing, Vec<NfPoly>) {
        let p = 101;
        let ring = HkyRing::new(p, Vec::new(), 1);
        let gens = vec![
            ring.from_raw(&[(vec![2], 1), (vec![0], p - 1)]),
            ring.from_raw(&[(vec![2], 1), (vec![0], p - 4)]),
            ring.from_raw(&[(vec![1], 1), (vec![0], p - 1)]),
        ];
        (ring, gens)
    }

    #[test]
    fn fall_degrees_with_reports_every_closure_before_returning() {
        let (ring, gens) = inconsistent();
        let limits = HkyLimits::default();
        let plain = fall_degrees(&ring, &gens, 1, 6, &limits);
        let mut seen: Vec<(u32, bool, Vec<u64>)> = Vec::new();
        let with = fall_degrees_with(&ring, &gens, 1, 6, &limits, |cl| {
            seen.push((cl.c, cl.contains_one, cl.dims.clone()));
        });
        assert!(plain.last_fall.is_some(), "the system is inconsistent");
        assert_eq!(with.last_fall, plain.last_fall);
        assert_eq!(with.falls, plain.falls);
        let expect: Vec<(u32, bool, Vec<u64>)> = plain
            .closures
            .iter()
            .map(|cl| (cl.c, cl.contains_one, cl.dims.clone()))
            .collect();
        assert_eq!(seen, expect);
    }

    #[test]
    fn resuming_at_a_later_degree_gives_the_same_last_fall() {
        let (ring, gens) = inconsistent();
        let limits = HkyLimits::default();
        let full = fall_degrees(&ring, &gens, 1, 6, &limits);
        let d = full.last_fall.expect("inconsistent");
        let resumed = fall_degrees(&ring, &gens, d, 6, &limits);
        assert_eq!(resumed.last_fall, Some(d));
    }
}
