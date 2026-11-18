//! Linearized root finding for point decomposition on y^2 + xy = x^3 + a x^2 + b over F_(2^n).
//!
//! P1 + P2 = T iff P1, P2, -T lie on a line y = lam (x + xT) + yT + xT.  Then
//!     e1 = x1 + x2 = lam^2 + lam + a + xT,      e2 = x1 x2 = xT lam^2 + ((yT + xT)^2 + b) / xT,
//! both F_2-affine in lam.  With e1 = sum t_i v_i in V, lam = Q(e1 + c1) + eps (Q a linear right inverse of
//! z -> z^2 + z, valid iff Tr(e1 + c1) = 0) and W = Q + id, the condition e2 in V2 = span(V V) is
//!     pi(xT W(e1)) + eps pi(xT) = pi(xT W(c1) + c2),      sum t_i Tr(v_i) = Tr(c1),
//! a linear system in (t, eps): k + 1 unknowns, n - dim V2 + 1 equations (pi: F_2^n -> F_2^n / V2).
//! Each solution gives x^2 + e1 x + e2, kept when its roots are in V.
//!
//! Three summands: for every abscissa x3 in V (walked in Gray-code order by index, lifted on the fly, both
//! signs) solve the two-summand problem for T = R - P3.  Nothing is stored per factor-base point, so a range
//! of indices of a factor base far too large to list (k = 44) can be scanned.

use crate::field::{Fe, Field};

pub type Pt = (Fe, Fe);

pub struct LMap<T> {
    t: Vec<[T; 256]>,
}

impl<T: Copy + Default + std::ops::BitXor<Output = T>> LMap<T> {
    pub fn new(images: &[T], n: usize) -> Self {
        let nb = (n + 7) / 8;
        let mut t = vec![[T::default(); 256]; nb];
        for j in 0..nb {
            for x in 1usize..256 {
                let bit = 8 * j + x.trailing_zeros() as usize;
                let im = if bit < n { images[bit] } else { T::default() };
                t[j][x] = t[j][x & (x - 1)] ^ im;
            }
        }
        LMap { t }
    }
    #[inline(always)]
    pub fn at(&self, x: Fe) -> T {
        let mut r = T::default();
        for (j, tab) in self.t.iter().enumerate() {
            r = r ^ tab[((x.0[j >> 3] >> ((j & 7) * 8)) & 255) as usize];
        }
        r
    }
}

/// Reduced xor basis as (pivot bit, vector).
fn echelon(vectors: impl IntoIterator<Item = Fe>) -> Vec<(usize, Fe)> {
    let mut piv: Vec<(usize, Fe)> = Vec::new();
    for mut v in vectors {
        for &(pb, pv) in &piv {
            if v.bit(pb) {
                v ^= pv;
            }
        }
        if let Some(pb) = v.top() {
            for q in piv.iter_mut() {
                if q.1.bit(pb) {
                    q.1 ^= v;
                }
            }
            piv.push((pb, v));
        }
    }
    piv
}

pub fn independent(basis: &[Fe]) -> bool {
    echelon(basis.iter().copied()).len() == basis.len()
}

/// Images of the unit vectors under a surjection F_2^n -> F_2^(n - dim) with kernel span(basis).
fn quotient_map(basis: &[Fe], n: usize) -> (Vec<u128>, usize) {
    let piv = echelon(basis.iter().copied());
    let mut pos = vec![usize::MAX; n];
    let mut nfree = 0;
    for b in 0..n {
        if !piv.iter().any(|&(pb, _)| pb == b) {
            pos[b] = nfree;
            nfree += 1;
        }
    }
    assert!(nfree <= 128);
    let mut imgs = vec![0u128; n];
    for b in 0..n {
        if pos[b] != usize::MAX {
            imgs[b] = 1u128 << pos[b];
        } else {
            let pv = piv.iter().find(|&&(pb, _)| pb == b).unwrap().1;
            for c in 0..n {
                if c != b && pv.bit(c) {
                    imgs[b] ^= 1u128 << pos[c];
                }
            }
        }
    }
    (imgs, nfree)
}

#[derive(Default, Clone, Debug)]
pub struct Stats {
    pub abscissas: u64,
    pub targets: u64,
    pub lin_solvable: u64,
    pub lin_candidates: u64,
    pub quad_split: u64,
    pub found: u64,
    pub degenerate: u64,
    pub max_kernel_dim: u64,
}

impl Stats {
    pub fn merge(&mut self, o: &Stats) {
        self.abscissas += o.abscissas;
        self.targets += o.targets;
        self.lin_solvable += o.lin_solvable;
        self.lin_candidates += o.lin_candidates;
        self.quad_split += o.quad_split;
        self.found += o.found;
        self.degenerate += o.degenerate;
        self.max_kernel_dim = self.max_kernel_dim.max(o.max_kernel_dim);
    }
}

pub struct Solver {
    pub f: Field,
    pub n: usize,
    pub k: usize,
    pub a: Fe,
    pub b: Fe,
    pub sqrtb: Fe,
    pub basis: Vec<Fe>,
    pub dim_v2: usize,
    /// width of the linear system: n - dim V2 + 1
    pub w: usize,
    wv: Vec<Fe>,
    trv: u64,
    trmask: Fe,
    q: LMap<Fe>,
    sqrt: LMap<Fe>,
    chk: LMap<u128>,
    pi: LMap<u64>,
    /// Row form of the linear system (see `d2_rows`); None when V2's annihilator is not delta <1, w, ..>.
    rows: Option<RowForm>,
    /// Use the row form when it exists (default); false forces the generic column form.
    pub use_rows: bool,
}

/// For V = <1, w, .., w^(k-1)> and a modulus x^n + (terms of degree <= 2k - 1), the annihilator of V2 under the
/// trace form is delta <1, w, .., w^(n-2k)> with delta = 1/f'(w).  Equation j of the system is then the trace
/// form against z_j = w^j delta xT, so its whole coefficient row is one table lookup L(z_j), and z_(j+1) = w z_j
/// is a shift: the k + 1 field multiplications of the column form become two.
struct RowForm {
    delta: Fe,
    nq: usize,
    /// z -> (Tr(z W(v_i)))_(i<k), Tr(z) at bit k
    l: LMap<u64>,
    /// u -> (Tr(w^j u))_(j<nq)
    g: LMap<u64>,
}

struct Cand {
    tx: Fe,
    ym: Fe,
    e1: Fe,
    e2: Fe,
    lam: Fe,
    p3: Pt,
}

const BLK: usize = 128;

impl Solver {
    pub fn new(f: Field, a: Fe, b: Fe, basis: Vec<Fe>) -> Solver {
        let n = f.n;
        let k = basis.len();
        assert!(k <= 62 && echelon(basis.iter().copied()).len() == k, "basis must be independent, k <= 62");
        // Artin-Schreier inverse: echelonize the images z^2 + z of the unit vectors, tracking preimages
        let mut piv: Vec<(usize, Fe, Fe)> = Vec::new(); // (pivot, image, preimage)
        for i in 0..n {
            let (mut im, mut pre) = (f.sq(Fe::unit(i)) ^ Fe::unit(i), Fe::unit(i));
            for &(pb, pim, ppre) in &piv {
                if im.bit(pb) {
                    im ^= pim;
                    pre ^= ppre;
                }
            }
            if let Some(pb) = im.top() {
                for q in piv.iter_mut() {
                    if q.1.bit(pb) {
                        q.1 ^= im;
                        q.2 ^= pre;
                    }
                }
                piv.push((pb, im, pre));
            }
        }
        let qimg: Vec<Fe> = (0..n)
            .map(|bit| {
                let (mut v, mut pre) = (Fe::unit(bit), Fe::ZERO);
                let mut order: Vec<&(usize, Fe, Fe)> = piv.iter().collect();
                order.sort_by(|x, y| y.0.cmp(&x.0));
                for &&(pb, pim, ppre) in &order {
                    if v.bit(pb) {
                        v ^= pim;
                        pre ^= ppre;
                    }
                }
                pre
            })
            .collect();
        let q = LMap::new(&qimg, n);
        let sqrtimg: Vec<Fe> = (0..n)
            .map(|i| {
                let mut t = Fe::unit(i);
                for _ in 0..n - 1 {
                    t = f.sq(t);
                }
                t
            })
            .collect();
        let sqrt = LMap::new(&sqrtimg, n);
        // Tr(e_i): Tr(x) = x + x^2 + ... + x^(2^(n-1)) is 0 or 1
        let mut trmask = Fe::ZERO;
        for i in 0..n {
            let (mut t, mut s) = (Fe::unit(i), Fe::ZERO);
            for _ in 0..n {
                s ^= t;
                t = f.sq(t);
            }
            assert!(s == Fe::ZERO || s == Fe::ONE);
            if s == Fe::ONE {
                trmask ^= Fe::unit(i);
            }
        }
        let mut prods = Vec::new();
        for i in 0..k {
            for j in i..k {
                prods.push(f.mul(basis[i], basis[j]));
            }
        }
        let v2: Vec<Fe> = echelon(prods).into_iter().map(|x| x.1).collect();
        let (pi128, wq) = quotient_map(&v2, n);
        assert!(wq + 1 <= 63, "n - dim V2 must be at most 62");
        let pi64: Vec<u64> = pi128.iter().map(|&x| x as u64).collect();
        let (chkimg, _) = quotient_map(&basis, n);
        let wv: Vec<Fe> = basis.iter().map(|&v| q.at(v) ^ v).collect();
        let mut trv = 0u64;
        for (i, v) in basis.iter().enumerate() {
            if v.and_parity(&trmask) {
                trv |= 1 << i;
            }
        }
        let sqrtb = sqrt.at(b);
        let tr = |x: Fe| x.and_parity(&trmask);
        let nq = wq;
        let rows = {
            // f'(w): odd exponents of x^n + sum x^e
            let pow = |e: usize| (0..e).fold(Fe::ONE, |acc, _| f.mulx(acc));
            let mut d = if n % 2 == 1 { pow(n - 1) } else { Fe::ZERO };
            for &e in &f.low {
                if e % 2 == 1 {
                    d ^= pow(e as usize - 1);
                }
            }
            if d.is_zero() || nq == 0 {
                None
            } else {
                let delta = f.inv(d);
                let mut ok = true;
                let mut eta = delta;
                for _ in 0..nq {
                    ok &= v2.iter().all(|&v| !tr(f.mul(eta, v)));
                    eta = f.mulx(eta);
                }
                if !ok {
                    None
                } else {
                    let limg: Vec<u64> = (0..n)
                        .map(|bit| {
                            let e = Fe::unit(bit);
                            let mut r = if tr(e) { 1u64 << k } else { 0 };
                            for i in 0..k {
                                if tr(f.mul(e, wv[i])) {
                                    r |= 1 << i;
                                }
                            }
                            r
                        })
                        .collect();
                    let mut trw = Vec::with_capacity(n + nq);
                    let mut t = Fe::ONE;
                    for _ in 0..n + nq {
                        trw.push(tr(t));
                        t = f.mulx(t);
                    }
                    let gimg: Vec<u64> =
                        (0..n).map(|bit| (0..nq).fold(0u64, |r, j| r | ((trw[j + bit] as u64) << j))).collect();
                    Some(RowForm { delta, nq, l: LMap::new(&limg, n), g: LMap::new(&gimg, n) })
                }
            }
        };
        Solver {
            rows, use_rows: true,
            n, k, a, b, sqrtb, dim_v2: v2.len(), w: wq + 1, wv, trv, trmask, q, sqrt,
            chk: LMap::new(&chkimg, n), pi: LMap::new(&pi64, n), basis, f,
        }
    }

    #[inline(always)]
    pub fn tr(&self, x: Fe) -> bool {
        x.and_parity(&self.trmask)
    }
    #[inline(always)]
    pub fn in_v(&self, x: Fe) -> bool {
        self.chk.at(x) == 0
    }

    /// The element of V with Gray-code index g.
    pub fn gray(&self, g: u64) -> Fe {
        let mut bits = g ^ (g >> 1);
        let mut x = Fe::ZERO;
        while bits != 0 {
            x ^= self.basis[bits.trailing_zeros() as usize];
            bits &= bits - 1;
        }
        x
    }

    pub fn on_curve(&self, p: Pt) -> bool {
        let f = &self.f;
        let x2 = f.sq(p.0);
        f.sq(p.1) ^ f.mul(p.0, p.1) == f.mul(x2, p.0) ^ f.mul(self.a, x2) ^ self.b
    }
    pub fn neg(&self, p: Pt) -> Pt {
        (p.0, p.0 ^ p.1)
    }
    /// Group law (None = the point at infinity); one inversion per call.  The reference implementation.
    pub fn add(&self, p: Option<Pt>, q: Option<Pt>) -> Option<Pt> {
        let f = &self.f;
        let (p, q) = match (p, q) {
            (None, q) => return q,
            (p, None) => return p,
            (Some(p), Some(q)) => (p, q),
        };
        if p.0 == q.0 {
            if q.1 == (p.0 ^ p.1) {
                return None;
            }
            let lam = p.0 ^ f.mul(p.1, f.inv(p.0));
            let x3 = f.sq(lam) ^ lam ^ self.a;
            return Some((x3, f.sq(p.0) ^ f.mul(lam ^ Fe::ONE, x3)));
        }
        let lam = f.mul(p.1 ^ q.1, f.inv(p.0 ^ q.0));
        let x3 = f.sq(lam) ^ lam ^ p.0 ^ q.0 ^ self.a;
        Some((x3, f.mul(lam, p.0 ^ x3) ^ x3 ^ p.1))
    }
    /// A point with abscissa x, if there is one.
    pub fn lift(&self, x: Fe) -> Option<Pt> {
        let f = &self.f;
        if x.is_zero() {
            return Some((x, self.sqrtb));
        }
        let c = x ^ self.a ^ f.mul(self.b, f.sq(f.inv(x)));
        if self.tr(c) {
            return None;
        }
        Some((x, f.mul(x, self.q.at(c))))
    }

    /// The linear system of one target (1/xT supplied); solutions are queued as candidates.
    #[inline(always)]
    fn d2_core(&self, t: Pt, ixt: Fe, p3: Pt, queue: &mut Vec<Cand>, st: &mut Stats) {
        match &self.rows {
            Some(rf) if self.use_rows => self.d2_rows(rf, t, ixt, p3, queue, st),
            _ => self.d2_cols(t, ixt, p3, queue, st),
        }
    }

    pub fn has_row_form(&self) -> bool {
        self.rows.is_some()
    }

    #[inline(always)]
    fn push_candidate(&self, sol: u64, t: Pt, ym: Fe, c1: Fe, c2: Fe, p3: Pt, queue: &mut Vec<Cand>) {
        let k = self.k;
        let mut e1 = Fe::ZERO;
        let mut tm = sol & ((1u64 << k) - 1);
        while tm != 0 {
            e1 ^= self.basis[tm.trailing_zeros() as usize];
            tm &= tm - 1;
        }
        let mut lam = self.q.at(e1 ^ c1);
        if (sol >> k) & 1 == 1 {
            lam ^= Fe::ONE;
        }
        let e2 = self.f.mul(t.0, lam ^ e1 ^ c1) ^ c2; // lam^2 = lam + e1 + c1
        queue.push(Cand { tx: t.0, ym, e1, e2, lam, p3 });
    }

    /// The same system by rows: equation j is the trace form against z_j = w^j delta xT.
    #[inline(always)]
    fn d2_rows(&self, rf: &RowForm, t: Pt, ixt: Fe, p3: Pt, queue: &mut Vec<Cand>, st: &mut Stats) {
        let f = &self.f;
        let k = self.k;
        st.targets += 1;
        let ym = t.1 ^ t.0;
        let c1 = self.a ^ t.0;
        let c2 = f.mul(f.sq(ym) ^ self.b, ixt);
        let mut z = f.mul(rf.delta, t.0);
        let rhs = rf.g.at(f.mul(z, self.q.at(c1) ^ c1) ^ f.mul(rf.delta, c2));
        let varmask = (1u64 << (k + 1)) - 1;
        // all rows first (the right-hand side in bit 63), then Gauss-Jordan: each pivot row is cleared from
        // every other row in one pass with no dependency between rows, which the compiler vectorizes
        let mut rows = [0u64; 64];
        let nrows = rf.nq + 1;
        let nr8 = (nrows + 7) & !7;
        rows[0] = self.trv | ((self.tr(c1) as u64) << 63);
        for j in 0..rf.nq {
            rows[j + 1] = rf.l.at(z) | (((rhs >> j) & 1) << 63);
            z = f.mulx(z);
        }
        let mut pivrow = [0u8; 64];
        let mut piv = [0u32; 64];
        let mut nb = 0usize;
        let mut pivmask = 0u64;
        for p in 0..nrows {
            let r = rows[p];
            let low = r & varmask;
            if low == 0 {
                if r >> 63 != 0 {
                    return;
                }
                continue;
            }
            let c = low.trailing_zeros();
            for row in rows[..nr8].iter_mut() {
                let mk = ((*row >> c) & 1).wrapping_neg();
                *row ^= r & mk;
            }
            rows[p] = r;
            pivrow[nb] = p as u8;
            piv[nb] = c;
            pivmask |= 1 << c;
            nb += 1;
        }
        st.lin_solvable += 1;
        let mut part = 0u64;
        for i in 0..nb {
            part |= (rows[pivrow[i] as usize] >> 63) << piv[i];
        }
        let mut free = varmask & !pivmask;
        let mut kern = [0u64; 64];
        let mut nk = 0usize;
        while free != 0 {
            let c = free.trailing_zeros();
            let mut kv = 1u64 << c;
            for i in 0..nb {
                kv |= ((rows[pivrow[i] as usize] >> c) & 1) << piv[i];
            }
            kern[nk] = kv;
            nk += 1;
            free &= free - 1;
        }
        st.max_kernel_dim = st.max_kernel_dim.max(nk as u64);
        let mut sol = part;
        for g in 0u64..(1u64 << nk) {
            if g != 0 {
                sol ^= kern[g.trailing_zeros() as usize];
            }
            st.lin_candidates += 1;
            self.push_candidate(sol, t, ym, c1, c2, p3, queue);
        }
    }

    /// The system by columns (any V): k + 1 multiplications by xT and the projection pi.
    #[inline(always)]
    fn d2_cols(&self, t: Pt, ixt: Fe, p3: Pt, queue: &mut Vec<Cand>, st: &mut Stats) {
        let f = &self.f;
        let k = self.k;
        st.targets += 1;
        let ym = t.1 ^ t.0;
        let c1 = self.a ^ t.0;
        let c2 = f.mul(f.sq(ym) ^ self.b, ixt);
        let trbit = 1u64 << (self.w - 1);
        let mut rhs = self.pi.at(f.mul(t.0, self.q.at(c1) ^ c1) ^ c2) | if self.tr(c1) { trbit } else { 0 };
        let mut bvec = [0u64; 64];
        let mut bcomb = [0u64; 64];
        let mut piv = [0u32; 64];
        let mut kern = [0u64; 64];
        let (mut nb, mut nk) = (0usize, 0usize);
        for i in 0..=k {
            let mut col = if i < k {
                self.pi.at(f.mul(t.0, self.wv[i])) | if (self.trv >> i) & 1 == 1 { trbit } else { 0 }
            } else {
                self.pi.at(t.0)
            };
            let mut comb = 1u64 << i;
            for j in 0..nb {
                let mk = ((col >> piv[j]) & 1).wrapping_neg();
                col ^= bvec[j] & mk;
                comb ^= bcomb[j] & mk;
            }
            if col != 0 {
                piv[nb] = 63 - col.leading_zeros();
                bvec[nb] = col;
                bcomb[nb] = comb;
                nb += 1;
            } else {
                kern[nk] = comb;
                nk += 1;
            }
        }
        let mut part = 0u64;
        for j in 0..nb {
            let mk = ((rhs >> piv[j]) & 1).wrapping_neg();
            rhs ^= bvec[j] & mk;
            part ^= bcomb[j] & mk;
        }
        if rhs != 0 {
            return;
        }
        st.lin_solvable += 1;
        st.max_kernel_dim = st.max_kernel_dim.max(nk as u64);
        let mut sol = part;
        for g in 0u64..(1u64 << nk) {
            if g != 0 {
                sol ^= kern[g.trailing_zeros() as usize];
            }
            st.lin_candidates += 1;
            self.push_candidate(sol, t, ym, c1, c2, p3, queue);
        }
    }

    /// The quadratics of the queued candidates, with shared inversions of e1.
    fn flush(&self, queue: &mut Vec<Cand>, out: &mut Vec<[Pt; 3]>, st: &mut Stats) {
        let f = &self.f;
        let e1: Vec<Fe> = queue.iter().map(|c| c.e1).collect();
        let mut ie = vec![Fe::ZERO; e1.len()];
        f.batch_inv(&e1, &mut ie);
        for (c, ie) in queue.iter().zip(ie) {
            let (x1, x2);
            if c.e1.is_zero() {
                x1 = self.sqrt.at(c.e2);
                x2 = x1;
            } else {
                let d = f.mul(c.e2, f.sq(ie));
                if self.tr(d) {
                    continue;
                }
                x1 = f.mul(c.e1, self.q.at(d));
                x2 = x1 ^ c.e1;
            }
            st.quad_split += 1;
            if !self.in_v(x1) {
                continue;
            }
            let mu = f.mul(c.lam, c.tx) ^ c.ym;
            st.found += 1;
            let mut dec = [c.p3, (x1, f.mul(c.lam, x1) ^ mu), (x2, f.mul(c.lam, x2) ^ mu)];
            dec.sort();
            out.push(dec);
        }
        queue.clear();
    }

    /// Two-summand decompositions of T by enumeration over V (k small): the fallback for x(T) = 0.
    fn d2_enum(&self, t: Pt, p3: Pt, out: &mut Vec<[Pt; 3]>) {
        assert!(self.k <= 24, "x(T) = 0 with a factor base too large to enumerate");
        for g in 0..(1u64 << self.k) {
            if let Some(p) = self.lift(self.gray(g)) {
                let signs = if p.0.is_zero() { vec![p] } else { vec![p, self.neg(p)] };
                for p in signs {
                    if let Some(q) = self.add(Some(t), Some(self.neg(p))) {
                        if self.in_v(q.0) {
                            let mut dec = [p3, p, q];
                            dec.sort();
                            out.push(dec);
                        }
                    }
                }
            }
        }
    }

    /// Two-summand decompositions of T (pairs sorted, deduplicated).
    pub fn decompose2(&self, t: Pt) -> (Vec<[Pt; 2]>, Stats) {
        let mut st = Stats::default();
        let mut out3 = Vec::new();
        let dummy = (Fe::ZERO, Fe::ZERO);
        if t.0.is_zero() {
            st.degenerate += 1;
            self.d2_enum(t, dummy, &mut out3);
        } else {
            let mut queue = Vec::new();
            self.d2_core(t, self.f.inv(t.0), dummy, &mut queue, &mut st);
            self.flush(&mut queue, &mut out3, &mut st);
        }
        let mut out: Vec<[Pt; 2]> = out3
            .into_iter()
            .map(|d| {
                let mut v: Vec<Pt> = Vec::new();
                let mut skipped = false;
                for p in d {
                    if p == dummy && !skipped {
                        skipped = true;
                    } else {
                        v.push(p);
                    }
                }
                [v[0], v[1]]
            })
            .collect();
        out.sort();
        out.dedup();
        (out, st)
    }

    /// Three-summand decompositions R = P1 + P2 + P3 with x(P3) the abscissas of Gray index in [lo, hi).
    /// Rows are sorted triples; the same triple can appear more than once across the whole index range.
    pub fn lrf3_range(&self, r: Pt, lo: u64, hi: u64) -> (Vec<[Pt; 3]>, Stats) {
        let f = &self.f;
        let mut st = Stats::default();
        let mut out = Vec::new();
        let mut queue: Vec<Cand> = Vec::with_capacity(512);
        let mut x = self.gray(lo);
        let mut xs = [Fe::ZERO; BLK];
        let mut dens = [Fe::ZERO; 2 * BLK];
        let mut invs = [Fe::ZERO; 2 * BLK];
        let mut ts: Vec<(Pt, Pt)> = Vec::with_capacity(2 * BLK);
        let mut tx = [Fe::ZERO; 2 * BLK];
        let mut itx = [Fe::ZERO; 2 * BLK];
        let mut g = lo;
        while g < hi {
            let cnt = ((hi - g) as usize).min(BLK);
            for j in 0..cnt {
                if g + j as u64 != lo {
                    x ^= self.basis[(g + j as u64).trailing_zeros() as usize];
                }
                xs[j] = x;
                dens[2 * j] = x;
                dens[2 * j + 1] = x ^ r.0;
            }
            st.abscissas += cnt as u64;
            f.batch_inv(&dens[..2 * cnt], &mut invs[..2 * cnt]);
            ts.clear();
            for j in 0..cnt {
                let xj = xs[j];
                if xj.is_zero() || dens[2 * j + 1].is_zero() {
                    // the 2-torsion point, or P3 = +-R: generic group law
                    if let Some(p) = self.lift(xj) {
                        let signs = if xj.is_zero() { vec![p] } else { vec![p, self.neg(p)] };
                        for p in signs {
                            if let Some(t) = self.add(Some(r), Some(self.neg(p))) {
                                ts.push((t, p));
                            }
                        }
                    }
                    continue;
                }
                let c = xj ^ self.a ^ f.mul(self.b, f.sq(invs[2 * j]));
                if self.tr(c) {
                    continue;
                }
                let y = f.mul(xj, self.q.at(c));
                let d = dens[2 * j + 1];
                for py in [y, y ^ xj] {
                    // T = R - P = R + (xj, py + xj)
                    let lam = f.mul(r.1 ^ py ^ xj, invs[2 * j + 1]);
                    let t0 = f.sq(lam) ^ lam ^ d ^ self.a;
                    let t1 = f.mul(lam, r.0 ^ t0) ^ t0 ^ r.1;
                    ts.push(((t0, t1), (xj, py)));
                }
            }
            for (i, t) in ts.iter().enumerate() {
                tx[i] = t.0 .0;
            }
            f.batch_inv(&tx[..ts.len()], &mut itx[..ts.len()]);
            for (i, &(t, p3)) in ts.iter().enumerate() {
                if t.0.is_zero() {
                    st.degenerate += 1;
                    self.d2_enum(t, p3, &mut out);
                } else {
                    self.d2_core(t, itx[i], p3, &mut queue, &mut st);
                }
                if queue.len() >= 256 {
                    self.flush(&mut queue, &mut out, &mut st);
                }
            }
            g += cnt as u64;
        }
        self.flush(&mut queue, &mut out, &mut st);
        (out, st)
    }

    /// The same search by enumeration: for P3 over [lo, hi) and every factor-base P2, test x(R - P3 - P2) in V.
    /// Additions share inversions in blocks, as in `lrf3_range`.  `fb` is the listed factor base.
    /// Returns the rows and the number of curve additions done.
    pub fn brute3_range(&self, r: Pt, fb: &[Pt], lo: u64, hi: u64) -> (Vec<[Pt; 3]>, u64) {
        let f = &self.f;
        let mut out = Vec::new();
        let mut adds = 0u64;
        let mut dens = [Fe::ZERO; BLK];
        let mut invs = [Fe::ZERO; BLK];
        for g in lo..hi {
            let Some(p) = self.lift(self.gray(g)) else { continue };
            let signs = if p.0.is_zero() { vec![p] } else { vec![p, self.neg(p)] };
            for p3 in signs {
                adds += 1;
                let Some(s) = self.add(Some(r), Some(self.neg(p3))) else { continue };
                for chunk in fb.chunks(BLK) {
                    for (j, p2) in chunk.iter().enumerate() {
                        dens[j] = s.0 ^ p2.0;
                    }
                    f.batch_inv(&dens[..chunk.len()], &mut invs[..chunk.len()]);
                    adds += chunk.len() as u64;
                    for (j, &p2) in chunk.iter().enumerate() {
                        let t = if dens[j].is_zero() {
                            match self.add(Some(s), Some(self.neg(p2))) {
                                Some(t) => t,
                                None => continue,
                            }
                        } else {
                            let lam = f.mul(s.1 ^ p2.1 ^ p2.0, invs[j]);
                            let t0 = f.sq(lam) ^ lam ^ dens[j] ^ self.a;
                            (t0, f.mul(lam, s.0 ^ t0) ^ t0 ^ s.1)
                        };
                        if self.in_v(t.0) {
                            let mut dec = [p3, p2, t];
                            dec.sort();
                            out.push(dec);
                        }
                    }
                }
            }
        }
        (out, adds)
    }

    /// Every point with abscissa in V (k small).
    pub fn factor_base(&self) -> Vec<Pt> {
        assert!(self.k <= 26);
        let mut fb = Vec::new();
        for g in 0..(1u64 << self.k) {
            if let Some(p) = self.lift(self.gray(g)) {
                fb.push(p);
                if !p.0.is_zero() {
                    fb.push(self.neg(p));
                }
            }
        }
        fb
    }
}
