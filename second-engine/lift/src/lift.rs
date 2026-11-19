//! The degree lift of experiments/ic-bounds-k14/LIFT.md, in ordinary monomial coordinates.
//!
//! Third implementation, written from the theorem (sections 1-2 of LIFT.md), not transliterated:
//! no sigma split at all (the theorem with the trivial symmetry group; Lemma 3 is not used), and V_D is
//! taken in its reduced echelon form in the file's own (degree-non-increasing) column order, so
//! that the decomposition of the independent implementation, V_D = V_P (+) V_N, has V_N = 0.
//!
//! Notation (theorem step -> code):
//!  * L = span(M_<=D); columns = the file's column list (`Layout`).  V_D = rows e_x + B[x] (Rref).
//!  * Pd = pivot columns ("regular" weight-D monomials are the pivots of weight D); N = free columns.
//!  * top = weight-(D+1) monomial; special  <=>  t = u n for a free n of weight D  (equivalently: some
//!    weight-D divisor of t is free, since u n has weight D + 1 iff no carry, i.e. n = t minus one bit);
//!    regular otherwise (all D + 1 divisors are pivots).                                [`Lift::new`]
//!  * T = F^Sp (+) F^N, Sp columns first; red(e_m) = Reg(m) | unit of T | -B[c] on N.  [`Tgt`]
//!  * Lemma 4: gamma(x, u) = red(u l_x) = red(e_{ux}) + B[x] R_u, R_u: e_n -> red(e_{un}) on N
//!    ("follower"; never has a regular coordinate: asserted).                         [`Lift::follower`]
//!  * Lemma 5: generators of one regular top t are differenced (gamma_i - gamma_1); generators with no
//!    regular coordinate go to K directly.                                            [`Lift::run`]
//!  * Lemma 2: only rows of pivot weight D are shifted (`include_low` shifts all rows).
//!  * Lemma 1 / fixed point: falls = rows of K (RREF, Sp first) with pivot in the N block; every fall
//!    row is shifted once at its state (all 2k shifts) until none is pending; `verify_falls` then
//!    shifts ALL final fall rows again and requires zero residuals (explicit Phi R_u c K).
//!  * Statements: dim V_{D+1} = dim V_D + |Reg| + dim K; V_{D+1} n L = V_D (+) F; 1 in V_{D+1} iff the
//!    monomial 1 is a pivot of V_D or its N column is a pivot of K (it is the last column, so its K
//!    row is the unit vector).                                                        [`Lift::finish`]
//!
//! Every vector is kept reduced modulo the current K on the live (non-pivot) T columns; when K grows,
//! stored vectors are synced: v <- v[live'] - v[new pivots] @ K[new pivots, live'] (exact because K
//! is in reduced echelon form and v was zero on the old pivots).

use crate::basis::{Rref, NONE};
use crate::field::{gemm, madds, rank_of, rref_extend, Fp};
use crate::layout::{FxMap, Layout};
use rayon::prelude::*;
use serde_json::json;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct LiftOpts {
    /// regular tops per batch
    pub batch_tops: usize,
    /// a fall round after every this many batches (0: only at the end)
    pub fall_every: usize,
    pub include_low: bool,
    pub verify_falls: bool,
    pub follower_cache_gb: f64,
    /// budget for one batch of generator products
    pub work_gb: f64,
    pub hash: bool,
    pub parts: bool,
    /// negative-control hooks
    pub test_no_fall_shifts: bool,
    pub test_drop_shift: Option<u32>,
}

impl Default for LiftOpts {
    fn default() -> Self {
        LiftOpts {
            batch_tops: 4096,
            fall_every: 1,
            include_low: false,
            verify_falls: true,
            follower_cache_gb: 4.0,
            work_gb: 2.0,
            hash: true,
            parts: false,
            test_no_fall_shifts: false,
            test_drop_shift: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Tgt {
    Reg,
    T(u32),
    Piv(u32),
}

struct KSpace {
    nt: usize,
    rows: Vec<Vec<u32>>,
    pivs: Vec<u32>,
    rowof: Vec<u32>,
    live: Vec<u32>,
    lpos: Vec<u32>,
    version: u64,
}

impl KSpace {
    fn new(nt: usize) -> KSpace {
        KSpace {
            nt,
            rows: vec![],
            pivs: vec![],
            rowof: vec![NONE; nt],
            live: (0..nt as u32).collect(),
            lpos: (0..nt as u32).collect(),
            version: 0,
        }
    }
    fn rank(&self) -> usize {
        self.rows.len()
    }
    /// rho(e_t) on the live columns
    fn rho_unit(&self, fp: &Fp, t: u32, out: &mut [u32]) {
        let r = self.rowof[t as usize];
        if r == NONE {
            out[self.lpos[t as usize] as usize] = fp.add(out[self.lpos[t as usize] as usize], 1);
        } else {
            let row = &self.rows[r as usize];
            for (o, &c) in out.iter_mut().zip(&self.live) {
                *o = fp.sub(*o, row[c as usize]);
            }
        }
    }
    /// add candidates (rows on the live columns, reduced mod K); returns number of new pivots
    fn insert(&mut self, fp: &Fp, cands: Vec<Vec<u32>>) -> usize {
        let c = self.live.len();
        let cands: Vec<Vec<u32>> = cands.into_iter().filter(|v| v.iter().any(|&x| x != 0)).collect();
        if cands.is_empty() || c == 0 {
            return 0;
        }
        let mut nr = Vec::new();
        let mut np = Vec::new();
        rref_extend(fp, &mut nr, &mut np, cands, c);
        if nr.is_empty() {
            return 0;
        }
        // old rows: r[L] -= r[newpiv] @ nr
        let nold = self.rows.len();
        if nold > 0 {
            let nn = np.len();
            let live = &self.live;
            let newcols: Vec<u32> = np.iter().map(|&q| live[q]).collect();
            let bm: Vec<u32> = nr.concat();
            self.rows.par_iter_mut().for_each(|row| {
                let a: Vec<u32> = newcols.iter().map(|&cc| row[cc as usize]).collect();
                if a.iter().all(|&x| x == 0) {
                    return;
                }
                let mut o: Vec<u32> = live.iter().map(|&cc| row[cc as usize]).collect();
                gemm(fp, &mut o, &a, &bm, 1, nn, c, true);
                for (j, &cc) in live.iter().enumerate() {
                    row[cc as usize] = o[j];
                }
            });
        }
        let added = nr.len();
        for (v, q) in nr.into_iter().zip(np) {
            let mut full = vec![0u32; self.nt];
            for (j, &cc) in self.live.iter().enumerate() {
                full[cc as usize] = v[j];
            }
            let col = self.live[q];
            debug_assert_eq!(full[col as usize], 1);
            self.rowof[col as usize] = self.rows.len() as u32;
            self.rows.push(full);
            self.pivs.push(col);
        }
        self.live = (0..self.nt as u32).filter(|&t| self.rowof[t as usize] == NONE).collect();
        self.lpos = vec![NONE; self.nt];
        for (j, &t) in self.live.iter().enumerate() {
            self.lpos[t as usize] = j as u32;
        }
        self.version += 1;
        added
    }
}

/// a matrix of rows reduced modulo K on a list of T columns
struct Fol {
    nrows: usize,
    cols: Vec<u32>,
    data: Vec<u32>,
    version: u64,
    last_use: u64,
}

impl Fol {
    fn bytes(&self) -> usize {
        self.data.len() * 4
    }
    fn sync(&mut self, fp: &Fp, k: &KSpace, counter: &mut u64) {
        if self.version == k.version {
            return;
        }
        let np: Vec<usize> = (0..self.cols.len()).filter(|&j| k.rowof[self.cols[j] as usize] != NONE).collect();
        if np.is_empty() {
            self.version = k.version;
            return;
        }
        let keep: Vec<usize> = (0..self.cols.len()).filter(|&j| k.rowof[self.cols[j] as usize] == NONE).collect();
        let c0 = self.cols.len();
        let n = self.nrows;
        let nk = keep.len();
        let mut out = vec![0u32; n * nk];
        let mut a = vec![0u32; n * np.len()];
        for i in 0..n {
            let row = &self.data[i * c0..(i + 1) * c0];
            for (jj, &j) in keep.iter().enumerate() {
                out[i * nk + jj] = row[j];
            }
            for (jj, &j) in np.iter().enumerate() {
                a[i * np.len() + jj] = row[j];
            }
        }
        let mut bm = vec![0u32; np.len() * nk];
        for (t, &j) in np.iter().enumerate() {
            let krow = &k.rows[k.rowof[self.cols[j] as usize] as usize];
            for (jj, &kk) in keep.iter().enumerate() {
                bm[t * nk + jj] = krow[self.cols[kk] as usize];
            }
        }
        gemm(fp, &mut out, &a, &bm, n, np.len(), nk, true);
        *counter += (n * np.len() * nk) as u64;
        self.cols = keep.iter().map(|&j| self.cols[j]).collect();
        self.data = out;
        self.version = k.version;
        debug_assert_eq!(self.cols, k.live);
    }
}

pub struct Lift<'a> {
    fp: Fp,
    lay: &'a Layout,
    vd: &'a Rref,
    opts: LiftOpts,
    d: u32,
    ns: u32,
    sp: Vec<u64>,
    sp_idx: FxMap<u64, u32>,
    nsp: usize,
    nt: usize,
    reg_tops: Vec<u64>,
    k: KSpace,
    fols: FxMap<u32, Fol>,
    tick: u64,
    fall_done: Vec<bool>,
    // counters
    pub c_gen_rows: u64,
    pub c_direct_rows: u64,
    pub c_fall_rows: u64,
    pub c_candidates: u64,
    pub c_follower_builds: u64,
    pub c_follower_evictions: u64,
    pub c_sync_madds: u64,
    pub c_rounds: u64,
    pub falls_new_per_round: Vec<u64>,
    pub peak_follower_bytes: usize,
    pub c_reg_in_follower: u64,
    pub t_fol: f64,
    pub t_ins: f64,
    pub t_prod: f64,
    pub t_falls: f64,
}

impl<'a> Lift<'a> {
    pub fn new(fp: Fp, lay: &'a Layout, vd: &'a Rref, opts: LiftOpts) -> Lift<'a> {
        let d = lay.d;
        let mon = lay.mon;
        let ns = mon.nshifts();
        // special tops: u n for free n of weight D (when that has weight D + 1)
        let mut sp: Vec<u64> = vd
            .free
            .par_iter()
            .filter(|&&c| lay.wt[c as usize] as u32 == d)
            .flat_map_iter(|&c| {
                let m = lay.masks[c as usize];
                (0..ns).map(move |s| mon.shift(m, s)).filter(|y| y.count_ones() == d + 1)
            })
            .collect();
        sp.par_sort_unstable_by(|a, b| b.cmp(a));
        sp.dedup();
        let mut sp_idx = FxMap::default();
        for (i, &m) in sp.iter().enumerate() {
            sp_idx.insert(m, i as u32);
        }
        let nsp = sp.len();
        let nt = nsp + vd.nf();
        let reg_tops: Vec<u64> = {
            let mut v: Vec<u64> = mon.of_weight(d + 1).into_par_iter().filter(|m| !sp_idx.contains_key(m)).collect();
            v.par_sort_unstable_by(|a, b| b.cmp(a));
            v
        };
        let nfree = vd.nf();
        Lift {
            fp,
            lay,
            vd,
            opts,
            d,
            ns,
            sp,
            sp_idx,
            nsp,
            nt,
            reg_tops,
            k: KSpace::new(nt),
            fols: FxMap::default(),
            tick: 0,
            fall_done: vec![false; nfree + nsp],
            c_gen_rows: 0,
            c_direct_rows: 0,
            c_fall_rows: 0,
            c_candidates: 0,
            c_follower_builds: 0,
            c_follower_evictions: 0,
            c_sync_madds: 0,
            c_rounds: 0,
            falls_new_per_round: vec![],
            peak_follower_bytes: 0,
            c_reg_in_follower: 0,
            t_fol: 0.0,
            t_ins: 0.0,
            t_prod: 0.0,
            t_falls: 0.0,
        }
    }

    #[inline]
    fn tgt(&self, y: u64) -> Tgt {
        if y.count_ones() == self.d + 1 {
            match self.sp_idx.get(&y) {
                Some(&i) => Tgt::T(i),
                None => Tgt::Reg,
            }
        } else {
            let c = self.lay.col(y).expect("shift left M_<=D") as usize;
            let r = self.vd.prow[c];
            if r != NONE {
                Tgt::Piv(r)
            } else {
                Tgt::T(self.nsp as u32 + self.vd.fpos[c])
            }
        }
    }

    fn gather_b(&self, rows: &[u32]) -> Vec<u32> {
        let nf = self.vd.nf();
        let mut a = vec![0u32; rows.len() * nf];
        a.par_chunks_mut(nf.max(1)).zip(rows.par_iter()).for_each(|(o, &r)| {
            if nf > 0 {
                o.copy_from_slice(self.vd.row(r as usize))
            }
        });
        a
    }

    /// R_s: row f = rho(red(e_{u_s f})) on the live columns
    fn follower(&mut self, key: u32) {
        let t = Instant::now();
        self.follower_inner(key);
        self.t_fol += t.elapsed().as_secs_f64();
    }

    fn follower_inner(&mut self, key: u32) {
        self.tick += 1;
        let tick = self.tick;
        if self.fols.contains_key(&key) {
            let mut f = self.fols.remove(&key).unwrap();
            f.sync(&self.fp, &self.k, &mut self.c_sync_madds);
            f.last_use = tick;
            self.fols.insert(key, f);
            return;
        }
        let nf = self.vd.nf();
        let nl = self.k.live.len();
        let mut data = vec![0u32; nf * nl];
        {
            let mon = self.lay.mon;
            let tg: Vec<Tgt> = self
                .vd
                .free
                .par_iter()
                .map(|&c| self.tgt(mon.shift(self.lay.masks[c as usize], key)))
                .collect();
            let nreg = tg.iter().filter(|t| **t == Tgt::Reg).count();
            assert_eq!(nreg, 0, "Lemma 4 violated: a follower row has a regular top coordinate");
            if nl > 0 {
                let hit: Vec<(usize, u32)> =
                    tg.iter().enumerate().filter_map(|(f, t)| if let Tgt::Piv(r) = t { Some((f, *r)) } else { None }).collect();
                let rows: Vec<u32> = hit.iter().map(|h| h.1).collect();
                let g = self.gather_b(&rows);
                let y = self.rho_neg_n(&g, hit.len());
                for (i, &(f, _)) in hit.iter().enumerate() {
                    data[f * nl..(f + 1) * nl].copy_from_slice(&y[i * nl..(i + 1) * nl]);
                }
                let k = &self.k;
                let fp = self.fp;
                data.par_chunks_mut(nl).zip(tg.par_iter()).for_each(|(o, t)| {
                    if let Tgt::T(tc) = t {
                        k.rho_unit(&fp, *tc, o)
                    }
                });
            }
        }
        self.c_follower_builds += 1;
        let f = Fol { nrows: nf, cols: self.k.live.clone(), data, version: self.k.version, last_use: tick };
        // evict least recently used (never Z) while over budget
        let budget = (self.opts.follower_cache_gb * 1e9) as usize;
        loop {
            let tot: usize = self.fols.values().map(|x| x.bytes()).sum::<usize>() + f.bytes();
            self.peak_follower_bytes = self.peak_follower_bytes.max(tot);
            if tot <= budget {
                break;
            }
            let victim = self.fols.iter().min_by_key(|(_, v)| v.last_use).map(|(kk, _)| *kk);
            match victim {
                Some(v) => {
                    self.fols.remove(&v);
                    self.c_follower_evictions += 1;
                }
                None => break,
            }
        }
        self.fols.insert(key, f);
    }

    /// rho(-G) on the live columns for G (n x #N) supported on the N block of T:
    /// -G @ Z with Z_f = e_f (f live) or -K_f (f a pivot of K); only the pivot part is a product.
    fn rho_neg_n(&mut self, g: &[u32], n: usize) -> Vec<u32> {
        let nl = self.k.live.len();
        let nf = self.vd.nf();
        let nsp = self.nsp;
        let mut out = vec![0u32; n * nl];
        if nl == 0 || n == 0 {
            return out;
        }
        let k = &self.k;
        let fp = self.fp;
        let pivf: Vec<usize> = (0..nf).filter(|&f| k.rowof[nsp + f] != NONE).collect();
        let livef: Vec<(usize, usize)> =
            (0..nf).filter(|&f| k.rowof[nsp + f] == NONE).map(|f| (f, k.lpos[nsp + f] as usize)).collect();
        out.par_chunks_mut(nl).zip(g.par_chunks(nf)).for_each(|(o, gr)| {
            for &(f, j) in &livef {
                o[j] = fp.neg(gr[f]);
            }
        });
        if !pivf.is_empty() {
            let np = pivf.len();
            let mut a = vec![0u32; n * np];
            a.par_chunks_mut(np).zip(g.par_chunks(nf)).for_each(|(o, gr)| {
                for (t, &f) in pivf.iter().enumerate() {
                    o[t] = gr[f];
                }
            });
            let mut km = vec![0u32; np * nl];
            km.par_chunks_mut(nl).zip(pivf.par_iter()).for_each(|(o, &f)| {
                let row = &k.rows[k.rowof[nsp + f] as usize];
                for (j, &c) in k.live.iter().enumerate() {
                    o[j] = row[c as usize];
                }
            });
            gemm(&fp, &mut out, &a, &km, n, np, nl, false);
        }
        out
    }

    /// Y = B[xs] @ R_s (+ eterms) on the live columns
    fn products(&mut self, s: u32, xs: &[u32], eterms: Option<&[Tgt]>) -> Vec<u32> {
        let t = Instant::now();
        let r = self.products_inner(s, xs, eterms);
        self.t_prod += t.elapsed().as_secs_f64();
        r
    }

    fn products_inner(&mut self, s: u32, xs: &[u32], eterms: Option<&[Tgt]>) -> Vec<u32> {
        let nl = self.k.live.len();
        let nf = self.vd.nf();
        let n = xs.len();
        let mut y = vec![0u32; n * nl];
        if nl == 0 || n == 0 {
            return y;
        }
        let a = self.gather_b(xs);
        self.follower(s);
        let fol = &self.fols[&s];
        gemm(&self.fp, &mut y, &a, &fol.data, n, nf, nl, false);
        if let Some(et) = eterms {
            let piv: Vec<(usize, u32)> =
                et.iter().enumerate().filter_map(|(i, t)| if let Tgt::Piv(r) = t { Some((i, *r)) } else { None }).collect();
            if !piv.is_empty() {
                let rows: Vec<u32> = piv.iter().map(|x| x.1).collect();
                let g = self.gather_b(&rows);
                let e = self.rho_neg_n(&g, piv.len()); // rho(-B[r])
                for (i, &(u, _)) in piv.iter().enumerate() {
                    for j in 0..nl {
                        y[u * nl + j] = self.fp.add(y[u * nl + j], e[i * nl + j]);
                    }
                }
            }
            for (u, t) in et.iter().enumerate() {
                if let Tgt::T(tc) = t {
                    self.k.rho_unit(&self.fp, *tc, &mut y[u * nl..(u + 1) * nl]);
                }
            }
        }
        y
    }

    fn insert(&mut self, cands: Vec<Vec<u32>>) -> usize {
        let t = Instant::now();
        self.c_candidates += cands.len() as u64;
        let r = self.k.insert(&self.fp, cands);
        self.t_ins += t.elapsed().as_secs_f64();
        r
    }

    /// one fall round: shift every pending fall row (its current state) by every u; returns #taken
    fn fall_round(&mut self) -> usize {
        let nsp = self.nsp;
        let nf = self.vd.nf();
        let pend: Vec<u32> =
            self.k.pivs.iter().copied().filter(|&t| t as usize >= nsp && !self.fall_done[t as usize]).collect();
        if pend.is_empty() {
            return 0;
        }
        self.c_rounds += 1;
        self.falls_new_per_round.push(pend.len() as u64);
        let mut phi = vec![0u32; pend.len() * nf];
        for (i, &t) in pend.iter().enumerate() {
            self.fall_done[t as usize] = true;
            let row = &self.k.rows[self.k.rowof[t as usize] as usize];
            phi[i * nf..(i + 1) * nf].copy_from_slice(&row[nsp..]);
        }
        let n = pend.len();
        let mut cands = Vec::new();
        for s in 0..self.ns {
            if Some(s) == self.opts.test_drop_shift {
                continue;
            }
            let nl = self.k.live.len();
            if nl == 0 {
                break;
            }
            // the states were captured above; insert per shift and per chunk of rows
            for c0 in (0..n).step_by(4096) {
                let c1 = (c0 + 4096).min(n);
                let nl = self.k.live.len();
                if nl == 0 {
                    break;
                }
                self.follower(s);
                let fol = &self.fols[&s];
                let mut y = vec![0u32; (c1 - c0) * nl];
                gemm(&self.fp, &mut y, &phi[c0 * nf..c1 * nf], &fol.data, c1 - c0, nf, nl, false);
                self.c_fall_rows += (c1 - c0) as u64;
                cands.extend(y.chunks(nl).map(|c| c.to_vec()));
                let cs = std::mem::take(&mut cands);
                self.insert(cs);
            }
        }
        n
    }

    fn fall_fixed_point(&mut self) {
        if self.opts.test_no_fall_shifts {
            return;
        }
        let t = Instant::now();
        while self.fall_round() > 0 {}
        self.t_falls += t.elapsed().as_secs_f64();
    }

    pub fn run(&mut self) {
        let d = self.d;
        let ns = self.ns;
        let mon = self.lay.mon;
        // coverage check (combinatorial): every divisor of a regular top is a pivot of weight D
        let lay = self.lay;
        let vd = self.vd;
        let bad = self
            .reg_tops
            .par_iter()
            .filter(|&&t| {
                let mut b = t;
                while b != 0 {
                    let bit = b & b.wrapping_neg();
                    b ^= bit;
                    let c = lay.col(t ^ bit).unwrap() as usize;
                    if vd.prow[c] == NONE {
                        return true;
                    }
                }
                false
            })
            .count();
        assert_eq!(bad, 0, "a regular top has a free divisor");
        // regular weight-D rows and their direct units
        let regx: Vec<u32> = (0..vd.rank() as u32).filter(|&r| lay.wt[vd.piv[r as usize] as usize] as u32 == d).collect();
        let lowx: Vec<u32> = if self.opts.include_low {
            (0..vd.rank() as u32).filter(|&r| (lay.wt[vd.piv[r as usize] as usize] as u32) < d).collect()
        } else {
            vec![]
        };
        // direct units: (row, s, eterm)
        let mut grouped_from_x = 0u64;
        let mut direct: Vec<Vec<(u32, Tgt)>> = vec![Vec::new(); ns as usize];
        for &r in regx.iter().chain(lowx.iter()) {
            let m = lay.masks[vd.piv[r as usize] as usize];
            for s in 0..ns {
                match self.tgt(mon.shift(m, s)) {
                    Tgt::Reg => grouped_from_x += 1,
                    t => direct[s as usize].push((r, t)),
                }
            }
        }
        assert_eq!(
            grouped_from_x,
            (d as u64 + 1) * self.reg_tops.len() as u64,
            "generator count of the regular tops is not D + 1 each"
        );
        let tops = std::mem::take(&mut self.reg_tops);
        let ntops = tops.len();
        let mut pos = 0usize;
        let mut bi = 0usize;
        let work = (self.opts.work_gb * 1e9) as usize;
        loop {
            if pos >= ntops && bi > 0 {
                break;
            }
            // batch size from the work budget: generator block (D + 1 rows per top) and its products
            let nl0 = self.k.live.len().max(1);
            let nbt = self.opts.batch_tops.max(1).min((work / ((d as usize + 1) * nl0 * 4 * 3)).max(1));
            let end = (pos + nbt).min(ntops);
            if self.k.live.is_empty() {
                // saturated: K = T, every further vector reduces to zero
                self.c_gen_rows += 0;
                break;
            }
            // --- grouped units of this batch of regular tops
            let tb = &tops[pos..end];
            // per shift: (row x, top index in batch, slot)
            let mut per_s: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); ns as usize];
            for (ti, &t) in tb.iter().enumerate() {
                let mut b = t;
                let mut slot = 0u32;
                while b != 0 {
                    let bit = b & b.wrapping_neg();
                    b ^= bit;
                    let s = bit.trailing_zeros();
                    let x = t ^ bit;
                    let r = vd.prow[lay.col(x).unwrap() as usize];
                    per_s[s as usize].push((r, ti as u32, slot));
                    slot += 1;
                }
            }
            let nl = self.k.live.len();
            let nper = d as usize + 1;
            let mut gens = vec![0u32; tb.len() * nper * nl];
            let mut have = vec![false; tb.len() * nper];
            for s in 0..ns {
                if Some(s) == self.opts.test_drop_shift {
                    continue;
                }
                let us = &per_s[s as usize];
                if us.is_empty() {
                    continue;
                }
                let xs: Vec<u32> = us.iter().map(|u| u.0).collect();
                let y = self.products(s, &xs, None);
                self.c_gen_rows += xs.len() as u64;
                for (i, &(_, ti, sl)) in us.iter().enumerate() {
                    let o = (ti as usize * nper + sl as usize) * nl;
                    gens[o..o + nl].copy_from_slice(&y[i * nl..(i + 1) * nl]);
                    have[ti as usize * nper + sl as usize] = true;
                }
            }
            // Lemma 5 differences
            let fp = self.fp;
            let cands: Vec<Vec<u32>> = (0..tb.len())
                .into_par_iter()
                .flat_map_iter(|ti| {
                    let base = ti * nper;
                    let first = (0..nper).find(|&s| have[base + s]);
                    let gens = &gens;
                    let have = &have;
                    (0..nper).filter_map(move |sl| {
                        let f = first?;
                        if sl == f || !have[base + sl] {
                            return None;
                        }
                        let a = &gens[(base + sl) * nl..(base + sl + 1) * nl];
                        let b = &gens[(base + f) * nl..(base + f + 1) * nl];
                        Some(a.iter().zip(b).map(|(&x, &y)| fp.sub(x, y)).collect::<Vec<u32>>())
                    })
                })
                .collect();
            drop(gens);
            self.insert(cands);
            // --- direct units: the share of this batch
            for s in 0..ns {
                if Some(s) == self.opts.test_drop_shift {
                    continue;
                }
                let all = &direct[s as usize];
                let lo = if ntops == 0 { 0 } else { all.len() * pos / ntops };
                let hi = if ntops == 0 { all.len() } else { all.len() * end / ntops };
                if lo == hi || self.k.live.is_empty() {
                    continue;
                }
                let part: Vec<(u32, Tgt)> = all[lo..hi].to_vec();
                for chunk in part.chunks(4096) {
                    let xs: Vec<u32> = chunk.iter().map(|u| u.0).collect();
                    let et: Vec<Tgt> = chunk.iter().map(|u| u.1).collect();
                    let y = self.products(s, &xs, Some(&et));
                    self.c_direct_rows += xs.len() as u64;
                    let nl = self.k.live.len();
                    let c: Vec<Vec<u32>> = if nl == 0 { vec![] } else { y.chunks(nl).map(|v| v.to_vec()).collect() };
                    self.insert(c);
                }
            }
            if self.opts.fall_every > 0 && (bi + 1) % self.opts.fall_every == 0 {
                self.fall_fixed_point();
            }
            pos = end;
            bi += 1;
        }
        self.fall_fixed_point();
    }

    /// shift ALL final fall rows again: every image must already be in K
    fn verify_falls(&mut self) -> (u64, u64) {
        let nsp = self.nsp;
        let nf = self.vd.nf();
        let falls: Vec<u32> = self.k.pivs.iter().copied().filter(|&t| t as usize >= nsp).collect();
        let mut nonzero = 0u64;
        let mut checked = 0u64;
        if self.k.live.is_empty() || falls.is_empty() {
            return (checked, nonzero);
        }
        for chunk in falls.chunks(4096) {
            let n = chunk.len();
            let mut phi = vec![0u32; n * nf];
            for (i, &t) in chunk.iter().enumerate() {
                let row = &self.k.rows[self.k.rowof[t as usize] as usize];
                phi[i * nf..(i + 1) * nf].copy_from_slice(&row[nsp..]);
            }
            for s in 0..self.ns {
                let nl = self.k.live.len();
                self.follower(s);
                let fol = &self.fols[&s];
                let mut y = vec![0u32; n * nl];
                gemm(&self.fp, &mut y, &phi, &fol.data, n, nf, nl, false);
                checked += n as u64;
                nonzero += y.chunks(nl).filter(|v| v.iter().any(|&x| x != 0)).count() as u64;
            }
        }
        (checked, nonzero)
    }

    /// final statements of the theorem
    pub fn finish(&mut self, t0: Instant) -> (serde_json::Value, Option<WSpace>) {
        let (vchecked, vnonzero) = if self.opts.verify_falls && !self.opts.test_no_fall_shifts {
            self.verify_falls()
        } else {
            (0, 0)
        };
        let lay = self.lay;
        let vd = self.vd;
        let n = lay.ncols();
        let nsp = self.nsp;
        let nf = vd.nf();
        let d = self.d;
        let rank_k = self.k.rank();
        let fall_t: Vec<u32> = self.k.pivs.iter().copied().filter(|&t| t as usize >= nsp).collect();
        let one_col = (n - 1) as u32;
        let one_in_vd = vd.prow[one_col as usize] != NONE;
        let contains_one =
            one_in_vd || { let t = nsp + vd.fpos[one_col as usize] as usize; self.k.rowof[t] != NONE };
        let n_reg = self.reg_tops_len();
        let dim_vd1 = vd.rank() as u64 + n_reg as u64 + rank_k as u64;
        // W = V_D (+) F
        let mut wpiv: Vec<u32> = vd.piv.clone();
        wpiv.extend(fall_t.iter().map(|&t| vd.free[t as usize - nsp]));
        wpiv.sort();
        let mut profile = crate::basis::profile_of(lay, &wpiv, d);
        profile.push(dim_vd1);
        let mut w = None;
        let mut sha = serde_json::Value::Null;
        if self.opts.hash || self.opts.parts {
            // RREF of W: fall rows (RREF on N, in N order = column order) and V_D rows reduced by them
            let fallf: Vec<usize> = fall_t.iter().map(|&t| t as usize - nsp).collect();
            let mut isfall = vec![false; nf];
            for &f in &fallf {
                isfall[f] = true;
            }
            let keepf: Vec<usize> = (0..nf).filter(|&f| !isfall[f]).collect();
            let nk = keepf.len();
            let nfl = fallf.len();
            // fall matrix restricted: (nfl x nk)
            let mut fm = vec![0u32; nfl * nk];
            for (i, &t) in fall_t.iter().enumerate() {
                let row = &self.k.rows[self.k.rowof[t as usize] as usize][nsp..];
                for (j, &f) in keepf.iter().enumerate() {
                    fm[i * nk + j] = row[f];
                }
            }
            let rank = vd.rank();
            // rows of W in pivot order, computed and hashed in chunks (no rank x #N copy)
            let mut rows_of: Vec<(u32, usize, bool)> = Vec::with_capacity(wpiv.len());
            for i in 0..rank {
                rows_of.push((vd.piv[i], i, false));
            }
            for (i, &f) in fallf.iter().enumerate() {
                rows_of.push((vd.free[f], i, true));
            }
            rows_of.sort();
            let wfree: Vec<u32> = keepf.iter().map(|&f| vd.free[f]).collect();
            let mut h = crate::basis::rref_sha_begin(lay, &wpiv, &wfree);
            let keep_w = self.opts.parts;
            let mut wb: Vec<u32> = Vec::new();
            for chunk in rows_of.chunks(4096) {
                let n = chunk.len();
                let mut o = vec![0u32; n * nk];
                let mut a = vec![0u32; n * nfl.max(1)];
                if nk > 0 {
                    o.par_chunks_mut(nk).zip(a.par_chunks_mut(nfl.max(1)))
                        .zip(chunk.par_iter())
                        .for_each(|((oo, aa), &(_, i, fall))| {
                            if fall {
                                oo.copy_from_slice(&fm[i * nk..(i + 1) * nk]);
                            } else {
                                let row = vd.row(i);
                                for (j, &f) in keepf.iter().enumerate() {
                                    oo[j] = row[f];
                                }
                                for (j, &f) in fallf.iter().enumerate() {
                                    aa[j] = row[f];
                                }
                            }
                        });
                    if nfl > 0 {
                        gemm(&self.fp, &mut o, &a, &fm, n, nfl, nk, true);
                    }
                }
                crate::basis::rref_sha_rows(&mut h, &o);
                if keep_w {
                    wb.extend_from_slice(&o);
                }
            }
            sha = json!(crate::basis::sha_finish(h));
            if keep_w {
                w = Some(WSpace { piv: wpiv.clone(), free: wfree, b: wb });
            }
        }
        let mut parts = serde_json::Value::Null;
        if self.opts.parts {
            if let Some(ws) = &w {
                parts = sigma_parts(&self.fp, lay, ws, d);
            }
        }
        let rec = json!({
            "lift": "hky_lift/0.1 (monomial coordinates, no sigma split)",
            "k": lay.k, "D": d, "D_target": d + 1, "p": self.fp.p,
            "n_cols_VD": n,
            "dim_VD": vd.rank(),
            "free_columns_VD": nf,
            "regular_weight_D_rows": (0..vd.rank()).filter(|&i| lay.wt[vd.piv[i] as usize] as u32 == d).count(),
            "free_weight_D_monomials": vd.free.iter().filter(|&&c| lay.wt[c as usize] as u32 == d).count(),
            "special_tops": nsp,
            "regular_tops": n_reg,
            "dim_T": self.nt,
            "dim_K": rank_k,
            "falls": fall_t.len(),
            "falls_new_per_round": self.falls_new_per_round,
            "rounds": self.c_rounds,
            "saturated": rank_k == self.nt,
            "fixed_point": !self.opts.test_no_fall_shifts && self.opts.test_drop_shift.is_none(),
            "verify_falls": {"rows_checked": vchecked, "nonzero": vnonzero, "ran": self.opts.verify_falls},
            "contains_one": contains_one,
            "one_in_VD": one_in_vd,
            "dim_VD1": dim_vd1,
            "dim_low": wpiv.len(),
            "profile": profile,
            "low_basis_sha256": sha,
            "sigma_parts": parts,
            "counts": {
                "generator_rows": self.c_gen_rows, "direct_rows": self.c_direct_rows,
                "fall_rows": self.c_fall_rows, "candidates": self.c_candidates,
                "follower_builds": self.c_follower_builds, "follower_evictions": self.c_follower_evictions,
                "follower_sync_madds": self.c_sync_madds,
                "seconds_followers": self.t_fol, "seconds_products_incl_followers": self.t_prod,
                "seconds_K_insert": self.t_ins, "seconds_fall_phase": self.t_falls,
                "madds_total": madds(),
            },
            "memory_bytes": {
                "B": vd.rank() * nf * 4,
                "K_full_width": rank_k * self.nt * 4,
                "followers_peak": self.peak_follower_bytes,
            },
            "settings": format!("{:?}", self.opts),
            "seconds": t0.elapsed().as_secs_f64(),
        });
        (rec, w)
    }

    fn reg_tops_len(&self) -> usize {
        // reg_tops was moved out in run(); recompute the count
        let total = crate::layout::binom(2 * self.lay.k as u64, self.d as u64 + 1) as usize;
        total - self.nsp
    }

    pub fn sp_masks(&self) -> &[u64] {
        &self.sp
    }
}

pub struct WSpace {
    pub piv: Vec<u32>,
    pub free: Vec<u32>,
    pub b: Vec<u32>,
}

/// dims of the sigma-invariant / anti-invariant parts of W and of W n M_<=j (small cases only)
pub fn sigma_parts(fp: &Fp, lay: &Layout, w: &WSpace, d: u32) -> serde_json::Value {
    let n = lay.ncols();
    let nk = w.free.len();
    let mon = lay.mon;
    let sig: Vec<usize> = lay.masks.iter().map(|&m| lay.col(mon.sigma(m)).unwrap() as usize).collect();
    let mut out = vec![];
    for sign in [1i64, -1] {
        let mut rows = Vec::new();
        let mut piv = Vec::new();
        let mut prof = Vec::new();
        for j in 0..=d {
            let vs: Vec<Vec<u32>> = (0..w.piv.len())
                .filter(|&i| lay.wt[w.piv[i] as usize] as u32 == j)
                .map(|i| {
                    let mut v = vec![0u32; n];
                    v[w.piv[i] as usize] = 1;
                    for f in 0..nk {
                        v[w.free[f] as usize] = w.b[i * nk + f];
                    }
                    let mut u = v.clone();
                    for c in 0..n {
                        let s = if sign == 1 { v[sig[c]] } else { fp.neg(v[sig[c]]) };
                        u[c] = fp.add(v[c], s);
                    }
                    u
                })
                .collect();
            rref_extend(fp, &mut rows, &mut piv, vs, n);
            prof.push(rows.len());
        }
        out.push(prof);
    }
    let total: usize = out[0].last().unwrap() + out[1].last().unwrap();
    json!({"part_profiles_low": out, "sigma_stable": total == w.piv.len()})
}

/// Random checks on the input (Lemma 2 needs V_D closed; the theorem needs g in V_D).
pub fn input_checks(fp: &Fp, lay: &Layout, vd: &Rref, g: Option<&[(u64, u32)]>, trials: usize, seed: u64) -> serde_json::Value {
    let n = lay.ncols();
    let d = lay.d;
    let mon = lay.mon;
    let rank = vd.rank();
    let nf = vd.nf();
    let mut s = seed | 1;
    let mut rnd = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    // closedness: random combination of low rows, every shift, residual must be zero
    let low: Vec<usize> = (0..rank).filter(|&i| (lay.wt[vd.piv[i] as usize] as u32) < d).collect();
    let mut bad = 0u64;
    let mut tested = 0u64;
    for _ in 0..trials {
        if low.is_empty() {
            break;
        }
        let mut coef = vec![0u32; rank];
        for &i in &low {
            coef[i] = (rnd() % fp.p) as u32;
        }
        let mut freev = vec![0u32; nf];
        gemm(fp, &mut freev, &coef, &vd.b, 1, rank, nf, false);
        let mut v = vec![0u32; n];
        for i in 0..rank {
            v[vd.piv[i] as usize] = coef[i];
        }
        for f in 0..nf {
            v[vd.free[f] as usize] = freev[f];
        }
        for sh in 0..mon.nshifts() {
            let mut u = vec![0u32; n];
            for c in 0..n {
                if v[c] != 0 {
                    assert!((lay.wt[c] as u32) < d);
                    u[lay.col(mon.shift(lay.masks[c], sh)).unwrap() as usize] = v[c];
                }
            }
            tested += 1;
            if vd.residual(fp, &u).iter().any(|&x| x != 0) {
                bad += 1;
            }
        }
    }
    let mut gres = serde_json::Value::Null;
    if let Some(g) = g {
        let gdeg = g.iter().map(|t| t.0.count_ones()).max().unwrap();
        let mut bad_g = 0;
        let mut tested_g = 0;
        // g itself, and random combinations of the m g
        let mut v = vec![0u32; n];
        for &(e, c) in g {
            v[lay.col(e).unwrap() as usize] = c;
        }
        tested_g += 1;
        if vd.residual(fp, &v).iter().any(|&x| x != 0) {
            bad_g += 1;
        }
        if gdeg <= d {
            for _ in 0..trials {
                let mut v = vec![0u32; n];
                for w in 0..=(d - gdeg) {
                    for m in mon.of_weight(w) {
                        let f = (rnd() % fp.p) as u32;
                        for (e, c) in crate::closure::translate(&mon, g, m) {
                            let col = lay.col(e).unwrap() as usize;
                            v[col] = fp.add(v[col], fp.mul(f, c));
                        }
                    }
                }
                tested_g += 1;
                if vd.residual(fp, &v).iter().any(|&x| x != 0) {
                    bad_g += 1;
                }
            }
        }
        gres = json!({"tested": tested_g, "not_in_VD": bad_g});
    }
    json!({"closed_random": {"shifted_combinations": tested, "leave_VD": bad}, "generators": gres})
}

pub fn rank_dense(fp: &Fp, vs: Vec<Vec<u32>>, c: usize) -> usize {
    rank_of(fp, vs, c)
}
