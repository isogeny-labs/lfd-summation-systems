//! Direct closure V_n in M_<=n (the oracle for validation; small k only).
//!
//! V_n = least subspace of span(M_<=n) containing NF(m g) for weight(m) <= n - deg g ("gens all",
//! as `csd/ic_bounds_big.py`) and closed under f -> u_{i,j} f for f in V_n supported on weight < n.
//! Written from that definition, independent of the lift: a reduced echelon basis kept on a
//! compacting list of non-pivot columns; every low pivot row is shifted at its state when taken
//! (closure invariant), and at the end an exhaustive sweep (every shift of every final low row
//! reduces to zero) certifies closedness.

use crate::basis::{Echelon, Rref, NONE};
use crate::field::{gemm, rref_extend, Fp};
use crate::layout::{Layout, Mon};
use rayon::prelude::*;

/// S_3(x1, x2, X) for y^2 = x^3 + a x + b (csd/ic_bounds_d5_verify.py s3_coeffs).
pub fn s3_coeffs(fp: &Fp, a: i64, b: i64, x: i64) -> Vec<((u32, u32), u32)> {
    let p = fp.p as i64;
    let mut c: Vec<((u32, u32), i128)> = Vec::new();
    let mut add = |e: (u32, u32), v: i128| {
        if let Some(t) = c.iter_mut().find(|t| t.0 == e) {
            t.1 += v;
        } else {
            c.push((e, v));
        }
    };
    let (a, b, x) = (a as i128, b as i128, x as i128);
    add((2, 0), x * x);
    add((1, 1), -2 * x * x);
    add((0, 2), x * x);
    add((2, 1), -2 * x);
    add((1, 2), -2 * x);
    add((1, 0), -2 * x * a);
    add((0, 1), -2 * x * a);
    add((0, 0), -4 * x * b);
    add((2, 2), 1);
    add((1, 1), -2 * a);
    add((0, 0), a * a);
    add((1, 0), -4 * b);
    add((0, 1), -4 * b);
    c.into_iter()
        .map(|(e, v)| (e, v.rem_euclid(p as i128) as u32))
        .filter(|t| t.1 != 0)
        .collect()
}

/// NF of S_3 modulo the chains as {mask: coeff}: x_i^e -> exponent e mod 2^k in block i.
pub fn chain_generator(fp: &Fp, k: u32, a: i64, b: i64, x: i64) -> Vec<(u64, u32)> {
    let n = 1u64 << k;
    let mut g: Vec<(u64, u32)> = Vec::new();
    for ((e1, e2), v) in s3_coeffs(fp, a, b, x) {
        let m = (e1 as u64 % n) | ((e2 as u64 % n) << k);
        if let Some(t) = g.iter_mut().find(|t| t.0 == m) {
            t.1 = fp.add(t.1, v);
        } else {
            g.push((m, v));
        }
    }
    g.retain(|t| t.1 != 0);
    g.sort();
    g
}

/// translate g by monomial m (group multiplication: blockwise addition mod 2^k)
pub fn translate(mon: &Mon, g: &[(u64, u32)], m: u64) -> Vec<(u64, u32)> {
    let k = mon.k;
    let nb = (1u64 << k) - 1;
    g.iter()
        .map(|&(e, v)| {
            let a1 = ((e & nb) + (m & nb)) & nb;
            let a2 = (((e >> k) & nb) + ((m >> k) & nb)) & nb;
            (a1 | (a2 << k), v)
        })
        .collect()
}

/// Reduced echelon basis on a compacting list of columns.
pub struct Basis {
    pub fp: Fp,
    pub n: usize,
    pub fcols: Vec<u32>,
    pub fpos: Vec<u32>,
    pub rows: Vec<Vec<u32>>,
    pub piv: Vec<u32>,
    pub prow: Vec<u32>,
    dead: usize,
}

impl Basis {
    pub fn new(fp: Fp, n: usize) -> Basis {
        Basis {
            fp,
            n,
            fcols: (0..n as u32).collect(),
            fpos: (0..n as u32).collect(),
            rows: Vec::new(),
            piv: Vec::new(),
            prow: vec![NONE; n],
            dead: 0,
        }
    }
    pub fn rank(&self) -> usize {
        self.rows.len()
    }
    /// residuals (on fcols) of full-width vectors
    pub fn reduce(&self, vs: &[Vec<u32>]) -> Vec<Vec<u32>> {
        let nv = vs.len();
        let c = self.fcols.len();
        let r = self.rows.len();
        let mut out = vec![0u32; nv * c];
        let mut a = vec![0u32; nv * r];
        for (i, v) in vs.iter().enumerate() {
            for (j, &fc) in self.fcols.iter().enumerate() {
                out[i * c + j] = v[fc as usize];
            }
            for (t, &pc) in self.piv.iter().enumerate() {
                a[i * r + t] = v[pc as usize];
            }
        }
        let bm: Vec<u32> = self.rows.concat();
        gemm(&self.fp, &mut out, &a, &bm, nv, r, c, true);
        out.chunks(c.max(1)).take(nv).map(|x| if c == 0 { vec![] } else { x.to_vec() }).collect()
    }
    /// add residuals (on fcols); returns the new pivot columns
    pub fn add(&mut self, res: Vec<Vec<u32>>) -> Vec<u32> {
        let c = self.fcols.len();
        let res: Vec<Vec<u32>> = res.into_iter().filter(|v| v.iter().any(|&x| x != 0)).collect();
        if res.is_empty() {
            return vec![];
        }
        let mut nr = Vec::new();
        let mut np = Vec::new();
        rref_extend(&self.fp, &mut nr, &mut np, res, c);
        // eliminate the new pivots from the old rows: old -= old[:, np] @ nr
        let r = self.rows.len();
        if r > 0 && !np.is_empty() {
            let mut a = vec![0u32; r * np.len()];
            for (i, row) in self.rows.iter().enumerate() {
                for (t, &q) in np.iter().enumerate() {
                    a[i * np.len() + t] = row[q];
                }
            }
            let mut flat: Vec<u32> = self.rows.concat();
            let bm: Vec<u32> = nr.concat();
            gemm(&self.fp, &mut flat, &a, &bm, r, np.len(), c, true);
            for (i, row) in self.rows.iter_mut().enumerate() {
                row.copy_from_slice(&flat[i * c..(i + 1) * c]);
            }
        }
        let mut out = Vec::new();
        for (row, q) in nr.into_iter().zip(np) {
            let col = self.fcols[q];
            self.prow[col as usize] = self.rows.len() as u32;
            self.rows.push(row);
            self.piv.push(col);
            out.push(col);
            self.dead += 1;
        }
        if self.dead * 4 > self.fcols.len() {
            self.compact();
        }
        out
    }
    fn compact(&mut self) {
        let keep: Vec<usize> = (0..self.fcols.len()).filter(|&j| self.prow[self.fcols[j] as usize] == NONE).collect();
        self.rows.par_iter_mut().for_each(|row| {
            let nrow: Vec<u32> = keep.iter().map(|&j| row[j]).collect();
            *row = nrow;
        });
        self.fcols = keep.iter().map(|&j| self.fcols[j]).collect();
        self.fpos = vec![NONE; self.n];
        for (j, &c) in self.fcols.iter().enumerate() {
            self.fpos[c as usize] = j as u32;
        }
        self.dead = 0;
    }
    /// full-width row
    pub fn full_row(&self, t: usize) -> Vec<u32> {
        let mut v = vec![0u32; self.n];
        for (j, &c) in self.fcols.iter().enumerate() {
            v[c as usize] = self.rows[t][j];
        }
        v[self.piv[t] as usize] = 1;
        v
    }
    /// convert to an Rref (pivots sorted, rows on the true free columns)
    pub fn to_rref(&self) -> Rref {
        let mut piv = self.piv.clone();
        piv.sort();
        let mut r = Rref::new_index(self.n, piv);
        let nf = r.nf();
        for (i, &pc) in r.piv.clone().iter().enumerate() {
            let t = self.prow[pc as usize] as usize;
            for (f, &c) in r.free.clone().iter().enumerate() {
                let j = self.fpos[c as usize];
                r.b[i * nf + f] = if j == NONE { 0 } else { self.rows[t][j as usize] };
            }
        }
        r
    }
}

pub struct ClosureOut {
    pub lay: Layout,
    pub rref: Rref,
    pub sweep_images: u64,
    pub closed: bool,
}

/// Direct closure V_n (order 0 columns of M_<=n).
pub fn closure(fp: Fp, k: u32, n: u32, g: &[(u64, u32)]) -> ClosureOut {
    let lay = Layout::order0(k, n);
    let mon = lay.mon;
    let nc = lay.ncols();
    let mut bs = Basis::new(fp, nc);
    let gdeg = g.iter().map(|t| t.0.count_ones()).max().unwrap();
    // generators m g, weight(m) <= n - gdeg
    let mut gens = Vec::new();
    if gdeg <= n {
        for w in 0..=(n - gdeg) {
            for m in mon.of_weight(w) {
                let mut v = vec![0u32; nc];
                for (e, c) in translate(&mon, g, m) {
                    let col = lay.col(e).expect("generator term outside M_<=n") as usize;
                    v[col] = fp.add(v[col], c);
                }
                gens.push(v);
            }
        }
    }
    for chunk in gens.chunks(512) {
        let r = bs.reduce(chunk);
        bs.add(r);
    }
    // shift tables on low columns (weight < n)
    let ns = mon.nshifts();
    let low0 = lay.wt.iter().position(|&w| (w as u32) < n).unwrap_or(nc);
    let shift_of = |c: usize, s: u32| -> usize { lay.col(mon.shift(lay.masks[c], s)).unwrap() as usize };
    let image = |v: &[u32], s: u32| -> Vec<u32> {
        let mut o = vec![0u32; nc];
        for c in low0..nc {
            if v[c] != 0 {
                o[shift_of(c, s)] = v[c];
            }
        }
        o
    };
    let mut processed = vec![false; nc];
    loop {
        let pend: Vec<u32> = bs
            .piv
            .iter()
            .copied()
            .filter(|&c| (c as usize) >= low0 && !processed[c as usize])
            .collect();
        if pend.is_empty() {
            break;
        }
        for chunk in pend.chunks(64) {
            let states: Vec<Vec<u32>> = chunk
                .iter()
                .map(|&c| {
                    processed[c as usize] = true;
                    bs.full_row(bs.prow[c as usize] as usize)
                })
                .collect();
            let imgs: Vec<Vec<u32>> = states
                .par_iter()
                .flat_map_iter(|st| (0..ns).map(move |s| (st, s)))
                .map(|(st, s)| image(st, s))
                .collect();
            let r = bs.reduce(&imgs);
            bs.add(r);
        }
    }
    // exhaustive certificate sweep
    let mut closed = false;
    let mut sweep_images = 0u64;
    for _round in 0..50 {
        let lows: Vec<u32> = bs.piv.iter().copied().filter(|&c| (c as usize) >= low0).collect();
        let mut bad = false;
        for chunk in lows.chunks(64) {
            let states: Vec<Vec<u32>> = chunk.iter().map(|&c| bs.full_row(bs.prow[c as usize] as usize)).collect();
            let imgs: Vec<Vec<u32>> = states
                .par_iter()
                .flat_map_iter(|st| (0..ns).map(move |s| (st, s)))
                .map(|(st, s)| image(st, s))
                .collect();
            sweep_images += imgs.len() as u64;
            let r = bs.reduce(&imgs);
            if r.iter().any(|v| v.iter().any(|&x| x != 0)) {
                bad = true;
                bs.add(r);
            }
        }
        if !bad {
            closed = true;
            break;
        }
    }
    let rref = bs.to_rref();
    ClosureOut { lay, rref, sweep_images, closed }
}

/// Export an RREF as an hky.basis/v1 echelon; with `scramble`, every row gets random multiples of
/// up to `scramble` rows with later pivots added (still echelon, no longer reduced).
pub fn to_echelon(fp: &Fp, lay: Layout, r: &Rref, d: u32, k: u32, scramble: usize, seed: u64) -> Echelon {
    let n = lay.ncols();
    let rank = r.rank();
    let nf = r.nf();
    let mut full: Vec<Vec<u32>> = (0..rank)
        .map(|i| {
            let mut v = vec![0u32; n];
            v[r.piv[i] as usize] = 1;
            for f in 0..nf {
                v[r.free[f] as usize] = r.b[i * nf + f];
            }
            v
        })
        .collect();
    if scramble > 0 {
        let mut s = seed | 1;
        let mut rnd = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        // process from the first row: adding later rows (original, reduced) keeps echelon
        let orig = full.clone();
        for i in 0..rank {
            for _ in 0..scramble {
                if i + 1 >= rank {
                    break;
                }
                let j = i + 1 + (rnd() as usize % (rank - i - 1));
                let f = (rnd() % (fp.p - 1) + 1) as u32;
                let nf_ = fp.neg(f);
                fp.axpy_sub(&mut full[i], nf_, &orig[j]);
            }
        }
    }
    let mut off = vec![0u64];
    let mut ents = Vec::new();
    for i in 0..rank {
        let pc = r.piv[i] as usize;
        assert_eq!(full[i][pc], 1);
        for c in pc + 1..n {
            if full[i][c] != 0 {
                ents.push((c as u32, full[i][c]));
            }
        }
        off.push(ents.len() as u64);
    }
    let profile = crate::basis::profile_of(&lay, &r.piv, d);
    let flags = if r.piv.last() == Some(&((n - 1) as u32)) { 1 } else { 0 };
    let nnz = ents.len() as u64;
    Echelon { p: fp.p, m: 2, k, d, column_order: 0, flags, profile, lay, piv: r.piv.clone(), off, ents, nnz, src: None }
}
