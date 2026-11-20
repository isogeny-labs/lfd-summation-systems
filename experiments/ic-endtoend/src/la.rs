//! Sparse linear algebra mod a prime N: 2-core peeling (first step of
//! structured Gaussian elimination), scalar Wiedemann, signed union-find.
use crate::field::*;

pub type Row = Vec<(u32, i64)>; // (column, small integer coefficient)

#[derive(Default, Clone, Copy, Debug, serde::Serialize)]
pub struct LaCost { pub muladds: u64, pub bookkeeping: u64, pub n_core: u64, pub rows_core: u64, pub nnz_core: u64, pub attempts: u64 }

/// Peel rows that hold a column of weight 1 until none is left; returns
/// the indices of the remaining rows (the 2-core).
pub fn peel(rows: &[Row], ncols: usize, cost: &mut LaCost) -> Vec<usize> {
    let mut w = vec![0u32; ncols];
    let mut col_rows: Vec<Vec<u32>> = vec![Vec::new(); ncols];
    for (i, r) in rows.iter().enumerate() {
        for &(c, _) in r { w[c as usize] += 1; col_rows[c as usize].push(i as u32); cost.bookkeeping += 1; }
    }
    let mut alive = vec![true; rows.len()];
    let mut stack: Vec<u32> = (0..ncols as u32).filter(|&c| w[c as usize] == 1).collect();
    while let Some(c) = stack.pop() {
        cost.bookkeeping += 1;
        if w[c as usize] != 1 { continue; }
        let r = *col_rows[c as usize].iter().find(|&&r| alive[r as usize]).unwrap();
        alive[r as usize] = false;
        for &(c2, _) in &rows[r as usize] {
            cost.bookkeeping += 1;
            w[c2 as usize] -= 1;
            if w[c2 as usize] == 1 { stack.push(c2); }
        }
    }
    (0..rows.len()).filter(|&i| alive[i]).collect()
}

fn coef(c: i64, n: u64) -> u64 { if c >= 0 { (c as u64) % n } else { n - ((-c) as u64 % n) } }

/// Berlekamp-Massey over F_N: returns connection polynomial C (C[0] = 1) and L.
fn bm(s: &[u64], n: u64, cost: &mut LaCost) -> (Vec<u64>, usize) {
    let mut c = vec![1u64]; let mut b = vec![1u64];
    let (mut l, mut m, mut bb) = (0usize, 1usize, 1u64);
    for i in 0..s.len() {
        let mut d = s[i];
        for j in 1..=l.min(c.len() - 1) { d = addmod(d, mulmod(c[j], s[i - j], n), n); cost.muladds += 1; }
        if d == 0 { m += 1; continue; }
        let coefv = mulmod(d, invmod(bb, n), n);
        let t = c.clone();
        if c.len() < b.len() + m { c.resize(b.len() + m, 0); }
        for j in 0..b.len() { c[j + m] = submod(c[j + m], mulmod(coefv, b[j], n), n); cost.muladds += 1; }
        if 2 * l <= i { l = i + 1 - l; b = t; bb = d; m = 1; } else { m += 1; }
    }
    c.resize(l + 1, 0);
    (c, l)
}

/// Kernel vector of B (B = M_core^T, square-padded) where `rows` are the
/// relations (columns of B). Returns lambda over the given rows.
pub fn kernel(rows: &[Row], ncols: usize, n: u64, seed: u64, cost: &mut LaCost) -> Option<Vec<u64>> {
    // B: ncols x nrows; apply B*v (v over relations) -> over columns
    let nr = rows.len();
    let dim = nr.max(ncols);
    let rows_m: Vec<Vec<(usize, u64)>> = rows.iter().map(|r| r.iter().map(|&(c, v)| (c as usize, coef(v, n))).collect()).collect();
    let nnz: u64 = rows_m.iter().map(|r| r.len() as u64).sum();
    cost.n_core = ncols as u64; cost.rows_core = nr as u64; cost.nnz_core = nnz;
    // A: dim x dim, A v = (B v padded with zeros)
    let apply = |v: &[u64], cost: &mut LaCost| -> Vec<u64> {
        let mut out = vec![0u64; dim];
        for (j, r) in rows_m.iter().enumerate() {
            if v[j] == 0 { continue; }
            for &(c, a) in r { out[c] = addmod(out[c], mulmod(a, v[j], n), n); cost.muladds += 1; }
        }
        out
    };
    let mut rng = Rng(seed);
    for _attempt in 0..4 {
        cost.attempts += 1;
        let u: Vec<u64> = (0..dim).map(|_| rng.below(n)).collect();
        let y: Vec<u64> = (0..dim).map(|_| if (0..dim).len() > 0 { rng.below(n) } else { 0 }).collect();
        let y: Vec<u64> = y.iter().enumerate().map(|(i, &v)| if i < nr { v } else { 0 }).collect();
        let mut seq = Vec::with_capacity(2 * dim + 2);
        let mut w = y.clone();
        for _ in 0..(2 * dim + 2) {
            let mut s = 0u64;
            for i in 0..dim { if w[i] != 0 && u[i] != 0 { s = addmod(s, mulmod(u[i], w[i], n), n); } }
            cost.muladds += dim as u64;
            seq.push(s);
            w = apply(&w, cost);
        }
        let (cpoly, l) = bm(&seq, n, cost);
        // minimal polynomial g(X) = X^l C(1/X): g_j = C[l - j]
        let g: Vec<u64> = (0..=l).map(|j| cpoly[l - j]).collect();
        let e = g.iter().take_while(|&&v| v == 0).count();
        if e == 0 || e > l { continue; }
        let h = &g[e..]; // h(X), h(0) != 0
        // z = h(A) y by Horner
        let mut z = vec![0u64; dim];
        for j in (0..h.len()).rev() {
            z = apply(&z, cost);
            for i in 0..dim { if y[i] != 0 { z[i] = addmod(z[i], mulmod(h[j], y[i], n), n); } }
            cost.muladds += dim as u64;
        }
        let mut cur = z;
        for _ in 0..=e {
            if cur.iter().all(|&v| v == 0) { break; }
            let nx = apply(&cur, cost);
            if nx.iter().all(|&v| v == 0) {
                return Some(cur[..nr].to_vec());
            }
            cur = nx;
        }
    }
    None
}

/// Signed union-find: rows with two entries are edges; reports when a row
/// does not raise the rank (a dependency exists).
pub struct SignedUf { parent: Vec<u32>, par: Vec<i8>, unbal: Vec<bool> }
impl SignedUf {
    pub fn new(n: usize) -> Self { SignedUf { parent: (0..n as u32).collect(), par: vec![1; n], unbal: vec![false; n] } }
    fn find(&mut self, x: u32, ops: &mut u64) -> (u32, i8) {
        let mut path = Vec::new();
        let mut cur = x;
        let mut s = 1i8;
        while self.parent[cur as usize] != cur { path.push(cur); s *= self.par[cur as usize]; cur = self.parent[cur as usize]; *ops += 1; }
        // path compression with parity
        let root = cur;
        let mut acc = s;
        for &v in &path { let pv = self.par[v as usize]; self.par[v as usize] = acc; self.parent[v as usize] = root; acc *= pv; }
        (root, s)
    }
    /// add a row; returns true when it closes a dependency.
    pub fn add_row(&mut self, row: &Row, ops: &mut u64) -> bool {
        *ops += 1;
        match row.len() {
            0 => true,
            1 => { // half-edge: unbalances its component
                let (r, _) = self.find(row[0].0, ops);
                if self.unbal[r as usize] { true } else { self.unbal[r as usize] = true; false }
            }
            2 => {
                let (i, ci) = row[0]; let (j, cj) = row[1];
                if ci.abs() != 1 || cj.abs() != 1 {
                    // treat general coefficients as an unbalanced edge in the union (conservative: never claims a false dependency
                    // unless the component is already unbalanced and connected)
                }
                let w: i8 = if (ci > 0) == (cj > 0) { -1 } else { 1 }; // phi_j = w phi_i
                let (ri, si) = self.find(i, ops); let (rj, sj) = self.find(j, ops);
                if ri != rj {
                    // phi_j = sj*phi_rj, phi_i = si*phi_ri; want sj*phi_rj = w si phi_ri
                    self.parent[rj as usize] = ri;
                    self.par[rj as usize] = w * si * sj;
                    let u = self.unbal[ri as usize] || self.unbal[rj as usize];
                    self.unbal[ri as usize] = u;
                    false
                } else if si * sj == w && ci.abs() == 1 && cj.abs() == 1 {
                    true // balanced cycle
                } else if self.unbal[ri as usize] {
                    true // second unbalanced cycle
                } else { self.unbal[ri as usize] = true; false }
            }
            _ => false,
        }
    }
}
