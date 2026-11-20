//! End-to-end PKM index calculus, m = 2 (pair-sum table, one lookup per
//! attempt) and m = 3 (pair-sum table + meet in the middle), with peeling and
//! Wiedemann, every operation counted.
use crate::curve::*;
use crate::field::*;
use crate::hash::HMap;
use crate::la::*;
use rayon::prelude::*;
use serde::Serialize;
use std::time::Instant;

#[derive(Serialize, Clone, Default)]
pub struct E2e {
    pub b: u32, pub p: u64, pub n: u64, pub m: u32, pub k: u32, pub f: usize, pub fb_points: usize,
    pub target_idx: u64, pub x_true: u64, pub x_found: u64, pub verified: bool,
    pub attempts: u64, pub relations: u64, pub dep_checks: u64,
    pub ga_setup: u64, pub ga_table: u64, pub ga_walk: u64, pub ga_search: u64, pub ga_resolve: u64,
    pub hash_inserts: u64, pub hash_lookups: u64, pub bookkeeping: u64,
    pub la: LaCost, pub la_bookkeeping: u64, pub table_entries: u64, pub table_bytes: u64,
    pub secs_table: f64, pub secs_search: f64, pub secs_la: f64,
    pub grob_cols_cert: f64, pub grob_d_cert: u32, pub grob_cols_ceiling: f64, pub grob_lower_paper: f64,
}
impl E2e {
    pub fn group_adds(&self) -> u64 { self.ga_setup + self.ga_table + self.ga_walk + self.ga_search + self.ga_resolve }
    pub fn native(&self) -> u64 { self.hash_inserts + self.hash_lookups + self.bookkeeping + self.la.muladds + self.la_bookkeeping }
}

pub fn binom(n: u64, k: u64) -> f64 { if k > n { return 0.0; } let mut r = 1.0; for i in 0..k { r = r * (n - i) as f64 / (i + 1) as f64; } r }
pub fn binom_le(n: u64, d: u64) -> f64 { (0..=d).map(|i| binom(n, i)).sum() }
/// certified LFD (m = 2) per the paper's table; conjecture floor((k+14)/4) beyond k = 14
pub fn lfd_cert(k: u32) -> u32 { match k { 0..=5 => 4, 6..=9 => 5, 10..=13 => 6, 14 => 7, _ => (k + 14) / 4 } }

pub struct Table { pub map: HMap, pub adds: u64, pub inserts: u64, pub secs: f64 }

/// x(P_i +- P_j) -> packed (i, j, sign)
pub fn build_table(c: &Curve, reps: &[Pt]) -> Table {
    let t0 = Instant::now();
    let parts: Vec<Vec<(u64, u64)>> = (0..reps.len()).into_par_iter().map_init(|| (Vec::new(), Vec::new(), Vec::new()), |(s, om, op), i| {
        batch_diff_sum(c, reps[i], &reps[i + 1..], s, om, op);
        let mut v = Vec::with_capacity(2 * om.len() + 1);
        let d = c.dbl(reps[i]);
        if !d.inf { v.push((d.x, pack(i, i, 1))); }
        for j in 0..om.len() { if om[j] != u64::MAX {
            v.push((op[j], pack(i, i + 1 + j, 1)));
            v.push((om[j], pack(i, i + 1 + j, 0)));
        } }
        v
    }).collect();
    let total: usize = parts.iter().map(|v| v.len()).sum();
    let mut map = HMap::with_capacity(total);
    for v in parts { for (x, val) in v { map.insert(x, val); } }
    Table { map, adds: total as u64, inserts: total as u64, secs: t0.elapsed().as_secs_f64() }
}
fn pack(i: usize, j: usize, s: u64) -> u64 { ((i as u64) << 32) | ((j as u64) << 1) | s }
fn unpack(v: u64) -> (usize, usize, i64) { ((v >> 32) as usize, ((v & 0xFFFF_FFFF) >> 1) as usize, if v & 1 == 1 { 1 } else { -1 }) }

/// The point P_i + s P_j and its row entries.
fn pair_point(c: &Curve, reps: &[Pt], v: u64) -> (Pt, Vec<(u32, i64)>) {
    let (i, j, s) = unpack(v);
    if i == j { (c.dbl(reps[i]), vec![(i as u32, 2)]) }
    else { let pj = if s > 0 { reps[j] } else { c.neg(reps[j]) }; (c.add(reps[i], pj), vec![(i as u32, 1), (j as u32, s)]) }
}
fn merge(mut row: Vec<(u32, i64)>) -> Row {
    row.sort_by_key(|e| e.0);
    let mut out: Row = Vec::new();
    for (cidx, v) in row { if let Some(l) = out.last_mut() { if l.0 == cidx { l.1 += v; continue; } } out.push((cidx, v)); }
    out.retain(|e| e.1 != 0);
    out
}

pub fn run(c: &Curve, b: u32, m: u32, k: u32, base: &Base, table: &Table, g: Pt, ti: u64, x_true: u64) -> E2e {
    let n = c.n;
    let reps = &base.reps;
    let f = reps.len();
    let q = c.mul(g, x_true).0;
    let mut e = E2e { b, p: c.f.p, n, m, k, f, fb_points: base.size(), target_idx: ti, x_true, ..Default::default() };
    e.ga_table = table.adds; e.hash_inserts = table.inserts; e.table_entries = table.map.len as u64; e.table_bytes = table.map.bytes() as u64; e.secs_table = table.secs;
    e.grob_d_cert = lfd_cert(k);
    e.grob_cols_cert = binom_le(2 * k as u64, lfd_cert(k) as u64);
    e.grob_cols_ceiling = binom_le(2 * k as u64, k as u64 + 2);
    let dlo = if k >= 10 { (k - 10) / 14 + 1 } else { 0 };
    e.grob_lower_paper = binom(2 * k as u64, dlo as u64);
    if m == 3 { e.grob_d_cert = 0; e.grob_cols_cert = f64::NAN; e.grob_cols_ceiling = binom_le(3 * k as u64, 3 * k as u64 + 1); e.grob_lower_paper = f64::NAN; }
    let mut rng = Rng(derive(0xE2E, (b as u64) << 40 | (m as u64) << 32 | (k as u64) << 16 | ti));
    // steps S_t = c_t G + d_t Q
    let mut steps = Vec::new();
    for _ in 0..32 {
        let (cc, dd) = (rng.below(n), rng.below(n));
        let (p1, o1) = c.mul(g, cc); let (p2, o2) = c.mul(q, dd);
        e.ga_setup += o1 + o2 + 1;
        steps.push((c.add(p1, p2), cc, dd));
    }
    let (a0, b0) = (1 + rng.below(n - 1), 1 + rng.below(n - 1));
    let (p1, o1) = c.mul(g, a0); let (p2, o2) = c.mul(q, b0);
    e.ga_setup += o1 + o2 + 1;
    let mut r = c.add(p1, p2);
    let (mut a, mut bb) = (a0, b0);
    let mut rows: Vec<Row> = Vec::new();
    let mut coefs: Vec<(u64, u64)> = Vec::new();
    let mut uf = SignedUf::new(f);
    let (mut s, mut om, mut op) = (Vec::new(), Vec::new(), Vec::new());
    let t0 = Instant::now();
    let delta = (f / 128).max(8);
    let mut next_check = (3 * f / 4) as u64;
    let mut la_secs = 0.0;
    loop {
        // advance the walk (random step index: no cycling)
        let t = rng.below(32) as usize;
        r = c.add(r, steps[t].0);
        a = addmod(a, steps[t].1, n); bb = addmod(bb, steps[t].2, n);
        e.ga_walk += 1; e.bookkeeping += 3; e.attempts += 1;
        if r.inf { continue; }
        let mut row: Option<Row> = None;
        if m == 2 {
            e.hash_lookups += 1;
            if let Some(v) = table.map.get(r.x) {
                let (sp, ent) = pair_point(c, reps, v);
                e.ga_resolve += 1;
                if sp == r { row = Some(merge(ent)); }
                else if sp == c.neg(r) { row = Some(merge(ent.into_iter().map(|(i, v)| (i, -v)).collect())); }
            }
        } else {
            batch_diff_sum(c, r, reps, &mut s, &mut om, &mut op);
            e.ga_search += 2 * f as u64;
            for l in 0..f {
                if om[l] == u64::MAX { continue; }
                for (xv, sl) in [(om[l], 1i64), (op[l], -1i64)] {
                    e.hash_lookups += 1;
                    if let Some(v) = table.map.get(xv) {
                        // D = R - sl*P_l = +-(P_i + s P_j)
                        let pl = if sl > 0 { reps[l] } else { c.neg(reps[l]) };
                        let d = c.sub(r, pl);
                        let (sp, ent) = pair_point(c, reps, v);
                        e.ga_resolve += 2;
                        let sign = if sp == d { 1 } else if sp == c.neg(d) { -1 } else { 0 };
                        if sign != 0 {
                            let mut ent2: Vec<(u32, i64)> = ent.into_iter().map(|(i, v)| (i, sign * v)).collect();
                            ent2.push((l as u32, sl));
                            row = Some(merge(ent2));
                        }
                    }
                    if row.is_some() { break; }
                }
                if row.is_some() { break; }
            }
        }
        let Some(row) = row else { continue };
        rows.push(row.clone()); coefs.push((a, bb));
        e.relations += 1;
        let try_solve = if m == 2 { uf.add_row(&row, &mut e.la_bookkeeping) } else if e.relations >= next_check {
            next_check += delta as u64; e.dep_checks += 1;
            let mut lc = LaCost::default();
            let core = peel(&rows, f, &mut lc);
            e.la_bookkeeping += lc.bookkeeping;
            let mut cols = std::collections::HashSet::new();
            for &i in &core { for &(cc, _) in &rows[i] { cols.insert(cc); } }
            core.len() > cols.len()
        } else { false };
        if !try_solve { continue; }
        let tl = Instant::now();
        let mut lc = LaCost::default();
        let core = peel(&rows, f, &mut lc);
        // reindex columns
        let mut idx = std::collections::HashMap::new();
        let core_rows: Vec<Row> = core.iter().map(|&i| rows[i].iter().map(|&(cc, v)| {
            let l = idx.len() as u32; let id = *idx.entry(cc).or_insert(l); (id, v) }).collect()).collect();
        let ker = kernel(&core_rows, idx.len(), n, derive(0x1A, e.relations), &mut lc);
        e.la.muladds += lc.muladds; e.la_bookkeeping += lc.bookkeeping; e.la.attempts += lc.attempts;
        e.la.n_core = lc.n_core; e.la.rows_core = lc.rows_core; e.la.nnz_core = lc.nnz_core;
        la_secs += tl.elapsed().as_secs_f64();
        if let Some(lam) = ker {
            let (mut num, mut den) = (0u64, 0u64);
            for (j, &ri) in core.iter().enumerate() {
                num = addmod(num, mulmod(lam[j], coefs[ri].0, n), n);
                den = addmod(den, mulmod(lam[j], coefs[ri].1, n), n);
            }
            e.la.muladds += 2 * core.len() as u64;
            if den != 0 {
                let x = mulmod((n - num) % n, invmod(den, n), n);
                if c.mul(g, x).0 == q {
                    e.x_found = x; e.verified = x == x_true;
                    break;
                }
            }
        }
        if m == 2 { /* keep collecting; the union-find stays in dependent state */ }
    }
    e.secs_search = t0.elapsed().as_secs_f64() - la_secs;
    e.secs_la = la_secs;
    e
}
