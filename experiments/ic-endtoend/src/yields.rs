//! Exact relation-yield measurement (census for m = 2, sampling with hash
//! lookups / meet in the middle for m = 2, 3). No Groebner bases.
use crate::curve::*;
use crate::field::*;
use crate::hash::HMap;
use rayon::prelude::*;
use serde::Serialize;
use std::time::Instant;

#[derive(Serialize, Clone)]
pub struct YieldCell {
    pub b: u32, pub p: u64, pub n: u64, pub curve: usize, pub k: u32, pub arm: String, pub m: u32,
    pub method: String, pub f: usize, pub fb_points: usize, pub lambda: f64, pub model_pr: f64,
    pub targets: u64, pub hits: u64, pub mean_reps: f64, pub p_hat: f64, pub p_lo: f64, pub p_hi: f64,
    pub ratio: f64, pub ratio_lo: f64, pub ratio_hi: f64, pub element_tests: u64, pub secs: f64,
}

pub fn wilson(h: u64, n: u64, z: f64) -> (f64, f64) {
    if n == 0 { return (0.0, 1.0); }
    let nf = n as f64; let ph = h as f64 / nf;
    let d = 1.0 + z * z / nf;
    let c = (ph + z * z / (2.0 * nf)) / d;
    let w = z * ((ph * (1.0 - ph) / nf + z * z / (4.0 * nf * nf)).sqrt()) / d;
    ((c - w).max(0.0), (c + w).min(1.0))
}

pub fn fact(m: u32) -> f64 { (1..=m).map(|i| i as f64).product() }

/// x-coordinates of all P_i +- P_j (i <= j), O excluded, as a Vec (with repeats).
pub fn pair_xs(c: &Curve, reps: &[Pt]) -> Vec<Vec<u64>> {
    (0..reps.len()).into_par_iter().map_init(|| (Vec::new(), Vec::new(), Vec::new()), |(s, om, op), i| {
        batch_diff_sum(c, reps[i], &reps[i + 1..], s, om, op);
        let mut v = Vec::with_capacity(2 * (reps.len() - i));
        let d = c.dbl(reps[i]);
        if !d.inf { v.push(d.x); }
        for j in 0..om.len() { if om[j] != u64::MAX { v.push(om[j]); v.push(op[j]); } }
        v
    }).collect()
}

fn finish(mut cell: YieldCell, z: f64) -> YieldCell {
    cell.p_hat = cell.hits as f64 / cell.targets.max(1) as f64;
    let (lo, hi) = wilson(cell.hits, cell.targets, z);
    cell.p_lo = lo; cell.p_hi = hi;
    cell.ratio = cell.p_hat / cell.model_pr; cell.ratio_lo = lo / cell.model_pr; cell.ratio_hi = hi / cell.model_pr;
    cell
}

pub fn census_m2(c: &Curve, base: &Base, hdr: YieldCell) -> YieldCell {
    let t0 = Instant::now();
    let parts = pair_xs(c, &base.reps);
    let mut all: Vec<u64> = parts.into_iter().flatten().collect();
    let tests = all.len() as u64;
    all.par_sort_unstable();
    all.dedup();
    let distinct = all.len() as u64;
    let mut cell = hdr;
    cell.method = "census".into();
    // exact: Pr[uniform nonzero R in sum set] = 2|S_x| / (N - 1)
    cell.p_hat = 2.0 * distinct as f64 / (c.n - 1) as f64;
    cell.p_lo = cell.p_hat; cell.p_hi = cell.p_hat;
    cell.targets = c.n - 1; cell.hits = 2 * distinct;
    cell.mean_reps = 2.0 * tests as f64 / (c.n - 1) as f64; // multiset count / (N-1)
    cell.ratio = cell.p_hat / cell.model_pr; cell.ratio_lo = cell.ratio; cell.ratio_hi = cell.ratio;
    cell.element_tests = tests;
    cell.secs = t0.elapsed().as_secs_f64();
    cell
}

/// Sampling. m = 2: x(R -+ Q) in base x-set. m = 3: x(R -+ Q) in T2.
pub fn sample(c: &Curve, base: &Base, t2: Option<&HMap>, hdr: YieldCell, seed: u64) -> YieldCell {
    let t0 = Instant::now();
    let mut cell = hdr;
    cell.method = "sample".into();
    let fb = base.size() as u64;
    let (mut targets, mut hits, mut repsum, mut tests) = (0u64, 0u64, 0u64, 0u64);
    let batch = 256u64;
    let set = t2.unwrap_or(&base.xset);
    while (hits < 400 && tests < (1u64 << 30) && targets < 1_000_000) || targets < 200 {
        let res: Vec<(bool, u64)> = (targets..targets + batch).into_par_iter()
            .map_init(|| (Vec::new(), Vec::new(), Vec::new()), |(s, om, op), t| {
                let mut rng = Rng(derive(seed, t));
                let r = c.random_point(&mut rng);
                batch_diff_sum(c, r, &base.reps, s, om, op);
                let mut h = 0u64;
                for i in 0..om.len() {
                    if om[i] != u64::MAX { if set.contains(om[i]) { h += 1; } if set.contains(op[i]) { h += 1; } }
                }
                (h > 0, h)
            }).collect();
        for (d, h) in res { if d { hits += 1; } repsum += h; }
        targets += batch; tests += batch * fb;
    }
    cell.targets = targets; cell.hits = hits; cell.element_tests = tests;
    cell.mean_reps = repsum as f64 / targets as f64 / cell.m as f64;
    cell.secs = t0.elapsed().as_secs_f64();
    finish(cell, 2.576)
}

pub fn build_t2(c: &Curve, reps: &[Pt]) -> HMap {
    let parts = pair_xs(c, reps);
    let n: usize = parts.iter().map(|v| v.len()).sum();
    let mut h = HMap::with_capacity(n);
    for v in parts { for x in v { h.insert(x, 0); } }
    h
}

pub fn header(b: u32, c: &Curve, ci: usize, k: u32, arm: &str, m: u32, base: &Base) -> YieldCell {
    let fbp = base.size();
    let lambda = (fbp as f64).powi(m as i32) / (fact(m) * c.n as f64);
    YieldCell { b, p: c.f.p, n: c.n, curve: ci, k, arm: arm.into(), m, method: String::new(), f: base.reps.len(),
        fb_points: fbp, lambda, model_pr: 1.0 - (-lambda).exp(), targets: 0, hits: 0, mean_reps: 0.0,
        p_hat: 0.0, p_lo: 0.0, p_hi: 0.0, ratio: 0.0, ratio_lo: 0.0, ratio_hi: 0.0, element_tests: 0, secs: 0.0 }
}

pub fn coset_t(f: &Fp, k: u32) -> u64 {
    let mut t = 3u64;
    loop { let tm = f.to_m(t); if f.pow(tm, 1u64 << k) != f.one { return tm; } t += 1; }
}

pub fn run_all(bs: &[u32], out: &mut dyn FnMut(&YieldCell)) {
    for &b in bs {
        let (_, kk) = find_prime(b);
        let curves = curves_for(b, 2);
        for m in [2u32, 3] {
            for k in 3..=kk {
                let nom = (1u64 << k) as f64;
                let lam = nom.powi(m as i32) / (fact(m) * curves[0].n as f64);
                if lam < 2f64.powi(-10) || lam > 2.0 { continue; }
                if m == 3 && 2.0 * (nom / 2.0) * (nom / 2.0) > (1u64 << 24) as f64 { continue; }
                for (ci, c) in curves.iter().enumerate() {
                    let mut arms: Vec<Base> = vec![Base::pkm(c, k, c.f.one)];
                    if ci == 0 {
                        arms.push(Base::pkm(c, k, coset_t(&c.f, k)));
                        let f = arms[0].reps.len();
                        let mut rng = Rng(derive(20261120, (b as u64) << 32 | (k as u64) << 8 | m as u64));
                        arms.push(Base::random(c, f, &mut rng));
                        let g = c.random_point(&mut rng);
                        if m == 2 { arms.push(Base::ap(c, f, g)); }
                    }
                    for base in &arms {
                        let hdr = header(b, c, ci, k, &base.label, m, base);
                        let seed = derive(0x5EED, (b as u64) << 40 | (ci as u64) << 32 | (k as u64) << 16 | (m as u64) << 8 | base.label.len() as u64);
                        if m == 2 {
                            let f = base.reps.len() as f64;
                            if 2.0 * f * f <= (1u64 << 26) as f64 { out(&census_m2(c, base, hdr.clone())); }
                            out(&sample(c, base, None, hdr, seed));
                        } else {
                            let t2 = build_t2(c, &base.reps);
                            out(&sample(c, base, Some(&t2), hdr, seed));
                        }
                    }
                }
            }
        }
    }
}
