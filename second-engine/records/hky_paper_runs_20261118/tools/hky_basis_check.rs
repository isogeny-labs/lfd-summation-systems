//! Round-trip checker for `hky.basis/v1` (BASIS_FORMAT.md), independent of the
//! engine's code: its own parser, its own chain product (exponent addition mod
//! 2^k, not the engine's square rules), its own reduction.
//!
//! `hky_basis_check --basis PATH [--sidecar PATH] [--closedness]`
//!
//! Checks: magic/version/encoding, trailing SHA-256, the sidecar's file and
//! content SHA-256, column list (count = Σ C(2k, d), distinct, degree ≤ D, the
//! declared order rule), rows (pivots strictly increasing, entries strictly
//! ascending after the pivot, coefficients in [1, p)), profile and
//! contains_one recomputed from pivots. With `--closedness` (needs the
//! sidecar's p, a, b, target abscissa): the generator NF(S_3(x1, x2, X)) lies
//! in the span, and u·h lies in the span for every basis row h of degree < D
//! and every chain variable u, i.e. the span W contains the generator and is
//! closed, so V_D ⊆ W. Prints one JSON line; exits 1 on any failure.

use rayon::prelude::*;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::HashMap;

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

struct Cur<'a> {
    b: &'a [u8],
    i: usize,
}

impl Cur<'_> {
    fn take(&mut self, n: usize) -> &[u8] {
        assert!(self.i + n <= self.b.len(), "truncated basis file");
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        s
    }
    fn u32(&mut self) -> u32 {
        u32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn u64(&mut self) -> u64 {
        u64::from_le_bytes(self.take(8).try_into().unwrap())
    }
}

fn binom(n: u64, k: u64) -> u64 {
    (0..k).fold(1u64, |acc, i| acc * (n - i) / (i + 1))
}

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    a * b % p
}

struct Basis {
    p: u64,
    k: u32,
    d: u32,
    cols: Vec<u64>,
    /// rows: pivot, then (col, coef) after it
    piv: Vec<u32>,
    rows: Vec<Vec<(u32, u32)>>,
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let get = |f: &str| argv.iter().position(|a| a == f).map(|i| argv[i + 1].clone());
    let path = get("--basis").expect("--basis");
    let side_path = get("--sidecar").unwrap_or(format!("{path}.json"));
    let closedness = argv.iter().any(|a| a == "--closedness");
    let bytes = std::fs::read(&path).expect("read basis");
    let mut fails: Vec<String> = Vec::new();
    let mut check = |ok: bool, what: &str| {
        if !ok {
            fails.push(what.to_string());
        }
    };
    let file_sha = hex(&Sha256::digest(&bytes));
    let n = bytes.len();
    assert!(n > 32 + 64);
    let content_sha = hex(&Sha256::digest(&bytes[..n - 32]));
    check(content_sha == hex(&bytes[n - 32..]), "trailing sha256");
    let side: Option<Value> = std::fs::read_to_string(&side_path)
        .ok()
        .map(|s| serde_json::from_str(&s).expect("sidecar JSON"));
    if let Some(s) = &side {
        check(s["file_sha256"] == json!(file_sha), "sidecar file_sha256");
        check(s["content_sha256"] == json!(content_sha), "sidecar content_sha256");
        check(s["bytes"] == json!(n as u64), "sidecar bytes");
    }
    let mut c = Cur { b: &bytes[..n - 32], i: 0 };
    check(c.take(8) == b"HKYBAS01", "magic");
    check(c.u32() == 1, "version");
    let m = c.u32();
    check(m == 2, "m = 2");
    let p = c.u64();
    let k = c.u32();
    let d = c.u32();
    let n_vars = c.u32();
    check(n_vars == 2 * k, "n_vars");
    check(c.u32() == 1, "mask_encoding");
    let order = c.u32();
    let flags = c.u32();
    let n_cols = c.u64() as usize;
    let rank = c.u64() as usize;
    let profile: Vec<u64> = (0..=d).map(|_| c.u64()).collect();
    let expect_cols: u64 = (0..=d as u64).map(|j| binom(2 * k as u64, j)).sum();
    check(n_cols as u64 == expect_cols, "n_cols = sum C(2k, d)");
    let cols: Vec<u64> = (0..n_cols).map(|_| c.u64()).collect();
    // column rule
    let kk = k as u64;
    let low = (1u64 << kk) - 1;
    let sigma = |m: u64| (m >> kk) | ((m & low) << kk);
    let ord_key = |m: u64| (m.count_ones(), m);
    let mut col_ok = cols.iter().all(|&m| m.count_ones() <= d && m >> (2 * kk) == 0);
    if order == 0 {
        col_ok &= cols.windows(2).all(|w| ord_key(w[0]) > ord_key(w[1]));
    } else if order == 1 {
        let mut i = 0;
        let mut last_rep: Option<u64> = None;
        while i < n_cols {
            let r = cols[i];
            let s = sigma(r);
            col_ok &= r >= s;
            if let Some(lr) = last_rep {
                col_ok &= ord_key(lr) > ord_key(r);
            }
            last_rep = Some(r);
            if s != r {
                col_ok &= i + 1 < n_cols && cols[i + 1] == s;
                i += 2;
            } else {
                i += 1;
            }
        }
    } else {
        col_ok = false;
    }
    let index: HashMap<u64, u32> = cols.iter().enumerate().map(|(i, &m)| (m, i as u32)).collect();
    col_ok &= index.len() == n_cols;
    check(col_ok, "column list and order rule");
    let mut piv = Vec::with_capacity(rank);
    let mut rows = Vec::with_capacity(rank);
    let mut rows_ok = true;
    for _ in 0..rank {
        let pv = c.u32();
        let t = c.u32() as usize;
        let mut r = Vec::with_capacity(t);
        let mut last = pv;
        for _ in 0..t {
            let cc = c.u32();
            let x = c.u32();
            rows_ok &= cc > last && (cc as usize) < n_cols && x != 0 && u64::from(x) < p;
            last = cc;
            r.push((cc, x));
        }
        if let Some(&prev) = piv.last() {
            rows_ok &= pv > prev;
        }
        rows_ok &= (pv as usize) < n_cols;
        piv.push(pv);
        rows.push(r);
    }
    check(rows_ok, "rows: pivots increasing, entries ascending after pivot, coefficients in [1,p)");
    check(c.i == c.b.len(), "no trailing bytes before the digest");
    let mut prof = vec![0u64; d as usize + 1];
    for &pv in &piv {
        prof[cols[pv as usize].count_ones() as usize] += 1;
    }
    for i in 1..prof.len() {
        prof[i] += prof[i - 1];
    }
    check(prof == profile, "profile recomputed from pivots");
    let one = piv.last().map(|&pv| cols[pv as usize] == 0).unwrap_or(false);
    check(one == (flags & 1 == 1), "contains_one flag");
    if let Some(s) = &side {
        check(s["profile"] == json!(prof), "sidecar profile");
        check(s["rank"] == json!(rank as u64), "sidecar rank");
        check(s["p"] == json!(p) && s["k"] == json!(k) && s["D"] == json!(d), "sidecar p,k,D");
    }
    // negative control of the closedness detector: `--drop-row I` removes row I (by pivot order)
    let (mut piv, mut rows) = (piv, rows);
    if let Some(i) = get("--drop-row") {
        let i: usize = i.parse().unwrap();
        piv.remove(i);
        rows.remove(i);
        eprintln!("negative control: dropped row {i}");
    }
    let b = Basis { p, k, d, cols, piv, rows };
    let mut closed = Value::Null;
    if closedness {
        let s = side.as_ref().expect("--closedness needs the sidecar");
        let a = s["instance"]["a"].as_u64().expect("sidecar instance.a");
        let bb = s["instance"]["b"].as_u64().expect("sidecar instance.b");
        let x = s["instance"]["target"][0].as_u64().expect("sidecar instance.target[0]");
        let (ok_gen, ok_prod, products) = check_closed(&b, &index, a, bb, x);
        check(ok_gen, "generator in span");
        check(ok_prod, "closed under chain-variable products");
        closed = json!({"generator_in_span": ok_gen, "closed": ok_prod, "products_checked": products,
            "mode": if argv.iter().any(|s| s == "--exact") { "exact" } else { "annihilator" }});
    }
    let out = json!({
        "schema": "hky.basis.check/v1", "basis": path, "file_sha256": file_sha,
        "content_sha256": content_sha, "p": p, "k": k, "D": d, "n_cols": n_cols, "rank": rank,
        "profile": prof, "contains_one": one, "column_order": order, "closedness": closed,
        "ok": fails.is_empty(), "failures": fails,
    });
    println!("{out}");
    if !fails.is_empty() {
        std::process::exit(1);
    }
}

/// True if `v` (sparse, any order) reduces to zero against the echelon.
fn in_span(b: &Basis, pivot_row: &[u32], v: &[(u32, u64)], acc: &mut Vec<u64>) -> bool {
    let p = b.p;
    if v.is_empty() {
        return true;
    }
    let mut lo = usize::MAX;
    for &(c, x) in v {
        acc[c as usize] = (acc[c as usize] + x) % p;
        lo = lo.min(c as usize);
    }
    let n = acc.len();
    let mut ok = true;
    for col in lo..n {
        let x = acc[col] % p;
        if x == 0 {
            continue;
        }
        let r = pivot_row[col];
        if r == u32::MAX {
            ok = false;
            break;
        }
        acc[col] = 0;
        let f = p - x;
        for &(cc, y) in &b.rows[r as usize] {
            let i = cc as usize;
            acc[i] = (acc[i] + mulm(f, u64::from(y), p)) % p;
        }
    }
    for a in acc.iter_mut().skip(lo) {
        *a = 0;
    }
    ok
}

/// A random element of the annihilator of the span: free columns random,
/// pivot columns solved from the last pivot backwards
/// (a[π_i] = −Σ_{c > π_i} h_i[c]·a[c]); cost = the basis's entries.
fn annihilator(b: &Basis, seed: u64) -> Vec<u64> {
    let p = b.p;
    let n = b.cols.len();
    let mut s = seed;
    let mut next = || {
        s = s.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = s;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        (z ^ (z >> 31)) % p
    };
    let mut a: Vec<u64> = (0..n).map(|_| next()).collect();
    for i in (0..b.piv.len()).rev() {
        let mut acc = 0u64;
        for &(c, y) in &b.rows[i] {
            acc = (acc + mulm(u64::from(y), a[c as usize], p)) % p;
        }
        a[b.piv[i] as usize] = (p - acc) % p;
    }
    a
}

/// Returns (generator in span [exact], closed, products checked). Closedness:
/// every product u·h (h a basis row of degree < D, u a chain variable) is
/// annihilated by `trials` random annihilators of the span — a product outside
/// the span passes one trial with probability 1/p, so a non-closed span passes
/// with probability ≤ p^-trials. `exact` reduces every product instead (small k).
fn check_closed(b: &Basis, index: &HashMap<u64, u32>, a: u64, bb: u64, x: u64) -> (bool, bool, u64) {
    let argv: Vec<String> = std::env::args().collect();
    let exact = argv.iter().any(|s| s == "--exact");
    let trials: u64 = argv
        .iter()
        .position(|s| s == "--trials")
        .map(|i| argv[i + 1].parse().unwrap())
        .unwrap_or(4);
    let p = b.p;
    let k = b.k as u64;
    let nn = 1u64 << k;
    let low = nn - 1;
    let n_cols = b.cols.len();
    let mut pivot_row = vec![u32::MAX; n_cols];
    for (i, &pv) in b.piv.iter().enumerate() {
        pivot_row[pv as usize] = i as u32;
    }
    let pi = p as i128;
    let (ai, bi, xi) = (a as i128, bb as i128, x as i128);
    let xx = xi * xi % pi;
    let terms: [((u64, u64), i128); 9] = [
        ((2, 0), xx),
        ((0, 2), xx),
        ((1, 1), -2 * xx - 2 * ai),
        ((2, 1), -2 * xi),
        ((1, 2), -2 * xi),
        ((1, 0), -2 * ai * xi - 4 * bi),
        ((0, 1), -2 * ai * xi - 4 * bi),
        ((0, 0), -4 * bi * xi + ai * ai),
        ((2, 2), 1),
    ];
    let mut g: Vec<(u32, u64)> = Vec::new();
    for ((e1, e2), v) in terms {
        let v = v.rem_euclid(pi) as u64;
        if v != 0 {
            let m = (e1 % nn) | ((e2 % nn) << k);
            g.push((index[&m], v));
        }
    }
    let mut acc = vec![0u64; n_cols];
    let gen_ok = in_span(b, &pivot_row, &g, &mut acc);
    let mul = |m: u64, v: usize| -> u64 {
        let (blk, j) = ((v as u64) / k, (v as u64) % k);
        let e = (m >> (blk * k)) & low;
        let e2 = (e + (1u64 << j)) & low;
        (m & !(low << (blk * k))) | (e2 << (blk * k))
    };
    let low_rows: Vec<usize> = (0..b.piv.len())
        .filter(|&i| b.cols[b.piv[i] as usize].count_ones() < b.d)
        .collect();
    let nv = 2 * b.k as usize;
    let seed0 = b.rows.len() as u64 ^ 0x5eed_0000;
    let anns: Vec<Vec<u64>> = if exact {
        Vec::new()
    } else {
        (0..trials).into_par_iter().map(|t| annihilator(b, seed0.wrapping_add(t * 7919))).collect()
    };
    // the annihilators really annihilate (sanity: every basis row)
    for an in &anns {
        let ok = (0..b.piv.len()).into_par_iter().all(|i| {
            let mut s = an[b.piv[i] as usize];
            for &(c, y) in &b.rows[i] {
                s = (s + mulm(u64::from(y), an[c as usize], p)) % p;
            }
            s == 0
        });
        assert!(ok, "annihilator construction");
    }
    let bad = std::sync::atomic::AtomicU64::new(0);
    low_rows.par_iter().for_each_init(
        || vec![0u64; if exact { n_cols } else { 0 }],
        |acc, &i| {
            let mut full: Vec<(u64, u64)> = vec![(b.cols[b.piv[i] as usize], 1)];
            full.extend(b.rows[i].iter().map(|&(c, y)| (b.cols[c as usize], u64::from(y))));
            for v in 0..nv {
                let prod: Vec<(u32, u64)> = full
                    .iter()
                    .map(|&(m, y)| (*index.get(&mul(m, v)).expect("product leaves M_<=D"), y))
                    .collect();
                let ok = if exact {
                    in_span(b, &pivot_row, &prod, acc)
                } else {
                    anns.iter().all(|an| {
                        prod.iter().fold(0u64, |s, &(c, y)| (s + mulm(y, an[c as usize], p)) % p) == 0
                    })
                };
                if !ok {
                    bad.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
            }
        },
    );
    let products = (low_rows.len() * nv) as u64;
    eprintln!(
        "closedness: mode={} products={products} failing={}",
        if exact { "exact".to_string() } else { format!("annihilator x{trials}") },
        bad.load(std::sync::atomic::Ordering::Relaxed)
    );
    (gen_ok, bad.into_inner() == 0, products)
}
