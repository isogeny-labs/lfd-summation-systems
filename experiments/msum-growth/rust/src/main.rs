//! `msum_direct`: closures V_D of the direct normalized m-summand prime chain
//! presentation `{NF(S_{m+1}(x_1..x_m, t))} ∪ G_{m,k}` on the extracted hky3
//! engine (`second-engine/engine`, unchanged). See `../PREREGISTRATION.md`.
//!
//! The generator is computed as in the paper (Section "Three to eight
//! summands"): the values of `S_{m+1}` on `H^m` (H = μ_{2^k} ⊂ F_p^*) from
//! Semaev's resultant recursion with numbers, then the inverse DFT over
//! `(Z/2^k)^m`. Chain variable `i*k + j` is `u_{i,j} = x_i^{2^j}`.
//!
//! Usage:
//!   msum_direct --m 3 --k 6 (--p P --a A --b B --t T [--y Y] | --p P --seed I)
//!               [--D d] [--verify on|off] [--sym on|off] [--ev-ranks D]
//!               [--budget S] [--max-entries N] [--out FILE] [--tag NAME]
//! Without `--D` only the generator, its checks, |Z| and (with `--ev-ranks`)
//! the evaluation ranks are computed.

use std::collections::BTreeMap;
use std::io::Write;
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::f4_fp_tower::SquareRule;
use crypto_lib::cryptanalysis::hky_fp_tower::{count_monomials, HkyLimits, HkyRing, NfPoly};
use crypto_lib::cryptanalysis::hky_fp_tower_fast::swap_first_two_blocks;
use crypto_lib::cryptanalysis::hky_fp_tower_fast2::{closure_fast2, Fast2Options};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

// ── arithmetic mod p ─────────────────────────────────────────────────────
fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn addm(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}
fn subm(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}
fn powm(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, a, p);
        }
        a = mulm(a, a, p);
        e >>= 1;
    }
    r
}
fn inv(a: u64, p: u64) -> u64 {
    assert!(a % p != 0, "inverse of 0");
    powm(a, p - 2, p)
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn sqrt_mod(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if powm(a, (p - 1) / 2, p) != 1 {
        return None;
    }
    // Tonelli–Shanks
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while powm(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, powm(z, q, p), powm(a, q, p), powm(a, (q + 1) / 2, p));
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mulm(tt, tt, p);
            i += 1;
        }
        let b = powm(c, 1 << (m - i - 1), p);
        m = i;
        c = mulm(b, b, p);
        t = mulm(t, c, p);
        r = mulm(r, b, p);
    }
    Some(r)
}

/// Order-2^k generator, the rule of `paper.rs` / `csd/ic_dreg.subgroup_2power`.
fn subgroup_generator(p: u64, k: usize) -> u64 {
    let n = 1u64 << k;
    assert!((p - 1) % n == 0, "2^{k} does not divide p - 1");
    for h in 2..p {
        let g = powm(h, (p - 1) / n, p);
        if powm(g, n / 2, p) != 1 {
            return g;
        }
    }
    unreachable!()
}

// ── summation polynomials by numbers ─────────────────────────────────────
#[derive(Clone, Copy)]
struct Curve {
    p: u64,
    a: u64,
    b: u64,
}

/// S_3(x1, x2, X) as [C, B, A] (coefficients of X^0, X^1, X^2).
fn s3_in_x(c: Curve, x1: u64, x2: u64) -> [u64; 3] {
    let p = c.p;
    let d = subm(x1, x2, p);
    let a2 = mulm(d, d, p);
    let s = addm(x1, x2, p);
    let pr = mulm(x1, x2, p);
    let b1 = addm(mulm(s, addm(pr, c.a, p), p), mulm(2, c.b, p), p);
    let b1 = subm(0, mulm(2, b1, p), p);
    let e = subm(pr, c.a, p);
    let c0 = subm(mulm(e, e, p), mulm(mulm(4, c.b, p), s, p), p);
    [c0, b1, a2]
}

fn det_mod(mut m: Vec<Vec<u64>>, p: u64) -> u64 {
    let n = m.len();
    let mut det = 1u64;
    for col in 0..n {
        let Some(piv) = (col..n).find(|&r| m[r][col] != 0) else {
            return 0;
        };
        if piv != col {
            m.swap(piv, col);
            det = subm(0, det, p);
        }
        det = mulm(det, m[col][col], p);
        let iv = inv(m[col][col], p);
        for r in col + 1..n {
            if m[r][col] != 0 {
                let f = mulm(m[r][col], iv, p);
                for c2 in col..n {
                    let v = mulm(f, m[col][c2], p);
                    m[r][c2] = subm(m[r][c2], v, p);
                }
            }
        }
    }
    det
}

/// Sylvester resultant of f (formal degree len-1) and g (formal degree len-1),
/// coefficients low to high.
fn resultant(f: &[u64], g: &[u64], p: u64) -> u64 {
    let (df, dg) = (f.len() - 1, g.len() - 1);
    let n = df + dg;
    let mut m = vec![vec![0u64; n]; n];
    for r in 0..dg {
        for i in 0..=df {
            m[r][r + i] = f[df - i];
        }
    }
    for r in 0..df {
        for i in 0..=dg {
            m[dg + r][r + i] = g[dg - i];
        }
    }
    det_mod(m, p)
}

/// Coefficients (low to high) of the polynomial of degree ≤ d through
/// (X, vals[X]) for X = 0..=d.
fn interpolate(vals: &[u64], p: u64) -> Vec<u64> {
    let n = vals.len();
    let mut m: Vec<Vec<u64>> = (0..n)
        .map(|x| {
            let mut row: Vec<u64> = (0..n).map(|j| powm(x as u64, j as u64, p)).collect();
            row.push(vals[x]);
            row
        })
        .collect();
    for col in 0..n {
        let piv = (col..n).find(|&r| m[r][col] != 0).unwrap();
        m.swap(piv, col);
        let iv = inv(m[col][col], p);
        for c2 in col..=n {
            m[col][c2] = mulm(m[col][c2], iv, p);
        }
        for r in 0..n {
            if r != col && m[r][col] != 0 {
                let f = m[r][col];
                for c2 in col..=n {
                    let v = mulm(f, m[col][c2], p);
                    m[r][c2] = subm(m[r][c2], v, p);
                }
            }
        }
    }
    (0..n).map(|r| m[r][n]).collect()
}

/// S_n(xs) for n = xs.len() >= 3, by
/// S_n(x1..xn) = Res_X(S_3(x1, x2, X), S_{n-1}(x3..xn, X)).
fn s_val(c: Curve, xs: &[u64]) -> u64 {
    let p = c.p;
    let n = xs.len();
    assert!(n >= 3);
    if n == 3 {
        let [c0, c1, c2] = s3_in_x(c, xs[0], xs[1]);
        let x = xs[2];
        return addm(addm(c0, mulm(c1, x, p), p), mulm(c2, mulm(x, x, p), p), p);
    }
    let f = s3_in_x(c, xs[0], xs[1]);
    let d = 1usize << (n - 3); // degree of S_{n-1} in its last variable
    let mut args: Vec<u64> = xs[2..].to_vec();
    args.push(0);
    let vals: Vec<u64> = (0..=d)
        .map(|x| {
            *args.last_mut().unwrap() = x as u64;
            s_val(c, &args)
        })
        .collect();
    let g = if n == 4 {
        // S_3(x3, x4, X) directly (degree 2 in X)
        s3_in_x(c, xs[2], xs[3]).to_vec()
    } else {
        interpolate(&vals, p)
    };
    let _ = vals;
    resultant(&f, &g, p)
}

// ── curve points (affine), for the generator check ───────────────────────
fn ec_add(c: Curve, p1: Option<(u64, u64)>, p2: Option<(u64, u64)>) -> Option<(u64, u64)> {
    let p = c.p;
    let (Some((x1, y1)), Some((x2, y2))) = (p1, p2) else {
        return p1.or(p2);
    };
    let l = if x1 == x2 {
        if addm(y1, y2, p) == 0 {
            return None;
        }
        mulm(
            addm(mulm(3, mulm(x1, x1, p), p), c.a, p),
            inv(mulm(2, y1, p), p),
            p,
        )
    } else {
        mulm(subm(y2, y1, p), inv(subm(x2, x1, p), p), p)
    };
    let x3 = subm(subm(mulm(l, l, p), x1, p), x2, p);
    let y3 = subm(mulm(l, subm(x1, x3, p), p), y1, p);
    Some((x3, y3))
}
fn random_point(c: Curve, rng: &mut StdRng) -> (u64, u64) {
    loop {
        let x = rng.gen_range(0..c.p);
        let r = addm(addm(mulm(mulm(x, x, c.p), x, c.p), mulm(c.a, x, c.p), c.p), c.b, c.p);
        if let Some(y) = sqrt_mod(r, c.p) {
            if rng.gen_bool(0.5) {
                return (x, y);
            }
            return (x, subm(0, y, c.p));
        }
    }
}

fn seed_rng(s: &str) -> StdRng {
    let h = Sha256::digest(s.as_bytes());
    StdRng::seed_from_u64(u64::from_le_bytes(h[..8].try_into().unwrap()))
}

// ── args ─────────────────────────────────────────────────────────────────
struct Args {
    m: usize,
    k: usize,
    p: u64,
    a: Option<u64>,
    b: Option<u64>,
    t: Option<u64>,
    y: Option<u64>,
    seed: Option<u64>,
    d: Option<u32>,
    verify: bool,
    sym: bool,
    ev_ranks: Option<u32>,
    budget: u64,
    max_entries: u64,
    out: Option<String>,
    tag: String,
}

fn parse() -> Args {
    let argv: Vec<String> = std::env::args().collect();
    let mut a = Args {
        m: 3,
        k: 0,
        p: 0,
        a: None,
        b: None,
        t: None,
        y: None,
        seed: None,
        d: None,
        verify: true,
        sym: true,
        ev_ranks: None,
        budget: 7200,
        max_entries: 1_200_000_000,
        out: None,
        tag: String::new(),
    };
    let mut i = 1;
    while i < argv.len() {
        let v = argv.get(i + 1).cloned().unwrap_or_default();
        match argv[i].as_str() {
            "--m" => a.m = v.parse().unwrap(),
            "--k" => a.k = v.parse().unwrap(),
            "--p" => a.p = v.parse().unwrap(),
            "--a" => a.a = Some(v.parse().unwrap()),
            "--b" => a.b = Some(v.parse().unwrap()),
            "--t" => a.t = Some(v.parse().unwrap()),
            "--y" => a.y = Some(v.parse().unwrap()),
            "--seed" => a.seed = Some(v.parse().unwrap()),
            "--D" => a.d = Some(v.parse().unwrap()),
            "--verify" => a.verify = v == "on",
            "--sym" => a.sym = v == "on",
            "--ev-ranks" => a.ev_ranks = Some(v.parse().unwrap()),
            "--budget" => a.budget = v.parse().unwrap(),
            "--max-entries" => a.max_entries = v.parse().unwrap(),
            "--out" => a.out = Some(v),
            "--tag" => a.tag = v,
            o => panic!("unknown flag {o}"),
        }
        i += 2;
    }
    assert!(a.m >= 3 && a.k >= 1 && a.p > 0);
    a
}

fn status_kb(field: &str) -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines()
        .find(|l| l.starts_with(field))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|x| x.parse().ok())
}

fn sha256_file(path: &std::path::Path) -> String {
    let bytes = std::fs::read(path).unwrap_or_default();
    Sha256::digest(&bytes).iter().map(|x| format!("{x:02x}")).collect()
}

fn main() {
    let args = parse();
    let started = Instant::now();
    let (m, k, p) = (args.m, args.k, args.p);
    assert!(is_prime(p) && p < 1 << 31, "p must be a prime below 2^31");
    let n = 1usize << k;
    assert!((p - 1) % n as u64 == 0 && p > 2 * n as u64, "need 2^k | p-1 and p > 2^(k+1)");
    // ── instance
    let (a, b, t, y, seed_label) = match args.seed {
        Some(i) => {
            let label = format!("msum-{p}-{i}");
            let mut rng = seed_rng(&label);
            loop {
                let a = rng.gen_range(0..p);
                let b = rng.gen_range(0..p);
                let disc = addm(mulm(4, powm(a, 3, p), p), mulm(27, mulm(b, b, p), p), p);
                if disc == 0 {
                    continue;
                }
                let c = Curve { p, a, b };
                let (t, y) = random_point(c, &mut rng);
                let carry = addm(addm(mulm(5, powm(t, 3, p), p), mulm(3, mulm(a, t, p), p), p), mulm(6, b, p), p);
                if carry == 0 {
                    continue;
                }
                break (a, b, t, Some(y), Some(label));
            }
        }
        None => (
            args.a.expect("--a") % p,
            args.b.expect("--b") % p,
            args.t.expect("--t") % p,
            args.y,
            None,
        ),
    };
    let curve = Curve { p, a, b };
    let disc = addm(mulm(4, powm(a, 3, p), p), mulm(27, mulm(b, b, p), p), p);
    assert!(disc != 0, "singular curve");
    if let Some(y) = y {
        let r = addm(addm(mulm(mulm(t, t, p), t, p), mulm(a, t, p), p), b, p);
        assert_eq!(mulm(y, y, p), r, "target not on the curve");
    }
    let carry = addm(addm(mulm(5, powm(t, 3, p), p), mulm(3, mulm(a, t, p), p), p), mulm(6, b, p), p);

    // ── check 1: the resultant recursion vanishes at sums of random points
    let mut rng = seed_rng(&format!("check-{p}-{a}-{b}-{t}"));
    let mut sum_checks = 0;
    for _ in 0..8 {
        let pts: Vec<(u64, u64)> = (0..m).map(|_| random_point(curve, &mut rng)).collect();
        let mut s = None;
        for q in &pts {
            s = ec_add(curve, s, Some(*q));
        }
        let Some((sx, _)) = s else { continue };
        let mut xs: Vec<u64> = pts.iter().map(|q| q.0).collect();
        xs.push(sx);
        assert_eq!(s_val(curve, &xs), 0, "S_{} does not vanish at a sum", m + 1);
        // and is not identically zero
        let mut xs2 = xs.clone();
        xs2[m] = addm(sx, 1, p);
        sum_checks += (s_val(curve, &xs2) != 0) as u32;
    }
    assert!(sum_checks > 0, "S_{} vanishes at non-sums", m + 1);

    // ── values on H^m, |Z|, inverse DFT
    let zeta = subgroup_generator(p, k);
    let hpow: Vec<u64> = (0..n).map(|e| powm(zeta, e as u64, p)).collect();
    let total = n.pow(m as u32);
    let tv = Instant::now();
    use rayon::prelude::*;
    let vals: Vec<u64> = (0..total)
        .into_par_iter()
        .map(|idx| {
            let mut xs = Vec::with_capacity(m + 1);
            let mut r = idx;
            for _ in 0..m {
                xs.push(hpow[r % n]);
                r /= n;
            }
            xs.push(t);
            s_val(curve, &xs)
        })
        .collect();
    let zeros: Vec<usize> = (0..total).filter(|&i| vals[i] == 0).collect();
    let values_s = tv.elapsed().as_secs_f64();
    // inverse DFT along each axis: c(e) = N^{-m} Σ_i f(i) ζ^{-e·i}
    let winv: Vec<u64> = (0..n).map(|e| powm(inv(zeta, p), e as u64, p)).collect();
    let mut coef = vals.clone();
    for axis in 0..m {
        let stride = n.pow(axis as u32);
        let lines: Vec<usize> = (0..total).filter(|&i| (i / stride) % n == 0).collect();
        let newv: Vec<(usize, Vec<u64>)> = lines
            .par_iter()
            .map(|&base| {
                let inl: Vec<u64> = (0..n).map(|j| coef[base + j * stride]).collect();
                let out: Vec<u64> = (0..n)
                    .map(|e| {
                        let mut acc: u128 = 0;
                        for (j, &v) in inl.iter().enumerate() {
                            acc += v as u128 * winv[(e * j) % n] as u128;
                        }
                        (acc % p as u128) as u64
                    })
                    .collect();
                (base, out)
            })
            .collect();
        for (base, out) in newv {
            for (e, v) in out.into_iter().enumerate() {
                coef[base + e * stride] = v;
            }
        }
    }
    let ninv = inv(powm(n as u64, m as u64, p), p);
    for v in coef.iter_mut() {
        *v = mulm(*v, ninv, p);
    }
    // ── the ring and the generator
    let nt = m * k;
    let mut rules = Vec::with_capacity(nt);
    for i in 0..m {
        for j in 0..k {
            let top = j == k - 1;
            rules.push(SquareRule {
                next: if top { None } else { Some((i * k + j + 1) as u8) },
                a: 0,
                b: if top { 0 } else { 1 },
                c: 0,
                d: if top { 1 } else { 0 },
            });
        }
    }
    let ring = HkyRing::new(p, rules, 0);
    let mut raw: Vec<(Vec<u32>, u64)> = Vec::new();
    let mut masks: BTreeMap<u64, u64> = BTreeMap::new();
    for (idx, &c) in coef.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let mut e = vec![0u32; nt];
        let mut r = idx;
        let mut mask = 0u64;
        for i in 0..m {
            let ei = r % n;
            r /= n;
            for j in 0..k {
                if (ei >> j) & 1 == 1 {
                    e[i * k + j] = 1;
                    mask |= 1 << (i * k + j);
                }
            }
        }
        raw.push((e, c));
        masks.insert(mask, c);
    }
    let g: NfPoly = ring.from_raw(&raw);
    assert_eq!(g.len(), masks.len(), "normal form merged distinct masks");
    let deg_g = HkyRing::degree(&g);
    // check 2: the normal form evaluates to S_{m+1} on H^m (64 random grid points + all zeros up to 64)
    let mut grid_checks = 0;
    let mut sample: Vec<usize> = (0..64).map(|_| rng.gen_range(0..total)).collect();
    sample.extend(zeros.iter().take(64));
    for &idx in &sample {
        let mut yv = vec![0u64; nt];
        let mut r = idx;
        for i in 0..m {
            let ei = r % n;
            r /= n;
            for j in 0..k {
                yv[i * k + j] = hpow[(ei << j) % n];
            }
        }
        assert_eq!(ring.eval(&g, &yv), vals[idx], "NF(S) differs from S on the grid");
        grid_checks += 1;
    }
    // check 3 (m = 3, k >= 3): coefficient of x1^3 x2^3 x3^3 against -8(5t^3+3at+6b)
    let c333 = if m == 3 && k >= 3 {
        let mask = 0b11u64 | (0b11 << k) | (0b11 << (2 * k));
        let c = masks.get(&mask).copied().unwrap_or(0);
        let expect = subm(0, mulm(8, carry, p), p);
        Some(json!({"coef": c, "expected_up_to_scalar": expect,
                    "ratio": if expect != 0 { mulm(c, inv(expect, p), p) } else { 0 }}))
    } else {
        None
    };
    // top-degree masks
    let top_masks: Vec<String> = masks
        .keys()
        .filter(|mk| mk.count_ones() == deg_g)
        .take(8)
        .map(|mk| format!("{mk:#x}"))
        .collect();
    // ── evaluation ranks on Z (proper ideals)
    let ev = args.ev_ranks.filter(|_| !zeros.is_empty()).map(|dmax| {
        let pts: Vec<Vec<u64>> = zeros
            .iter()
            .map(|&idx| {
                let mut yv = vec![0u64; nt];
                let mut r = idx;
                for i in 0..m {
                    let ei = r % n;
                    r /= n;
                    for j in 0..k {
                        yv[i * k + j] = hpow[(ei << j) % n];
                    }
                }
                yv
            })
            .collect();
        let nz = pts.len();
        let mut basis: Vec<(usize, Vec<u64>)> = Vec::new(); // (pivot, row normalised)
        let mut ranks = Vec::new();
        for d in 0..=dmax {
            if basis.len() < nz {
                // enumerate masks of popcount d in nt bits
                let mut masks_d: Vec<u64> = Vec::new();
                if d == 0 {
                    masks_d.push(0);
                } else if (d as usize) <= nt {
                    let mut x: u64 = (1u64 << d) - 1;
                    while x < (1u64 << nt) {
                        masks_d.push(x);
                        let c = x & x.wrapping_neg();
                        let r = x + c;
                        x = (((r ^ x) >> 2) / c) | r;
                    }
                }
                for mk in masks_d {
                    if basis.len() == nz {
                        break;
                    }
                    let mut v: Vec<u64> = pts
                        .iter()
                        .map(|yv| {
                            let mut acc = 1u64;
                            let mut bits = mk;
                            while bits != 0 {
                                let j = bits.trailing_zeros() as usize;
                                bits &= bits - 1;
                                acc = mulm(acc, yv[j], p);
                            }
                            acc
                        })
                        .collect();
                    for (pv, row) in &basis {
                        let f = v[*pv];
                        if f != 0 {
                            for (x, &r) in v.iter_mut().zip(row) {
                                *x = subm(*x, mulm(f, r, p), p);
                            }
                        }
                    }
                    if let Some(pv) = v.iter().position(|&x| x != 0) {
                        let iv = inv(v[pv], p);
                        for x in v.iter_mut() {
                            *x = mulm(*x, iv, p);
                        }
                        basis.push((pv, v));
                    }
                }
            }
            ranks.push(basis.len() as u64);
        }
        let ideal_dim: Vec<u64> = (0..=dmax)
            .map(|d| count_monomials(nt, 0, d) as u64 - ranks[d as usize])
            .collect();
        json!({"rank_ev_le_j": ranks, "ideal_dim_le_j": ideal_dim})
    });
    let exe = std::env::current_exe().ok();
    let bin_sha = exe.as_deref().map(sha256_file).unwrap_or_default();
    let mut rec = json!({
        "schema": "msum_direct.closure/v1",
        "tag": args.tag,
        "presentation": "direct normalized chain: {NF(S_{m+1})} + G_{m,k}",
        "m": m, "k": k, "p": p, "a": a, "b": b, "t": t, "y": y, "seed": seed_label,
        "carry_5t3_3at_6b": carry,
        "zeta": zeta,
        "n_chain_vars": nt,
        "generator_terms": g.len(),
        "generator_degree": deg_g,
        "generator_top_masks_first8": top_masks,
        "c333": c333,
        "Z": zeros.len(),
        "unit": zeros.is_empty(),
        "checks": {"sum_vanishing": 8, "sum_nonvanishing_hits": sum_checks, "grid_eval": grid_checks},
        "values_s": values_s,
        "ev": ev,
        "binary_sha256": bin_sha,
        "rayon_threads": rayon::current_num_threads(),
    });
    eprintln!(
        "instance m={m} k={k} p={p} a={a} b={b} t={t}: deg g={deg_g}, terms={}, |Z|={}, c333={:?}",
        g.len(),
        zeros.len(),
        rec["c333"]
    );
    // ── the closure
    if let Some(d) = args.d {
        let limits = HkyLimits {
            max_cols: u128::MAX,
            max_entries: args.max_entries,
            budget: Duration::from_secs(args.budget),
            ..HkyLimits::default()
        };
        let opts = Fast2Options {
            symmetry: if args.sym { swap_first_two_blocks(k, m) } else { None },
            verify_closed: args.verify,
            ..Fast2Options::default()
        };
        let tc = Instant::now();
        let (cl, st) = closure_fast2(&ring, &[g.clone()], d, &limits, &opts);
        let closure_s = tc.elapsed().as_secs_f64();
        let low_profile_ok = cl.censored.is_none();
        rec["D"] = json!(d);
        rec["columns"] = json!(count_monomials(nt, 0, d) as u64);
        rec["closure"] = json!({
            "cols": cl.cols as u64,
            "profile": cl.dims,
            "rank": cl.dims.last().copied(),
            "contains_one": cl.contains_one,
            "censored": cl.censored,
            "closed": low_profile_ok,
            "verify_closed": args.verify && !cl.contains_one && low_profile_ok,
            "verify_products": st.verify_products,
            "verify_nonzero": st.verify_nonzero,
            "verify_s": st.verify_s,
            "rows_reduced": cl.rows_reduced,
            "entries": cl.entries,
            "muladds": cl.muladds,
            "symmetry": st.symmetry,
            "part_cols": st.part_cols,
            "part_rank": st.part_rank,
            "peak_words": st.peak_words,
            "threads": st.threads,
            "simd": st.simd,
            "seconds": closure_s,
        });
        if args.verify && !cl.contains_one && cl.censored.is_none() {
            assert_eq!(st.verify_nonzero, 0, "closedness check failed");
        }
        eprintln!(
            "closure D={d}: cols={} profile={:?} one={} censored={:?} verify_nonzero={} {:.1}s",
            cl.cols, cl.dims, cl.contains_one, cl.censored, st.verify_nonzero, closure_s
        );
    }
    rec["wall_s"] = json!(started.elapsed().as_secs_f64());
    rec["maxrss_mb"] = json!(status_kb("VmHWM:").map(|kb| kb / 1024));
    let line = rec.to_string();
    println!("{line}");
    if let Some(o) = &args.out {
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(o).expect("--out");
        writeln!(f, "{line}").unwrap();
    }
    let _: Value = rec;
}
