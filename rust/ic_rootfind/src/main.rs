//! ic_rootfind: point decomposition on binary curves without a Groebner basis (see solver.rs).
//!
//!   ic_rootfind selftest
//!   ic_rootfind solve  --n N --k K --a A --b HEX --xr HEX --yr HEX [--m 2|3]     decompositions as JSON
//!   ic_rootfind bench  --n N --k K [--a A --b HEX] [--log2-range L] [--threads T] [--seed S] [--brute-log2 L]
//!
//! `bench` plants one decomposition R = P1 + P2 + P3, scans a range of 2^L abscissas x(P3) that contains a
//! planted one, checks that the planted decomposition is found and that every row re-adds to R, and reports
//! the throughput.  With --brute-log2 it also times the enumeration inner loop (batched additions and
//! membership tests) on 2^L additions.

mod field;
mod solver;

use field::{Fe, Field};
use solver::{Pt, Solver, Stats};
use std::time::Instant;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn fe(&mut self, n: usize) -> Fe {
        let mut v = [self.next(), self.next(), self.next()];
        for i in 0..3 {
            let base = 64 * i;
            if base >= n {
                v[i] = 0;
            } else if base + 64 > n {
                v[i] &= (1u64 << (n - base)) - 1;
            }
        }
        Fe(v)
    }
}

fn standard_basis(k: usize) -> Vec<Fe> {
    (0..k).map(Fe::unit).collect()
}

fn random_point(s: &Solver, rng: &mut Rng) -> Pt {
    loop {
        if let Some(p) = s.lift(rng.fe(s.n)) {
            return if rng.next() & 1 == 1 { s.neg(p) } else { p };
        }
    }
}

fn sum(s: &Solver, pts: &[Pt]) -> Option<Pt> {
    pts.iter().fold(None, |acc, &p| s.add(acc, Some(p)))
}

fn canon(mut rows: Vec<[Pt; 3]>) -> Vec<[Pt; 3]> {
    rows.sort();
    rows.dedup();
    rows
}

fn lrf3_threads(s: &Solver, r: Pt, lo: u64, hi: u64, threads: usize) -> (Vec<[Pt; 3]>, Stats) {
    if threads <= 1 {
        return s.lrf3_range(r, lo, hi);
    }
    let step = ((hi - lo) / (8 * threads as u64)).max(1);
    let next = std::sync::atomic::AtomicU64::new(lo);
    let parts: Vec<(Vec<[Pt; 3]>, Stats)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                sc.spawn(|| {
                    let (mut rows, mut st) = (Vec::new(), Stats::default());
                    loop {
                        let a = next.fetch_add(step, std::sync::atomic::Ordering::Relaxed);
                        if a >= hi {
                            break;
                        }
                        let (r2, s2) = s.lrf3_range(r, a, (a + step).min(hi));
                        rows.extend(r2);
                        st.merge(&s2);
                    }
                    (rows, st)
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let (mut rows, mut st) = (Vec::new(), Stats::default());
    for (r2, s2) in parts {
        rows.extend(r2);
        st.merge(&s2);
    }
    (rows, st)
}

fn selftest() {
    let mut rng = Rng(1);
    // field arithmetic against the bit-by-bit product, one-, two- and three-word paths
    for &n in &[7usize, 17, 31, 61, 63, 64, 65, 89, 127, 128, 129, 131, 163, 191] {
        let f = Field::new(n);
        for _ in 0..200 {
            let (a, b) = (rng.fe(n), rng.fe(n));
            assert_eq!(f.mul(a, b), f.mul_slow(a, b), "mul n={n}");
            assert_eq!(f.sq(a), f.mul_slow(a, a), "sq n={n}");
            if !a.is_zero() {
                assert_eq!(f.mul(a, f.inv(a)), Fe::ONE, "inv n={n}");
            }
        }
        let xs: Vec<Fe> = (0..37).map(|i| if i % 9 == 4 { Fe::ZERO } else { rng.fe(n) }).collect();
        let mut inv = vec![Fe::ZERO; xs.len()];
        f.batch_inv(&xs, &mut inv);
        for (x, i) in xs.iter().zip(&inv) {
            assert!(if x.is_zero() { i.is_zero() } else { f.mul(*x, *i) == Fe::ONE });
        }
    }
    assert_eq!(Field::new(131).low, vec![8, 3, 2, 0], "the ECC2K-130 pentanomial");
    println!("field arithmetic: ok (14 degrees)");

    // random instances: linear solver against enumeration, complete index range
    let (mut instances, mut nonempty, mut total, mut rowform) = (0, 0, 0, 0);
    for &n in &[9usize, 11, 13, 16, 17, 19] {
        for trial in 0..6 {
            for k in 3..=6 {
                let f = Field::new(n);
                let a = if rng.next() & 1 == 1 { Fe::ONE } else { Fe::ZERO };
                let b = if trial % 2 == 0 { Fe::ONE } else { loop { let b = rng.fe(n); if !b.is_zero() { break b } } };
                let basis = if trial % 3 == 0 { standard_basis(k) } else {
                    loop {
                        let c: Vec<Fe> = (0..k).map(|_| rng.fe(n)).collect();
                        if solver::independent(&c) { break c }
                    }
                };
                let mut s = Solver::new(f, a, b, basis);
                if trial % 3 == 0 { assert!(s.has_row_form() || n - s.dim_v2 == 0, "row form expected for the standard basis"); }
                let fb = s.factor_base();
                for p in &fb { assert!(s.on_curve(*p) && s.in_v(p.0)); }
                let r = random_point(&s, &mut rng);
                let (lrf, _) = s.lrf3_range(r, 0, 1 << k);
                let (bru, _) = s.brute3_range(r, &fb, 0, 1 << k);
                let (lrf, bru) = (canon(lrf), canon(bru));
                assert_eq!(lrf, bru, "n={n} k={k} trial={trial}");
                for d in &lrf { assert_eq!(sum(&s, d), Some(r)); }
                assert_eq!(canon(lrf3_threads(&s, r, 0, 1 << k, 3).0), bru);
                if s.has_row_form() {
                    rowform += 1;
                    s.use_rows = false;       // the generic column form on the same instance
                    assert_eq!(canon(s.lrf3_range(r, 0, 1 << k).0), bru);
                    s.use_rows = true;
                }
                // two summands: T = sum of two factor-base points always has a decomposition
                let t = s.add(Some(fb[(rng.next() % fb.len() as u64) as usize]), Some(fb[(rng.next() % fb.len() as u64) as usize]));
                if let Some(t) = t {
                    let (d2, _) = s.decompose2(t);
                    assert!(!d2.is_empty());
                    for d in &d2 { assert_eq!(sum(&s, d), Some(t)); assert!(s.in_v(d[0].0) && s.in_v(d[1].0)); }
                    let mut want: Vec<[Pt; 2]> = Vec::new();
                    for &p in &fb {
                        if let Some(q) = s.add(Some(t), Some(s.neg(p))) {
                            if s.in_v(q.0) { let mut d = [p, q]; d.sort(); want.push(d); }
                        }
                    }
                    want.sort();
                    want.dedup();
                    assert_eq!(d2, want, "m=2 n={n} k={k}");
                }
                instances += 1;
                nonempty += (!lrf.is_empty()) as u32;
                total += lrf.len();
            }
        }
    }
    println!("m = 3 and m = 2 against enumeration: ok ({instances} instances, {nonempty} with decompositions, {total} decompositions; row and column forms compared on {rowform})");

    // planted decompositions in multi-word fields, where random targets have none
    for &(n, k) in &[(67usize, 8usize), (71, 10), (97, 18), (131, 44), (163, 52)] {
        let mut s = Solver::new(Field::new(n), Fe::ZERO, Fe::ONE, standard_basis(k));
        assert!(s.has_row_form());
        for it in 0..4 {
            s.use_rows = it % 2 == 0;
            let (r, planted, g) = plant(&s, &mut rng);
            let lo = g & !0x3ff;
            let hi = (lo + 1024).min(1 << k);
            let (rows, _) = s.lrf3_range(r, lo, hi);
            assert!(canon(rows.clone()).contains(&planted), "planted n={n} k={k}");
            for d in &rows { assert_eq!(sum(&s, d), Some(r)); for p in d { assert!(s.on_curve(*p) && s.in_v(p.0)); } }
            if k <= 10 {
                let fb = s.factor_base();
                assert_eq!(canon(rows), canon(s.brute3_range(r, &fb, lo, hi).0));
            }
        }
    }
    println!("planted decompositions, n = 67, 71, 97, 131, 163: ok");
    println!("SELFTEST PASSED");
}

/// R = P1 + P2 + P3 with random factor-base points; returns R, the sorted triple and the Gray index of x(P3).
fn plant(s: &Solver, rng: &mut Rng) -> (Pt, [Pt; 3], u64) {
    loop {
        let mut pts = Vec::new();
        let mut g3 = 0;
        while pts.len() < 3 {
            let g = rng.next() & ((1u64 << s.k) - 1);
            if let Some(p) = s.lift(s.gray(g)) {
                pts.push(if rng.next() & 1 == 1 { s.neg(p) } else { p });
                g3 = g;
            }
        }
        if let Some(r) = sum(s, &pts) {
            let mut d = [pts[0], pts[1], pts[2]];
            d.sort();
            return (r, d, g3);
        }
    }
}

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1).cloned())
}

fn pt_json(p: Pt) -> String {
    format!("[\"{}\",\"{}\"]", p.0.hex(), p.1.hex())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(String::as_str).unwrap_or("");
    let num = |name: &str, default: u64| arg(&args, name).map(|v| v.parse().unwrap()).unwrap_or(default);
    match cmd {
        "selftest" => selftest(),
        "mulbench" => {
            let n = num("--n", 131) as usize;
            let f = Field::new(n);
            let mut rng = Rng(3);
            let (mut a, b) = (rng.fe(n), rng.fe(n));
            let t0 = Instant::now();
            for _ in 0..20_000_000u64 { a = f.mul(a, b); }
            println!("n={n} mul {:.2} ns  ({})", t0.elapsed().as_secs_f64() / 2e7 * 1e9, a.hex());
            let t0 = Instant::now();
            for _ in 0..20_000_000u64 { a = f.sq(a); }
            println!("n={n} sq  {:.2} ns  ({})", t0.elapsed().as_secs_f64() / 2e7 * 1e9, a.hex());
        }
        "solve" => {
            let (n, k, m) = (num("--n", 0) as usize, num("--k", 0) as usize, num("--m", 3));
            let a = if num("--a", 0) == 1 { Fe::ONE } else { Fe::ZERO };
            let b = Fe::from_hex(&arg(&args, "--b").unwrap());
            let r = (Fe::from_hex(&arg(&args, "--xr").unwrap()), Fe::from_hex(&arg(&args, "--yr").unwrap()));
            let s = Solver::new(Field::new(n), a, b, standard_basis(k));
            assert!(s.on_curve(r));
            let rows: Vec<String> = if m == 2 {
                s.decompose2(r).0.iter().map(|d| format!("[{},{}]", pt_json(d[0]), pt_json(d[1]))).collect()
            } else {
                canon(s.lrf3_range(r, 0, 1 << k).0).iter()
                    .map(|d| format!("[{},{},{}]", pt_json(d[0]), pt_json(d[1]), pt_json(d[2]))).collect()
            };
            println!("[{}]", rows.join(","));
        }
        "bench" => {
            let (n, k) = (num("--n", 131) as usize, num("--k", 44) as usize);
            let a = if num("--a", 0) == 1 { Fe::ONE } else { Fe::ZERO };
            let b = arg(&args, "--b").map(|h| Fe::from_hex(&h)).unwrap_or(Fe::ONE);
            let (l, threads, seed) = (num("--log2-range", 20), num("--threads", 1) as usize, num("--seed", 1));
            let t0 = Instant::now();
            let mut s = Solver::new(Field::new(n), a, b, standard_basis(k));
            let t_setup = t0.elapsed().as_secs_f64();
            if args.iter().any(|x| x == "--columns") { s.use_rows = false; }
            let form = if s.has_row_form() && s.use_rows { "rows" } else { "columns" };
            let mut rng = Rng(seed);
            let (r, planted, g) = plant(&s, &mut rng);
            let size = 1u64 << l.min(k as u64);
            let lo = g & !(size - 1);
            let t0 = Instant::now();
            let (rows, st) = lrf3_threads(&s, r, lo, lo + size, threads);
            let t = t0.elapsed().as_secs_f64();
            let rows = canon(rows);
            for d in &rows {
                assert_eq!(sum(&s, d), Some(r), "a row does not sum to R");
                for p in d { assert!(s.on_curve(*p) && s.in_v(p.0)); }
            }
            let found = rows.contains(&planted);
            assert!(found, "planted decomposition not found");
            let mut extra = String::new();
            if let Some(bl) = arg(&args, "--brute-log2") {
                // the enumeration inner loop: 2^bl additions S - P2 with membership tests, P2 over a listed
                // slice of the factor base (2^10 abscissas), repeated over P3
                let bl: u32 = bl.parse().unwrap();
                let mut fb = Vec::new();
                for gg in 0..1024u64 {
                    if let Some(p) = s.lift(s.gray(lo + gg)) { fb.push(p); fb.push(s.neg(p)); }
                }
                let t0 = Instant::now();
                let (mut adds, mut gg, mut hits) = (0u64, lo, 0usize);
                while adds < (1u64 << bl) {
                    let (rws, a2) = s.brute3_range(r, &fb, gg, gg + 16);
                    adds += a2;
                    hits += rws.len();
                    gg += 16;
                }
                let tb = t0.elapsed().as_secs_f64();
                extra = format!(",\"brute_adds\":{adds},\"t_brute\":{tb:.6},\"brute_ns_per_add\":{:.2},\"brute_hits\":{hits}", tb / adds as f64 * 1e9);
            }
            println!(
                "{{\"form\":\"{form}\",\"n\":{n},\"k\":{k},\"a\":{},\"b\":\"{}\",\"modulus_low_terms\":{:?},\"dimV2\":{},\"unknowns\":{},\"equations\":{},\"seed\":{seed},\"threads\":{threads},\"log2_range\":{},\"range_lo\":{lo},\"xR\":\"{}\",\"yR\":\"{}\",\"t_setup\":{t_setup:.6},\"t_lrf\":{t:.6},\"abscissas\":{},\"targets\":{},\"lin_solvable\":{},\"lin_candidates\":{},\"quad_split\":{},\"max_kernel_dim\":{},\"degenerate\":{},\"rows\":{},\"planted_found\":{found},\"ns_per_abscissa\":{:.2},\"ns_per_target\":{:.2}{extra}}}",
                if a == Fe::ONE { 1 } else { 0 }, b.hex(), s.f.low, s.dim_v2, k + 1, s.w, size.trailing_zeros(), r.0.hex(), r.1.hex(),
                st.abscissas, st.targets, st.lin_solvable, st.lin_candidates, st.quad_split, st.max_kernel_dim, st.degenerate, rows.len(),
                t / st.abscissas as f64 * 1e9, t / st.targets as f64 * 1e9
            );
        }
        _ => {
            eprintln!("usage: ic_rootfind selftest | solve ... | bench ...   (see src/main.rs)");
            std::process::exit(2);
        }
    }
}
