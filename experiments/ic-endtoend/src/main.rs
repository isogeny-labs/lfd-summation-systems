mod curve; mod e2e; mod field; mod hash; mod la; mod report; mod rho; mod yields;
use curve::*;
use field::*;
use rayon::prelude::*;
use std::io::Write;

fn bits(args: &[String], def: &[u32]) -> Vec<u32> {
    args.iter().position(|a| a == "--bits").map(|i| args[i + 1].split(',').map(|s| s.parse().unwrap()).collect()).unwrap_or(def.to_vec())
}
fn opt(args: &[String], name: &str, def: u64) -> u64 { args.iter().position(|a| a == name).map(|i| args[i + 1].parse().unwrap()).unwrap_or(def) }
fn appender(path: &str) -> std::fs::File { std::fs::OpenOptions::new().create(true).append(true).open(path).unwrap() }

pub fn generator(c: &Curve, b: u32) -> Pt { c.random_point(&mut Rng(derive(0xC0FFEE, b as u64))) }
pub fn target(c: &Curve, b: u32, i: u64) -> u64 { 1 + derive(0x7A6, (b as u64) << 8 | i) % (c.n - 1) }

/// e2e cells: k*-1..k*+1 around the predicted optimum, under the table caps
pub fn e2e_ks(n: u64, kk: u32, m: u32) -> Vec<u32> {
    let cap = if m == 2 { (1u64 << 25) as f64 } else { (1u64 << 24) as f64 };
    let cost = |k: u32| { let f = (1u64 << (k - 1)) as f64; if m == 2 { 2.0 * f * f + n as f64 / (4.0 * f) } else { 2.0 * f * f + 1.4 * n as f64 / f } };
    let ok = |k: u32| { let f = (1u64 << (k - 1)) as f64; k >= 3 && k <= kk && 2.0 * f * f <= cap };
    let ks: Vec<u32> = (3..=kk).filter(|&k| ok(k)).collect();
    let best = *ks.iter().min_by(|&&x, &&y| cost(x).partial_cmp(&cost(y)).unwrap()).unwrap();
    [best - 1, best, best + 1].into_iter().filter(|&k| ok(k)).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = std::env::var("ICEE_OUT").unwrap_or("results".into());
    std::fs::create_dir_all(&out_dir).unwrap();
    match args.get(1).map(|s| s.as_str()) {
        Some("curves") => {
            for b in bits(&args, &[16, 20, 24, 28, 32, 36, 40]) {
                let (p, kk) = find_prime(b);
                for (i, c) in curves_for(b, 2).iter().enumerate() {
                    println!("{{\"b\":{b},\"p\":{p},\"K\":{kk},\"curve\":{i},\"a\":{},\"bcoef\":{},\"N\":{}}}", c.a_raw, c.b_raw, c.n);
                }
            }
        }
        Some("yield") => {
            let mut fh = appender(&format!("{out_dir}/yield.jsonl"));
            yields::run_all(&bits(&args, &[16, 20, 24, 28, 32, 36, 40]), &mut |cell| {
                let s = serde_json::to_string(cell).unwrap();
                eprintln!("b={} m={} k={} c={} {} {} f={} lam={:.4} hits={}/{} r={:.3} [{:.3},{:.3}] {:.1}s", cell.b, cell.m, cell.k, cell.curve, cell.arm, cell.method, cell.f, cell.lambda, cell.hits, cell.targets, cell.ratio, cell.ratio_lo, cell.ratio_hi, cell.secs);
                writeln!(fh, "{s}").unwrap();
            });
        }
        Some("e2e") => {
            let ntargets = opt(&args, "--targets", 3);
            let ms: Vec<u32> = args.iter().position(|a| a == "--m").map(|i| args[i + 1].split(',').map(|s| s.parse().unwrap()).collect()).unwrap_or(vec![2, 3]);
            let mut fh = appender(&format!("{out_dir}/e2e.jsonl"));
            for b in bits(&args, &[20, 24, 28, 32, 36, 40]) {
                let (_, kk) = find_prime(b);
                let c = &curves_for(b, 1)[0];
                let g = generator(c, b);
                for &m in &ms {
                    for k in e2e_ks(c.n, kk, m) {
                        let base = Base::pkm(c, k, c.f.one);
                        let table = e2e::build_table(c, &base.reps);
                        let res: Vec<e2e::E2e> = (0..ntargets).into_par_iter().map(|ti| e2e::run(c, b, m, k, &base, &table, g, ti, target(c, b, ti))).collect();
                        for r in res {
                            eprintln!("b={} m={} k={} f={} t={} ok={} att={} rel={} GA={} native={} LA={} core={} {:.1}s+{:.1}s", b, m, k, r.f, r.target_idx, r.verified, r.attempts, r.relations, r.group_adds(), r.native(), r.la.muladds, r.la.n_core, r.secs_search, r.secs_la);
                            writeln!(fh, "{}", serde_json::to_string(&r).unwrap()).unwrap();
                        }
                    }
                }
            }
        }
        Some("rho") => {
            let trials = opt(&args, "--trials", 400);
            let mut fh = appender(&format!("{out_dir}/rho.jsonl"));
            for b in bits(&args, &[20, 24, 28, 32, 36, 40]) {
                let c = &curves_for(b, 1)[0];
                let g = generator(c, b);
                for neg in [true, false] {
                    let t0 = std::time::Instant::now();
                    let res: Vec<rho::RhoTrial> = (0..trials).into_par_iter().map(|i| {
                        let x = target(c, b, i);
                        let q = c.mul(g, x).0;
                        rho::rho(c, g, q, x, neg, derive(0x2B0, (b as u64) << 32 | (neg as u64) << 31 | i), b)
                    }).collect();
                    let mean = res.iter().map(|t| t.group_adds() as f64).sum::<f64>() / trials as f64;
                    let ok = res.iter().filter(|t| t.ok).count();
                    let exp = if neg { (std::f64::consts::PI * c.n as f64 / 4.0).sqrt() } else { (std::f64::consts::PI * c.n as f64 / 2.0).sqrt() };
                    eprintln!("b={b} neg={neg} ok={ok}/{trials} mean GA={mean:.0} ratio={:.3} {:.1}s", mean / exp, t0.elapsed().as_secs_f64());
                    for r in res { writeln!(fh, "{}", serde_json::to_string(&r).unwrap()).unwrap(); }
                }
            }
        }
        Some("selftest") => {
            for c in curves_for(16, 2) {
                let f = &c.f;
                let mut cnt = 1u64;
                for x in 0..f.p { let r = c.rhs(f.to_m(x)); cnt += if r == 0 { 1 } else if f.is_square(r) { 2 } else { 0 }; }
                println!("naive #E = {cnt}, bsgs N = {}", c.n);
                assert_eq!(cnt, c.n);
            }
            for b in [24u32, 40] {
                let c = &curves_for(b, 1)[0];
                let mut rng = Rng(7);
                let reps: Vec<Pt> = (0..50).map(|_| c.random_point(&mut rng)).collect();
                let r = c.random_point(&mut rng);
                let (mut s, mut om, mut op) = (vec![], vec![], vec![]);
                batch_diff_sum(c, r, &reps, &mut s, &mut om, &mut op);
                for i in 0..50 { assert_eq!(om[i], c.sub(r, reps[i]).x); assert_eq!(op[i], c.add(r, reps[i]).x); }
                let g = generator(c, b);
                assert!(c.mul(g, c.n).0.inf);
                println!("b={b}: batch ok, N*G = O");
            }
            println!("selftest ok");
        }
        Some("report") => report::run(&out_dir),
        _ => eprintln!("usage: ic-endtoend curves|yield|e2e|rho|report [--bits 20,24] [--trials n] [--targets n] [--m 2,3]"),
    }
}
