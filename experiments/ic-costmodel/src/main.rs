//! Fully charged cost model of PKM prime-field summation-polynomial index
//! calculus against Pollard rho (see PREREGISTRATION.md). Every quantity is
//! carried as log2. Commands: `kd`, `validate`, `sweep` (writes the tables and
//! the figure). Every number this program prints is [pcc:unverified].
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;

const NEG: f64 = f64::NEG_INFINITY;

// ---------- log2 helpers ----------
struct Lf { t: Vec<f64> } // log2 n!
impl Lf {
    fn new(max: usize) -> Lf {
        let mut t = vec![0.0; max + 1];
        for i in 1..=max { t[i] = t[i - 1] + (i as f64).log2(); }
        Lf { t }
    }
    fn fact(&self, n: usize) -> f64 { self.t[n] }
    fn binom(&self, n: usize, k: usize) -> f64 { if k > n { NEG } else { self.t[n] - self.t[k] - self.t[n - k] } }
    fn binom_le(&self, n: usize, d: usize) -> f64 { let mut s = NEG; for i in 0..=d.min(n) { s = lse(s, self.binom(n, i)); } s }
}
fn lse(a: f64, b: f64) -> f64 {
    if a == NEG { return b; } if b == NEG { return a; }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (1.0 + (lo - hi).exp2()).log2()
}
fn l2(x: f64) -> f64 { if x <= 0.0 { NEG } else { x.log2() } }

// ---------- constants (PREREGISTRATION.md) ----------
fn rho_m(m: usize) -> f64 { [0.0, 0.0, 0.53, 0.92, 0.977, 0.992, 0.997, 0.999, 0.9996][m] }
const SETUP_PER_BIT: f64 = 94.0;
const BOOK_PER_ATTEMPT: f64 = 3.0;
const LA2_PER_REL: f64 = 15.0;
const LA_C: f64 = 4.73;
const CORE_FRAC: f64 = 0.65;
const PEEL_PER_CLASS: f64 = 115.0;
const RHO_NATIVE_PER_GA: f64 = 2.25;
const DENSE_TO_END: f64 = 7.5e-3;
const DENSE_TO_ONE: f64 = 4.5e-3;
const PRICES: [f64; 3] = [0.0, 0.1, 1.0];

// ---------- degrees ----------
/// Lambda(w): largest l with 2^(l-1) < 7776 * C(2l, <= w)^2.
fn lambda_w(lf: &Lf, w: usize) -> usize {
    let mut best = 0;
    for l in 1..1000 { if (l as f64 - 1.0) < 7776f64.log2() + 2.0 * lf.binom_le(2 * l, w) { best = l; } }
    best
}
/// m = 2 proved lower bound on failed attempts: max(4 [thm:three], D+1 for every D with k >= k_D [thm:growth]).
fn d_lb2(kd: &[usize], k: usize) -> usize {
    let mut d = 4;
    for (dd, &kdd) in kd.iter().enumerate() { if dd >= 2 && k >= kdd { d = d.max(dd + 1); } }
    let paper = if k >= 10 { (k - 10) / 14 + 1 } else { 0 };
    d.max(paper)
}
fn law2(k: usize) -> usize { (k + 14) / 4 }
fn d_law(m: usize, k: usize) -> usize {
    let d = if m == 2 { law2(k) } else if k + 1 < m { m * k + 1 } else { m * (m - 1) + 1 + law2(k) - 4 };
    d.min(d_ceil(m, k))
}
fn d_ceil(m: usize, k: usize) -> usize { if m == 2 { k + 2 } else { m * k + 1 } }
fn d_lb(m: usize, k: usize, kd: &[usize]) -> usize {
    match m { 2 => d_lb2(kd, k), 3 => if k >= 3 { 7 } else { 2 }, 4 => if k >= 4 { 13 } else { 2 }, _ => 2 }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Sc { Lb, Law, Ceil, LawDense, EnumOnly, LbOnset, Fixed(usize) }
impl Sc {
    fn name(self) -> &'static str { match self { Sc::Lb => "LB", Sc::Law => "LAW", Sc::Ceil => "CEIL", Sc::LawDense => "LAW-dense", Sc::EnumOnly => "enum-only", Sc::LbOnset => "LB+onset", Sc::Fixed(_) => "fixed-D" } }
}

#[derive(Clone, Debug)]
struct Best { total: f64, ga: f64, nat: f64, m: usize, k: usize, route: String, d: usize, mem: f64, la: f64, attempts: f64 }

/// Cost of IC at (n = log2 N, m, k) for one route; returns (log2 GA, log2 native, route label, D, log2 memory bytes, log2 LA, log2 attempts).
/// `f_override`: log2 f (validation with measured f); `only_a`: restrict route (a) to that a.
fn routes(lf: &Lf, kd: &[usize], n: f64, m: usize, k: usize, sc: Sc, lfc: Option<f64>, only_a: Option<usize>) -> Vec<(f64, f64, String, usize, f64, f64, f64)> {
    let lfc = lfc.unwrap_or(k as f64 - 1.0); // log2 f
    let lbig = lfc + 1.0; // log2 |F|
    let llam = m as f64 * lbig - lf.fact(m) - n;
    let lam = if llam > 60.0 { f64::INFINITY } else { llam.exp2() };
    let lpr = if llam < -40.0 { llam } else { l2(-(-lam).exp_m1()) };
    let q = (-lam).exp(); // failure probability
    let lr = rho_m(m).log2() + lfc;
    let la = lr - lpr; // attempts
    let setup = (SETUP_PER_BIT * n).log2();
    let resolve = (m as f64).log2() + lr;
    let la_cost = if m == 2 { LA2_PER_REL.log2() + lr } else {
        lse((LA_C * m as f64 * CORE_FRAC * CORE_FRAC).log2() + 2.0 * lfc, PEEL_PER_CLASS.log2() + lfc)
    };
    let ga0 = lse(lse(setup, la), resolve);
    let nat0 = lse(BOOK_PER_ATTEMPT.log2() + la, la_cost);
    let mut out = Vec::new();
    // route (a)
    let a_lo = (m + 1) / 2;
    for a in a_lo..=m {
        if let Some(oa) = only_a { if a != oa { continue; } }
        let b = m - a;
        let lt = a as f64 * lfc + (a as f64 - 1.0) - lf.fact(a);
        let (t_ga, t_ins) = if a == 1 { (NEG, lfc) } else { (lt, lt) };
        let ls = if b == 0 { NEG } else { b as f64 * lbig - lf.fact(b) };
        let lookups = if b == 0 { 0.0 } else { ls };
        let ga = lse(lse(ga0, t_ga), la + ls);
        let nat = lse(lse(nat0, t_ins), la + lookups);
        out.push((ga, nat, format!("enum a={a}"), 0, 4.0 + lt.max(lfc), la_cost, la));
    }
    if sc == Sc::EnumOnly || only_a.is_some() { return out; }
    // route (b)
    let vars = m * k;
    let (d, per) = match sc {
        Sc::Lb | Sc::LbOnset | Sc::Fixed(_) => {
            let d = match sc {
                Sc::Lb => d_lb(m, k, kd),
                // sensitivity: m >= 5 with the measured onset lfd >= deg g + 1 (empirical, not proved)
                Sc::LbOnset => if m <= 4 { d_lb(m, k, kd) } else if k + 1 < m { m * k + 1 } else { m * (m - 1) + 1 },
                Sc::Fixed(d) => if m == 2 { d.max(d_lb(m, k, kd)) } else { d },
                _ => unreachable!(),
            }.min(vars);
            // preregistered floor: the top layer C(mk, D); the onset sensitivity uses the largest layer <= D (also a floor on |M_<=D|), since C(mk, D) -> 1 as D -> mk
            let c = if sc == Sc::LbOnset { lf.binom(vars, d.min(vars / 2)) } else { lf.binom(vars, d) };
            let per = if m == 2 { lse(l2(q) + c, l2(1.0 - q)) } else { c };
            (d, per)
        }
        Sc::Law => { let d = d_law(m, k); (d, 2.0 * lf.binom_le(vars, d)) }
        Sc::LawDense => { let d = d_law(m, k); (d, DENSE_TO_ONE.log2() + 3.0 * lf.binom_le(vars, d)) }
        Sc::Ceil => { let d = d_ceil(m, k); (d, DENSE_TO_END.log2() + 3.0 * lf.binom_le(vars, d)) }
        Sc::EnumOnly => unreachable!(),
    };
    out.push((ga0, lse(nat0, la + per), "groebner".to_string(), d, NEG, la_cost, la));
    out
}

fn total(ga: f64, nat: f64, pi: f64) -> f64 { if pi == 0.0 { ga } else { lse(ga, pi.log2() + nat) } }

fn best_at(lf: &Lf, kd: &[usize], n: usize, m_set: &[usize], sc: Sc, pi: f64) -> Best {
    let mut best: Option<Best> = None;
    for &m in m_set {
        for k in 2..=(n - 2) {
            for (ga, nat, route, d, mem, la, att) in routes(lf, kd, n as f64, m, k, sc, None, None) {
                let t = total(ga, nat, pi);
                if best.as_ref().map_or(true, |b| t < b.total) {
                    best = Some(Best { total: t, ga, nat, m, k, route, d, mem, la, attempts: att });
                }
            }
        }
    }
    best.unwrap()
}
fn rho_cost(n: f64, pi: f64) -> f64 {
    let ga = 0.5 * (n + (std::f64::consts::PI / 4.0).log2());
    if pi == 0.0 { ga } else { ga + (1.0 + RHO_NATIVE_PER_GA * pi).log2() }
}

// ---------- validation against Study A ----------
fn num(v: &serde_json::Value, path: &[&str]) -> f64 {
    let mut x = v; for p in path { x = &x[*p]; } x.as_f64().unwrap_or(0.0)
}
fn validate(lf: &Lf, kd: &[usize], dir: &str) -> (String, bool) {
    let e2e = fs::read_to_string(format!("{dir}/e2e.jsonl")).expect("e2e.jsonl");
    // (b, m, k) -> (n, sum f, sum GA, sum native, count)
    let mut cells: BTreeMap<(u64, u64, u64), (f64, f64, f64, f64, f64)> = BTreeMap::new();
    for line in e2e.lines().filter(|l| !l.trim().is_empty()) {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        if !v["verified"].as_bool().unwrap_or(false) { continue; }
        let key = (num(&v, &["b"]) as u64, num(&v, &["m"]) as u64, num(&v, &["k"]) as u64);
        let ga: f64 = ["ga_setup", "ga_table", "ga_walk", "ga_search", "ga_resolve"].iter().map(|f| num(&v, &[f])).sum();
        let nat = num(&v, &["hash_inserts"]) + num(&v, &["hash_lookups"]) + num(&v, &["bookkeeping"]) + num(&v, &["la", "muladds"]) + num(&v, &["la_bookkeeping"]);
        let e = cells.entry(key).or_insert((num(&v, &["n"]), 0.0, 0.0, 0.0, 0.0));
        e.1 += num(&v, &["f"]); e.2 += ga; e.3 += nat; e.4 += 1.0;
    }
    let rho = fs::read_to_string(format!("{dir}/rho.jsonl")).expect("rho.jsonl");
    let mut rcells: BTreeMap<u64, (f64, f64, f64, f64)> = BTreeMap::new();
    for line in rho.lines().filter(|l| !l.trim().is_empty()) {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        if !v["neg"].as_bool().unwrap_or(false) { continue; }
        let ga: f64 = ["walk_adds", "start_adds", "replay_adds", "setup_adds", "escape_adds"].iter().map(|f| num(&v, &[f])).sum();
        let nat = num(&v, &["hash_ops"]) + num(&v, &["bookkeeping"]);
        let e = rcells.entry(num(&v, &["b"]) as u64).or_insert((num(&v, &["n"]), 0.0, 0.0, 0.0));
        e.1 += ga; e.2 += nat; e.3 += 1.0;
    }
    let mut s = String::new();
    writeln!(s, "<!-- generated by `ic-costmodel validate`; every number [pcc:unverified] -->\n").unwrap();
    writeln!(s, "## Validation against Study A (end to end, `m = 2, 3`, route (a) with `a = 2`)\n").unwrap();
    writeln!(s, "`pred/meas` = model total / measured mean total (3 verified targets per cell). Model with nominal `f = 2^(k-1)`; in brackets with the measured `f`. Scored: `log2 p >= 28`, every price within `[1/2, 2]`, median `|log2 ratio|` <= 0.5.\n").unwrap();
    writeln!(s, "| log2 p | m | k | f meas | measured @0 | @0.1 | @1 | predicted @0 | @0.1 | @1 | pred/meas @0 | @0.1 | @1 | (meas. f) @0 / @0.1 / @1 | scored |").unwrap();
    writeln!(s, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|").unwrap();
    let mut scored: Vec<f64> = Vec::new();
    let mut all_in = true;
    let mut worst: (f64, String) = (0.0, String::new());
    // best cell per (b, m): measured and predicted ratio to measured rho
    let mut best_meas: BTreeMap<(u64, u64), [(f64, u64); 3]> = BTreeMap::new();
    let mut best_pred: BTreeMap<(u64, u64), [(f64, u64); 3]> = BTreeMap::new();
    for (&(b, m, k), &(nn, fs_, ga, nat, c)) in &cells {
        let n = nn.log2();
        let fmeas = fs_ / c;
        let meas: Vec<f64> = PRICES.iter().map(|&pi| ga / c + pi * nat / c).collect();
        let r = routes(lf, kd, n, m as usize, k as usize, Sc::EnumOnly, None, Some(2));
        let rm = routes(lf, kd, n, m as usize, k as usize, Sc::EnumOnly, Some(fmeas.log2()), Some(2));
        let pred: Vec<f64> = PRICES.iter().map(|&pi| total(r[0].0, r[0].1, pi).exp2()).collect();
        let predm: Vec<f64> = PRICES.iter().map(|&pi| total(rm[0].0, rm[0].1, pi).exp2()).collect();
        let is_scored = b >= 28;
        let mut row_ok = true;
        for i in 0..3 {
            let lr = (pred[i] / meas[i]).log2();
            if is_scored {
                scored.push(lr.abs());
                if lr.abs() > 1.0 { row_ok = false; all_in = false; }
                if lr.abs() > worst.0 { worst = (lr.abs(), format!("log2 p = {b}, m = {m}, k = {k}, price {}", PRICES[i])); }
            }
            let bm = best_meas.entry((b, m)).or_insert([(f64::INFINITY, 0); 3]);
            if meas[i] < bm[i].0 { bm[i] = (meas[i], k); }
            let bp = best_pred.entry((b, m)).or_insert([(f64::INFINITY, 0); 3]);
            if pred[i] < bp[i].0 { bp[i] = (pred[i], k); }
        }
        writeln!(s, "| {b} | {m} | {k} | {:.0} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.2} | {:.2} | {:.2} | {:.2} / {:.2} / {:.2} | {} |",
            fmeas, meas[0], meas[1], meas[2], pred[0], pred[1], pred[2], pred[0] / meas[0], pred[1] / meas[1], pred[2] / meas[2],
            predm[0] / meas[0], predm[1] / meas[1], predm[2] / meas[2],
            if !is_scored { "no".into() } else if row_ok { "in".to_string() } else { "**OUT**".to_string() }).unwrap();
    }
    scored.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if scored.is_empty() { 0.0 } else { scored[scored.len() / 2] };
    // rho
    writeln!(s, "\n### Rho model against the measured rho (negation map, 400 trials per size)\n").unwrap();
    writeln!(s, "| log2 p | measured mean GA | model `sqrt(pi N/4)` | model/meas @0 | @0.1 | @1 | scored |").unwrap();
    writeln!(s, "|---|---|---|---|---|---|---|").unwrap();
    let mut rho_ok = true;
    let mut rho_meas: BTreeMap<u64, [f64; 3]> = BTreeMap::new();
    for (&b, &(nn, ga, nat, c)) in &rcells {
        let n = nn.log2();
        let meas: Vec<f64> = PRICES.iter().map(|&pi| ga / c + pi * nat / c).collect();
        rho_meas.insert(b, [meas[0], meas[1], meas[2]]);
        let pred: Vec<f64> = PRICES.iter().map(|&pi| rho_cost(n, pi).exp2()).collect();
        let sc = b >= 32;
        let mut ok = true;
        for i in 0..3 { let r = pred[i] / meas[i]; if sc && !(0.75..=1.33).contains(&r) { ok = false; rho_ok = false; } }
        writeln!(s, "| {b} | {:.3e} | {:.3e} | {:.2} | {:.2} | {:.2} | {} |", meas[0], pred[0], pred[0] / meas[0], pred[1] / meas[1], pred[2] / meas[2],
            if !sc { "no" } else if ok { "in" } else { "**OUT**" }).unwrap();
    }
    writeln!(s, "\n### Best cell per size: IC/rho, measured against predicted (both divided by the measured rho)\n").unwrap();
    writeln!(s, "| log2 p | m | measured best k, IC/rho @0 / @0.1 / @1 | predicted best k (over the measured k), IC/rho @0 / @0.1 / @1 |").unwrap();
    writeln!(s, "|---|---|---|---|").unwrap();
    for (&(b, m), bm) in &best_meas {
        let bp = &best_pred[&(b, m)];
        let rm = rho_meas[&b];
        writeln!(s, "| {b} | {m} | {}, {:.1} / {:.1} / {:.1} | {}, {:.1} / {:.1} / {:.1} |", bm[0].1, bm[0].0 / rm[0], bm[1].0 / rm[1], bm[2].0 / rm[2],
            bp[0].1, bp[0].0 / rm[0], bp[1].0 / rm[1], bp[2].0 / rm[2]).unwrap();
    }
    let pass = all_in && med <= 0.5 && rho_ok;
    writeln!(s, "\n**Verdict:** scored cell-prices: {}; all within `[1/2, 2]`: {}; median `|log2(pred/meas)|` = {:.3} (bound 0.5); worst {:.3} ({}); rho within `[0.75, 1.33]` at `log2 p >= 32`: {}. **Validation {}.**",
        scored.len(), all_in, med, worst.0, worst.1, rho_ok, if pass { "PASSES" } else { "FAILS" }).unwrap();
    (s, pass)
}

// ---------- sweep, tables and figure ----------
const SCEN: [Sc; 3] = [Sc::Lb, Sc::Law, Sc::Ceil];
const HEAD_N: [usize; 6] = [40, 64, 80, 128, 160, 256];

fn sweep(lf: &Lf, kd: &[usize], pass: bool) -> (String, String, BTreeMap<(usize, usize, usize), f64>) {
    let all_m: Vec<usize> = (2..=8).collect();
    let mut jl = String::new();
    let mut curve: BTreeMap<(usize, usize, usize), f64> = BTreeMap::new(); // (scen idx, price idx, n) -> log2
    let scen_all = [Sc::Lb, Sc::Law, Sc::Ceil, Sc::LawDense, Sc::EnumOnly, Sc::LbOnset];
    for n in 40..=256usize {
        for (si, &sc) in scen_all.iter().enumerate() {
            for (pi_i, &pi) in PRICES.iter().enumerate() {
                let b = best_at(lf, kd, n, &all_m, sc, pi);
                curve.insert((si, pi_i, n), b.total);
                let mut per_m = String::new();
                for &m in &all_m {
                    let bm = best_at(lf, kd, n, &[m], sc, pi);
                    write!(per_m, "{}\"{m}\":{{\"log2\":{:.3},\"k\":{},\"route\":\"{}\"}}", if m > 2 { "," } else { "" }, bm.total, bm.k, bm.route).unwrap();
                }
                writeln!(jl, "{{\"log2N\":{n},\"scenario\":\"{}\",\"price\":{pi},\"log2_ic\":{:.3},\"log2_rho\":{:.3},\"m\":{},\"k\":{},\"route\":\"{}\",\"D\":{},\"log2_ga\":{:.3},\"log2_native\":{:.3},\"log2_attempts\":{:.3},\"log2_la\":{:.3},\"log2_mem_bytes\":{:.3},\"per_m\":{{{}}}}}",
                    sc.name(), b.total, rho_cost(n as f64, pi), b.m, b.k, b.route, b.d, b.ga, b.nat, b.attempts, b.la, b.mem, per_m).unwrap();
            }
        }
    }
    // headline
    let mut s = String::new();
    writeln!(s, "<!-- generated by `ic-costmodel sweep`; every number [pcc:unverified] -->\n").unwrap();
    if !pass { writeln!(s, "**FALSIFIED MODEL: the preregistered validation failed; numbers below are reported as such.**\n").unwrap(); }
    writeln!(s, "## Minimum IC cost against rho, `m = 2..8`, `k` free\n").unwrap();
    writeln!(s, "Entry: `log2(IC / rho)` (negative = IC cheaper) and the minimizer `m, k, route`. rho = `sqrt(pi N/4)(1 + 2.25 price)`.\n").unwrap();
    for (si, &sc) in scen_all.iter().enumerate() {
        writeln!(s, "### {}\n", sc.name()).unwrap();
        writeln!(s, "| log2 N | log2 rho @0 | price 0 | price 0.1 | price 1 | IC beats rho (all three prices) |").unwrap();
        writeln!(s, "|---|---|---|---|---|---|").unwrap();
        for &n in &HEAD_N {
            let mut cells = Vec::new(); let mut all = true;
            for (pi_i, &pi) in PRICES.iter().enumerate() {
                let b = best_at(lf, kd, n, &all_m, sc, pi);
                let d = b.total - rho_cost(n as f64, pi);
                if d >= 0.0 { all = false; }
                let _ = pi_i;
                cells.push(format!("{:+.1} (m={}, k={}, {}{})", d, b.m, b.k, b.route, if b.route == "groebner" { format!(", D={}", b.d) } else { String::new() }));
            }
            writeln!(s, "| {n} | {:.1} | {} | {} | {} | {} |", rho_cost(n as f64, 0.0), cells[0], cells[1], cells[2], if all { "**yes**" } else { "no" }).unwrap();
        }
        writeln!(s).unwrap();
        let _ = si;
    }
    // per-m at price 0.1 and 1 for the three headline scenarios
    for &pi in &[0.1, 1.0] {
        writeln!(s, "### Per `m`, price {pi}: `log2(IC/rho)` of the best `k` (route)\n").unwrap();
        writeln!(s, "| scenario | log2 N | {} |", all_m.iter().map(|m| format!("m={m}")).collect::<Vec<_>>().join(" | ")).unwrap();
        writeln!(s, "|---|---|{}", "---|".repeat(all_m.len())).unwrap();
        for &sc in &SCEN {
            for &n in &HEAD_N {
                let cells: Vec<String> = all_m.iter().map(|&m| {
                    let b = best_at(lf, kd, n, &[m], sc, pi);
                    format!("{:+.1} (k={}, {})", b.total - rho_cost(n as f64, pi), b.k, if b.route == "groebner" { "G".to_string() } else { b.route.replace("enum ", "") })
                }).collect();
                writeln!(s, "| {} | {n} | {} |", sc.name(), cells.join(" | ")).unwrap();
            }
        }
        writeln!(s).unwrap();
    }
    // first n where IC beats rho at all prices, per scenario
    writeln!(s, "### First `log2 N` in 40..256 at which IC beats rho at all three prices\n").unwrap();
    writeln!(s, "| scenario | first log2 N | last log2 N | number of sizes |").unwrap();
    writeln!(s, "|---|---|---|---|").unwrap();
    for (si, &sc) in scen_all.iter().enumerate() {
        let wins: Vec<usize> = (40..=256).filter(|&n| (0..3).all(|p| curve[&(si, p, n)] < rho_cost(n as f64, PRICES[p]))).collect();
        writeln!(s, "| {} | {} | {} | {} |", sc.name(), wins.first().map_or("none".into(), |x| x.to_string()), wins.last().map_or("none".into(), |x| x.to_string()), wins.len()).unwrap();
    }
    // least constant degree D (top-layer floor C(mk, D) on every attempt) that keeps IC at or above rho at price 0.1 and at price 1
    writeln!(s, "\n### Least constant degree a lower bound would need to keep a given `m` at or above rho\n").unwrap();
    writeln!(s, "Route (b) charged `C(mk, D)` per attempt for a constant `D` (for `m = 2` never below the proved floor), route (a) available; entry: least `D` with `IC >= rho` at price 0.1 / at price 1, for every `k`. `-` = no route-(b) degree is needed (route (a) and the floor already lose); `>60` = none up to 60.\n").unwrap();
    writeln!(s, "| log2 N | {} |", (2..=8).map(|m| format!("m={m}")).collect::<Vec<_>>().join(" | ")).unwrap();
    writeln!(s, "|---|{}", "---|".repeat(7)).unwrap();
    for &n in &HEAD_N {
        let mut cells = Vec::new();
        for m in 2..=8usize {
            let mut parts = Vec::new();
            for &pi in &[0.1, 1.0] {
                let rho = rho_cost(n as f64, pi);
                let mut need: Option<usize> = None;
                for d in 1..=60usize { if best_at(lf, kd, n, &[m], Sc::Fixed(d), pi).total >= rho { need = Some(d); break; } }
                parts.push(match need { Some(1) => "-".to_string(), Some(d) => d.to_string(), None => ">60".to_string() });
            }
            cells.push(parts.join(" / "));
        }
        writeln!(s, "| {n} | {} |", cells.join(" | ")).unwrap();
    }
    (jl, s, curve)
}

fn figure(curve: &BTreeMap<(usize, usize, usize), f64>) -> String {
    let (w, h) = (1140.0, 420.0);
    let (pw, ph, ml, mt) = (300.0, 300.0, 60.0, 50.0);
    let gap = 50.0;
    let ymax = 160.0; let ymin = 0.0;
    let xs = |n: f64| (n - 40.0) / 216.0 * pw;
    let ys = |y: f64| ph - (y.clamp(ymin, ymax) - ymin) / (ymax - ymin) * ph;
    let colors = ["#2b6cb0", "#dd6b20", "#c53030", "#805ad5", "#718096"];
    let names = ["LB (proved floor)", "LAW (measured, extrapolated)", "CEIL (proved ceiling)", "LAW, dense solver", "enumeration only"];
    let mut s = String::new();
    writeln!(s, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\" font-family=\"sans-serif\" font-size=\"11\">").unwrap();
    writeln!(s, "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>").unwrap();
    writeln!(s, "<text x=\"{}\" y=\"18\" font-size=\"13\" text-anchor=\"middle\">PKM index calculus (min over m = 2..8, k) against Pollard rho: log2 cost in GAE [pcc:unverified]</text>", w / 2.0).unwrap();
    for (pi_i, &pi) in PRICES.iter().enumerate() {
        let ox = ml + pi_i as f64 * (pw + gap);
        writeln!(s, "<g transform=\"translate({ox},{mt})\">").unwrap();
        writeln!(s, "<rect width=\"{pw}\" height=\"{ph}\" fill=\"none\" stroke=\"#999\"/>").unwrap();
        writeln!(s, "<text x=\"{}\" y=\"-8\" text-anchor=\"middle\">native operations at {pi} GAE</text>", pw / 2.0).unwrap();
        for t in (0..=160).step_by(20) { let y = ys(t as f64);
            writeln!(s, "<line x1=\"0\" x2=\"{pw}\" y1=\"{y}\" y2=\"{y}\" stroke=\"#eee\"/><text x=\"-6\" y=\"{}\" text-anchor=\"end\">{t}</text>", y + 4.0).unwrap(); }
        for t in [40, 64, 80, 128, 160, 256] { let x = xs(t as f64);
            writeln!(s, "<line x1=\"{x}\" x2=\"{x}\" y1=\"0\" y2=\"{ph}\" stroke=\"#eee\"/><text x=\"{x}\" y=\"{}\" text-anchor=\"middle\">{t}</text>", ph + 14.0).unwrap(); }
        writeln!(s, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\">log2 N</text>", pw / 2.0, ph + 30.0).unwrap();
        if pi_i == 0 { writeln!(s, "<text transform=\"translate(-40,{}) rotate(-90)\" text-anchor=\"middle\">log2 cost (GAE)</text>", ph / 2.0).unwrap(); }
        for si in 0..5 {
            let mut d = String::new();
            for n in 40..=256usize { let v = curve[&(si, pi_i, n)]; write!(d, "{}{:.1},{:.1} ", if n == 40 { "M" } else { "L" }, xs(n as f64), ys(v)).unwrap(); }
            let dash = if si >= 3 { " stroke-dasharray=\"5,3\"" } else { "" };
            writeln!(s, "<path d=\"{d}\" fill=\"none\" stroke=\"{}\" stroke-width=\"{}\"{dash}/>", colors[si], if si >= 3 { 1.2 } else { 2.0 }).unwrap();
        }
        let mut d = String::new();
        for n in 40..=256usize { write!(d, "{}{:.1},{:.1} ", if n == 40 { "M" } else { "L" }, xs(n as f64), ys(rho_cost(n as f64, pi))).unwrap(); }
        writeln!(s, "<path d=\"{d}\" fill=\"none\" stroke=\"black\" stroke-width=\"2\"/>").unwrap();
        writeln!(s, "</g>").unwrap();
    }
    let ly = h - 22.0;
    let mut x = ml;
    for si in 0..5 { writeln!(s, "<line x1=\"{x}\" x2=\"{}\" y1=\"{ly}\" y2=\"{ly}\" stroke=\"{}\" stroke-width=\"2\"/><text x=\"{}\" y=\"{}\">{}</text>", x + 20.0, colors[si], x + 24.0, ly + 4.0, names[si]).unwrap(); x += 185.0; }
    writeln!(s, "<line x1=\"{x}\" x2=\"{}\" y1=\"{ly}\" y2=\"{ly}\" stroke=\"black\" stroke-width=\"2\"/><text x=\"{}\" y=\"{}\">rho sqrt(pi N/4)</text>", x + 20.0, x + 24.0, ly + 4.0).unwrap();
    writeln!(s, "</svg>").unwrap();
    s
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("all");
    let lf = Lf::new(4096);
    let kd: Vec<usize> = (0..12).map(|d| if d < 2 { usize::MAX } else { d + lambda_w(&lf, d - 1) }).collect();
    let here = env!("CARGO_MANIFEST_DIR");
    let adir = format!("{here}/../ic-endtoend/results");
    match cmd {
        "kd" => {
            for d in 2..12 { println!("D = {d}: Lambda(D-1) = {}, k_D = {}, 14D+10 = {}", lambda_w(&lf, d - 1), kd[d], 14 * d + 10); }
            for k in [2, 6, 10, 14, 20, 27, 39, 52, 64, 77, 100, 127, 200, 254] {
                println!("k = {k}: D_LB(m=2) = {}, paper floor((k-10)/14)+1 = {}, LAW(m=2) = {}, LAW(m=3) = {}, LAW(m=4) = {}, LAW(m=8) = {}", d_lb2(&kd, k), if k >= 10 { (k - 10) / 14 + 1 } else { 0 }, d_law(2, k), d_law(3, k), d_law(4, k), d_law(8, k));
            }
        }
        "validate" | "all" | "sweep" => {
            let (v, pass) = validate(&lf, &kd, &adir);
            fs::write(format!("{here}/results/validation.md"), &v).unwrap();
            println!("{}", v.lines().last().unwrap_or(""));
            if cmd == "validate" { return; }
            let (jl, head, curve) = sweep(&lf, &kd, pass);
            fs::write(format!("{here}/results/sweep.jsonl"), jl).unwrap();
            fs::write(format!("{here}/results/headline.md"), &head).unwrap();
            fs::write(format!("{here}/results/fig-cost.svg"), figure(&curve)).unwrap();
            println!("{head}");
        }
        _ => eprintln!("usage: ic-costmodel [kd|validate|sweep]"),
    }
}
