//! hky_lift: degree lift V_D -> V_{D+1} n M_<=D from an hky.basis/v1 export.
//!
//! Subcommands:
//!   lift     --basis F [--a A --b B --X X] [--include-low] [--batch-tops N] [--fall-every N]
//!            [--no-verify-falls] [--no-hash] [--parts] [--cache-gb G] [--trials N] [--no-refuse]
//!            [--test-no-fall-shifts] [--test-drop-shift S] [--out J] [--expect J]
//!   closure  --p P --a A --b B --X X --k K --D D [--order 0|1] [--export F] [--scramble N] [--seed S] [--out J]
//!   validate --p P --a A --b B --X X --k K --D D [--order 0|1] [--scramble N] [--dir DIR] [--parts]
//!            [--expect-dim N --expect-one 0|1 --expect-profile a,b,..]  (lift vs direct closure)
//!   verify   --got J --expected J        (field-by-field comparison with a recorded lift JSON)
//!   corrupt  --in F --out G --mode drop-row:I | perturb:I | flip-byte:OFF
//!
//! Exit codes: 0 ok; 2 input refused; 3 mismatch against an expectation; 1 usage/other error.

use hky_lift::basis::{read_hky_basis, write_hky_basis, Rref};
use hky_lift::closure::{chain_generator, closure, to_echelon};
use hky_lift::field::{madds, Fp};
use hky_lift::layout::Layout;
use hky_lift::lift::{input_checks, Lift, LiftOpts};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::time::Instant;

struct Args {
    kv: HashMap<String, String>,
    flags: Vec<String>,
}

impl Args {
    fn parse(v: &[String]) -> Args {
        let mut kv = HashMap::new();
        let mut flags = Vec::new();
        let mut i = 0;
        while i < v.len() {
            let a = &v[i];
            if let Some(k) = a.strip_prefix("--") {
                if i + 1 < v.len() && !v[i + 1].starts_with("--") {
                    kv.insert(k.to_string(), v[i + 1].clone());
                    i += 2;
                    continue;
                }
                flags.push(k.to_string());
            }
            i += 1;
        }
        Args { kv, flags }
    }
    fn get(&self, k: &str) -> Option<&str> {
        self.kv.get(k).map(|s| s.as_str())
    }
    fn num<T: std::str::FromStr>(&self, k: &str) -> Option<T> {
        self.get(k).map(|s| s.parse().unwrap_or_else(|_| die(&format!("--{k}: bad number {s}"))))
    }
    fn req<T: std::str::FromStr>(&self, k: &str) -> T {
        self.num(k).unwrap_or_else(|| die(&format!("--{k} is required")))
    }
    fn flag(&self, k: &str) -> bool {
        self.flags.iter().any(|f| f == k)
    }
}

fn die(m: &str) -> ! {
    eprintln!("error: {m}");
    std::process::exit(1)
}

fn maxrss_mb() -> f64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for l in s.lines() {
        if let Some(r) = l.strip_prefix("VmHWM:") {
            let kb: f64 = r.trim().trim_end_matches("kB").trim().parse().unwrap_or(0.0);
            return kb / 1024.0;
        }
    }
    0.0
}

fn opts_from(a: &Args) -> LiftOpts {
    let mut o = LiftOpts::default();
    if let Some(n) = a.num("batch-tops") {
        o.batch_tops = n;
    }
    if let Some(n) = a.num("fall-every") {
        o.fall_every = n;
    }
    if let Some(g) = a.num("cache-gb") {
        o.follower_cache_gb = g;
    }
    if let Some(g) = a.num("work-gb") {
        o.work_gb = g;
    }
    o.include_low = a.flag("include-low");
    o.verify_falls = !a.flag("no-verify-falls");
    o.hash = !a.flag("no-hash");
    o.parts = a.flag("parts");
    o.test_no_fall_shifts = a.flag("test-no-fall-shifts");
    o.test_drop_shift = a.num("test-drop-shift");
    o
}

/// run the lift on a file; returns (record, refused)
fn run_lift(path: &str, a: &Args, g_override: Option<Vec<(u64, u32)>>) -> (Value, bool) {
    let t0 = Instant::now();
    let (ech, info) = match read_hky_basis(path) {
        Ok(x) => x,
        Err(e) => {
            return (json!({"status": "refused", "reason": e, "basis": path}), true);
        }
    };
    let fp = Fp::new(ech.p);
    let t_read = t0.elapsed().as_secs_f64();
    let vd = Rref::from_echelon(&fp, &ech);
    let t_rref = t0.elapsed().as_secs_f64();
    // instance for the generator check: flags, else sidecar
    let mut inst: Option<(i64, i64, i64)> = None;
    if let (Some(aa), Some(bb), Some(xx)) = (a.num::<i64>("a"), a.num::<i64>("b"), a.num::<i64>("X")) {
        inst = Some((aa, bb, xx));
    } else if let Ok(s) = std::fs::read_to_string(format!("{path}.json")) {
        if let Ok(v) = serde_json::from_str::<Value>(&s) {
            let i = if v["instance"].is_null() { &v["extra"]["instance"] } else { &v["instance"] };
            if let (Some(aa), Some(bb), Some(xx)) = (i["a"].as_i64(), i["b"].as_i64(), i["target"][0].as_i64()) {
                inst = Some((aa, bb, xx));
            }
        }
    }
    let g = g_override.or_else(|| inst.map(|(aa, bb, xx)| chain_generator(&fp, ech.k, aa, bb, xx)));
    let trials = a.num("trials").unwrap_or(3usize);
    let checks = input_checks(&fp, &ech.lay, &vd, g.as_deref(), trials, 0x9e3779b97f4a7c15);
    let t_checks = t0.elapsed().as_secs_f64();
    let leave = checks["closed_random"]["leave_VD"].as_u64().unwrap_or(0);
    let gbad = checks["generators"]["not_in_VD"].as_u64().unwrap_or(0);
    let refuse = (leave > 0 || gbad > 0) && !a.flag("no-refuse");
    if refuse {
        return (
            json!({"status": "refused",
                   "reason": if gbad > 0 {"the generator g (or some m g) is not in the stored space"} else {"the stored space is not closed: shifted random low elements leave it"},
                   "input_checks": checks, "basis": path, "file_sha256": info.file_sha256}),
            true,
        );
    }
    let opts = opts_from(a);
    let mut lf = Lift::new(fp, &ech.lay, &vd, opts);
    lf.run();
    let (mut rec, _w) = lf.finish(t0);
    rec["status"] = json!("done");
    rec["input"] = json!({
        "basis": path, "file_sha256": info.file_sha256, "content_sha256": info.content_sha256,
        "bytes": info.bytes, "column_order": ech.column_order, "profile_VD": ech.profile,
        "instance_a_b_X": inst, "nnz_stored": ech.nnz,
    });
    rec["input_checks"] = checks;
    rec["timing_s"] = json!({"read": t_read, "rref": t_rref - t_read, "checks": t_checks - t_rref,
                            "lift": t0.elapsed().as_secs_f64() - t_checks, "total": t0.elapsed().as_secs_f64()});
    rec["maxrss_mb"] = json!(maxrss_mb());
    rec["threads"] = json!(rayon::current_num_threads());
    rec["binary_sha256"] = json!(self_sha());
    (rec, false)
}

fn self_sha() -> String {
    std::env::current_exe()
        .ok()
        .and_then(|p| std::fs::read(p).ok())
        .map(|b| hky_lift::basis::sha256_hex(&b))
        .unwrap_or_default()
}

/// field-by-field comparison of our record with a recorded lift JSON (either implementation)
fn verify(got: &Value, exp: &Value) -> Value {
    let mut rows = vec![];
    let mut ok = true;
    let mut cmp = |name: &str, g: Value, e: Value| {
        if e.is_null() {
            rows.push(json!({"field": name, "expected": null, "got": g, "result": "not in expected record"}));
            return;
        }
        let same = g == e;
        ok &= same;
        rows.push(json!({"field": name, "expected": e, "got": g, "result": if same {"equal"} else {"DIFFERENT"}}));
    };
    cmp("k", got["k"].clone(), exp["k"].clone());
    cmp("D", got["D"].clone(), exp["D"].clone());
    cmp("p", got["p"].clone(), exp["p"].clone());
    cmp("dim_VD", got["dim_VD"].clone(), exp["dim_VD"].clone());
    cmp("contains_one", got["contains_one"].clone(), exp["contains_one"].clone());
    cmp("fixed_point", got["fixed_point"].clone(), exp["fixed_point"].clone());
    cmp("dim_VD1", got["dim_VD1"].clone(), exp["dim_VD1"].clone());
    cmp("profile", got["profile"].clone(), exp["profile"].clone());
    // dim F = sum of the per-part falls in the split implementations
    let ef = match &exp["falls"] {
        Value::Array(v) => json!(v.iter().map(|x| x.as_u64().unwrap_or(0)).sum::<u64>()),
        x => x.clone(),
    };
    cmp("falls (dim F = dim(V_{D+1} n L) - dim V_D)", got["falls"].clone(), ef);
    let gl = got["dim_low"].clone();
    let el = match &exp["low_part_ranks"] {
        Value::Array(v) => json!(v.iter().map(|x| x.as_u64().unwrap_or(0)).sum::<u64>()),
        _ => Value::Null,
    };
    cmp("dim(V_{D+1} n M_<=D)", gl, el);
    if !got["sigma_parts"].is_null() {
        cmp(
            "part_profiles_low",
            got["sigma_parts"]["part_profiles_low"].clone(),
            exp["part_profiles_low"].clone(),
        );
    }
    let regsum = match &exp["regular_top_orbits"] {
        Value::Array(v) => json!(v.iter().map(|x| x.as_u64().unwrap_or(0)).sum::<u64>()),
        _ => Value::Null,
    };
    rows.push(json!({"field": "regular_tops (orbit count summed over parts; not the same notion)", "expected": regsum, "got": got["regular_tops"], "result": "informational"}));
    json!({"all_equal": ok, "fields": rows})
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 {
        die("usage: hky_lift lift|closure|validate|verify|corrupt ... (see src/main.rs)");
    }
    let a = Args::parse(&argv[2..]);
    match argv[1].as_str() {
        "lift" => {
            let path = a.get("basis").unwrap_or_else(|| die("--basis required")).to_string();
            let (mut rec, refused) = run_lift(&path, &a, None);
            let mut code = if refused { 2 } else { 0 };
            if let Some(e) = a.get("expect") {
                let exp: Value = serde_json::from_str(&std::fs::read_to_string(e).unwrap()).unwrap();
                let v = verify(&rec, &exp);
                if !v["all_equal"].as_bool().unwrap() {
                    code = 3;
                }
                rec["verify"] = v;
                rec["verify_expected_file"] = json!(e);
            }
            let s = serde_json::to_string_pretty(&rec).unwrap();
            if let Some(o) = a.get("out") {
                std::fs::write(o, &s).unwrap();
            }
            println!("{s}");
            std::process::exit(code);
        }
        "closure" => {
            let t0 = Instant::now();
            let p: u64 = a.req("p");
            let fp = Fp::new(p);
            let k: u32 = a.req("k");
            let d: u32 = a.req("D");
            let g = chain_generator(&fp, k, a.req("a"), a.req("b"), a.req("X"));
            let out = closure(fp, k, d, &g);
            let rref = &out.rref;
            let prof = hky_lift::basis::profile_of(&out.lay, &rref.piv, d);
            let one = rref.piv.last() == Some(&((out.lay.ncols() - 1) as u32));
            let mut rec = json!({"k": k, "D": d, "p": p, "dim": rref.rank(), "n_cols": out.lay.ncols(),
                "profile": prof, "contains_one": one, "closed": out.closed, "sweep_images": out.sweep_images,
                "rref_sha": rref.sha(&out.lay), "seconds": t0.elapsed().as_secs_f64(), "madds": madds(), "maxrss_mb": maxrss_mb()});
            if let Some(path) = a.get("export") {
                let order: u32 = a.num("order").unwrap_or(0);
                let ech = if order == 0 {
                    to_echelon(&fp, out.lay, rref, d, k, a.num("scramble").unwrap_or(0), a.num("seed").unwrap_or(7))
                } else {
                    let lay1 = Layout::order1(k, d);
                    let r1 = reorder(&fp, &out.lay, rref, &lay1);
                    let mut e = to_echelon(&fp, lay1, &r1, d, k, a.num("scramble").unwrap_or(0), a.num("seed").unwrap_or(7));
                    e.column_order = 1;
                    e
                };
                let sha = write_hky_basis(
                    path,
                    &ech,
                    json!({"instance": {"prime": p, "a": a.req::<i64>("a"), "b": a.req::<i64>("b"), "target": [a.req::<i64>("X")], "k": k}}),
                )
                .unwrap();
                rec["export"] = json!({"path": path, "file_sha256": sha, "order": order});
            }
            let s = serde_json::to_string_pretty(&rec).unwrap();
            if let Some(o) = a.get("out") {
                std::fs::write(o, &s).unwrap();
            }
            println!("{s}");
        }
        "validate" => {
            let code = validate(&a);
            std::process::exit(code);
        }
        "verify" => {
            let g: Value = serde_json::from_str(&std::fs::read_to_string(a.get("got").unwrap()).unwrap()).unwrap();
            let e: Value = serde_json::from_str(&std::fs::read_to_string(a.get("expected").unwrap()).unwrap()).unwrap();
            let v = verify(&g, &e);
            println!("{}", serde_json::to_string_pretty(&v).unwrap());
            std::process::exit(if v["all_equal"].as_bool().unwrap() { 0 } else { 3 });
        }
        "corrupt" => corrupt(&a),
        "bench" => {
            let v = bench_gemm(a.num("n").unwrap_or(2048), a.num("m").unwrap_or(3000), a.num("c").unwrap_or(4000));
            println!("{v}");
        }
        x => die(&format!("unknown subcommand {x}")),
    }
}

/// RREF of the same space in another column order (small cases)
fn reorder(fp: &Fp, from: &Layout, r: &Rref, to: &Layout) -> Rref {
    let n = to.ncols();
    let nf = r.nf();
    let vs: Vec<Vec<u32>> = (0..r.rank())
        .map(|i| {
            let mut v = vec![0u32; n];
            v[to.col(from.masks[r.piv[i] as usize]).unwrap() as usize] = 1;
            for f in 0..nf {
                v[to.col(from.masks[r.free[f] as usize]).unwrap() as usize] = r.b[i * nf + f];
            }
            v
        })
        .collect();
    let mut rows = Vec::new();
    let mut piv = Vec::new();
    hky_lift::field::rref_extend(fp, &mut rows, &mut piv, vs, n);
    let mut out = Rref::new_index(n, piv.iter().map(|&x| x as u32).collect());
    let nf2 = out.nf();
    let free = out.free.clone();
    for (i, row) in rows.iter().enumerate() {
        for (f, &c) in free.iter().enumerate() {
            out.b[i * nf2 + f] = row[c as usize];
        }
    }
    out
}

fn validate(a: &Args) -> i32 {
    let t0 = Instant::now();
    let p: u64 = a.req("p");
    let fp = Fp::new(p);
    let (k, d): (u32, u32) = (a.req("k"), a.req("D"));
    let (ca, cb, cx): (i64, i64, i64) = (a.req("a"), a.req("b"), a.req("X"));
    let g = chain_generator(&fp, k, ca, cb, cx);
    let dir = a.get("dir").unwrap_or("/tmp").to_string();
    let order: u32 = a.num("order").unwrap_or(0);
    let scr: usize = a.num("scramble").unwrap_or(3);
    let tag = format!("p{p}-a{ca}-b{cb}-X{cx}-k{k}-D{d}-o{order}");
    // 1. V_D by direct closure, exported (scrambled echelon) and read back by the lift
    let vd = closure(fp, k, d, &g);
    assert!(vd.closed);
    let t_vd = t0.elapsed().as_secs_f64();
    let path = format!("{dir}/vd-{tag}.hkybasis");
    let (ech, lay_o) = if order == 0 {
        (to_echelon(&fp, Layout::order0(k, d), &vd.rref, d, k, scr, 11), Layout::order0(k, d))
    } else {
        let l1 = Layout::order1(k, d);
        let r1 = reorder(&fp, &vd.lay, &vd.rref, &l1);
        let mut e = to_echelon(&fp, Layout::order1(k, d), &r1, d, k, scr, 11);
        e.column_order = 1;
        (e, l1)
    };
    // with --vd-basis F the lift reads F (e.g. the engine's export) instead of our own export, and
    // the stored space is compared with our direct closure V_D (canonical hash in F's column order)
    let ext = a.get("vd-basis").map(|s| s.to_string());
    let mut vd_cmp = Value::Null;
    let path = match &ext {
        None => {
            write_hky_basis(&path, &ech, json!({"instance": {"prime": p, "a": ca, "b": cb, "target": [cx], "k": k}})).unwrap();
            path
        }
        Some(f) => {
            let (e2, info) = read_hky_basis(f).unwrap_or_else(|x| die(&x));
            let r2 = Rref::from_echelon(&fp, &e2);
            let ours = Rref::from_echelon(&fp, &ech);
            vd_cmp = json!({"file": f, "file_sha256": info.file_sha256, "column_order": e2.column_order,
                "same_order_as_ours": e2.column_order == ech.column_order,
                "rank_file": r2.rank(), "rank_ours": ours.rank(),
                "sha_file": r2.sha(&e2.lay), "sha_ours": ours.sha(&ech.lay),
                "equal": e2.column_order == ech.column_order && r2.sha(&e2.lay) == ours.sha(&ech.lay)});
            f.clone()
        }
    };
    let mut la = Args { kv: a.kv.clone(), flags: a.flags.clone() };
    la.kv.remove("a");
    let t1 = Instant::now();
    let (lrec, refused) = run_lift(&path, &la, Some(g.clone()));
    let t_lift = t1.elapsed().as_secs_f64();
    if refused {
        println!("{}", serde_json::to_string_pretty(&lrec).unwrap());
        return 2;
    }
    // 2. direct closure V_{D+1} (oracle)
    let t2 = Instant::now();
    let v1 = closure(fp, k, d + 1, &g);
    let t_direct = t2.elapsed().as_secs_f64();
    let n1 = v1.lay.ncols();
    let prof1: Vec<u64> = hky_lift::basis::profile_of(&v1.lay, &v1.rref.piv, d + 1);
    let one1 = v1.rref.piv.last() == Some(&((n1 - 1) as u32));
    // low rows of the direct RREF (pivot weight <= D) restricted to M_<=D, re-expressed in the
    // lift's column order (order 0: same RREF; order 1: re-echelonised)
    let lay_d0 = Layout::order0(k, d);
    let off = n1 - lay_d0.ncols();
    let lowpiv: Vec<u32> = v1.rref.piv.iter().filter(|&&c| c as usize >= off).map(|&c| c - off as u32).collect();
    let mut low = Rref::new_index(lay_d0.ncols(), lowpiv);
    let nfl = low.nf();
    let nf1 = v1.rref.nf();
    for (i, &pc) in low.piv.clone().iter().enumerate() {
        let r = v1.rref.prow[pc as usize + off] as usize;
        for (f, &c) in low.free.clone().iter().enumerate() {
            let f1 = v1.rref.fpos[c as usize + off];
            low.b[i * nfl + f] = if f1 == u32::MAX { 0 } else { v1.rref.b[r * nf1 + f1 as usize] };
        }
        // entries of a low row on top columns must be zero (they are pivots-free top columns)
        for f in 0..nf1 {
            if (v1.rref.free[f] as usize) < off {
                assert_eq!(v1.rref.b[r * nf1 + f], 0, "low row with a top entry");
            }
        }
    }
    let low_o = if order == 0 { low } else { reorder(&fp, &lay_d0, &low, &lay_o) };
    let direct_sha = low_o.sha(&lay_o);
    let mut parts = Value::Null;
    if a.flag("parts") {
        let w = hky_lift::lift::WSpace { piv: low_o.piv.clone(), free: low_o.free.clone(), b: low_o.b.clone() };
        parts = hky_lift::lift::sigma_parts(&fp, &lay_o, &w, d);
    }
    let checks = [
        ("dim_VD1", json!(v1.rref.rank()), lrec["dim_VD1"].clone()),
        ("profile", json!(prof1), lrec["profile"].clone()),
        ("contains_one", json!(one1), lrec["contains_one"].clone()),
        ("low_basis_sha256", json!(direct_sha), lrec["low_basis_sha256"].clone()),
    ];
    let mut ok = true;
    if !vd_cmp.is_null() {
        ok &= vd_cmp["equal"].as_bool().unwrap_or(false);
    }
    let mut rows = vec![];
    for (n, e, gv) in checks.iter() {
        let same = e == gv;
        ok &= same;
        rows.push(json!({"field": n, "direct": e, "lift": gv, "equal": same}));
    }
    if !parts.is_null() {
        let same = parts == lrec["sigma_parts"];
        ok &= same;
        rows.push(json!({"field": "sigma_parts", "direct": parts, "lift": lrec["sigma_parts"], "equal": same}));
    }
    // optional recorded expectations (e.g. from validation.json of the paper repository)
    if let Some(e) = a.num::<u64>("expect-dim") {
        let same = json!(e) == lrec["dim_VD1"];
        ok &= same;
        rows.push(json!({"field": "dim_VD1 vs recorded", "recorded": e, "lift": lrec["dim_VD1"], "equal": same}));
    }
    if let Some(e) = a.num::<u32>("expect-one") {
        let same = json!(e == 1) == lrec["contains_one"];
        ok &= same;
        rows.push(json!({"field": "contains_one vs recorded", "recorded": e == 1, "lift": lrec["contains_one"], "equal": same}));
    }
    if let Some(e) = a.get("expect-profile") {
        let ev: Vec<u64> = e.split(',').map(|x| x.parse().unwrap()).collect();
        let same = json!(ev) == lrec["profile"];
        ok &= same;
        rows.push(json!({"field": "profile vs recorded", "recorded": ev, "lift": lrec["profile"], "equal": same}));
    }
    let rec = json!({
        "case": tag, "k": k, "D": d, "p": p, "order": order, "scramble": scr,
        "dim_VD": vd.rref.rank(), "direct_VD1_closed": v1.closed, "vd_basis_vs_our_closure": vd_cmp,
        "all_equal": ok, "checks": rows,
        "lift": {"dim_T": lrec["dim_T"], "dim_K": lrec["dim_K"], "special_tops": lrec["special_tops"],
                 "regular_tops": lrec["regular_tops"], "falls": lrec["falls"], "rounds": lrec["rounds"],
                 "verify_falls": lrec["verify_falls"], "input_checks": lrec["input_checks"],
                 "counts": lrec["counts"], "settings": lrec["settings"]},
        "seconds": {"closure_VD": t_vd, "lift_incl_read": t_lift, "direct_VD1": t_direct},
        "maxrss_mb": maxrss_mb(), "threads": rayon::current_num_threads(),
    });
    println!("{}", serde_json::to_string(&rec).unwrap());
    if ext.is_none() {
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(format!("{path}.json"));
    }
    if ok {
        0
    } else {
        3
    }
}

fn corrupt(a: &Args) {
    let inp = a.get("in").unwrap();
    let outp = a.get("out").unwrap();
    let mode = a.get("mode").unwrap();
    if let Some(off) = mode.strip_prefix("flip-byte:") {
        let mut b = std::fs::read(inp).unwrap();
        let o: usize = off.parse().unwrap();
        b[o] ^= 1;
        std::fs::write(outp, b).unwrap();
        let _ = std::fs::copy(format!("{inp}.json"), format!("{outp}.json"));
        println!("{{\"mode\": \"{mode}\", \"note\": \"raw byte flipped, trailer not updated\"}}");
        return;
    }
    let (mut e, _) = hky_lift::basis::read_hky_basis_opt(inp, true).unwrap_or_else(|x| die(&x));
    let fp = Fp::new(e.p);
    let note;
    if let Some(i) = mode.strip_prefix("drop-row:") {
        let i: usize = i.parse().unwrap();
        let i = i.min(e.piv.len() - 1);
        let (s, t) = (e.off[i] as usize, e.off[i + 1] as usize);
        let removed_piv = e.piv[i];
        e.ents.drain(s..t);
        e.piv.remove(i);
        let w = (t - s) as u64;
        e.off.remove(i + 1);
        for x in e.off.iter_mut().skip(i + 1) {
            *x -= w;
        }
        e.profile = hky_lift::basis::profile_of(&e.lay, &e.piv, e.d);
        let n = e.lay.ncols() as u32;
        e.flags = if e.piv.last() == Some(&(n - 1)) { 1 } else { 0 };
        note = format!("row {i} (pivot column {removed_piv}, weight {}) removed; header and trailer rewritten", e.lay.wt[removed_piv as usize]);
    } else if let Some(i) = mode.strip_prefix("perturb:") {
        let i: usize = i.parse().unwrap();
        let i = i.min(e.piv.len() - 1);
        let (s, t) = (e.off[i] as usize, e.off[i + 1] as usize);
        if s < t {
            e.ents[s].1 = fp.add(e.ents[s].1, 1);
            if e.ents[s].1 == 0 {
                e.ents[s].1 = 1;
            }
            note = format!("row {i}: coefficient at column {} + 1", e.ents[s].0);
        } else {
            // empty row: add an entry on the last column
            let n = e.lay.ncols() as u32;
            e.ents.insert(s, (n - 1, 1));
            for x in e.off.iter_mut().skip(i + 1) {
                *x += 1;
            }
            note = format!("row {i}: entry 1 added at the last column");
        }
    } else {
        die("mode: drop-row:I | perturb:I | flip-byte:OFF");
    }
    let side: Value = std::fs::read_to_string(format!("{inp}.json"))
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or(Value::Null);
    let sha = write_hky_basis(outp, &e, json!({"corrupted_from": inp, "note": note, "instance": side["extra"]["instance"]})).unwrap();
    // keep the instance where run_lift looks for it
    let mut s2: Value = serde_json::from_str(&std::fs::read_to_string(format!("{outp}.json")).unwrap()).unwrap();
    s2["instance"] = side["instance"].clone();
    if s2["instance"].is_null() {
        s2["instance"] = side["extra"]["instance"].clone();
    }
    std::fs::write(format!("{outp}.json"), serde_json::to_string_pretty(&s2).unwrap()).unwrap();
    println!("{}", json!({"mode": mode, "note": note, "file_sha256": sha}));
}

/// gemm rate on this host: hky_lift-bench N M C
pub fn bench_gemm(n: usize, m: usize, c: usize) -> Value {
    let fp = Fp::new(65537);
    let a: Vec<u32> = (0..n * m).map(|i| (i as u64 * 2654435761 % 65537) as u32).collect();
    let b: Vec<u32> = (0..m * c).map(|i| (i as u64 * 40503 % 65537) as u32).collect();
    let mut o = vec![0u32; n * c];
    let t = Instant::now();
    hky_lift::field::gemm(&fp, &mut o, &a, &b, n, m, c, false);
    let s = t.elapsed().as_secs_f64();
    json!({"n": n, "m": m, "c": c, "seconds": s, "gmadds_per_s": (n * m * c) as f64 / s / 1e9, "threads": rayon::current_num_threads()})
}
