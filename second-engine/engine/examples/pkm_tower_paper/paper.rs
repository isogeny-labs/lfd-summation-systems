//! `pkm_tower_pilot --instance <paper instance JSON> --D <degree> …`
//!
//! Builds exactly the paper's normalized two-summand prime chain system
//! (`lfd-summation-systems/paper/curve-lfd-bounds.tex`; `csd/ic_dreg.py`
//! `pkm_generator`, `csd/ic_bounds_d5.py`) on the pilot's Kummer tower ring and
//! computes one closure V_D. Correspondence (derived in
//! `research/hky_paper_runs_20261118/README.md`): Kummer tower with scale 1 and
//! top constant 1, so the engine's y_{i,j} is the paper's chain variable
//! u_{i,j} = X_i^{2^j} with no rescaling; x_i = y_{i,0}; the generator is
//! NF(S_3(x_1, x_2, X_R)) with X_R the target abscissa. The generator the
//! pilot derives is compared, coefficient by coefficient, with S_3's
//! coefficients rewritten directly, and with the instance file's `generator`.
//!
//! Flags: `--instance PATH` `--D N` `[--engine hky3|hky2]` `[--sym on|off]`
//! `[--out PATH]` `[--export-basis PATH]` `[--checkpoint-dir DIR]`
//! `[--checkpoint-every MINUTES]` `[--checkpoint-every-batches N]`
//! `[--resume DIR]` `[--max-entries N]` `[--budget SECONDS]`
//! `[--hky3-verify-closed on]` (engine-opt: re-check closedness at the end; doubles memory)
//! `[--hky3-skip on]` (engine-skip: skip the products PROOF.md proves to lie in the span; same space)
//! `[--hky3-fault N]` (negative controls only)
//! `[--override-a A]` `[--override-x X]` (negative controls: change the curve
//! coefficient or the target abscissa; the cross-check against the file's
//! generator is then skipped and the record says so).

use super::*;
use crypto_lib::cryptanalysis::hky_closure_io::{self, IoOptions};
use crypto_lib::cryptanalysis::hky_fp_tower::NfPoly;
use sha2::{Digest, Sha256};

fn sha256_file(path: &str) -> String {
    let mut h = Sha256::new();
    let mut f = File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = std::io::Read::read(&mut f, &mut buf).expect("read");
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    h.finalize().iter().map(|x| format!("{x:02x}")).collect()
}

fn status_kb(field: &str) -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines()
        .find(|l| l.starts_with(field))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|x| x.parse().ok())
}

/// The order-2^k subgroup generator of `csd/ic_dreg.subgroup_2power`.
pub fn subgroup_generator(p: u64, k: usize) -> u64 {
    let n = 1u64 << k;
    assert!((p - 1).is_multiple_of(n), "2^{k} does not divide p − 1");
    for h in 2..p {
        let g = powm(h, (p - 1) / n, p);
        if powm(g, n / 2, p) != 1 {
            return g;
        }
    }
    unreachable!()
}

/// `ic_dreg.s3_prime_coeffs` then `pkm_generator`: S_3(x1, x2, X) as
/// {chain mask: coefficient}, x_i^e ↦ bits of e mod 2^k in block i.
pub fn paper_generator(p: u64, a: u64, b: u64, x: u64, k: usize) -> BTreeMap<u64, u64> {
    let (a, b, x) = (a % p, b % p, x % p);
    let m = |v: i128| -> u64 { (v.rem_euclid(p as i128)) as u64 };
    let (ai, bi, xi) = (a as i128, b as i128, x as i128);
    let pi = p as i128;
    let xx = (xi * xi) % pi;
    let c: [((u64, u64), u64); 9] = [
        ((2, 0), m(xx)),
        ((0, 2), m(xx)),
        ((1, 1), m(-2 * xx - 2 * ai)),
        ((2, 1), m(-2 * xi)),
        ((1, 2), m(-2 * xi)),
        ((1, 0), m(-2 * ai * xi - 4 * bi)),
        ((0, 1), m(-2 * ai * xi - 4 * bi)),
        ((0, 0), m(-4 * bi * xi + ai * ai)),
        ((2, 2), 1),
    ];
    let n = 1u64 << k;
    let mut g = BTreeMap::new();
    for ((e1, e2), v) in c {
        if v == 0 {
            continue;
        }
        let key = (e1 % n) | ((e2 % n) << k);
        let e = g.entry(key).or_insert(0u64);
        *e = (*e + v) % p;
    }
    g.retain(|_, v| *v != 0);
    g
}

struct PArgs {
    instance: String,
    d: u32,
    engine: Engine,
    sym: bool,
    out: Option<String>,
    io: IoOptions,
    max_entries: u64,
    budget: u64,
    override_a: Option<u64>,
    override_x: Option<u64>,
    verify_closed: bool,
    skip: bool,
    skip_order: u8,
    skip_witness: u8,
    skip_witness_max: u16,
    fault: u8,
}

fn parse(argv: &[String]) -> PArgs {
    let mut a = PArgs {
        instance: String::new(),
        d: 0,
        engine: Engine::Hky3,
        sym: true,
        out: None,
        io: IoOptions::default(),
        max_entries: 40_000_000_000,
        budget: 60 * 86_400,
        override_a: None,
        override_x: None,
        verify_closed: false,
        skip: false,
        skip_order: 1,
        skip_witness: 1,
        skip_witness_max: 0,
        fault: 0,
    };
    let mut i = 1;
    while i < argv.len() {
        let v = argv.get(i + 1).cloned().unwrap_or_default();
        match argv[i].as_str() {
            "--instance" => a.instance = v,
            "--D" => a.d = v.parse().expect("--D"),
            "--engine" => {
                a.engine = match v.as_str() {
                    "hky3" => Engine::Hky3,
                    "hky2" => Engine::Hky2,
                    o => panic!("--engine hky3|hky2, not {o}"),
                }
            }
            "--sym" => a.sym = v == "on",
            "--out" => a.out = Some(v),
            "--export-basis" => a.io.export = Some(v),
            "--checkpoint-dir" => a.io.checkpoint_dir = Some(v),
            "--checkpoint-every" => {
                a.io.checkpoint_every_s = 60.0 * v.parse::<f64>().expect("--checkpoint-every")
            }
            "--checkpoint-every-batches" => {
                a.io.checkpoint_every_batches = v.parse().expect("--checkpoint-every-batches")
            }
            "--resume" => {
                a.io.checkpoint_dir = Some(v);
                a.io.resume = true;
            }
            "--max-entries" => a.max_entries = v.parse().expect("--max-entries"),
            "--budget" => a.budget = v.parse().expect("--budget"),
            "--override-a" => a.override_a = Some(v.parse().expect("--override-a")),
            "--override-x" => a.override_x = Some(v.parse().expect("--override-x")),
            "--hky3-verify-closed" => a.verify_closed = v == "on",
            "--hky3-skip" => a.skip = v == "on",
            "--hky3-skip-order" => a.skip_order = v.parse().expect("--hky3-skip-order"),
            "--hky3-skip-witness" => {
                a.skip_witness = match v.as_str() {
                    "off" => 0,
                    "leads" => 1,
                    "span" | "on" => 2,
                    o => o.parse().expect("--hky3-skip-witness off|leads|span"),
                }
            }
            "--hky3-skip-witness-max" => {
                a.skip_witness_max = v.parse().expect("--hky3-skip-witness-max")
            }
            "--hky3-fault" => a.fault = v.parse().expect("--hky3-fault"),
            o => panic!("--instance mode: unknown flag `{o}`"),
        }
        i += 2;
    }
    assert!(!a.instance.is_empty() && a.d > 0, "--instance PATH --D N");
    if a.engine == Engine::Hky2 {
        assert!(
            a.io.checkpoint_dir.is_none() && a.io.export.is_none(),
            "checkpoint/resume and export are implemented for --engine hky3"
        );
    }
    a
}

pub fn run(argv: &[String]) {
    let pa = parse(argv);
    let started = Instant::now();
    let text = std::fs::read_to_string(&pa.instance)
        .unwrap_or_else(|e| panic!("--instance {}: {e}", pa.instance));
    let inst: Value = serde_json::from_str(&text).expect("instance JSON");
    let num = |k: &str| -> u64 {
        inst[k]
            .as_u64()
            .or_else(|| inst[k].as_i64().map(|x| x as u64))
            .or_else(|| inst[k].to_string().parse().ok())
            .unwrap_or_else(|| panic!("instance field `{k}`"))
    };
    let p = num("prime");
    let k = num("k") as usize;
    let a_file = num("a");
    let b_file = num("b");
    let tx: u64 = inst["target"][0].to_string().parse().expect("target[0]");
    let ty: u64 = inst["target"][1].to_string().parse().expect("target[1]");
    assert!(is_prime_u64(p), "p = {p} is not prime");
    let a_c = pa.override_a.unwrap_or(a_file) % p;
    let x_r = pa.override_x.unwrap_or(tx) % p;
    let b_c = b_file % p;
    let control = pa.override_a.is_some() || pa.override_x.is_some();
    if !control {
        // the target is a point of the curve
        let rhs = (mulm(mulm(tx, tx, p), tx, p) + mulm(a_c, tx, p) + b_c) % p;
        assert_eq!(mulm(ty, ty, p), rhs, "target is not on the curve");
    }
    assert!(nonsingular(a_c, b_c, p), "singular curve");
    // ── the tower: Kummer, scale 1 (x = y_0), top constant 1 (X^{2^k} = 1)
    let zeta = subgroup_generator(p, k);
    let tower = Tower {
        kind: Kind::Kummer,
        t: k,
        v: (0..1u64 << k).map(|e| powm(zeta, e, p)).collect(),
        scale: 1,
        steps: vec![Step::Square; k],
        c: 1,
        info: json!({"g": 1, "zeta": zeta, "presentation": "paper: subgroup, c1 = c2 = 1"}),
        cover: None,
    };
    tower.check(p);
    let curve = PrimeCurve { p, a: a_c, b: b_c };
    let target = Target {
        planted: false,
        x_r,
        solution: None,
        xs: Vec::new(),
    };
    let mut rng = StdRng::seed_from_u64(0);
    let (eqs, n, _, _) = system(
        Control::Tower,
        &tower,
        &curve,
        2,
        &target,
        1,
        None,
        &mut rng,
    );
    assert_eq!(
        n,
        2 * k,
        "two blocks of k chain variables, no free variable"
    );
    let ring = HkyRing::new(p, tower_ring(&tower, 2, p).rules, n - 2 * k);
    let gens: Vec<NfPoly> = eqs
        .iter()
        .map(|e| ring.from_raw(&raw_terms(e)))
        .filter(|g| !g.is_empty())
        .collect();
    assert_eq!(gens.len(), 1, "one generator after the chain normal form");
    // ── cross-checks of the generator
    let mine: BTreeMap<u64, u64> = gens[0]
        .iter()
        .map(|(m, f, x)| {
            assert!(f.iter().all(|&e| e == 0));
            (*m, *x % p)
        })
        .filter(|(_, x)| *x != 0)
        .collect();
    let direct = paper_generator(p, a_c, b_c, x_r, k);
    assert_eq!(
        mine, direct,
        "pilot NF(S3) differs from the paper's pkm_generator"
    );
    let mut file_gen_match: Option<bool> = None;
    if let Some(obj) = inst["generator"].as_object() {
        let mut fg = BTreeMap::new();
        for (key, v) in obj {
            let (e1, e2) = key.split_once(',').expect("generator key e1,e2");
            let (e1, e2): (u64, u64) = (e1.parse().unwrap(), e2.parse().unwrap());
            let n2 = 1u64 << k;
            let v: u64 = v.to_string().parse().expect("generator coefficient");
            let e = fg.entry((e1 % n2) | ((e2 % n2) << k)).or_insert(0u64);
            *e = (*e + v) % p;
        }
        fg.retain(|_, v| *v != 0);
        let ok = fg == mine;
        if !control {
            assert!(ok, "the instance file's generator differs from NF(S3)");
        }
        file_gen_match = Some(ok);
    }
    eprintln!(
        "paper instance: p={p} k={k} a={a_c} b={b_c} X_R={x_r} D={} generator terms={} file_generator_match={file_gen_match:?} control={control}",
        pa.d,
        mine.len()
    );
    // ── the closure
    let limits = HkyLimits {
        max_cols: u128::MAX,
        max_entries: pa.max_entries,
        budget: Duration::from_secs(pa.budget),
        ..HkyLimits::default()
    };
    let perm = if pa.sym {
        hky_fp_tower_fast::swap_first_two_blocks(k, 2)
    } else {
        None
    };
    hky_closure_io::install_stop_handlers();
    let (cl, engine_json, io_json) = match pa.engine {
        Engine::Hky3 => {
            let opts = Fast2Options {
                symmetry: perm,
                io: pa.io.clone(),
                verify_closed: pa.verify_closed,
                skip: pa.skip,
                skip_order: pa.skip_order,
                skip_witness: pa.skip_witness,
                skip_witness_max: pa.skip_witness_max,
                fault: pa.fault,
                ..Fast2Options::default()
            };
            let (cl, st) = hky_fp_tower_fast2::closure_fast2(&ring, &gens, pa.d, &limits, &opts);
            let ex = st.io.export.as_ref().map(|e| {
                json!({"path": e.path, "content_sha256": e.content_sha256, "file_sha256": e.file_sha256,
                       "bytes": e.bytes, "n_cols": e.n_cols, "rank": e.rank, "profile": e.profile,
                       "column_order": e.column_order, "contains_one": e.contains_one, "seconds": e.seconds,
                       "error": e.error})
            });
            let io = json!({"checkpoints": st.io.checkpoints, "last_checkpoint_bytes": st.io.last_checkpoint_bytes,
                "checkpoint_s": st.io.checkpoint_s, "resumed_from": st.io.resumed_from, "resumes": st.io.resumes,
                "resumed_at_batch": st.io.resumed_at_batch, "elapsed_total_s": st.io.elapsed_total_s,
                "checkpoint_failures": st.io.checkpoint_failures, "last_checkpoint_error": st.io.last_checkpoint_error,
                "export": ex});
            let o = &opts;
            let ej = json!({"name": "hky_fp_tower_fast2", "symmetry": st.symmetry, "part_cols": st.part_cols,
                "part_rank": st.part_rank, "threads": st.threads, "simd": st.simd, "batches": st.batches,
                "interreductions": st.interreductions, "products": st.products, "zero": st.zero,
                "core_muladds": st.core_muladds, "zero_muladds": st.zero_muladds, "swept_batches": st.swept_batches,
                "verify_products": st.verify_products, "verify_nonzero": st.verify_nonzero, "verify_s": st.verify_s,
                "final_row_bytes": st.final_row_bytes, "sweep_bytes": st.sweep_bytes,
                "skip": st.skip_on, "queued_items": st.queued_items, "skipped_items": st.skipped_items,
                "covered_cols": st.covered_cols, "witness_muladds": st.witness_muladds, "witness_s": st.witness_s, "cover_hist": st.cover_hist, "rank_by_deg": st.rank_by_deg, "witness_seq_s": st.witness_seq_s, "witness_survivors": st.witness_survivors, "witness_rows": st.witness_rows, "coverable_items": st.coverable_items, "coverable_ops": st.coverable_ops,
                "ops_by_src": st.ops_by_src, "prods_by_src": st.prods_by_src, "muladds": st.muladds,
                "peak_words": st.peak_words, "stored_words": st.stored_words,
                "parallel_s": st.parallel_s, "second_s": st.second_s, "interreduce_s": st.interreduce_s,
                "options": {"batch_min": o.batch_min, "batch_max": o.batch_max, "target_pivots": o.target_pivots,
                  "interreduce": o.interreduce, "interreduce_min": o.interreduce_min,
                  "interreduce_block": o.interreduce_block, "rebucket": o.rebucket, "par_second": o.par_second,
                  "bitmap": o.bitmap, "deg_weight": o.deg_weight, "dense": [o.dense_num, o.dense_den],
                  "np_epochs": o.np_epochs, "sweep": o.sweep, "sweep_ratio": o.sweep_ratio,
                  "sweep_tile_bytes": o.sweep_tile_bytes, "sweep_max_words": o.sweep_max_words,
                  "verify_closed": o.verify_closed, "skip": o.skip, "skip_order": o.skip_order, "skip_witness": o.skip_witness, "skip_witness_max": o.skip_witness_max, "fault": o.fault}});
            (cl, ej, io)
        }
        _ => {
            let opts = FastOptions {
                symmetry: perm,
                ..FastOptions::default()
            };
            let (cl, st) = hky_fp_tower_fast::closure_fast(&ring, &gens, pa.d, &limits, &opts);
            let ej = json!({"name": "hky_fp_tower_fast", "symmetry": st.symmetry, "part_cols": st.part_cols,
                "part_rank": st.part_rank, "threads": st.threads, "simd": st.simd, "batches": st.batches});
            (cl, ej, Value::Null)
        }
    };
    let wall = started.elapsed().as_secs_f64();
    let exe = std::env::current_exe().ok();
    let binary_sha256 = exe
        .as_ref()
        .and_then(|p| p.to_str())
        .map(sha256_file)
        .unwrap_or_default();
    let commit = option_env!("HKY_ENGINE_COMMIT").unwrap_or("unknown");
    let inst_sha = sha256_file(&pa.instance);
    let generator: BTreeMap<String, u64> = mine.iter().map(|(m, x)| (format!("{m}"), *x)).collect();
    let rec = json!({
        "schema": "hky_paper_runs.closure/v1",
        "instance": {"path": pa.instance, "sha256": inst_sha, "prime": p, "k": k, "a": a_file, "b": b_file,
                     "target": [tx, ty], "seed": inst["seed"], "unit_nonresidue": inst["unit_nonresidue"]},
        "control": if control { json!({"override_a": pa.override_a, "override_x": pa.override_x}) } else { Value::Null },
        "system": {"p": p, "k": k, "a": a_c, "b": b_c, "x_r": x_r, "zeta": zeta,
                   "generator_masks": generator, "generator_matches_pkm_generator": true,
                   "generator_matches_file": file_gen_match},
        "D": pa.d,
        "columns": cl.cols as u64,
        "rank": cl.dims.last().copied(),
        "profile": cl.dims,
        "contains_one": cl.contains_one,
        "closed": cl.censored.is_none(),
        "censored": cl.censored,
        "rows_reduced": cl.rows_reduced,
        "entries": cl.entries,
        "muladds": cl.muladds,
        "closure_ms_this_process": cl.ms,
        "wall_s_this_process": wall,
        "maxrss_mb": status_kb("VmHWM:").map(|kb| kb / 1024),
        "engine": engine_json,
        "io": io_json,
        "binary_sha256": binary_sha256,
        "engine_commit": commit,
        "rayon_threads": rayon::current_num_threads(),
    });
    // sidecar of an exported basis
    if let Some(ex) = rec["io"]["export"]
        .as_object()
        .filter(|e| e["error"].is_null())
    {
        let path = ex["path"].as_str().unwrap().to_string();
        let side = json!({
            "schema": "hky.basis.sidecar/v1",
            "basis_file": std::path::Path::new(&path).file_name().and_then(|s| s.to_str()),
            "file_sha256": ex["file_sha256"], "content_sha256": ex["content_sha256"], "bytes": ex["bytes"],
            "p": p, "m": 2, "k": k, "D": pa.d, "n_vars": 2 * k,
            "n_cols": ex["n_cols"], "rank": ex["rank"], "contains_one": ex["contains_one"],
            "profile": ex["profile"], "column_order": ex["column_order"], "mask_encoding": 1,
            "instance": {"path": pa.instance, "sha256": rec["instance"]["sha256"], "prime": p, "a": a_c, "b": b_c,
                         "target": [x_r, if control { Value::Null } else { json!(ty) }], "k": k,
                         "control": rec["control"]},
            "engine": {"name": rec["engine"]["name"], "commit": commit, "binary_sha256": binary_sha256,
                       "options": rec["engine"]["options"]},
            "closure": {"muladds": cl.muladds, "rows_reduced": cl.rows_reduced,
                        "seconds": rec["io"]["elapsed_total_s"], "resumed": rec["io"]["resumes"].as_u64().unwrap_or(0) > 0},
            "created_utc_unix_s": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0),
        });
        std::fs::write(format!("{path}.json"), format!("{side:#}\n")).expect("write basis sidecar");
    }
    println!("{rec}");
    if let Some(o) = &pa.out {
        let tmp = format!("{o}.tmp");
        std::fs::write(&tmp, format!("{rec:#}\n")).expect("write --out");
        std::fs::rename(&tmp, o).expect("rename --out");
    }
    eprintln!(
        "paper closure: k={k} D={} cols={} profile={:?} one={} censored={:?} muladds={} wall={wall:.1}s",
        pa.d, cl.cols, cl.dims, cl.contains_one, cl.censored, cl.muladds
    );
    if cl.censored.is_some() {
        std::process::exit(3);
    }
}
