//! Bounded independent F4 reading of the frozen split-chain text input.
//! Build as an example of isogeny-labs/crypto. The Python/heap closure source
//! is not imported. A positive reading still needs an independent certificate.
use std::env;
use std::fs;
use std::collections::HashSet;
use std::sync::Arc;
use std::time::Duration;

// The copied F4 module has an optional callback type for a signature engine.
// This adapter never sets that callback; the type completes its interface.
mod sig_fp_tower { pub struct SigStep; }
mod f4_fp_tower;
use f4_fp_tower::{
    f4_tower, F4Trace, SquareRule, StepTrace, TowerF4Options, TowerRing, TraceStep,
};
use serde_json::json;
use sha2::{Digest, Sha256};

fn ints(line: &str) -> Vec<u64> {
    line.split_whitespace().map(|s| s.parse().expect("integer")).collect()
}

fn step(number: usize, s: &StepTrace) {
    eprintln!("step={} degree={} pairs={} cols={} fresh={} elapsed_ms={:.0}",
              number, s.degree, s.critical_pairs+s.tower_pairs, s.cols, s.fresh, s.ms);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert!((4..=5).contains(&args.len()) || (args.len()==6 && args[4]=="--replay"),
            "usage: lfd_split_f4_20261124 INPUT SECONDS MAX_DENSE_ENTRIES [TRACE_OUT | --replay TRACE_IN]");
    let path = &args[1];
    let seconds: u64 = args[2].parse().expect("seconds");
    let max_dense: u64 = args[3].parse().expect("max_dense");
    let data = fs::read_to_string(path).expect("input");
    let input_sha256 = format!("{:x}", Sha256::digest(data.as_bytes()));
    let lines: Vec<&str> = data.lines().collect();
    let header = ints(lines[0]);
    assert!(header.len() >= 5);
    let (p, m, k, aux, bound) = (header[0], header[1] as usize,
                                 header[2] as usize, header[3] as usize,
                                 header[4] as u32);
    assert!(p > 3 && p < (1 << 24) && m >= 3 && k >= 2 && aux == m - 2);
    let n = m * k;
    assert!(n + 6 * aux < 64 && bound < 63);
    assert_eq!(ints(lines[1]), vec![(m - 1) as u64]);
    assert_eq!(lines.len(), m + 1);
    let mut rules = Vec::with_capacity(n);
    for block in 0..m {
        for level in 0..k {
            rules.push(if level + 1 < k {
                SquareRule { next: Some((block*k+level+1) as u8),
                             a: 0, b: 1, c: 0, d: 0 }
            } else {
                SquareRule { next: None, a: 0, b: 0, c: 0, d: 1 }
            });
        }
    }
    let ring = TowerRing::new(p, rules, aux);
    let mut input = Vec::new();
    for line in &lines[2..] {
        let vals = ints(line);
        assert_eq!(vals.len(), 1 + 2 * vals[0] as usize);
        let terms: Vec<(Vec<u32>, u64)> = vals[1..].chunks_exact(2).map(|pair| {
            let (packed, coef) = (pair[0], pair[1]);
            assert!(coef > 0 && coef < p);
            let mask = packed & ((1u64 << n) - 1);
            let mut e = vec![0u32; n + aux];
            for (j, ej) in e.iter_mut().take(n).enumerate() {
                *ej = ((mask >> j) & 1) as u32;
            }
            let mut rebuilt = mask;
            for j in 0..aux {
                e[n+j] = ((packed >> (n+6*j)) & 63) as u32;
                rebuilt += (e[n+j] as u64) << (n+6*j);
            }
            assert_eq!(packed, rebuilt);
            (e, coef)
        }).collect();
        let poly = ring.from_raw(&terms);
        assert!(!poly.is_empty());
        input.push(poly);
    }
    assert_eq!(input.len(), m-1);
    let mut opts = TowerF4Options::new(bound)
        .with_budget(Duration::from_secs(seconds))
        .with_max_dense(max_dense)
        .on_step(step);
    let replaying = args.len() == 6;
    let mut replayed_trace_sha256 = None;
    if replaying {
        let bytes = fs::read(&args[5]).expect("trace input");
        replayed_trace_sha256 = Some(format!("{:x}", Sha256::digest(&bytes)));
        let saved: serde_json::Value = serde_json::from_slice(&bytes).expect("trace JSON");
        assert_eq!(saved["m"].as_u64(), Some(m as u64));
        assert_eq!(saved["k"].as_u64(), Some(k as u64));
        assert_eq!(saved["inconsistent"].as_bool(), Some(true));
        let mut trace = F4Trace::default();
        trace.steps = saved["steps"].as_array().expect("steps").iter().map(|s| {
            let kept: HashSet<(u128,u32)> = s["kept"].as_array().expect("kept").iter().map(|pair| {
                let a = pair.as_array().expect("pair");
                (a[0].as_str().expect("mono").parse().expect("mono integer"),
                 a[1].as_u64().expect("element") as u32)
            }).collect();
            let fresh_lms = s["fresh_lms"].as_array().expect("fresh_lms").iter()
                .map(|v| v.as_str().expect("lead").parse().expect("lead integer")).collect();
            TraceStep { degree: s["degree"].as_u64().expect("degree") as u32,
                        kept, fresh_lms }
        }).collect();
        opts.replay = Some(Arc::new(trace));
    } else {
        opts.record_trace = true;
    }
    let rep = f4_tower(&input, &ring, &opts);
    let mut trace_sha256 = None;
    if let (Some(dest), Some(trace)) = (if replaying {None} else {args.get(4)}, rep.learned.as_ref()) {
        let steps: Vec<_> = trace.steps.iter().map(|s| {
            let mut kept: Vec<_> = s.kept.iter().map(|&(mono, element)|
                (mono.to_string(), element)).collect();
            kept.sort();
            json!({"degree": s.degree, "kept": kept,
                   "fresh_lms": s.fresh_lms.iter().map(|m| m.to_string()).collect::<Vec<_>>()})
        }).collect();
        let bytes = serde_json::to_vec(&json!({
            "input_sha256": input_sha256, "p": p, "m": m, "k": k,
            "inconsistent": rep.inconsistent, "solving_degree_max": rep.solving_degree_max,
            "steps": steps
        })).expect("trace JSON");
        fs::write(dest, &bytes).expect("trace output");
        trace_sha256 = Some(format!("{:x}", Sha256::digest(&bytes)));
    }
    println!("{}", json!({
        "input": path, "input_sha256": input_sha256,
        "p": p, "m": m, "k": k, "aux": aux, "degree_cap": bound,
        "seconds_cap": seconds, "max_dense_entries": max_dense,
        "inconsistent": rep.inconsistent,
        "degree_reached": rep.degree_reached,
        "solving_degree_max": rep.solving_degree_max,
        "last_productive_degree": rep.last_productive_degree,
        "steps": rep.steps, "max_rows": rep.max_rows, "max_cols": rep.max_cols,
        "max_dense_entries_seen": rep.max_dense_entries,
        "timed_out": rep.timed_out, "oversize": rep.oversize,
        "pairs_above_bound": rep.pairs_above_bound,
        "learned_steps": rep.learned.as_ref().map(|t| t.steps.len()),
        "learned_rows": rep.learned.as_ref().map(|t| t.kept_rows()),
        "trace_sha256": trace_sha256,
        "replay_mode": replaying, "replayed_trace_sha256": replayed_trace_sha256,
        "replay_failure": rep.replay_failure,
        "elapsed_ms": rep.ms,
    }));
}
