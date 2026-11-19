//! The paper's unit-instance recipe (`lfd-summation-systems/scripts/ic_bounds_big_run.py`
//! `cmd_instance` with `--unit-nonresidue`, no `--want`), reimplemented in Rust.
//!
//! The target is drawn exactly as Python's `random.Random(seed)` would draw it
//! (MT19937 seeded by `init_by_array([seed])`, `randrange(1, p)` by
//! rejection on `getrandbits(bit_length(p − 1))`), lifted with the paper's
//! Tonelli–Shanks `sqrt_mod`, and kept at the first draw whose abscissa is a
//! quadratic non-residue. The instance JSON has the schema of the paper's
//! `experiments/ic-bounds-k14/instances/*.json`, plus a `unitness` block that
//! checks Lemma `lem:unit`'s hypotheses explicitly and enumerates the scalar
//! solution set.
//!
//! Usage: `hky_paper_instance --prime P --a A --b 0 --seed S --k K --out PATH [--D 6]`
//!        `hky_paper_instance --selftest`   (reproduces the paper's two unit targets)

use std::collections::BTreeMap;

struct Mt {
    mt: [u32; 624],
    i: usize,
}

impl Mt {
    fn init_genrand(s: u32) -> Mt {
        let mut mt = [0u32; 624];
        mt[0] = s;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt { mt, i: 624 }
    }
    /// CPython `random.seed(n)` for an int n ≥ 0.
    fn python_seed(n: u64) -> Mt {
        let mut key: Vec<u32> = Vec::new();
        let mut x = n;
        while x > 0 {
            key.push(x as u32);
            x >>= 32;
        }
        if key.is_empty() {
            key.push(0);
        }
        let mut m = Mt::init_genrand(19650218);
        let mt = &mut m.mt;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..624.max(key.len()) {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..623 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1566083941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
        m
    }
    fn next_u32(&mut self) -> u32 {
        if self.i >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
                let mut v = self.mt[(k + 397) % 624] ^ (y >> 1);
                if y & 1 == 1 {
                    v ^= 0x9908_b0df;
                }
                self.mt[k] = v;
            }
            self.i = 0;
        }
        let mut y = self.mt[self.i];
        self.i += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
    fn getrandbits(&mut self, k: u32) -> u64 {
        assert!((1..=32).contains(&k));
        u64::from(self.next_u32() >> (32 - k))
    }
    /// `randrange(start, stop)` (CPython ≥ 3.2, `_randbelow_with_getrandbits`).
    fn randrange(&mut self, start: u64, stop: u64) -> u64 {
        let n = stop - start;
        let k = 64 - n.leading_zeros();
        loop {
            let r = self.getrandbits(k);
            if r < n {
                return start + r;
            }
        }
    }
}

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
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

fn legendre(a: u64, p: u64) -> i32 {
    let a = a % p;
    if a == 0 {
        return 0;
    }
    if powm(a, (p - 1) / 2, p) == 1 {
        1
    } else {
        -1
    }
}

/// `csd/field.py sqrt_mod`, step for step (the same root).
fn sqrt_mod(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if legendre(a, p) != 1 {
        return None;
    }
    if p % 4 == 3 {
        return Some(powm(a, (p + 1) / 4, p));
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while legendre(z, p) != -1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, powm(z, q, p), powm(a, q, p), powm(a, (q + 1) / 2, p));
    while t != 1 {
        let (mut i, mut t2) = (0u32, t);
        while t2 != 1 {
            t2 = mulm(t2, t2, p);
            i += 1;
        }
        let b = powm(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mulm(b, b, p);
        t = mulm(mulm(t, b, p), b, p);
        r = mulm(r, b, p);
    }
    Some(r)
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

fn rhs(x: u64, a: u64, b: u64, p: u64) -> u64 {
    (mulm(mulm(x, x, p), x, p) + mulm(a, x, p) + b) % p
}

/// `ic_curve_sweep.random_rational_target` + `PrimeCurve.lift_x`.
fn random_rational_target(rng: &mut Mt, a: u64, b: u64, p: u64) -> (u64, u64) {
    for _ in 0..5000 {
        let x = rng.randrange(1, p);
        let r = rhs(x, a, b, p);
        if r == 0 {
            return (x, 0);
        }
        if legendre(r, p) != 1 {
            continue;
        }
        return (x, sqrt_mod(r, p).unwrap());
    }
    panic!("No rational target found.");
}

/// The recipe: first draw whose abscissa is a non-residue (want = false).
fn draw(p: u64, a: u64, b: u64, seed: u64) -> ((u64, u64), u64) {
    let mut rng = Mt::python_seed(seed);
    for attempt in 0..5000u64 {
        let r = random_rational_target(&mut rng, a, b, p);
        if powm(r.0, (p - 1) / 2, p) != p - 1 {
            continue;
        }
        return (r, attempt + 1);
    }
    panic!("no non-residue target in 5000 attempts");
}

/// `ic_dreg.s3_prime_coeffs`: {(e1, e2): coefficient}.
fn s3_coeffs(p: u64, a: u64, b: u64, x: u64) -> BTreeMap<(u64, u64), u64> {
    let pi = p as i128;
    let (a, b, x) = (a as i128, b as i128, x as i128);
    let xx = x * x % pi;
    let raw: [((u64, u64), i128); 9] = [
        ((2, 0), xx),
        ((0, 2), xx),
        ((1, 1), -2 * xx - 2 * a),
        ((2, 1), -2 * x),
        ((1, 2), -2 * x),
        ((1, 0), -2 * a * x - 4 * b),
        ((0, 1), -2 * a * x - 4 * b),
        ((0, 0), -4 * b * x + a * a),
        ((2, 2), 1),
    ];
    raw.iter()
        .map(|&(e, v)| (e, v.rem_euclid(pi) as u64))
        .filter(|&(_, v)| v != 0)
        .collect()
}

fn subgroup_generator(p: u64, k: u32) -> u64 {
    let n = 1u64 << k;
    for h in 2..p {
        let g = powm(h, (p - 1) / n, p);
        if powm(g, n / 2, p) != 1 {
            return g;
        }
    }
    unreachable!()
}

/// Every (E1, E2) with S_3(g^E1, g^E2, X) = 0, by solving the quadratic in
/// x2 for each x1 of the subgroup and testing membership (x2^{2^k} = 1).
fn solutions(p: u64, a: u64, b: u64, x: u64, k: u32) -> u64 {
    let c = s3_coeffs(p, a, b, x);
    let co = |e1, e2| c.get(&(e1, e2)).copied().unwrap_or(0);
    let n = 1u64 << k;
    let g = subgroup_generator(p, k);
    let inv2 = (p + 1) / 2;
    let mut count = 0u64;
    let mut x1 = 1u64;
    for _ in 0..n {
        let pw = [1, x1, mulm(x1, x1, p)];
        let coef = |e2: u64| (0..3).fold(0u64, |s, e1| (s + mulm(co(e1 as u64, e2), pw[e1], p)) % p);
        let (qa, qb, qc) = (coef(2), coef(1), coef(0));
        let mut roots: Vec<u64> = Vec::new();
        if qa == 0 {
            if qb != 0 {
                roots.push(mulm(p - qc % p, powm(qb, p - 2, p), p) % p);
            } else {
                assert!(qc != 0, "S3 vanishes identically in x2");
            }
        } else {
            let disc = (mulm(qb, qb, p) + p - mulm(4, mulm(qa, qc, p), p)) % p;
            if let Some(s) = sqrt_mod(disc, p) {
                let ia = powm(mulm(2, qa, p), p - 2, p);
                for r in [(p - qb + s) % p, (2 * p - qb - s) % p] {
                    roots.push(mulm(r, ia, p));
                }
                roots.sort_unstable();
                roots.dedup();
            }
        }
        for r in roots {
            if r != 0 && powm(r, n, p) == 1 {
                count += 1;
            }
        }
        let _ = inv2;
        x1 = mulm(x1, g, p);
    }
    count
}

fn binom(n: u64, k: u64) -> u64 {
    (0..k).fold(1u64, |acc, i| acc * (n - i) / (i + 1))
}

fn columns(k: u64, d: u64) -> u64 {
    (0..=d.min(2 * k)).map(|j| binom(2 * k, j)).sum()
}

fn instance_json(p: u64, a: u64, b: u64, seed: u64, k: u32, d_eval: u64) -> String {
    assert!(is_prime(p), "{p} is not prime");
    assert!((p - 1) % (1u64 << k) == 0, "2^k must divide p − 1 (the paper's make_instance assert)");
    assert!((4 * mulm(mulm(a, a, p), a, p) + 27 * mulm(b, b, p)) % p != 0, "singular curve");
    let ((x, y), attempts) = draw(p, a % p, b % p, seed);
    let z = solutions(p, a % p, b % p, x, k);
    // Lemma lem:unit, hypotheses checked one by one
    let two_k1 = (p - 1) % (1u64 << (k + 1)) == 0;
    let x_nonres = legendre(x, p) == -1;
    let on_curve = mulm(y, y, p) == rhs(x, a % p, b % p, p);
    let y_nonzero = y != 0;
    let j1728 = b % p == 0;
    // a primitive 2^{k+1}-th root of unity exists, and every element of H is a square
    let gk = subgroup_generator(p, k);
    let h_in_squares = legendre(gk, p) == 1;
    let unit = two_k1 && x_nonres && on_curve && y_nonzero && j1728 && h_in_squares && z == 0;
    let gen: Vec<String> = s3_coeffs(p, a % p, b % p, x)
        .iter()
        .map(|((e1, e2), v)| format!("  \"{e1},{e2}\": {v}"))
        .collect();
    let src = include_str!("hky_paper_instance.rs");
    let src_hash = {
        use sha2::{Digest, Sha256};
        Sha256::digest(src.as_bytes())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>()
    };
    format!(
        "{{\n \"k\": {k},\n \"prime\": {p},\n \"a\": {a},\n \"b\": {b},\n \"seed\": {seed},\n \"want\": false,\n \"unit_nonresidue\": true,\n \"target\": [\n  {x},\n  {y}\n ],\n \"target_x_is_square\": {},\n \"attempts\": {attempts},\n \"points\": {z},\n \"Z\": [],\n \"columns\": {{\n  \"4\": {},\n  \"5\": {},\n  \"6\": {}\n }},\n \"Z_bound_2_2k\": {},\n \"generator\": {{\n{}\n }},\n \"eval_rank_below_D\": 0,\n \"onto\": true,\n \"eval_rank_D\": {d_eval},\n \"source_sha256\": {{\n  \"generator\": \"{src_hash}\",\n  \"recipe\": \"lfd-summation-systems scripts/ic_bounds_big_run.py cmd_instance --unit-nonresidue (reimplemented in Rust: research/hky_paper_runs_20261118/tools/hky_paper_instance.rs)\"\n }},\n \"unitness\": {{\n  \"lemma\": \"curve-lfd-bounds.tex lem:unit\",\n  \"two_pow_k_plus_1_divides_p_minus_1\": {two_k1},\n  \"subgroup_generator\": {gk},\n  \"subgroup_in_squares\": {h_in_squares},\n  \"x_nonresidue\": {x_nonres},\n  \"target_on_curve\": {on_curve},\n  \"y_nonzero\": {y_nonzero},\n  \"j_1728\": {j1728},\n  \"scalar_solutions_enumerated\": {z},\n  \"unit_ideal\": {unit}\n }}\n}}\n",
        legendre(x, p) == 1,
        columns(k as u64, 4),
        columns(k as u64, 5),
        columns(k as u64, 6),
        2u64 << k,
        gen.join(",\n"),
    )
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    if argv.iter().any(|a| a == "--selftest") {
        // the paper's two unit targets (instances unit-p65537-a1-s10, unit-p163841-a3-s21)
        let (r1, at1) = draw(65537, 1, 0, 10);
        let (r2, at2) = draw(163841, 3, 0, 21);
        println!("p=65537 a=1 seed=10 -> {:?} attempts {at1} (paper: (4271, 11720))", r1);
        println!("p=163841 a=3 seed=21 -> {:?} attempts {at2} (paper: (73730, 28291))", r2);
        assert_eq!(r1, (4271, 11720));
        assert_eq!(r2, (73730, 28291));
        for (p, a, x, k) in [(65537u64, 1u64, 4271u64, 14u32), (163841, 3, 73730, 14)] {
            let z = solutions(p, a, 0, x, k);
            println!("p={p} k={k}: scalar solutions {z} (paper: 0)");
            assert_eq!(z, 0);
        }
        // the enumeration finds solutions where they exist (positive control):
        // a residue abscissa that is a sum of two subgroup points
        let (p, a, k) = (65537u64, 1u64, 6u32);
        let g = subgroup_generator(p, k);
        let mut found = 0;
        'o: for e1 in 0..(1u64 << k) {
            for e2 in 0..(1u64 << k) {
                let (x1, x2) = (powm(g, e1, p), powm(g, e2, p));
                let (r1, r2) = (rhs(x1, a, 0, p), rhs(x2, a, 0, p));
                if x1 == x2 || legendre(r1, p) != 1 || legendre(r2, p) != 1 {
                    continue;
                }
                let (y1, y2) = (sqrt_mod(r1, p).unwrap(), sqrt_mod(r2, p).unwrap());
                let lam = mulm((y2 + p - y1) % p, powm((x2 + p - x1) % p, p - 2, p), p);
                let x3 = (mulm(lam, lam, p) + 2 * p - x1 - x2) % p;
                found = solutions(p, a, 0, x3, k);
                println!("positive control: X = x(P1+P2) = {x3}: {found} scalar solutions (> 0 expected)");
                break 'o;
            }
        }
        assert!(found > 0);
        println!("selftest ok");
        return;
    }
    if argv.iter().any(|a| a == "--planted") {
        // negative-control target: X = x(P1 + P2) with x(P1), x(P2) in the order-2^k subgroup,
        // so the system has scalar solutions and its ideal is not the unit ideal.
        let get = |f: &str| argv.iter().position(|a| a == f).map(|i| argv[i + 1].clone());
        let p: u64 = get("--prime").expect("--prime").parse().unwrap();
        let a: u64 = get("--a").expect("--a").parse().unwrap();
        let k: u32 = get("--k").expect("--k").parse().unwrap();
        let g = subgroup_generator(p, k);
        for e1 in 1..(1u64 << k) {
            for e2 in (e1 + 1)..(1u64 << k) {
                let (x1, x2) = (powm(g, e1, p), powm(g, e2, p));
                let (r1, r2) = (rhs(x1, a, 0, p), rhs(x2, a, 0, p));
                if legendre(r1, p) != 1 || legendre(r2, p) != 1 {
                    continue;
                }
                let (y1, y2) = (sqrt_mod(r1, p).unwrap(), sqrt_mod(r2, p).unwrap());
                let lam = mulm((y2 + p - y1) % p, powm((x2 + p - x1) % p, p - 2, p), p);
                let x3 = (mulm(lam, lam, p) + 2 * p - x1 - x2) % p;
                let y3 = mulm(lam, (x1 + p - x3) % p, p);
                let y3 = (y3 + p - y1) % p;
                println!("{{\"x\": {x3}, \"y\": {y3}, \"e1\": {e1}, \"e2\": {e2}, \"scalar_solutions\": {}}}", solutions(p, a, 0, x3, k));
                return;
            }
        }
    }
    let get = |f: &str| -> Option<String> {
        argv.iter().position(|a| a == f).map(|i| argv[i + 1].clone())
    };
    let p: u64 = get("--prime").expect("--prime").parse().unwrap();
    let a: u64 = get("--a").expect("--a").parse().unwrap();
    let b: u64 = get("--b").unwrap_or("0".into()).parse().unwrap();
    let seed: u64 = get("--seed").expect("--seed").parse().unwrap();
    let k: u32 = get("--k").expect("--k").parse().unwrap();
    let d: u64 = get("--D").unwrap_or("6".into()).parse().unwrap();
    let js = instance_json(p, a, b, seed, k, d);
    match get("--out") {
        Some(o) => std::fs::write(&o, &js).expect("write --out"),
        None => print!("{js}"),
    }
    eprint!("{js}");
}
