//! **Paper-instance runner of the second closure engine** (standalone extract).
//!
//! This file is the `--instance` path of isogeny-labs/crypto
//! `examples/pkm_tower_pilot.rs` at commit 7014d5dd2a, and nothing else: the
//! items that `pkm_tower_paper/paper.rs` (unchanged) uses through
//! `use super::*`, copied verbatim from the pilot, except where a comment
//! marked `EXTRACT:` says what was cut and why (the branches of the pilot's
//! other controls, towers, F4 engines and field backends, which the paper
//! instances never reach). `../PROVENANCE.md` lists every change.
//!
//! Usage: `pkm_tower_pilot --instance PATH --D N [flags of paper.rs]`.

#![recursion_limit = "256"]
// EXTRACT: the items below carry the variants and fields the pilot's other
// modes use; the instance path does not read them all.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::fs::File;
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::f4_fp_tower::{SquareRule, TowerRing};
use crypto_lib::cryptanalysis::hky_fp_tower::{HkyLimits, HkyRing};
use crypto_lib::cryptanalysis::hky_fp_tower_fast::{self, FastOptions};
use crypto_lib::cryptanalysis::hky_fp_tower_fast2::{self, Fast2Options};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};

#[path = "pkm_tower_paper/paper.rs"]
mod paper;

/// EXTRACT: `ic_boundary::PrimeCurve` of the crypto crate, fields only (the
/// loader builds it and `summation` reads `p`, `a`, `b`).
#[derive(Clone, Debug)]
struct PrimeCurve {
    p: u64,
    a: u64,
    b: u64,
}

/// EXTRACT: the pilot's `--cover-map` type. The paper instances never set a
/// cover (`cover: None`), so the type is kept only as the field's type and
/// has no values.
#[derive(Clone, Copy, Debug)]
enum CoverMap {}

// ── Field arithmetic (pilot lines 180–212) ─────────────────────────

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn addm(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}

fn subm(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn powm(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
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

fn invm(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "inverse of zero");
    powm(a, p - 2, p)
}

// ── Polynomials while building a system ────────────────────────────


/// A polynomial over `F_p` in `n` variables, as a map from exponent
/// vectors to non-zero coefficients.
#[derive(Clone, Debug)]
struct Pol {
    n: usize,
    terms: BTreeMap<Vec<u32>, u64>,
}

impl Pol {
    fn zero(n: usize) -> Self {
        Pol {
            n,
            terms: BTreeMap::new(),
        }
    }

    fn constant(n: usize, c: u64, p: u64) -> Self {
        let mut f = Pol::zero(n);
        if !c.is_multiple_of(p) {
            f.terms.insert(vec![0; n], c % p);
        }
        f
    }

    fn var(n: usize, i: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        let mut f = Pol::zero(n);
        f.terms.insert(e, 1);
        f
    }

    fn add(&self, o: &Pol, p: u64) -> Pol {
        let mut f = self.clone();
        for (e, c) in &o.terms {
            let v = f.terms.entry(e.clone()).or_insert(0);
            *v = addm(*v, *c, p);
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn scale(&self, k: u64, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e, c) in &self.terms {
            let v = mulm(*c, k, p);
            if v != 0 {
                f.terms.insert(e.clone(), v);
            }
        }
        f
    }

    fn sub(&self, o: &Pol, p: u64) -> Pol {
        self.add(&o.scale(p - 1, p), p)
    }

    fn mul(&self, o: &Pol, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e1, c1) in &self.terms {
            for (e2, c2) in &o.terms {
                let e: Vec<u32> = e1.iter().zip(e2).map(|(a, b)| a + b).collect();
                let v = f.terms.entry(e).or_insert(0);
                *v = addm(*v, mulm(*c1, *c2, p), p);
            }
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn eval(&self, x: &[u64], p: u64) -> u64 {
        let mut acc = 0u64;
        for (e, c) in &self.terms {
            let mut t = *c;
            for (i, k) in e.iter().enumerate() {
                if *k > 0 {
                    t = mulm(t, powm(x[i], *k as u64, p), p);
                }
            }
            acc = addm(acc, t, p);
        }
        acc
    }

    // EXTRACT: `degree` and `to_poly` (F4 input) are not used by the instance path.
}

/// Semaev's third summation polynomial for `y² = x³ + ax + b`, on three
/// polynomial arguments:
/// `(X1 − X2)²X3² − 2((X1 + X2)(X1X2 + a) + 2b)X3 + (X1X2 − a)² − 4b(X1 + X2)`.
fn s3(x1: &Pol, x2: &Pol, x3: &Pol, a: u64, b: u64, p: u64) -> Pol {
    let n = x1.n;
    let d = x1.sub(x2, p);
    let t1 = d.mul(&d, p).mul(x3, p).mul(x3, p);
    let s = x1.add(x2, p);
    let pr = x1.mul(x2, p);
    let inner = s
        .mul(&pr.add(&Pol::constant(n, a, p), p), p)
        .add(&Pol::constant(n, mulm(2, b, p), p), p);
    let t2 = inner.mul(x3, p).scale(p - 2, p);
    let q = pr.sub(&Pol::constant(n, a, p), p);
    let t3 = q.mul(&q, p).sub(&s.scale(mulm(4, b, p), p), p);
    t1.add(&t2, p).add(&t3, p)
}

// ── Towers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Kummer,
    Dickson,
    Isogeny,
    /// The Dickson tower in its torus case (`2^t | p + 1`). Last, so the
    /// cell seeds of the other kinds (`kind as u64`) are unchanged.
    DicksonTorus,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Kummer => "kummer",
            Kind::Dickson => "dickson",
            Kind::Isogeny => "isogeny",
            Kind::DicksonTorus => "dickson-torus",
        }
    }
}

/// One level of a tower: how `y_{j+1}` (or the final constant) follows
/// from `y_j`.
#[derive(Clone, Copy, Debug)]
enum Step {
    /// `y_{j+1} = y_j²` (Kummer).
    Square,
    /// `y_{j+1} = y_j² − 2` (Dickson `D_2`).
    Dickson,
    /// `y_{j+1} = y_j + τ/(y_j − ξ)`, the x-map of a 2-isogeny with
    /// kernel `(ξ, 0)`.
    Velu { xi: u64, tau: u64 },
}

struct Tower {
    kind: Kind,
    t: usize,
    /// The factor-base abscissa set, `2^t` distinct elements of `F_p`.
    v: Vec<u64>,
    /// `x = scale · y_0` (Kummer: the coset representative `g`).
    scale: u64,
    steps: Vec<Step>,
    /// The value the last step must reach.
    c: u64,
    info: Value,
    /// `--cover-map`: the factor base is `x∘φ(V)`, not `V` (see [`CoverMap`]).
    /// `None` for every tower unless the flag is given.
    cover: Option<CoverMap>,
}

impl Tower {
    fn apply(&self, j: usize, y: u64, p: u64) -> u64 {
        match self.steps[j] {
            Step::Square => mulm(y, y, p),
            Step::Dickson => subm(mulm(y, y, p), 2, p),
            Step::Velu { xi, tau } => addm(y, mulm(tau, invm(subm(y, xi, p), p), p), p),
        }
    }

    /// `y_0 … y_{t−1}` for an abscissa of `V`.
    fn values(&self, x: u64, p: u64) -> Vec<u64> {
        let mut y = vec![0u64; self.t];
        y[0] = mulm(x, invm(self.scale, p), p);
        for j in 0..self.t - 1 {
            y[j + 1] = self.apply(j, y[j], p);
        }
        y
    }

    /// The `t` equations of one block whose variables start at `base`.
    fn equations(&self, n: usize, base: usize, p: u64) -> Vec<Pol> {
        (0..self.t)
            .map(|j| {
                let yj = Pol::var(n, base + j);
                let next = if j + 1 < self.t {
                    Pol::var(n, base + j + 1)
                } else {
                    Pol::constant(n, self.c, p)
                };
                let sq = yj.mul(&yj, p);
                match self.steps[j] {
                    Step::Square => next.sub(&sq, p),
                    Step::Dickson => next.sub(&sq.sub(&Pol::constant(n, 2, p), p), p),
                    Step::Velu { xi, tau } => {
                        let lhs = next.mul(&yj.sub(&Pol::constant(n, xi, p), p), p);
                        let rhs = sq
                            .sub(&yj.scale(xi, p), p)
                            .add(&Pol::constant(n, tau, p), p);
                        lhs.sub(&rhs, p)
                    }
                }
            })
            .collect()
    }

    /// The tower equations of one block (variables from `base`) as the
    /// square rules of [`f4_fp_tower`]: each equation solved for `y_j²`.
    fn square_rules(&self, base: usize, p: u64) -> Vec<SquareRule> {
        (0..self.t)
            .map(|j| {
                let top = j + 1 == self.t;
                let next = (!top).then_some((base + j + 1) as u8);
                match (self.steps[j], top) {
                    // y_{j+1} = y_j²; on top, y_j² = c
                    (Step::Square, false) => SquareRule {
                        next,
                        a: 0,
                        b: 1,
                        c: 0,
                        d: 0,
                    },
                    (Step::Square, true) => SquareRule {
                        next,
                        a: 0,
                        b: 0,
                        c: 0,
                        d: self.c % p,
                    },
                    // y_{j+1} = y_j² − 2; on top, y_j² = c + 2
                    (Step::Dickson, false) => SquareRule {
                        next,
                        a: 0,
                        b: 1,
                        c: 0,
                        d: 2 % p,
                    },
                    (Step::Dickson, true) => SquareRule {
                        next,
                        a: 0,
                        b: 0,
                        c: 0,
                        d: addm(self.c, 2, p),
                    },
                    // y_{j+1}(y_j − ξ) = y_j² − ξ·y_j + τ, so
                    // y_j² = y_j·y_{j+1} − ξ·y_{j+1} + ξ·y_j − τ
                    (Step::Velu { xi, tau }, false) => SquareRule {
                        next,
                        a: 1,
                        b: subm(0, xi, p),
                        c: xi,
                        d: subm(0, tau, p),
                    },
                    // c·(y_j − ξ) = y_j² − ξ·y_j + τ, so
                    // y_j² = (c + ξ)·y_j − (c·ξ + τ)
                    (Step::Velu { xi, tau }, true) => SquareRule {
                        next,
                        a: 0,
                        b: 0,
                        c: addm(self.c, xi, p),
                        d: subm(0, addm(mulm(self.c, xi, p), tau, p), p),
                    },
                }
            })
            .collect()
    }

    /// Every element of `V` satisfies the tower, `|V| = 2^t` and the
    /// elements are distinct.  Panics otherwise: a wrong tower would make
    /// every degree below meaningless.
    fn check(&self, p: u64) {
        let mut sorted = self.v.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 1 << self.t, "{:?}: V has repeats", self.kind);
        let n = self.t;
        let eqs = self.equations(n, 0, p);
        for &x in &self.v {
            let y = self.values(x, p);
            for (j, e) in eqs.iter().enumerate() {
                assert_eq!(
                    e.eval(&y, p),
                    0,
                    "{:?}: level {j} fails at x = {x}",
                    self.kind
                );
            }
        }
    }
}


fn nonsingular(a: u64, b: u64, p: u64) -> bool {
    let d = addm(
        mulm(4, mulm(mulm(a, a, p), a, p), p),
        mulm(27, mulm(b, b, p), p),
        p,
    );
    d != 0
}

// ── Targets ────────────────────────────────────────────────────────

/// A target and, when planted, the full solution of its system.
struct Target {
    planted: bool,
    x_r: u64,
    /// Values of every variable at the planted solution.
    solution: Option<Vec<u64>>,
    /// The planted abscissae, for the naive control.
    xs: Vec<u64>,
}

// ── Systems ────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Control {
    Tower,
    Ladder,
    Naive,
    Null,
    Multiples,
    /// A random factor base of the same size: two unknowns `x_1, x_2`,
    /// membership `∏_{f ∈ F}(x_i − f)`, `F` random abscissae on the curve
    /// (with the planted ones), `|F|` that of the cell's base (m = 2 only).
    RandBase,
}

impl Control {
    fn name(self) -> &'static str {
        match self {
            Control::Tower => "tower",
            Control::Ladder => "ladder",
            Control::Naive => "naive",
            Control::Null => "null",
            Control::Multiples => "multiples",
            Control::RandBase => "randbase",
        }
    }
}

/// The summation part: `S₃(x_1, x_2, x_R)` at `m = 2`; for `m ≥ 3` the
/// chain `S₃(x_1, x_2, u_1)`, `S₃(u_k, x_{k+2}, u_{k+1})`,
/// `S₃(u_{m−2}, x_m, x_R)`, whose free unknowns `u_k` are the abscissae of
/// the partial sums.
fn summation(tower: &Tower, curve: &PrimeCurve, m: usize, n: usize, x_r: u64) -> Vec<Pol> {
    let p = curve.p;
    let (a, b) = (curve.a, curve.b);
    let x = |i: usize| Pol::var(n, i * tower.t).scale(tower.scale, p);
    let xr = Pol::constant(n, x_r, p);
    // EXTRACT: the pilot's `--cover-map` branch (S₃ at x∘φ(u_i)) is cut: the
    // paper instances have no cover map (`cover: None`).
    if m == 2 {
        return vec![s3(&x(0), &x(1), &xr, a, b, p)];
    }
    let u = |k: usize| Pol::var(n, m * tower.t + k);
    let mut eqs = vec![s3(&x(0), &x(1), &u(0), a, b, p)];
    for k in 1..m - 2 {
        eqs.push(s3(&u(k - 1), &x(k + 1), &u(k), a, b, p));
    }
    eqs.push(s3(&u(m - 3), &x(m - 1), &xr, a, b, p));
    eqs
}

fn system(
    control: Control,
    tower: &Tower,
    curve: &PrimeCurve,
    m: usize,
    target: &Target,
    g: usize,
    null_width: Option<usize>,
    rng: &mut StdRng,
) -> (Vec<Pol>, usize, Option<Vec<u64>>, std::ops::Range<usize>) {
    let p = curve.p;
    let t = tower.t;
    // EXTRACT: the pilot's `naive` and `randbase` controls (early returns),
    // the `null` coefficient draws, the `multiples` equations and the
    // `ladder` quadrics are cut; the paper instances run the `tower`
    // control, for which the pilot's code is the lines below, unchanged
    // except that the `multiples` term of `n` is 0 and `sol` is the target's.
    assert_eq!(control, Control::Tower, "the extract runs the tower control only");
    let _ = (g, null_width, &mut *rng);
    let n = m * t + (m - 2);
    let mut eqs = summation(tower, curve, m, n, target.x_r);
    let towers_from = eqs.len();
    for i in 0..m {
        eqs.extend(tower.equations(n, i * t, p));
    }
    let towers = towers_from..eqs.len();
    let sol = target.solution.clone();
    (eqs, n, sol, towers)
}

/// The tower quotient ring of an `m`-summand system: `m` blocks of the
/// tower's rules, then the chain's `m − 2` free unknowns.
fn tower_ring(tower: &Tower, m: usize, p: u64) -> TowerRing {
    let mut rules = Vec::with_capacity(m * tower.t);
    for i in 0..m {
        rules.extend(tower.square_rules(i * tower.t, p));
    }
    TowerRing {
        p,
        rules,
        n_free: m - 2,
    }
}

fn raw_terms(e: &Pol) -> Vec<(Vec<u32>, u64)> {
    e.terms.iter().map(|(k, v)| (k.clone(), *v)).collect()
}


#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Engine {
    F4,
    Tower,
    Sig,
    /// `f4_tower_generic`, over the backend `--backend` picks.
    Generic,
    /// `hky_fp_tower`: the Huang–Kosters–Yeo closures `V_{F,c}` for
    /// `c = --hky-cmin … --hky-cmax`, and the first and last fall degree
    /// read off them (`research/prime_fall_degree_20261027/`). No F4 runs;
    /// the naive control runs too (free variables only).
    Hky,
    /// `hky_fp_tower_fast`: the same closures on the faster instrument
    /// (`research/hky_engine_speed_20261117/`), `--hky2-level`,
    /// `--hky2-sym`, `--hky2-batch`, `--hky2-order`.
    Hky2,
    /// `hky_fp_tower_fast2`: the second speed pass (`--hky3-*` options;
    /// `--hky2-sym` applies to it too).
    Hky3,
}


/// Deterministic Miller–Rabin for `n < 2⁶⁴` (the first twelve prime bases).
fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for q in BASES {
        if n.is_multiple_of(q) {
            return n == q;
        }
    }
    let mulm = |a: u64, b: u64| ((u128::from(a) * u128::from(b)) % u128::from(n)) as u64;
    let powm = |mut a: u64, mut e: u64| {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = mulm(r, a);
            }
            a = mulm(a, a);
            e >>= 1;
        }
        r
    };
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'bases: for a in BASES {
        let mut x = powm(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// `proth:K`: the smallest prime `c·2^K + 1` with `c` odd (so

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    // EXTRACT: the pilot's `main` dispatches `--instance` to `paper::run`
    // (unchanged here); its other modes are not part of the extract.
    if argv.iter().any(|x| x == "--instance") {
        paper::run(&argv);
        return;
    }
    eprintln!("usage: pkm_tower_pilot --instance PATH --D N [flags; see examples/pkm_tower_paper/paper.rs]");
    std::process::exit(2);
}
