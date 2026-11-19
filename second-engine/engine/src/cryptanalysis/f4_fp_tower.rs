//! The two items of isogeny-labs/crypto `src/cryptanalysis/f4_fp_tower.rs`
//! (commit 7014d5dd2a) that the extracted engine and instance loader use,
//! copied verbatim: the square rule of one tower variable (used by
//! `hky_fp_tower`) and the fields of the tower ring that the loader reads.
//! The F4 engine of that file is not part of this extract.

/// How one tower variable's square rewrites:
/// `y_j² = a·y_j·y_next + b·y_next + c·y_j + d`, all coefficients in
/// `[0, p)`. On a tower's top level `next` is `None` and `a = b = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SquareRule {
    pub next: Option<u8>,
    pub a: u64,
    pub b: u64,
    pub c: u64,
    pub d: u64,
}

/// The tower quotient ring (fields only; the F4 methods are not extracted).
#[derive(Clone, Debug)]
pub struct TowerRing {
    pub p: u64,
    /// One rule per tower variable.
    pub rules: Vec<SquareRule>,
    pub n_free: usize,
}
