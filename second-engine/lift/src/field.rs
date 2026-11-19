//! Exact arithmetic in F_p, p < 2^31, and the dense modular products the lift is built from.
//!
//! Entries are stored as `u32` in [0, p).  Products accumulate in `u64`; the accumulator is
//! reduced after at most `lim` products, where `lim (p-1)^2 + (p-1) <= u64::MAX`, so no partial
//! sum overflows (for p = 65537, lim is about 4.3e9: in practice one reduction per output entry).

use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

pub static MADDS: AtomicU64 = AtomicU64::new(0);

pub fn madds() -> u64 {
    MADDS.load(Ordering::Relaxed)
}

#[derive(Clone, Copy, Debug)]
pub struct Fp {
    pub p: u64,
    pub lim: u64,
}

impl Fp {
    pub fn new(p: u64) -> Fp {
        assert!(p >= 3 && p < (1u64 << 31), "p must be an odd prime below 2^31");
        let q = (p - 1) * (p - 1);
        let lim = ((u64::MAX - (p - 1)) / q).max(1);
        Fp { p, lim }
    }
    #[inline]
    pub fn mul(&self, a: u32, b: u32) -> u32 {
        ((a as u64 * b as u64) % self.p) as u32
    }
    #[inline]
    pub fn add(&self, a: u32, b: u32) -> u32 {
        let s = a as u64 + b as u64;
        (if s >= self.p { s - self.p } else { s }) as u32
    }
    #[inline]
    pub fn sub(&self, a: u32, b: u32) -> u32 {
        if a >= b {
            a - b
        } else {
            (a as u64 + self.p - b as u64) as u32
        }
    }
    #[inline]
    pub fn neg(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (self.p - a as u64) as u32
        }
    }
    pub fn pow(&self, mut a: u64, mut e: u64) -> u32 {
        let mut r = 1u64;
        a %= self.p;
        while e > 0 {
            if e & 1 == 1 {
                r = r * a % self.p;
            }
            a = a * a % self.p;
            e >>= 1;
        }
        r as u32
    }
    pub fn inv(&self, a: u32) -> u32 {
        assert!(a != 0, "inverse of 0");
        self.pow(a as u64, self.p - 2)
    }
    /// reduce a signed integer into [0, p)
    pub fn from_i64(&self, v: i64) -> u32 {
        v.rem_euclid(self.p as i64) as u32
    }

    /// v <- v - f * row (elementwise, mod p)
    pub fn axpy_sub(&self, v: &mut [u32], f: u32, row: &[u32]) {
        if f == 0 {
            return;
        }
        let nf = (self.p - f as u64) as u64;
        let p = self.p;
        for (x, &r) in v.iter_mut().zip(row) {
            *x = ((*x as u64 + nf * r as u64) % p) as u32;
        }
    }

    pub fn scale(&self, v: &mut [u32], f: u32) {
        let p = self.p;
        for x in v.iter_mut() {
            *x = (*x as u64 * f as u64 % p) as u32;
        }
    }
}

const CB: usize = 512;

static SIMPLE: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| std::env::var("HKY_LIFT_MICRO_KERNEL").is_err());

/// out (n x c) <- out - A (n x m) @ B (m x c)   if `sub`, else out + A @ B.  All mod p.
/// Row-parallel (rayon), column-blocked, zero entries of A skipped.
pub fn gemm(fp: &Fp, out: &mut [u32], a: &[u32], b: &[u32], n: usize, m: usize, c: usize, sub: bool) {
    assert_eq!(out.len(), n * c);
    assert_eq!(a.len(), n * m);
    assert_eq!(b.len(), m * c);
    if n == 0 || c == 0 || m == 0 {
        return;
    }
    let p = fp.p;
    let lim = fp.lim;
    if n < rayon::current_num_threads() && c > CB {
        // few rows: parallel over column blocks instead
        let nb = c.div_ceil(CB);
        let blocks: Vec<(usize, Vec<u32>, u64)> = (0..nb)
            .into_par_iter()
            .map(|bi| {
                let j0 = bi * CB;
                let j1 = (j0 + CB).min(c);
                let w = j1 - j0;
                let mut res = vec![0u32; n * w];
                let mut acc = vec![0u64; w];
                let mut work = 0u64;
                for i in 0..n {
                    for v in acc.iter_mut() {
                        *v = 0;
                    }
                    let mut cnt = 0u64;
                    for l in 0..m {
                        let x = a[i * m + l] as u64;
                        if x == 0 {
                            continue;
                        }
                        work += w as u64;
                        let brow = &b[l * c + j0..l * c + j1];
                        for (s, &y) in acc.iter_mut().zip(brow) {
                            *s += x * y as u64;
                        }
                        cnt += 1;
                        if cnt >= lim {
                            for s in acc.iter_mut() {
                                *s %= p;
                            }
                            cnt = 0;
                        }
                    }
                    for (j, s) in acc.iter().enumerate() {
                        let y = out[i * c + j0 + j] as u64;
                        let r = *s % p;
                        res[i * w + j] = (if sub { (y + p - r) % p } else { (y + r) % p }) as u32;
                    }
                }
                (j0, res, work)
            })
            .collect();
        let mut total = 0;
        for (j0, res, work) in blocks {
            let w = res.len() / n;
            for i in 0..n {
                out[i * c + j0..i * c + j0 + w].copy_from_slice(&res[i * w..(i + 1) * w]);
            }
            total += work;
        }
        MADDS.fetch_add(total, Ordering::Relaxed);
        return;
    }
    // row tiles of RT rows share every streamed block of B (RT-fold less B traffic); column blocks of
    // TB keep the RT x TB u64 accumulators in L1/L2
    const RT: usize = 8;
    const TB: usize = 512;
    let total: u64 = out
        .par_chunks_mut(RT * c)
        .zip(a.par_chunks(RT * m))
        .map(|(otile, atile)| {
            let r = otile.len() / c;
            let mut acc = vec![0u64; RT * TB];
            let mut work = 0u64;
            // inner indices where at least one row of the tile is nonzero
            let ls: Vec<usize> = (0..m).filter(|&l| (0..r).any(|i| atile[i * m + l] != 0)).collect();
            for i in 0..r {
                work += atile[i * m..(i + 1) * m].iter().filter(|&&x| x != 0).count() as u64 * c as u64;
            }
            if ls.is_empty() {
                return 0u64;
            }
            let mut j0 = 0;
            while j0 < c {
                let j1 = (j0 + TB).min(c);
                let w = j1 - j0;
                for v in acc[..r * TB].iter_mut() {
                    *v = 0;
                }
                tile_kernel(&mut acc, atile, m, r, &ls, b, c, j0, w, p, lim);
                for i in 0..r {
                    let o = &mut otile[i * c + j0..i * c + j1];
                    let ac = &acc[i * TB..i * TB + w];
                    if sub {
                        for (y, s) in o.iter_mut().zip(ac) {
                            let rr = *s % p;
                            *y = ((*y as u64 + p - rr) % p) as u32;
                        }
                    } else {
                        for (y, s) in o.iter_mut().zip(ac) {
                            let rr = *s % p;
                            *y = ((*y as u64 + rr) % p) as u32;
                        }
                    }
                }
                j0 = j1;
            }
            work
        })
        .sum();
    MADDS.fetch_add(total, Ordering::Relaxed);
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn tile_body(acc: &mut [u64], atile: &[u32], m: usize, r: usize, ls: &[usize], b: &[u32], c: usize, j0: usize, w: usize, p: u64, lim: u64) {
    const TB: usize = 512;
    if (ls.len() as u64) < lim && !*SIMPLE {
        // register-blocked micro-kernel: 4 rows x 8 columns of u64 accumulators stay in registers
        // across an inner block of LB indices; no reduction is needed (fewer than lim products)
        const LB: usize = 256;
        let w8 = w / 8 * 8;
        for l0 in (0..ls.len()).step_by(LB) {
            let lb = &ls[l0..(l0 + LB).min(ls.len())];
            let mut i0 = 0;
            while i0 < r {
                let q = (r - i0).min(4);
                let mut av = [[0u32; LB]; 4];
                for qq in 0..q {
                    for (t, &l) in lb.iter().enumerate() {
                        av[qq][t] = atile[(i0 + qq) * m + l];
                    }
                }
                let mut jj = 0;
                while jj < w8 {
                    let mut s = [[0u64; 8]; 4];
                    for qq in 0..4 {
                        if qq < q {
                            s[qq].copy_from_slice(&acc[(i0 + qq) * TB + jj..(i0 + qq) * TB + jj + 8]);
                        }
                    }
                    for (t, &l) in lb.iter().enumerate() {
                        let o = l * c + j0 + jj;
                        let bs: &[u32; 8] = b[o..o + 8].try_into().unwrap();
                        let bb: [u32; 8] = *bs;
                        for qq in 0..4 {
                            let x = av[qq][t];
                            for u in 0..8 {
                                s[qq][u] += x as u64 * bb[u] as u64;
                            }
                        }
                    }
                    for qq in 0..q {
                        acc[(i0 + qq) * TB + jj..(i0 + qq) * TB + jj + 8].copy_from_slice(&s[qq]);
                    }
                    jj += 8;
                }
                // remaining columns
                for qq in 0..q {
                    for (t, &l) in lb.iter().enumerate() {
                        let x = av[qq][t];
                        if x == 0 {
                            continue;
                        }
                        for u in w8..w {
                            acc[(i0 + qq) * TB + u] += x as u64 * b[l * c + j0 + u] as u64;
                        }
                    }
                }
                i0 += 4;
            }
        }
        return;
    }
    let mut cnt = 0u64;
    for &l in ls {
        let brow = &b[l * c + j0..l * c + j0 + w];
        for i in 0..r {
            let x = atile[i * m + l];
            if x == 0 {
                continue;
            }
            let ac = &mut acc[i * TB..i * TB + w];
            for (s, &y) in ac.iter_mut().zip(brow) {
                *s += x as u64 * y as u64;
            }
        }
        cnt += 1;
        if cnt >= lim {
            for s in acc[..r * TB].iter_mut() {
                *s %= p;
            }
            cnt = 0;
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn tile_avx2(acc: &mut [u64], atile: &[u32], m: usize, r: usize, ls: &[usize], b: &[u32], c: usize, j0: usize, w: usize, p: u64, lim: u64) {
    tile_body(acc, atile, m, r, ls, b, c, j0, w, p, lim)
}

#[allow(clippy::too_many_arguments)]
fn tile_kernel(acc: &mut [u64], atile: &[u32], m: usize, r: usize, ls: &[usize], b: &[u32], c: usize, j0: usize, w: usize, p: u64, lim: u64) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2 (checked at run time)
            unsafe { tile_avx2(acc, atile, m, r, ls, b, c, j0, w, p, lim) };
            return;
        }
    }
    tile_body(acc, atile, m, r, ls, b, c, j0, w, p, lim)
}

/// Gauss-Jordan: reduced row echelon form of the row space of `cands` (each of length c).
/// Returns rows (sorted by pivot) and pivots.  `rows`/`piv` may be pre-seeded with an RREF set
/// (zero on each other's pivots) that the candidates are added to; seeded rows are updated in place.
pub fn rref_extend(fp: &Fp, rows: &mut Vec<Vec<u32>>, piv: &mut Vec<usize>, cands: Vec<Vec<u32>>, c: usize) {
    const BLK: usize = 128;
    let mut it = cands.into_iter().peekable();
    while it.peek().is_some() {
        let mut blk: Vec<Vec<u32>> = it.by_ref().take(BLK).collect();
        // reduce the block against the current rows (one pass: rows are zero on other pivots)
        if !rows.is_empty() {
            let r = rows.len();
            let mut a = vec![0u32; blk.len() * r];
            for (i, v) in blk.iter().enumerate() {
                for (t, &q) in piv.iter().enumerate() {
                    a[i * r + t] = v[q];
                }
            }
            let mut flat: Vec<u32> = blk.concat();
            let bm: Vec<u32> = rows.concat();
            gemm(fp, &mut flat, &a, &bm, blk.len(), r, c, true);
            for (i, v) in blk.iter_mut().enumerate() {
                v.copy_from_slice(&flat[i * c..(i + 1) * c]);
            }
        }
        // sequential insertion within the block
        let base = rows.len();
        for mut v in blk {
            for t in base..rows.len() {
                let q = piv[t];
                let f = v[q];
                if f != 0 {
                    fp.axpy_sub(&mut v, f, &rows[t]);
                }
            }
            let lead = match v.iter().position(|&x| x != 0) {
                Some(l) => l,
                None => continue,
            };
            let inv = fp.inv(v[lead]);
            fp.scale(&mut v, inv);
            rows.par_iter_mut().for_each(|r| {
                let f = r[lead];
                if f != 0 {
                    fp.axpy_sub(r, f, &v);
                }
            });
            rows.push(v);
            piv.push(lead);
        }
    }
    // sort by pivot
    let mut idx: Vec<usize> = (0..rows.len()).collect();
    idx.sort_by_key(|&i| piv[i]);
    let mut nr = Vec::with_capacity(rows.len());
    let mut np = Vec::with_capacity(rows.len());
    let mut old: Vec<Option<Vec<u32>>> = std::mem::take(rows).into_iter().map(Some).collect();
    for i in idx {
        nr.push(old[i].take().unwrap());
        np.push(piv[i]);
    }
    *rows = nr;
    *piv = np;
}

pub fn rank_of(fp: &Fp, cands: Vec<Vec<u32>>, c: usize) -> usize {
    let mut rows = Vec::new();
    let mut piv = Vec::new();
    rref_extend(fp, &mut rows, &mut piv, cands, c);
    rows.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gemm_matches_naive() {
        for &p in &[65537u64, 2147483647u64, 12289] {
            let fp = Fp::new(p);
            let (n, m, c) = (7, 9, 4100);
            let mut s = 12345u64;
            let mut rnd = || {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                (s % p) as u32
            };
            let a: Vec<u32> = (0..n * m).map(|_| rnd()).collect();
            let b: Vec<u32> = (0..m * c).map(|_| rnd()).collect();
            let o0: Vec<u32> = (0..n * c).map(|_| rnd()).collect();
            let mut o = o0.clone();
            gemm(&fp, &mut o, &a, &b, n, m, c, true);
            for i in 0..n {
                for j in 0..c {
                    let mut acc: i128 = o0[i * c + j] as i128;
                    for l in 0..m {
                        acc -= a[i * m + l] as i128 * b[l * c + j] as i128;
                    }
                    assert_eq!(o[i * c + j] as i128, acc.rem_euclid(p as i128));
                }
            }
        }
    }
}
