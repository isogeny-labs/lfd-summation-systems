//! `hky.basis/v1` reader and writer (spec: research/hky_paper_runs_20261118/BASIS_FORMAT.adopted.md,
//! copied from wt/paper-engine, sha256 05c5d4e5...) and the reduced echelon form the lift needs.

use crate::field::{gemm, Fp};
use crate::layout::{binom, Layout};
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::io::Write;

pub const MAGIC: &[u8; 8] = b"HKYBAS01";

/// A pivoted row echelon form in a degree-non-increasing column order (not necessarily reduced).
pub struct Echelon {
    pub p: u64,
    pub m: u32,
    pub k: u32,
    pub d: u32,
    pub column_order: u32,
    pub flags: u32,
    pub profile: Vec<u64>,
    pub lay: Layout,
    /// pivot column of row i (strictly increasing)
    pub piv: Vec<u32>,
    /// row i's further entries are ents[off[i]..off[i+1]] (column ascending, > pivot, value in [1,p))
    pub off: Vec<u64>,
    pub ents: Vec<(u32, u32)>,
    /// number of stored further entries
    pub nnz: u64,
    /// file-backed rows: (file, byte offset of the entries of row i)
    pub src: Option<(std::fs::File, Vec<u64>)>,
}

pub struct ReadInfo {
    pub file_sha256: String,
    pub content_sha256: String,
    pub bytes: u64,
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

pub fn sha256_hex(b: &[u8]) -> String {
    hex(&Sha256::digest(b))
}


pub fn profile_of(lay: &Layout, piv: &[u32], d: u32) -> Vec<u64> {
    let mut cnt = vec![0u64; d as usize + 1];
    for &c in piv {
        cnt[lay.wt[c as usize] as usize] += 1;
    }
    let mut acc = 0;
    cnt.iter()
        .map(|&x| {
            acc += x;
            acc
        })
        .collect()
}

/// Streaming reader: SHA-256 over the file in one buffered pass, every check of the spec, and only the
/// byte offset of each row kept in memory (rows are re-read on demand with pread).  With
/// `in_memory`, the row entries are also loaded (used by the corruption tool and small tests).
pub fn read_hky_basis_opt(path: &str, in_memory: bool) -> Result<(Echelon, ReadInfo), String> {
    use std::io::{BufReader, Read};
    let file = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let len = file.metadata().map_err(|e| e.to_string())?.len();
    if len < 64 + 32 {
        return Err("file too short".into());
    }
    let body_len = len - 32;
    struct S {
        r: BufReader<std::fs::File>,
        h: Sha256,
        hf: Sha256,
        pos: u64,
        lim: u64,
    }
    impl S {
        fn take(&mut self, n: usize) -> Result<Vec<u8>, String> {
            if self.pos + n as u64 > self.lim {
                return Err("truncated file".into());
            }
            let mut b = vec![0u8; n];
            self.r.read_exact(&mut b).map_err(|e| e.to_string())?;
            self.h.update(&b);
            self.hf.update(&b);
            self.pos += n as u64;
            Ok(b)
        }
        fn u32(&mut self) -> Result<u32, String> {
            Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
        }
        fn u64(&mut self) -> Result<u64, String> {
            Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
        }
    }
    let mut r = S {
        r: BufReader::with_capacity(1 << 22, file.try_clone().map_err(|e| e.to_string())?),
        h: Sha256::new(),
        hf: Sha256::new(),
        pos: 0,
        lim: body_len,
    };
    if r.take(8)? != MAGIC {
        return Err("bad magic".into());
    }
    let version = r.u32()?;
    if version != 1 {
        return Err(format!("version {version} != 1"));
    }
    let m = r.u32()?;
    let p = r.u64()?;
    let k = r.u32()?;
    let d = r.u32()?;
    let n_vars = r.u32()?;
    let enc = r.u32()?;
    let column_order = r.u32()?;
    let flags = r.u32()?;
    let n_cols = r.u64()?;
    let rank = r.u64()?;
    if m != 2 {
        return Err(format!("m = {m}: the lift is for two summands"));
    }
    if n_vars != m * k || enc != 1 {
        return Err("n_vars or mask_encoding".into());
    }
    if !(3..1 << 31).contains(&p) {
        return Err("p out of range".into());
    }
    if flags & !1 != 0 {
        return Err("unknown flag bits".into());
    }
    let want: u64 = (0..=d as u64).map(|w| binom(n_vars as u64, w)).sum();
    if n_cols != want {
        return Err(format!("n_cols {n_cols} != binomial sum {want}"));
    }
    if rank > n_cols {
        return Err("rank > n_cols".into());
    }
    let mut profile = Vec::with_capacity(d as usize + 1);
    for _ in 0..=d {
        profile.push(r.u64()?);
    }
    let mut masks = Vec::with_capacity(n_cols as usize);
    let mb = r.take(8 * n_cols as usize)?;
    for c in mb.chunks(8) {
        masks.push(u64::from_le_bytes(c.try_into().unwrap()));
    }
    drop(mb);
    let lay = Layout::from_masks(k, d, masks)?;
    if column_order > 1 || !lay.is_order(column_order) {
        return Err(format!("column list does not follow column_order {column_order}"));
    }
    let mut piv = Vec::with_capacity(rank as usize);
    let mut boff = Vec::with_capacity(rank as usize);
    let mut off = vec![0u64];
    let mut ents = Vec::new();
    let mut nnz = 0u64;
    let mut prev: i64 = -1;
    for _ in 0..rank {
        let pc = r.u32()?;
        if (pc as i64) <= prev || pc as u64 >= n_cols {
            return Err("pivots not strictly increasing".into());
        }
        prev = pc as i64;
        let t = r.u32()?;
        boff.push(r.pos);
        let eb = r.take(8 * t as usize)?;
        let mut last = pc;
        for e in eb.chunks(8) {
            let c = u32::from_le_bytes(e[..4].try_into().unwrap());
            let v = u32::from_le_bytes(e[4..].try_into().unwrap());
            if c <= last || c as u64 >= n_cols {
                return Err("row entries not strictly ascending after the pivot".into());
            }
            if v == 0 || v as u64 >= p {
                return Err("coefficient not in [1, p)".into());
            }
            last = c;
            if in_memory {
                ents.push((c, v));
            }
        }
        nnz += t as u64;
        piv.push(pc);
        if in_memory {
            off.push(ents.len() as u64);
        }
    }
    if r.pos != body_len {
        return Err("trailing bytes after the rows".into());
    }
    let mut trailer = [0u8; 32];
    r.r.read_exact(&mut trailer).map_err(|e| e.to_string())?;
    r.hf.update(trailer);
    let content = r.h.finalize();
    if content.as_slice() != trailer {
        return Err("trailing SHA-256 does not match the content (corrupted file)".into());
    }
    let pr = profile_of(&lay, &piv, d);
    if pr != profile {
        return Err(format!("profile {profile:?} does not match the pivots {pr:?}"));
    }
    let one = (n_cols - 1) as u32;
    let has_one = piv.last() == Some(&one);
    if has_one != (flags & 1 == 1) {
        return Err("contains_one flag does not match the pivots".into());
    }
    let info = ReadInfo { file_sha256: hex(&r.hf.finalize()), content_sha256: hex(&trailer), bytes: len };
    // row byte offsets (start of the entries) and counts for pread
    let mut fileoff = Vec::with_capacity(boff.len() + 1);
    fileoff.extend_from_slice(&boff);
    let src = if in_memory { None } else { Some((file, fileoff)) };
    if in_memory {
        Ok((Echelon { p, m, k, d, column_order, flags, profile, lay, piv, off, ents, nnz, src }, info))
    } else {
        Ok((Echelon { p, m, k, d, column_order, flags, profile, lay, piv, off: vec![], ents: vec![], nnz, src }, info))
    }
}

pub fn read_hky_basis(path: &str) -> Result<(Echelon, ReadInfo), String> {
    read_hky_basis_opt(path, false)
}

impl Echelon {
    /// further entries of row i (column, coefficient)
    pub fn row_entries(&self, i: usize) -> Vec<(u32, u32)> {
        match &self.src {
            None => self.ents[self.off[i] as usize..self.off[i + 1] as usize].to_vec(),
            Some((f, offs)) => {
                use std::os::unix::fs::FileExt;
                // the u32 count precedes the entries
                let mut cnt = [0u8; 4];
                f.read_exact_at(&mut cnt, offs[i] - 4).expect("pread");
                let t = u32::from_le_bytes(cnt) as usize;
                let mut b = vec![0u8; 8 * t];
                f.read_exact_at(&mut b, offs[i]).expect("pread");
                b.chunks(8)
                    .map(|e| {
                        (u32::from_le_bytes(e[..4].try_into().unwrap()), u32::from_le_bytes(e[4..].try_into().unwrap()))
                    })
                    .collect()
            }
        }
    }
}

/// Write hky.basis/v1 (+ a minimal sidecar).  Returns the file sha256.
pub fn write_hky_basis(path: &str, e: &Echelon, sidecar_extra: serde_json::Value) -> Result<String, String> {
    let mut b: Vec<u8> = Vec::new();
    b.extend_from_slice(MAGIC);
    b.extend_from_slice(&1u32.to_le_bytes());
    b.extend_from_slice(&e.m.to_le_bytes());
    b.extend_from_slice(&e.p.to_le_bytes());
    b.extend_from_slice(&e.k.to_le_bytes());
    b.extend_from_slice(&e.d.to_le_bytes());
    b.extend_from_slice(&(e.m * e.k).to_le_bytes());
    b.extend_from_slice(&1u32.to_le_bytes());
    b.extend_from_slice(&e.column_order.to_le_bytes());
    b.extend_from_slice(&e.flags.to_le_bytes());
    b.extend_from_slice(&(e.lay.ncols() as u64).to_le_bytes());
    b.extend_from_slice(&(e.piv.len() as u64).to_le_bytes());
    for &x in &e.profile {
        b.extend_from_slice(&x.to_le_bytes());
    }
    for &mk in &e.lay.masks {
        b.extend_from_slice(&mk.to_le_bytes());
    }
    for i in 0..e.piv.len() {
        b.extend_from_slice(&e.piv[i].to_le_bytes());
        let s = &e.ents[e.off[i] as usize..e.off[i + 1] as usize];
        b.extend_from_slice(&(s.len() as u32).to_le_bytes());
        for &(c, v) in s {
            b.extend_from_slice(&c.to_le_bytes());
            b.extend_from_slice(&v.to_le_bytes());
        }
    }
    let content = Sha256::digest(&b);
    b.extend_from_slice(&content);
    let fsha = sha256_hex(&b);
    let mut f = std::fs::File::create(path).map_err(|e| e.to_string())?;
    f.write_all(&b).map_err(|e| e.to_string())?;
    let side = serde_json::json!({
        "schema": "hky.basis.sidecar/v1",
        "basis_file": std::path::Path::new(path).file_name().unwrap().to_string_lossy(),
        "file_sha256": fsha,
        "content_sha256": hex(&content),
        "bytes": b.len(),
        "p": e.p, "m": e.m, "k": e.k, "D": e.d, "n_vars": e.m * e.k,
        "n_cols": e.lay.ncols(), "rank": e.piv.len(), "contains_one": e.flags & 1 == 1,
        "profile": e.profile, "column_order": e.column_order, "mask_encoding": 1,
        "writer": "hky_lift (fixture writer, not the engine)",
        "extra": sidecar_extra,
    });
    std::fs::write(format!("{path}.json"), serde_json::to_string_pretty(&side).unwrap()).map_err(|e| e.to_string())?;
    Ok(fsha)
}

pub const NONE: u32 = u32::MAX;

/// Reduced row echelon form: row i = e_{piv[i]} + sum_f b[i][f] e_{free[f]}.
pub struct Rref {
    pub ncols: usize,
    pub piv: Vec<u32>,
    pub free: Vec<u32>,
    /// col -> row index (if pivot) or NONE
    pub prow: Vec<u32>,
    /// col -> index in free (if free) or NONE
    pub fpos: Vec<u32>,
    pub b: Vec<u32>,
}

impl Rref {
    pub fn rank(&self) -> usize {
        self.piv.len()
    }
    pub fn nf(&self) -> usize {
        self.free.len()
    }
    pub fn row(&self, i: usize) -> &[u32] {
        let nf = self.free.len();
        &self.b[i * nf..(i + 1) * nf]
    }
    pub fn new_index(ncols: usize, piv: Vec<u32>) -> Rref {
        let mut prow = vec![NONE; ncols];
        for (i, &c) in piv.iter().enumerate() {
            prow[c as usize] = i as u32;
        }
        let free: Vec<u32> = (0..ncols as u32).filter(|&c| prow[c as usize] == NONE).collect();
        let mut fpos = vec![NONE; ncols];
        for (f, &c) in free.iter().enumerate() {
            fpos[c as usize] = f as u32;
        }
        let b = vec![0u32; piv.len() * free.len()];
        Rref { ncols, piv, free, prow, fpos, b }
    }

    /// Fully reduce a pivoted echelon form by blocked back-substitution (rows processed from the
    /// last pivot to the first; rows of one block in parallel against the finished rows below it,
    /// then sequentially inside the block).  Memory: rank x #free u32.
    pub fn from_echelon(fp: &Fp, e: &Echelon) -> Rref {
        let n = e.lay.ncols();
        let mut r = Rref::new_index(n, e.piv.clone());
        let nf = r.nf();
        let rank = r.rank();
        const BS: usize = 2048;
        let mut hi = rank;
        while hi > 0 {
            let lo = hi.saturating_sub(BS);
            let (done_lo, done_hi) = r.b.split_at_mut(hi * nf);
            let blk = &mut done_lo[lo * nf..];
            let below: &[u32] = done_hi; // rows hi.. (finished), row t at (t - hi) * nf
            let prow = &r.prow;
            let fpos = &r.fpos;
            let p = fp.p;
            let lim = fp.lim;
            blk.par_chunks_mut(nf.max(1)).enumerate().for_each(|(ii, out)| {
                if nf == 0 {
                    return;
                }
                let i = lo + ii;
                let s = e.row_entries(i);
                let mut acc = vec![0u64; nf];
                let mut cnt = 0u64;
                for &(c, v) in &s {
                    let f = fpos[c as usize];
                    if f != NONE {
                        acc[f as usize] += v as u64;
                        continue;
                    }
                    let t = prow[c as usize] as usize;
                    if t >= hi {
                        let nv = p - v as u64;
                        let row = &below[(t - hi) * nf..(t - hi + 1) * nf];
                        for (a, &y) in acc.iter_mut().zip(row) {
                            *a += nv * y as u64;
                        }
                        cnt += 1;
                        if cnt + 1 >= lim {
                            for a in acc.iter_mut() {
                                *a %= p;
                            }
                            cnt = 0;
                        }
                    }
                }
                for (o, a) in out.iter_mut().zip(&acc) {
                    *o = (*a % p) as u32;
                }
            });
            // inside the block, from its last row up
            if nf > 0 {
                for i in (lo..hi).rev() {
                    let s = e.row_entries(i);
                    for &(c, v) in &s {
                        let t = r.prow[c as usize];
                        if t != NONE && (t as usize) < hi && (t as usize) > i {
                            let t = t as usize;
                            let (a, bb) = r.b.split_at_mut(t * nf);
                            let row_t = &bb[..nf];
                            fp.axpy_sub(&mut a[i * nf..(i + 1) * nf], v, row_t);
                        }
                    }
                }
            }
            hi = lo;
        }
        r
    }

    /// Canonical hash of the subspace: depends on the space and the column order only.
    pub fn sha(&self, lay: &Layout) -> String {
        rref_sha(lay, &self.piv, &self.free, &self.b)
    }

    /// residual of a full-width vector (on the free columns); zero iff v is in the span
    pub fn residual(&self, fp: &Fp, v: &[u32]) -> Vec<u32> {
        let nf = self.nf();
        let mut out: Vec<u32> = self.free.iter().map(|&c| v[c as usize]).collect();
        let a: Vec<u32> = self.piv.iter().map(|&c| v[c as usize]).collect();
        gemm(fp, &mut out, &a, &self.b, 1, self.rank(), nf, true);
        out
    }
}

pub fn rref_sha(lay: &Layout, piv: &[u32], free: &[u32], b: &[u32]) -> String {
    let mut h = rref_sha_begin(lay, piv, free);
    rref_sha_rows(&mut h, b);
    sha_finish(h)
}

/// streaming form of [`rref_sha`]: header, then the rows in pivot order with [`rref_sha_rows`]
pub fn rref_sha_begin(lay: &Layout, piv: &[u32], free: &[u32]) -> Sha256 {
    let mut h = Sha256::new();
    h.update(b"hky.lift.rref/v1");
    h.update((lay.ncols() as u64).to_le_bytes());
    h.update((piv.len() as u64).to_le_bytes());
    // column masks in order (the hash is relative to the column order)
    for &m in &lay.masks {
        h.update(m.to_le_bytes());
    }
    for &c in piv {
        h.update(c.to_le_bytes());
    }
    for &c in free {
        h.update(c.to_le_bytes());
    }
    h
}

pub fn rref_sha_rows(h: &mut Sha256, b: &[u32]) {
    for chunk in b.chunks(1 << 16) {
        let mut tmp = Vec::with_capacity(chunk.len() * 4);
        for &x in chunk {
            tmp.extend_from_slice(&x.to_le_bytes());
        }
        h.update(&tmp);
    }
}

pub fn sha_finish(h: Sha256) -> String {
    hex(&h.finalize())
}
