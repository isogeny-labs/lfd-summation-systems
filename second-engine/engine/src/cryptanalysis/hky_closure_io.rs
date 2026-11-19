//! # Checkpoint/resume and basis export for `hky_fp_tower_fast2` (`hky3`).
//!
//! Kept apart from the reduction internals: `closure_fast2` calls this module
//! at three delimited hooks (resume after setup, a periodic checkpoint at a
//! batch boundary, export at the end) and hands it one struct,
//! [`ClosureState`], which holds every piece of state the main loop carries
//! between batches. Nothing here changes the schedule: a resumed closure
//! performs the same batches in the same order as an uninterrupted one, so
//! its echelon, dims and exported basis are identical (the checkpoint test in
//! `research/hky_paper_runs_20261118/README.md`).
//!
//! * Checkpoint: `<dir>/state.ckpt`, written as `state.ckpt.tmp`, fsynced,
//!   renamed (atomic on POSIX), directory fsynced; trailing SHA-256 of the
//!   content; a fingerprint of the ring, generators, degree, column split and
//!   every schedule option, checked on resume. `<dir>/state.json` is a small
//!   human/bash-readable summary of the last checkpoint.
//! * Export: `hky.basis/v1`, specified in
//!   `research/hky_paper_runs_20261118/BASIS_FORMAT.md`.
//! * Stop: SIGTERM/SIGINT/SIGUSR1 (after [`install_stop_handlers`]) make the
//!   closure write a checkpoint at the next batch boundary and return a
//!   closure censored `"signal"` (the pod memory guard uses this).

use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use sha2::{Digest, Sha256};

use super::hky_fp_tower::{HkyRing, NfPoly};
use super::hky_fp_tower_fast::{Row, RowS, Space, NONE};
use super::hky_fp_tower_fast2::{ClosureState, Echelon, Fast2Options, Fast2Stats, Item};

/// What the io hooks do. All off by default.
#[derive(Clone, Debug, Default)]
pub struct IoOptions {
    /// Directory of `state.ckpt`. Without it nothing is checkpointed.
    pub checkpoint_dir: Option<String>,
    /// Write a checkpoint when this many seconds passed since the last one
    /// (0: no time trigger).
    pub checkpoint_every_s: f64,
    /// Test trigger: every this many batches (0: off).
    pub checkpoint_every_batches: u64,
    /// Continue from `<checkpoint_dir>/state.ckpt` when it exists.
    pub resume: bool,
    /// Write the closed basis (`hky.basis/v1`) here at the end.
    pub export: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct ExportReport {
    pub path: String,
    pub content_sha256: String,
    pub file_sha256: String,
    pub bytes: u64,
    pub n_cols: u64,
    pub rank: u64,
    pub profile: Vec<u64>,
    pub column_order: u32,
    pub contains_one: bool,
    pub seconds: f64,
    /// The write error, if the export failed (no file then).
    pub error: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct IoReport {
    pub checkpoints: u64,
    pub last_checkpoint_bytes: u64,
    pub checkpoint_s: f64,
    pub resumed_from: Option<String>,
    pub resumes: u64,
    pub resumed_at_batch: u64,
    pub export: Option<ExportReport>,
    pub elapsed_total_s: f64,
    /// Checkpoints that could not be written (disk full, I/O error).
    pub checkpoint_failures: u64,
    pub last_checkpoint_error: Option<String>,
}

// ── stop requests ──────────────────────────────────────────────────────

static STOP: AtomicBool = AtomicBool::new(false);

pub fn request_stop() {
    STOP.store(true, Ordering::SeqCst);
}

pub fn stop_requested() -> bool {
    STOP.load(Ordering::Relaxed)
}

extern "C" fn on_signal(_: libc::c_int) {
    STOP.store(true, Ordering::SeqCst);
}

/// SIGTERM, SIGINT and SIGUSR1 request a checkpoint and a clean stop.
pub fn install_stop_handlers() {
    let h = on_signal as extern "C" fn(libc::c_int) as libc::sighandler_t;
    unsafe {
        libc::signal(libc::SIGTERM, h);
        libc::signal(libc::SIGINT, h);
        libc::signal(libc::SIGUSR1, h);
    }
}

// ── the checkpoint clock ───────────────────────────────────────────────

pub struct Clock {
    on: bool,
    every_s: f64,
    every_b: u64,
    last: Instant,
    last_b: u64,
}

impl Clock {
    pub fn new(io: &IoOptions) -> Self {
        Clock {
            on: io.checkpoint_dir.is_some()
                && (io.checkpoint_every_s > 0.0 || io.checkpoint_every_batches > 0),
            every_s: io.checkpoint_every_s,
            every_b: io.checkpoint_every_batches,
            last: Instant::now(),
            last_b: 0,
        }
    }
    pub fn due(&mut self, batches: u64) -> bool {
        if !self.on {
            return false;
        }
        if self.every_b > 0 && batches >= self.last_b + self.every_b {
            self.last_b = batches;
            return true;
        }
        self.every_s > 0.0 && self.last.elapsed().as_secs_f64() >= self.every_s
    }
    pub fn reset(&mut self) {
        self.last = Instant::now();
    }
}

// ── fingerprint ────────────────────────────────────────────────────────

/// SHA-256 of everything that determines the closure's schedule and
/// result: ring, generators, degree, part sizes, schedule options.
pub fn fingerprint(
    ring: &HkyRing,
    gens: &[NfPoly],
    c: u32,
    o: &Fast2Options,
    part_cols: &[u64],
) -> [u8; 32] {
    let mut h = Sha256::new();
    h.update(b"hky3-closure-state/v2");
    h.update(ring.p.to_le_bytes());
    h.update((ring.n_free as u64).to_le_bytes());
    for r in &ring.rules {
        h.update([r.next.unwrap_or(255u8)]);
        for x in [r.a, r.b, r.c, r.d] {
            h.update(x.to_le_bytes());
        }
    }
    h.update((gens.len() as u64).to_le_bytes());
    for g in gens {
        h.update((g.len() as u64).to_le_bytes());
        for (m, f, x) in g {
            h.update(m.to_le_bytes());
            h.update(f);
            h.update(x.to_le_bytes());
        }
    }
    h.update(c.to_le_bytes());
    for &pc in part_cols {
        h.update(pc.to_le_bytes());
    }
    for x in [
        o.batch_min,
        o.batch_max,
        o.target_pivots,
        o.interreduce_min,
        o.interreduce_block,
        o.dense_num,
        o.dense_den,
        o.deg_weight,
    ] {
        h.update((x as u64).to_le_bytes());
    }
    h.update(o.interreduce.to_bits().to_le_bytes());
    h.update([
        o.rebucket as u8,
        o.par_second as u8,
        o.bitmap as u8,
        o.fault,
    ]);
    // engine-opt: the row storage (epoch coordinates) and the sweep schedule
    h.update([o.np_epochs as u8, o.sweep as u8]);
    h.update(o.sweep_ratio.to_le_bytes());
    h.update((o.sweep_tile_bytes as u64).to_le_bytes());
    h.update((o.sweep_max_words as u64).to_le_bytes());
    // engine-skip: another schedule (nothing added when off, so earlier
    // checkpoints keep their fingerprint)
    if o.skip {
        h.update(b"skip/v1");
        h.update([o.skip_order, o.skip_witness]);
        h.update(o.skip_witness_max.to_le_bytes());
    }
    match &o.symmetry {
        None => h.update([0u8]),
        Some(p) => {
            h.update([1u8]);
            for &x in p {
                h.update((x as u64).to_le_bytes());
            }
        }
    }
    h.finalize().into()
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

// ── hashing streams ────────────────────────────────────────────────────

struct HW<W: Write> {
    w: W,
    h: Sha256,
    n: u64,
    /// The first write error; later writes are skipped (disk full must not
    /// kill a multi-day closure: the caller discards the file and goes on).
    err: Option<std::io::Error>,
}

impl<W: Write> HW<W> {
    fn new(w: W) -> Self {
        HW {
            w,
            h: Sha256::new(),
            n: 0,
            err: None,
        }
    }
    fn put(&mut self, b: &[u8]) {
        self.h.update(b);
        self.n += b.len() as u64;
        if self.err.is_none() {
            if let Err(e) = self.w.write_all(b) {
                self.err = Some(e);
            }
        }
    }
    fn u8(&mut self, x: u8) {
        self.put(&[x]);
    }
    fn u16(&mut self, x: u16) {
        self.put(&x.to_le_bytes());
    }
    fn u32(&mut self, x: u32) {
        self.put(&x.to_le_bytes());
    }
    fn u64(&mut self, x: u64) {
        self.put(&x.to_le_bytes());
    }
    fn f64(&mut self, x: f64) {
        self.u64(x.to_bits());
    }
    fn pairs(&mut self, r: &[(u32, u32)]) {
        let mut buf = Vec::with_capacity(r.len() * 8);
        for &(c, x) in r {
            buf.extend_from_slice(&c.to_le_bytes());
            buf.extend_from_slice(&x.to_le_bytes());
        }
        self.put(&buf);
    }
    fn u32s(&mut self, v: &[u32]) {
        let mut buf = Vec::with_capacity(v.len() * 4);
        for &x in v {
            buf.extend_from_slice(&x.to_le_bytes());
        }
        self.put(&buf);
    }
    /// Appends the digest (not hashed) and returns it.
    fn finish(mut self) -> std::io::Result<([u8; 32], W, u64)> {
        let d: [u8; 32] = self.h.clone().finalize().into();
        if let Some(e) = self.err.take() {
            return Err(e);
        }
        self.w.write_all(&d)?;
        let n = self.n + 32;
        Ok((d, self.w, n))
    }
}

/// Flush, fsync and rename a finished stream; any error removes `tmp`.
fn seal(
    res: std::io::Result<([u8; 32], BufWriter<File>, u64)>,
    tmp: &Path,
    fin: &Path,
) -> std::io::Result<([u8; 32], u64)> {
    let out = (|| {
        let (d, bw, n) = res?;
        let f = bw.into_inner().map_err(|e| e.into_error())?;
        f.sync_all()?;
        drop(f);
        fs::rename(tmp, fin)?;
        Ok((d, n))
    })();
    if out.is_err() {
        let _ = fs::remove_file(tmp);
    }
    out
}

struct HR<R: Read> {
    r: R,
    h: Sha256,
}

impl<R: Read> HR<R> {
    fn get(&mut self, n: usize) -> Vec<u8> {
        let mut b = vec![0u8; n];
        self.r
            .read_exact(&mut b)
            .expect("hky_closure_io: truncated checkpoint");
        self.h.update(&b);
        b
    }
    fn u8(&mut self) -> u8 {
        self.get(1)[0]
    }
    fn u16(&mut self) -> u16 {
        u16::from_le_bytes(self.get(2).try_into().unwrap())
    }
    fn u32(&mut self) -> u32 {
        u32::from_le_bytes(self.get(4).try_into().unwrap())
    }
    fn u64(&mut self) -> u64 {
        u64::from_le_bytes(self.get(8).try_into().unwrap())
    }
    fn f64(&mut self) -> f64 {
        f64::from_bits(self.u64())
    }
    fn pairs(&mut self, n: usize) -> Row {
        let b = self.get(n * 8);
        b.chunks_exact(8)
            .map(|c| {
                (
                    u32::from_le_bytes(c[0..4].try_into().unwrap()),
                    u32::from_le_bytes(c[4..8].try_into().unwrap()),
                )
            })
            .collect()
    }
    fn u32s(&mut self, n: usize) -> Vec<u32> {
        let b = self.get(n * 4);
        b.chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect()
    }
}

fn fsync_dir(dir: &Path) {
    if let Ok(d) = File::open(dir) {
        let _ = d.sync_all();
    }
}

// ── checkpoint ─────────────────────────────────────────────────────────

const CKPT_MAGIC: &[u8; 8] = b"HKYCKP01";

fn stats_fields(s: &Fast2Stats) -> ([f64; 4], [u64; 14]) {
    (
        [s.parallel_s, s.second_s, s.serial_s, s.interreduce_s],
        [
            s.batches,
            s.rounds,
            s.interreductions,
            s.muladds,
            s.interreduce_muladds,
            s.second_muladds,
            s.products,
            s.zero,
            s.reducers,
            s.dense_rows,
            s.sparse_rows,
            s.stored_words,
            s.peak_words,
            s.io.checkpoints,
        ],
    )
}

/// Bytes needed for `x` (1..=4).
fn width(x: u32) -> usize {
    (((32 - x.leading_zeros()) as usize).div_ceil(8)).max(1)
}

fn pack(buf: &mut Vec<u8>, x: u32, w: usize) {
    buf.extend_from_slice(&x.to_le_bytes()[..w]);
}

/// A row in the checkpoint. engine-skip round: rows are written packed
/// (tags 2 and 3: each column and coefficient in the fewest bytes its row
/// needs; 3 + 3 bytes per sparse term and 3 per dense word for p < 2^24 and
/// fewer than 2^24 columns, against 8 and 4); tags 0 and 1 (unpacked, the
/// earlier format) are still read.
fn write_row<W: Write>(w: &mut HW<W>, r: &RowS) {
    match r {
        RowS::Sparse(r) => {
            let cw = width(r.iter().map(|t| t.0).max().unwrap_or(0));
            let vw = width(r.iter().map(|t| t.1).max().unwrap_or(0));
            w.u8(2);
            w.u32(r.len() as u32);
            w.u8((cw | vw << 4) as u8);
            let mut buf = Vec::with_capacity(r.len() * (cw + vw));
            for &(c, x) in r {
                pack(&mut buf, c, cw);
                pack(&mut buf, x, vw);
            }
            w.put(&buf);
        }
        RowS::Dense { lead, v, nnz } => {
            let vw = width(v.iter().copied().max().unwrap_or(0));
            w.u8(3);
            w.u32(*lead);
            w.u32(v.len() as u32);
            w.u32(*nnz);
            w.u8(vw as u8);
            let mut buf = Vec::with_capacity(v.len() * vw);
            for &x in v {
                pack(&mut buf, x, vw);
            }
            w.put(&buf);
        }
    }
}

fn unpack(b: &[u8]) -> u32 {
    let mut x = [0u8; 4];
    x[..b.len()].copy_from_slice(b);
    u32::from_le_bytes(x)
}

fn read_row<R: Read>(r: &mut HR<R>) -> RowS {
    match r.u8() {
        0 => {
            let n = r.u32() as usize;
            RowS::Sparse(r.pairs(n))
        }
        1 => {
            let lead = r.u32();
            let n = r.u32() as usize;
            let nnz = r.u32();
            RowS::Dense {
                lead,
                v: r.u32s(n),
                nnz,
            }
        }
        2 => {
            let n = r.u32() as usize;
            let wd = r.u8() as usize;
            let (cw, vw) = (wd & 15, wd >> 4);
            assert!(
                (1..=4).contains(&cw) && (1..=4).contains(&vw),
                "hky_closure_io: bad widths"
            );
            let b = r.get(n * (cw + vw));
            RowS::Sparse(
                b.chunks_exact(cw + vw)
                    .map(|c| (unpack(&c[..cw]), unpack(&c[cw..])))
                    .collect(),
            )
        }
        3 => {
            let lead = r.u32();
            let n = r.u32() as usize;
            let nnz = r.u32();
            let vw = r.u8() as usize;
            assert!((1..=4).contains(&vw), "hky_closure_io: bad width");
            let b = r.get(n * vw);
            RowS::Dense {
                lead,
                v: b.chunks_exact(vw).map(unpack).collect(),
                nnz,
            }
        }
        t => panic!("hky_closure_io: bad row tag {t}"),
    }
}

fn lead_of(r: &RowS) -> u32 {
    match r {
        RowS::Sparse(r) => r[0].0,
        RowS::Dense { lead, .. } => *lead,
    }
}

/// Writes `st` atomically to `<dir>/state.ckpt` (no-op without a dir).
pub(super) fn checkpoint_now(io: &IoOptions, fp: &[u8; 32], st: &ClosureState, rep: &mut IoReport) {
    let Some(dir) = io.checkpoint_dir.as_deref() else {
        return;
    };
    let t0 = Instant::now();
    let dir = PathBuf::from(dir);
    fs::create_dir_all(&dir).expect("hky_closure_io: create checkpoint dir");
    let tmp = dir.join("state.ckpt.tmp");
    let fin = dir.join("state.ckpt");
    let f = match File::create(&tmp) {
        Ok(f) => f,
        Err(e) => {
            rep.checkpoint_failures += 1;
            rep.last_checkpoint_error = Some(e.to_string());
            eprintln!("hky3 checkpoint FAILED to create {} ({e})", tmp.display());
            return;
        }
    };
    let mut w = HW::new(BufWriter::with_capacity(8 << 20, f));
    w.put(CKPT_MAGIC);
    w.put(fp);
    w.u64(st.rows_reduced);
    w.u8(st.contains_one as u8);
    w.u64(st.yield_q16);
    w.u64(st.fresh_since_ir as u64);
    w.u64(st.rank_at_ir as u64);
    w.f64(st.elapsed_before_s);
    let (fl, us) = stats_fields(&st.stats);
    for x in fl {
        w.f64(x);
    }
    for (i, x) in us.into_iter().enumerate() {
        // the checkpoint being written counts itself
        w.u64(if i == 13 { x + 1 } else { x });
    }
    w.u64(st.stats.io.resumes);
    if st.stats.skip_on {
        // engine-skip: counters and the cover table (only in checkpoints of
        // skip runs; the fingerprint tells a reader whether they are here)
        w.u64(st.stats.queued_items);
        w.u64(st.stats.skipped_items);
        w.u32(st.cover.len() as u32);
        for cv in &st.cover {
            w.u64(cv.len() as u64);
            let wide: Vec<u32> = cv.iter().map(|&x| u32::from(x)).collect();
            w.u32s(&wide);
        }
    }
    w.u32(st.ech.len() as u32);
    for e in &st.ech {
        w.u64(e.pivot.len() as u64);
        w.u64(e.rows.len() as u64);
        for r in &e.rows {
            write_row(&mut w, r);
        }
        // the epoch columns: stored tails are in these coordinates
        w.u64(e.np_cols.len() as u64);
        w.u32s(&e.np_cols);
    }
    w.u32(st.queue.len() as u32);
    for b in &st.queue {
        w.u64(b.len() as u64);
        for it in b {
            match it {
                Item::Row(q, r) => {
                    w.u8(0);
                    w.u8(*q);
                    w.u32(r.len() as u32);
                    w.pairs(r);
                }
                Item::Prod { part, row, var } => {
                    w.u8(1);
                    w.u8(*part);
                    w.u32(*row);
                    w.u16(*var);
                }
            }
        }
    }
    let (digest, bytes) = match seal(w.finish(), &tmp, &fin) {
        Ok(x) => x,
        Err(e) => {
            // a failed checkpoint (disk full, I/O error) is reported, not fatal:
            // the previous checkpoint stays, the closure goes on
            rep.checkpoint_failures += 1;
            rep.last_checkpoint_error = Some(e.to_string());
            eprintln!("hky3 checkpoint FAILED ({e}); tmp removed, previous checkpoint kept");
            return;
        }
    };
    fsync_dir(&dir);
    let rank: usize = st.ech.iter().map(|e| e.rows.len()).sum();
    let queued: usize = st.queue.iter().map(|b| b.len()).sum();
    let secs = t0.elapsed().as_secs_f64();
    rep.checkpoints = st.stats.io.checkpoints + 1;
    rep.last_checkpoint_bytes = bytes;
    rep.checkpoint_s += secs;
    let meta = format!(
        "{{\"schema\":\"hky3.checkpoint.meta/v1\",\"sha256\":\"{}\",\"bytes\":{bytes},\"batches\":{},\"rank\":{rank},\"queued\":{queued},\"muladds\":{},\"elapsed_s\":{:.1},\"write_s\":{secs:.2},\"fingerprint\":\"{}\",\"unix_s\":{}}}\n",
        hex(&digest),
        st.stats.batches,
        st.stats.muladds,
        st.elapsed_before_s,
        hex(fp),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0)
    );
    let mt = dir.join("state.json.tmp");
    fs::write(&mt, meta.as_bytes()).expect("hky_closure_io: write state.json");
    fs::rename(&mt, dir.join("state.json")).expect("hky_closure_io: rename state.json");
    fsync_dir(&dir);
    eprintln!(
        "hky3 checkpoint: {} bytes, batch {}, rank {rank}, queued {queued}, {secs:.1} s",
        bytes, st.stats.batches
    );
}

/// The checkpoint of `io`, if resuming is asked for and one exists.
/// Panics on a corrupt checkpoint or one of another closure (fingerprint).
pub(super) fn try_resume(
    io: &IoOptions,
    fp: &[u8; 32],
    space: &Space,
    skip: bool,
) -> Option<ClosureState> {
    if !io.resume {
        return None;
    }
    let dir = io.checkpoint_dir.as_deref()?;
    let path = Path::new(dir).join("state.ckpt");
    if !path.exists() {
        eprintln!(
            "hky3 resume: no checkpoint at {}, starting fresh",
            path.display()
        );
        return None;
    }
    let t0 = Instant::now();
    let f = File::open(&path).expect("hky_closure_io: open checkpoint");
    let mut r = HR {
        r: BufReader::with_capacity(8 << 20, f),
        h: Sha256::new(),
    };
    assert_eq!(
        &r.get(8)[..],
        CKPT_MAGIC,
        "hky_closure_io: not a checkpoint"
    );
    let got = r.get(32);
    assert_eq!(
        &got[..],
        fp,
        "hky_closure_io: the checkpoint is of another closure (fingerprint {} != {})",
        hex(&got),
        hex(fp)
    );
    let rows_reduced = r.u64();
    let contains_one = r.u8() != 0;
    let yield_q16 = r.u64();
    let fresh_since_ir = r.u64() as usize;
    let rank_at_ir = r.u64() as usize;
    let elapsed_before_s = r.f64();
    let fl: Vec<f64> = (0..4).map(|_| r.f64()).collect();
    let us: Vec<u64> = (0..14).map(|_| r.u64()).collect();
    let resumes = r.u64();
    let mut stats = Fast2Stats {
        parallel_s: fl[0],
        second_s: fl[1],
        serial_s: fl[2],
        interreduce_s: fl[3],
        batches: us[0],
        rounds: us[1],
        interreductions: us[2],
        muladds: us[3],
        interreduce_muladds: us[4],
        second_muladds: us[5],
        products: us[6],
        zero: us[7],
        reducers: us[8],
        dense_rows: us[9],
        sparse_rows: us[10],
        stored_words: us[11],
        peak_words: us[12],
        ..Fast2Stats::default()
    };
    stats.io.checkpoints = us[13];
    stats.io.resumes = resumes;
    stats.io.resumed_at_batch = us[0];
    let mut cover: Vec<Vec<u16>> = Vec::new();
    if skip {
        stats.skip_on = true;
        stats.queued_items = r.u64();
        stats.skipped_items = r.u64();
        let nc = r.u32() as usize;
        for q in 0..nc {
            let n = r.u64() as usize;
            assert_eq!(
                n,
                space.parts[q].keys.len(),
                "hky_closure_io: cover of part {q}"
            );
            cover.push(
                r.u32s(n)
                    .into_iter()
                    .map(|x| u16::try_from(x).expect("hky_closure_io: cover entry"))
                    .collect(),
            );
        }
    }
    let np = r.u32() as usize;
    assert_eq!(np, space.parts.len(), "hky_closure_io: part count");
    let mut ech = Vec::with_capacity(np);
    for q in 0..np {
        let n = r.u64() as usize;
        assert_eq!(
            n,
            space.parts[q].keys.len(),
            "hky_closure_io: part {q} columns"
        );
        let nr = r.u64() as usize;
        let mut rows = Vec::with_capacity(nr);
        let mut pivot = vec![NONE; n];
        for i in 0..nr {
            let row = read_row(&mut r);
            let l = lead_of(&row) as usize;
            assert!(l < n && pivot[l] == NONE, "hky_closure_io: bad pivot");
            pivot[l] = i as u32;
            rows.push(row);
        }
        let nnp = r.u64() as usize;
        let np_cols = r.u32s(nnp);
        assert!(
            np_cols.windows(2).all(|w| w[0] < w[1])
                && np_cols.last().is_none_or(|&c| (c as usize) < n),
            "hky_closure_io: bad epoch columns"
        );
        ech.push(Echelon::with_epoch(pivot, rows, np_cols));
    }
    let nb = r.u32() as usize;
    let mut queue = Vec::with_capacity(nb);
    for _ in 0..nb {
        let n = r.u64() as usize;
        let mut b = Vec::with_capacity(n);
        for _ in 0..n {
            b.push(match r.u8() {
                0 => {
                    let q = r.u8();
                    let len = r.u32() as usize;
                    Item::Row(q, r.pairs(len))
                }
                1 => {
                    let part = r.u8();
                    let row = r.u32();
                    let var = r.u16();
                    Item::Prod { part, row, var }
                }
                t => panic!("hky_closure_io: bad item tag {t}"),
            });
        }
        queue.push(b);
    }
    let digest: [u8; 32] = r.h.clone().finalize().into();
    let mut trailer = [0u8; 32];
    r.r.read_exact(&mut trailer)
        .expect("hky_closure_io: checkpoint without trailer");
    assert_eq!(
        digest, trailer,
        "hky_closure_io: checkpoint SHA-256 mismatch"
    );
    let mut extra = [0u8; 1];
    assert_eq!(
        r.r.read(&mut extra).unwrap_or(0),
        0,
        "hky_closure_io: trailing bytes after the checkpoint"
    );
    let rank: usize = ech.iter().map(|e| e.rows.len()).sum();
    eprintln!(
        "hky3 resume: {} batch {} rank {rank} sha256 {} ({:.1} s)",
        path.display(),
        stats.batches,
        hex(&digest),
        t0.elapsed().as_secs_f64()
    );
    Some(ClosureState {
        ech,
        queue,
        rows_reduced,
        contains_one,
        yield_q16,
        fresh_since_ir,
        rank_at_ir,
        stats,
        elapsed_before_s,
        cover,
    })
}

// ── basis export (hky.basis/v1) ────────────────────────────────────────

const BASIS_MAGIC: &[u8; 8] = b"HKYBAS01";
const MASK_BITS: u64 = (1u64 << 56) - 1;

fn key_mask(k: u128) -> u64 {
    (k as u64) & MASK_BITS
}

fn sigma(mask: u64, perm: &[usize]) -> u64 {
    let mut out = 0u64;
    let mut b = mask;
    while b != 0 {
        let j = b.trailing_zeros() as usize;
        b &= b - 1;
        out |= 1u64 << perm[j];
    }
    out
}

/// Sorted merge `a + s·b` of two sorted sparse vectors, times `scale`.
fn merge(a: &[(u32, u64)], b: &[(u32, u64)], s_neg: bool, scale: u64, p: u64) -> Vec<(u32, u64)> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    let bv = |x: u64| if s_neg { (p - x) % p } else { x };
    while i < a.len() || j < b.len() {
        let (c, x) = if j >= b.len() || (i < a.len() && a[i].0 < b[j].0) {
            i += 1;
            (a[i - 1].0, a[i - 1].1)
        } else if i >= a.len() || b[j].0 < a[i].0 {
            j += 1;
            (b[j - 1].0, bv(b[j - 1].1))
        } else {
            i += 1;
            j += 1;
            (a[i - 1].0, (a[i - 1].1 + bv(b[j - 1].1)) % p)
        };
        let x = x * scale % p;
        if x != 0 {
            out.push((c, x));
        }
    }
    out
}

/// Writes the closed basis of `ech` as `hky.basis/v1` (BASIS_FORMAT.md).
pub(super) fn export_basis(
    path: &str,
    ring: &HkyRing,
    space: &Space,
    ech: &[Echelon],
    perm: Option<&[usize]>,
    contains_one: bool,
    c: u32,
) -> ExportReport {
    let t0 = Instant::now();
    let p = ring.p;
    assert_eq!(
        ring.n_free, 0,
        "hky.basis/v1: free variables are not supported"
    );
    let n_vars = ring.n_tower() as u32;
    assert!(
        n_vars.is_multiple_of(2),
        "hky.basis/v1: m = 2 blocks expected"
    );
    let k = n_vars / 2;
    let split = space.split();
    let keys0 = &space.parts[0].keys;
    // columns in file order; pos0[i] = file column of part-0 column i
    let mut cols: Vec<u64> = Vec::with_capacity(keys0.len() * 2);
    let mut pos0: Vec<u32> = Vec::with_capacity(keys0.len());
    let mut paired: Vec<bool> = Vec::with_capacity(keys0.len());
    // part-1 column -> part-0 column, and its inverse
    let mut m1_to_0: Vec<u32> = Vec::new();
    let mut m0_to_1: Vec<u32> = vec![NONE; keys0.len()];
    let perm_v: Vec<usize> = if split {
        perm.expect("split space without a permutation").to_vec()
    } else {
        Vec::new()
    };
    for (i, &key) in keys0.iter().enumerate() {
        let r = key_mask(key);
        assert_eq!(
            u128::from(r) | (u128::from(r.count_ones()) << 120),
            key,
            "hky.basis/v1: key layout"
        );
        pos0.push(cols.len() as u32);
        cols.push(r);
        if split {
            let s = sigma(r, &perm_v);
            assert!(s <= r, "representative is the larger mask");
            if s != r {
                cols.push(s);
                paired.push(true);
            } else {
                paired.push(false);
            }
        } else {
            paired.push(false);
            let _ = i;
        }
    }
    if split {
        let keys1 = &space.parts[1].keys;
        let mut j = 0usize;
        for (i, &key) in keys0.iter().enumerate() {
            if j < keys1.len() && keys1[j] == key {
                m1_to_0.push(i as u32);
                m0_to_1[i] = j as u32;
                j += 1;
            }
        }
        assert_eq!(
            j,
            keys1.len(),
            "hky.basis/v1: part-1 columns are not a subsequence"
        );
    }
    let n_cols = cols.len();
    let deg_of = |fc: usize| cols[fc].count_ones() as usize;
    // pivots and profile (first pass)
    let mut counts = vec![0u64; c as usize + 1];
    let mut rank = 0u64;
    if contains_one {
        for &m in &cols {
            counts[m.count_ones() as usize] += 1;
        }
        rank = n_cols as u64;
    } else {
        for i in 0..keys0.len() {
            let has0 = ech[0].pivot[i] != NONE;
            let has1 = split && m0_to_1[i] != NONE && ech[1].pivot[m0_to_1[i] as usize] != NONE;
            let n = has0 as u64 + has1 as u64;
            counts[deg_of(pos0[i] as usize)] += n;
            rank += n;
        }
    }
    let mut profile = counts.clone();
    for d in 1..profile.len() {
        profile[d] += profile[d - 1];
    }
    // write
    let tmp = format!("{path}.tmp");
    let f = File::create(&tmp).unwrap_or_else(|e| panic!("hky.basis/v1: create {tmp}: {e}"));
    let mut w = HW::new(BufWriter::with_capacity(8 << 20, f));
    w.put(BASIS_MAGIC);
    w.u32(1);
    w.u32(2);
    w.u64(p);
    w.u32(k);
    w.u32(c);
    w.u32(n_vars);
    w.u32(1);
    w.u32(if split { 1 } else { 0 });
    w.u32(contains_one as u32);
    w.u64(n_cols as u64);
    w.u64(rank);
    for &x in &profile {
        w.u64(x);
    }
    {
        let mut buf = Vec::with_capacity(n_cols * 8);
        for &m in &cols {
            buf.extend_from_slice(&m.to_le_bytes());
        }
        w.put(&buf);
    }
    let mut written = 0u64;
    let mut emit = |w: &mut HW<BufWriter<File>>, piv: u32, rest: &[(u32, u64)]| {
        w.u32(piv);
        w.u32(rest.len() as u32);
        let mut buf = Vec::with_capacity(rest.len() * 8);
        let mut last = piv;
        for &(cc, x) in rest {
            assert!(cc > last && x != 0 && x < p, "hky.basis/v1: row order");
            last = cc;
            buf.extend_from_slice(&cc.to_le_bytes());
            buf.extend_from_slice(&(x as u32).to_le_bytes());
        }
        w.put(&buf);
        written += 1;
    };
    if contains_one {
        for fc in 0..n_cols {
            emit(&mut w, fc as u32, &[]);
        }
    } else {
        let inv2 = p.div_ceil(2);
        // rows are read through `Echelon::row_full`: stored tails are in
        // epoch coordinates (engine-opt)
        let lift = |q: usize, ri: usize| -> Vec<(u32, u64)> {
            let full = ech[q].row_full(ri);
            let mut v = Vec::with_capacity(full.len() * 2);
            for &(col, x) in &full {
                let x = u64::from(x);
                let i = if q == 0 {
                    col as usize
                } else {
                    m1_to_0[col as usize] as usize
                };
                let at = pos0[i];
                v.push((at, x));
                if paired[i] {
                    v.push((at + 1, if q == 0 { x } else { (p - x) % p }));
                }
            }
            v
        };
        for i in 0..keys0.len() {
            let r0 = (ech[0].pivot[i] != NONE).then(|| ech[0].pivot[i] as usize);
            let r1 = (split && m0_to_1[i] != NONE && ech[1].pivot[m0_to_1[i] as usize] != NONE)
                .then(|| ech[1].pivot[m0_to_1[i] as usize] as usize);
            let at = pos0[i];
            let one = |v: &[(u32, u64)], piv: u32| {
                assert!(
                    !v.is_empty() && v[0] == (piv, 1),
                    "hky.basis/v1: row is not monic at its pivot"
                );
            };
            match (r0, r1) {
                (None, None) => {}
                (Some(a), None) => {
                    let v = lift(0, a);
                    one(&v, at);
                    emit(&mut w, at, &v[1..]);
                }
                (None, Some(b)) => {
                    let v = lift(1, b);
                    one(&v, at);
                    emit(&mut w, at, &v[1..]);
                }
                (Some(a), Some(b)) => {
                    let (va, vb) = (lift(0, a), lift(1, b));
                    let s = merge(&va, &vb, false, inv2, p);
                    let d = merge(&va, &vb, true, inv2, p);
                    one(&s, at);
                    one(&d, at + 1);
                    emit(&mut w, at, &s[1..]);
                    emit(&mut w, at + 1, &d[1..]);
                }
            }
        }
    }
    assert_eq!(written, rank, "hky.basis/v1: rank");
    let (digest, bytes) = match seal(w.finish(), Path::new(&tmp), Path::new(path)) {
        Ok(x) => x,
        Err(e) => {
            eprintln!("hky.basis/v1: export FAILED ({e}); {tmp} removed");
            return ExportReport {
                path: path.to_string(),
                n_cols: n_cols as u64,
                rank,
                profile,
                column_order: if split { 1 } else { 0 },
                contains_one,
                seconds: t0.elapsed().as_secs_f64(),
                error: Some(e.to_string()),
                ..ExportReport::default()
            };
        }
    };
    // whole-file hash, read back from disk
    let mut h = Sha256::new();
    let mut rd = BufReader::with_capacity(8 << 20, File::open(path).expect("hky.basis/v1: reopen"));
    let mut buf = vec![0u8; 8 << 20];
    loop {
        let n = rd.read(&mut buf).expect("hky.basis/v1: read back");
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    ExportReport {
        path: path.to_string(),
        content_sha256: hex(&digest),
        file_sha256: hex(&h.finalize()),
        bytes,
        n_cols: n_cols as u64,
        rank,
        profile,
        column_order: if split { 1 } else { 0 },
        contains_one,
        seconds: t0.elapsed().as_secs_f64(),
        error: None,
    }
}
