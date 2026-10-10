//! EXPERIMENT (LOSAT_X_DPCAPTURE=<dir>, LOSAT_X_DPCAPTURE_EVERY=<n>, default 25): diagnostic
//! capture of the X-drop DP problems that reach `xdrop_align`, for the replay benchmark of the
//! vector kernels (`LOSAT_X_DPREPLAY`, test `replay_captured_problems` in `x_xdrop_avx2.rs`).
//!
//! No NCBI counterpart: a recording of the inputs and the result of calls of the port of
//! `Blast_SemiGappedAlign` / `ALIGN_EX` / `s_BlastAlignPackedNucl`
//! (c++/src/algo/blast/core/blast_gapalign.c at 598d8ae6); it does not change any value NCBI
//! computes. Off unless LOSAT_X_DPCAPTURE is set (read once).
//!
//! Every Nth call of a thread (the 1st, the N+1st, ...) is written as one self-contained record
//! to `<dir>/xdp-<pid>-t<k>.xdp` (one file per thread, opened for append): the traceback flag,
//! check_fence, the row residues `q.get(a)` for `a = 1..=len1`, the column residues `s.get(k)`
//! for `k = 0..=len2` (materialised, so a replay needs no packed or reversed form), the score
//! matrix as `n x n` i32 (`Tables` decoded from its 16-bit planes; entries below -32768 were
//! clamped there and act the same in every kernel), gap_open, gap_extend, x_drop, len1, len2, and
//! the result (offsets, score, fence_hit, cells, edit operations). Little-endian layout:
//!
//! ```text
//! u32 magic "XDP1", u32 body length
//! u8 version (1), u8 tb, u8 check_fence, u8 kind (0 Rows28, 1 Flat, 2 Func, 3 Tables)
//! u8 handled (the call returned Some), u8 fence_hit, u16 0
//! i32 gap_open, i32 gap_extend, i32 x_drop
//! u32 len1, u32 len2, u32 n, u32 a_offset, u32 b_offset, i32 score, u64 cells, u32 n_ops
//! u8 rows[len1], u8 cols[len2 + 1], i32 matrix[n * n], (u8 op, u32 count)[n_ops]
//! ```

use super::{ColSeq, RowSeq, Scores, XdropResult, TAB};
use std::cell::{Cell, RefCell};
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::OnceLock;

const MAGIC: [u8; 4] = *b"XDP1";
const VERSION: u8 = 1;
const KIND_ROWS28: u8 = 0;
const KIND_FLAT: u8 = 1;
const KIND_FUNC: u8 = 2;
const KIND_TABLES: u8 = 3;

struct Config {
    dir: PathBuf,
    every: u64,
}

fn config() -> Option<&'static Config> {
    static CONFIG: OnceLock<Option<Config>> = OnceLock::new();
    CONFIG
        .get_or_init(|| {
            let dir = std::env::var_os("LOSAT_X_DPCAPTURE")?;
            let every = std::env::var("LOSAT_X_DPCAPTURE_EVERY")
                .ok()
                .and_then(|v| v.parse::<u64>().ok())
                .filter(|&n| n >= 1)
                .unwrap_or(25);
            Some(Config {
                dir: PathBuf::from(dir),
                every,
            })
        })
        .as_ref()
}

/// LOSAT_X_DPCAPTURE is set.
pub(super) fn switch_on() -> bool {
    config().is_some()
}

thread_local! {
    static CALLS: Cell<u64> = const { Cell::new(0) };
    static INSIDE: Cell<bool> = const { Cell::new(false) };
    static FILE: RefCell<Option<std::fs::File>> = const { RefCell::new(None) };
}

/// Counts one call of this thread; true for every Nth one (the first included). Calls made
/// inside `without_capture` are not counted.
pub(super) fn take_sample() -> bool {
    let Some(cfg) = config() else {
        return false;
    };
    if INSIDE.with(|f| f.get()) {
        return false;
    }
    CALLS.with(|c| {
        let n = c.get();
        c.set(n + 1);
        n % cfg.every == 0
    })
}

/// Runs `f` (the call being recorded) without counting or recording the nested `xdrop_align`.
pub(super) fn without_capture<R>(f: impl FnOnce() -> R) -> R {
    INSIDE.with(|g| g.set(true));
    let r = f();
    INSIDE.with(|g| g.set(false));
    r
}

/// The score matrix as (kind, n, n x n values).
fn matrix(scores: &Scores<'_>) -> Option<(u8, usize, Vec<i32>)> {
    match *scores {
        Scores::Rows28(m) => Some((
            KIND_ROWS28,
            28,
            m.iter().flat_map(|row| row.iter().copied()).collect(),
        )),
        Scores::Flat { data, n } => Some((KIND_FLAT, n, data.get(..n * n)?.to_vec())),
        Scores::Func { f, n } => {
            if n > 64 {
                return None;
            }
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n {
                for j in 0..n {
                    v.push(f(i as u8, j as u8));
                }
            }
            Some((KIND_FUNC, n, v))
        }
        Scores::Tables { t, n } => {
            if n > 32 {
                return None;
            }
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n {
                for j in 0..n {
                    let lo = t.data[i * TAB + j];
                    let hi = t.data[i * TAB + 32 + j];
                    v.push(i16::from_le_bytes([lo, hi]) as i32);
                }
            }
            Some((KIND_TABLES, n, v))
        }
    }
}

/// One record (see the module header), or None when a row residue cannot be read or the matrix
/// is not representable (the call is then not recorded).
fn encode(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    tb: bool,
    check_fence: bool,
    result: Option<&XdropResult>,
    ops: &[(u8, u32)],
) -> Option<Vec<u8>> {
    let (kind, n, m) = matrix(scores)?;
    let mut rows = Vec::with_capacity(len1);
    for a in 1..=len1 {
        rows.push(q.get(a)?);
    }
    let ops: &[(u8, u32)] = if result.is_some() { ops } else { &[] };
    let res = result.copied().unwrap_or(XdropResult {
        a_offset: 0,
        b_offset: 0,
        score: 0,
        fence_hit: false,
        cells: 0,
    });
    let mut b = Vec::with_capacity(64 + len1 + len2 + 4 * m.len() + 5 * ops.len());
    b.extend_from_slice(&MAGIC);
    b.extend_from_slice(&0u32.to_le_bytes());
    b.extend_from_slice(&[VERSION, u8::from(tb), u8::from(check_fence), kind]);
    b.extend_from_slice(&[u8::from(result.is_some()), u8::from(res.fence_hit), 0, 0]);
    for v in [gap_open, gap_extend, x_drop] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    for v in [len1, len2, n, res.a_offset, res.b_offset] {
        b.extend_from_slice(&(v as u32).to_le_bytes());
    }
    b.extend_from_slice(&res.score.to_le_bytes());
    b.extend_from_slice(&res.cells.to_le_bytes());
    b.extend_from_slice(&(ops.len() as u32).to_le_bytes());
    b.extend_from_slice(&rows);
    b.extend((0..=len2).map(|k| s.get(k)));
    for v in &m {
        b.extend_from_slice(&v.to_le_bytes());
    }
    for &(op, count) in ops {
        b.push(op);
        b.extend_from_slice(&count.to_le_bytes());
    }
    let body = (b.len() - 8) as u32;
    b[4..8].copy_from_slice(&body.to_le_bytes());
    Some(b)
}

/// Appends one record to this thread's file.
pub(super) fn record(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    tb: bool,
    check_fence: bool,
    result: Option<&XdropResult>,
    ops: &[(u8, u32)],
) {
    static NEXT_FILE: AtomicUsize = AtomicUsize::new(0);
    static WARNED: AtomicBool = AtomicBool::new(false);
    let Some(cfg) = config() else {
        return;
    };
    let Some(buf) = encode(
        q,
        s,
        scores,
        len1,
        len2,
        gap_open,
        gap_extend,
        x_drop,
        tb,
        check_fence,
        result,
        ops,
    ) else {
        return;
    };
    FILE.with(|cell| {
        let mut file = cell.borrow_mut();
        if file.is_none() {
            #[cfg(not(target_family = "wasm"))]
            let pid = std::process::id();
            #[cfg(target_family = "wasm")]
            let pid = 0u32;
            let k = NEXT_FILE.fetch_add(1, Ordering::Relaxed);
            let path = cfg.dir.join(format!("xdp-{pid}-t{k:03}.xdp"));
            let opened = std::fs::create_dir_all(&cfg.dir).and_then(|()| {
                std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&path)
            });
            match opened {
                Ok(f) => *file = Some(f),
                Err(e) => {
                    if !WARNED.swap(true, Ordering::Relaxed) {
                        eprintln!("LOSAT_X_DPCAPTURE: cannot open {}: {e}", path.display());
                    }
                    return;
                }
            }
        }
        if let Some(f) = file.as_mut() {
            if let Err(e) = f.write_all(&buf) {
                if !WARNED.swap(true, Ordering::Relaxed) {
                    eprintln!("LOSAT_X_DPCAPTURE: write failed: {e}");
                }
            }
        }
    });
}

// ---------------------------------------------------------------------------
// Replay side (tests only)
// ---------------------------------------------------------------------------

/// One recorded call, read back from a capture file.
#[cfg(test)]
pub(super) struct Captured {
    pub(super) tb: bool,
    pub(super) check_fence: bool,
    pub(super) kind: u8,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) x_drop: i32,
    pub(super) len1: usize,
    pub(super) len2: usize,
    pub(super) n: usize,
    pub(super) rows: Vec<u8>,
    pub(super) cols: Vec<u8>,
    pub(super) matrix: Vec<i32>,
    /// What the recorded call returned (None: the vector kernels declined it).
    pub(super) result: Option<XdropResult>,
    pub(super) ops: Vec<(u8, u32)>,
    rows28: Option<std::rc::Rc<[[i32; 28]; 28]>>,
    tables: Option<std::rc::Rc<super::StaticTables>>,
}

#[cfg(test)]
impl Captured {
    /// Row residues in the forward `Bytes` form: `get(a) = rows[a - 1]`.
    pub(super) fn row_seq(&self) -> RowSeq<'_> {
        RowSeq::Bytes {
            data: &self.rows,
            base: -1,
            step: 1,
        }
    }

    /// Column residues in the forward form: `get(k) = cols[k]` (0 past the end).
    pub(super) fn col_seq(&self) -> ColSeq<'_> {
        ColSeq {
            data: &self.cols,
            base: 0,
            step: 1,
            zero_from: usize::MAX,
        }
    }

    /// The score source of the recorded call's kind (`Func` replays as `Flat`, `Tables` as tables
    /// rebuilt from the decoded matrix, which gives the same table bytes).
    pub(super) fn scores(&self) -> Scores<'_> {
        if let Some(m) = &self.rows28 {
            return Scores::Rows28(m);
        }
        if let Some(t) = &self.tables {
            return Scores::Tables { t, n: self.n };
        }
        Scores::Flat {
            data: &self.matrix,
            n: self.n,
        }
    }

    /// Replays the call through `xdrop_align` with the forward `Bytes` rows and columns.
    pub(super) fn replay(&self, sc: &mut super::XdropScratch) -> Option<XdropResult> {
        super::xdrop_align(
            &self.row_seq(),
            &self.col_seq(),
            &self.scores(),
            self.len1,
            self.len2,
            self.gap_open,
            self.gap_extend,
            self.x_drop,
            self.tb,
            self.check_fence,
            sc,
        )
    }
}

/// Reads every `*.xdp` file of `dir` (sorted by name), records in file order.
#[cfg(test)]
pub(super) fn load_dir(dir: &std::path::Path) -> std::io::Result<Vec<Captured>> {
    use std::collections::HashMap;
    use std::rc::Rc;
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "xdp"))
        .collect();
    files.sort();
    let bad = |what: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, what.to_string());
    let mut tables: HashMap<Vec<i32>, Rc<super::StaticTables>> = HashMap::new();
    let mut rows28: HashMap<Vec<i32>, Rc<[[i32; 28]; 28]>> = HashMap::new();
    let mut out = Vec::new();
    for path in files {
        let data = std::fs::read(&path)?;
        let mut pos = 0usize;
        while pos < data.len() {
            if data.len() < pos + 8 || data[pos..pos + 4] != MAGIC {
                return Err(bad("bad record header"));
            }
            let body = u32::from_le_bytes(data[pos + 4..pos + 8].try_into().unwrap()) as usize;
            let rec = data
                .get(pos + 8..pos + 8 + body)
                .ok_or_else(|| bad("truncated record"))?;
            pos += 8 + body;
            let mut r = Reader { b: rec, at: 0 };
            let version = r.u8()?;
            if version != VERSION {
                return Err(bad("unknown record version"));
            }
            let tb = r.u8()? != 0;
            let check_fence = r.u8()? != 0;
            let kind = r.u8()?;
            let handled = r.u8()? != 0;
            let fence_hit = r.u8()? != 0;
            r.take(2)?;
            let gap_open = r.i32()?;
            let gap_extend = r.i32()?;
            let x_drop = r.i32()?;
            let len1 = r.u32()? as usize;
            let len2 = r.u32()? as usize;
            let n = r.u32()? as usize;
            let a_offset = r.u32()? as usize;
            let b_offset = r.u32()? as usize;
            let score = r.i32()?;
            let cells = u64::from_le_bytes(r.take(8)?.try_into().unwrap());
            let n_ops = r.u32()? as usize;
            let rows_v = r.take(len1)?.to_vec();
            let cols = r.take(len2 + 1)?.to_vec();
            let mut matrix = Vec::with_capacity(n * n);
            for _ in 0..n * n {
                matrix.push(r.i32()?);
            }
            let mut ops = Vec::with_capacity(n_ops);
            for _ in 0..n_ops {
                let op = r.u8()?;
                ops.push((op, r.u32()?));
            }
            let mut c = Captured {
                tb,
                check_fence,
                kind,
                gap_open,
                gap_extend,
                x_drop,
                len1,
                len2,
                n,
                rows: rows_v,
                cols,
                matrix,
                result: handled.then_some(XdropResult {
                    a_offset,
                    b_offset,
                    score,
                    fence_hit,
                    cells,
                }),
                ops,
                rows28: None,
                tables: None,
            };
            if kind == KIND_ROWS28 && n == 28 {
                let m = rows28.entry(c.matrix.clone()).or_insert_with(|| {
                    let mut m = [[0i32; 28]; 28];
                    for (i, row) in m.iter_mut().enumerate() {
                        row.copy_from_slice(&c.matrix[i * 28..i * 28 + 28]);
                    }
                    Rc::new(m)
                });
                c.rows28 = Some(m.clone());
            } else if kind == KIND_TABLES {
                if let Some(t) = tables.get(&c.matrix) {
                    c.tables = Some(t.clone());
                } else {
                    let mm = c.matrix.clone();
                    let f = move |a: u8, b: u8| mm[a as usize * n + b as usize];
                    let t: Rc<super::StaticTables> = super::build_tables(&f, n)
                        .ok_or_else(|| bad("tables do not rebuild"))?
                        .into();
                    tables.insert(c.matrix.clone(), t.clone());
                    c.tables = Some(t);
                }
            }
            out.push(c);
        }
    }
    Ok(out)
}

#[cfg(test)]
struct Reader<'a> {
    b: &'a [u8],
    at: usize,
}

#[cfg(test)]
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> std::io::Result<&'a [u8]> {
        let s = self
            .b
            .get(self.at..self.at + n)
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "short record"))?;
        self.at += n;
        Ok(s)
    }
    fn u8(&mut self) -> std::io::Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u32(&mut self) -> std::io::Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> std::io::Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record written by `encode` reads back with the same fields.
    #[test]
    fn capture_record_round_trips() {
        let q = [3u8, 5, 7, 1];
        let s = [2u8, 4, 6, 8, 10];
        let m: Vec<i32> = (0..16 * 16).map(|i| (i % 7) - 3).collect();
        let rows = RowSeq::Bytes {
            data: &q,
            base: 4,
            step: -1,
        };
        let cols = ColSeq {
            data: &s,
            base: 4,
            step: -1,
            zero_from: 4,
        };
        let scores = Scores::Flat { data: &m, n: 16 };
        let res = XdropResult {
            a_offset: 3,
            b_offset: 2,
            score: 17,
            fence_hit: false,
            cells: 99,
        };
        let ops = [(3u8, 2u32), (0, 1)];
        let buf = encode(
            &rows,
            &cols,
            &scores,
            4,
            4,
            11,
            1,
            38,
            true,
            true,
            Some(&res),
            &ops,
        )
        .unwrap();
        let dir = std::env::temp_dir().join(format!("losat-xdp-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("xdp-test.xdp");
        std::fs::write(&path, [buf.clone(), buf].concat()).unwrap();
        let got = load_dir(&dir).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(got.len(), 2);
        let c = &got[1];
        assert!(c.tb && c.check_fence);
        assert_eq!((c.kind, c.n, c.len1, c.len2), (KIND_FLAT, 16, 4, 4));
        assert_eq!((c.gap_open, c.gap_extend, c.x_drop), (11, 1, 38));
        // rows read backwards from base 4: q[3], q[2], q[1], q[0]
        assert_eq!(c.rows, vec![1, 7, 5, 3]);
        // columns read backwards from base 4, zero from k = 4
        assert_eq!(c.cols, vec![10, 8, 6, 4, 0]);
        assert_eq!(c.matrix, m);
        assert_eq!(c.result, Some(res));
        assert_eq!(c.ops, ops.to_vec());
        for a in 1..=4 {
            assert_eq!(c.row_seq().get(a), rows.get(a));
        }
        for k in 0..=6 {
            assert_eq!(c.col_seq().get(k), cols.get(k));
        }
    }
}
