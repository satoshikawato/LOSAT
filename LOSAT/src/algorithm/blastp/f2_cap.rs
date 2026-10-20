//! F2 PROBE (exp/f2-bp-cellcap, measurement only, not for merging): with `F2CAP=<dir>` the default
//! seed stage of blastp writes what the subject scan reads and writes for the first query block's
//! first `F2CAP_N` subjects (default 2000): the lookup table (presence vector, backbone, overflow),
//! the subjects, the diagonal table's constants and initial state, and per candidate the values
//! the extension closure leaves in its diagonal entry. A replay outside LOSAT (a small crate on the
//! workstation) reads the files. Without `F2CAP` nothing here runs: the call site tests `enabled()`
//! once per subject and calls `x_scan_two_hit` as before; the capture is out of line and calls
//! `x_scan_two_hit` with a closure of its own type, so the default instance is unchanged.

use super::x_seed_fused::x_scan_two_hit;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::{BlastAaLookupTable, AA_HITS_PER_CELL};
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

struct Config {
    dir: PathBuf,
    max_subjects: usize,
}

fn config() -> &'static Option<Config> {
    static C: OnceLock<Option<Config>> = OnceLock::new();
    C.get_or_init(|| {
        let dir = std::env::var_os("F2CAP")?;
        let max_subjects = std::env::var("F2CAP_N")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(2000);
        Some(Config {
            dir: PathBuf::from(dir),
            max_subjects,
        })
    })
}

#[inline]
pub(super) fn enabled() -> bool {
    config().is_some()
}

#[derive(Default)]
struct State {
    started: bool,
    done: bool,
    lookup_ptr: usize,
    thread: Option<std::thread::ThreadId>,
    subjects: Vec<u8>,
    index: Vec<u8>,
    cands: Vec<u8>,
    nsub: u64,
    cells: u64,
    hits: u64,
    ncand: u64,
    discontinuities: u64,
    last_post_hash: u64,
    table_len: usize,
    diag_mask: i32,
    window: i32,
    word_size: i32,
}

static STATE: Mutex<Option<State>> = Mutex::new(None);

/// FxHash over the table's packed words (the replay computes the same).
fn table_hash(table: &[DiagStruct]) -> u64 {
    let mut h = 0u64;
    for d in table {
        h = (h.rotate_left(5) ^ u64::from(d.raw_bits())).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
    h
}

fn as_bytes<T: Copy>(v: &[T]) -> &[u8] {
    // SAFETY: plain-old-data slices (u32, i32, BackboneCell of i32s), read as bytes.
    unsafe { std::slice::from_raw_parts(v.as_ptr().cast::<u8>(), std::mem::size_of_val(v)) }
}

fn write_file(dir: &std::path::Path, name: &str, bytes: &[u8]) {
    let mut f = std::fs::File::create(dir.join(name)).expect("F2CAP: create file");
    f.write_all(bytes).expect("F2CAP: write file");
}

/// Independent count of the scan's cells (presence bit set) and hits, by a plain rolling index.
fn count_cells_hits(lookup: &BlastAaLookupTable, subject: &[u8], subject_length: usize) -> (u64, u64) {
    let wl = lookup.word_length as usize;
    let cs = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let s_last = subject_length.saturating_sub(wl);
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0) as usize;
    let mut index = 0usize;
    for i in 0..wl - 1 {
        index = (index << cs) | byte_at(i);
    }
    let (mut cells, mut hits) = (0u64, 0u64);
    for s in 0..=s_last {
        index = ((index << cs) | byte_at(s + wl - 1)) & mask;
        if lookup.pv[index >> 5] >> (index & 31) & 1 != 0 {
            cells += 1;
            hits += lookup.get_hits(index).len() as u64;
        }
    }
    (cells, hits)
}

#[cold]
#[inline(never)]
fn start(st: &mut State, cfg: &Config, lookup: &BlastAaLookupTable, diag_array: &[DiagStruct]) {
    std::fs::create_dir_all(&cfg.dir).expect("F2CAP: create dir");
    write_file(&cfg.dir, "pv.bin", as_bytes(&lookup.pv));
    write_file(&cfg.dir, "backbone.bin", as_bytes(&lookup.backbone));
    write_file(&cfg.dir, "overflow.bin", as_bytes(&lookup.overflow));
    let init: Vec<u32> = diag_array.iter().map(|d| d.raw_bits()).collect();
    write_file(&cfg.dir, "table_init.bin", as_bytes(&init));
    st.started = true;
    st.table_len = diag_array.len();
    st.last_post_hash = table_hash(diag_array);
    st.lookup_ptr = lookup as *const BlastAaLookupTable as usize;
    st.thread = Some(std::thread::current().id());
    let mut m = String::new();
    m += &format!("word_length={}\n", lookup.word_length);
    m += &format!("charsize={}\n", lookup.charsize);
    m += &format!("mask={}\n", lookup.mask);
    m += &format!("pv_len={}\n", lookup.pv.len());
    m += &format!("backbone_len={}\n", lookup.backbone.len());
    m += &format!("overflow_len={}\n", lookup.overflow.len());
    m += &format!("hits_per_cell={AA_HITS_PER_CELL}\n");
    m += &format!(
        "backbone_cell_bytes={}\n",
        std::mem::size_of::<crate::algorithm::tblastx::lookup::BackboneCell>()
    );
    m += &format!("table_len={}\n", diag_array.len());
    m += &format!("max_subjects={}\n", cfg.max_subjects);
    std::fs::write(cfg.dir.join("meta_start.txt"), m).expect("F2CAP: meta");
}

#[cold]
#[inline(never)]
fn finish(st: &mut State, cfg: &Config) {
    write_file(&cfg.dir, "subjects.bin", &st.subjects);
    write_file(&cfg.dir, "subjects.idx", &st.index);
    write_file(&cfg.dir, "cands.bin", &st.cands);
    let mut m = String::new();
    m += &format!("nsub={}\n", st.nsub);
    m += &format!("cells={}\n", st.cells);
    m += &format!("hits={}\n", st.hits);
    m += &format!("candidates={}\n", st.ncand);
    m += &format!("discontinuities={}\n", st.discontinuities);
    m += &format!("diag_mask={}\n", st.diag_mask);
    m += &format!("window={}\n", st.window);
    m += &format!("word_size={}\n", st.word_size);
    m += &format!("table_len={}\n", st.table_len);
    std::fs::write(cfg.dir.join("meta_end.txt"), m).expect("F2CAP: meta end");
    st.done = true;
    eprintln!(
        "[F2CAP] wrote {} subjects, {} cells, {} hits, {} candidates to {}",
        st.nsub,
        st.cells,
        st.hits,
        st.ncand,
        cfg.dir.display()
    );
}

/// `x_scan_two_hit` with the capture; the same arguments, the same effects.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn scan<F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let cfg = config().as_ref().expect("F2CAP config");
    let mut guard = STATE.lock().unwrap();
    let st = guard.get_or_insert_with(State::default);
    if !st.done {
        let ptr = lookup as *const BlastAaLookupTable as usize;
        if !st.started {
            start(st, cfg, lookup, diag_array);
        }
        if st.lookup_ptr != ptr || st.thread != Some(std::thread::current().id()) {
            // The next query block (or another thread): stop here.
            finish(st, cfg);
        }
    }
    if st.done {
        drop(guard);
        x_scan_two_hit(
            lookup,
            subject,
            subject_length,
            diag_array,
            diag_mask,
            diag_offset,
            window,
            word_size,
            candidate,
        );
        return;
    }
    st.diag_mask = diag_mask;
    st.window = window;
    st.word_size = word_size;
    let pre_hash = table_hash(diag_array);
    if pre_hash != st.last_post_hash {
        st.discontinuities += 1;
    }
    let (cells, hits) = count_cells_hits(lookup, subject, subject_length);
    st.cells += cells;
    st.hits += hits;
    let mut recs: Vec<u8> = Vec::new();
    let mut n = 0u32;
    x_scan_two_hit(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        |q: i32, s: i32, lh: i32, d: i32, e: &mut DiagStruct| {
            let before = e.raw_bits();
            candidate(q, s, lh, d, e);
            let after = e.raw_bits();
            for v in [q, s, lh, d] {
                recs.extend_from_slice(&v.to_le_bytes());
            }
            recs.extend_from_slice(&before.to_le_bytes());
            recs.extend_from_slice(&after.to_le_bytes());
            n += 1;
        },
    );
    let post_hash = table_hash(diag_array);
    st.last_post_hash = post_hash;
    // idx record (56 bytes): subject offset u64, slice len u32, aa_len u32, diag_offset i32,
    // n candidates u32, cells u64, hits u64, pre hash u64, post hash u64 -> 8+4+4+4+4+8+8+8+8 = 56.
    st.index.extend_from_slice(&(st.subjects.len() as u64).to_le_bytes());
    st.index.extend_from_slice(&(subject.len() as u32).to_le_bytes());
    st.index.extend_from_slice(&(subject_length as u32).to_le_bytes());
    st.index.extend_from_slice(&diag_offset.to_le_bytes());
    st.index.extend_from_slice(&n.to_le_bytes());
    st.index.extend_from_slice(&cells.to_le_bytes());
    st.index.extend_from_slice(&hits.to_le_bytes());
    st.index.extend_from_slice(&pre_hash.to_le_bytes());
    st.index.extend_from_slice(&post_hash.to_le_bytes());
    st.subjects.extend_from_slice(subject);
    st.cands.extend_from_slice(&recs);
    st.ncand += u64::from(n);
    st.nsub += 1;
    if st.nsub as usize >= cfg.max_subjects {
        finish(st, cfg);
    }
}
