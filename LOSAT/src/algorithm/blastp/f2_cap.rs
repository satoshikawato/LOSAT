//! F2 PROBE (exp/f2-bp-extcap, measurement only, not for merging): with `F2CAP=<dir>` the default
//! seed stage of blastp writes what the candidate closure (context check, `two_hit` extension,
//! diagonal entry update, init HSP) reads and writes, for the first query block's first
//! `F2CAP_N` subjects (default 2000): the query block (concatenated, sentinels as they are), the
//! contexts' frame bases / X-drops / cutoffs, the subjects (whole `aa_seq`), and per candidate the
//! arguments, the diagonal entry before and after, the context index and the extension's result.
//! A replay outside LOSAT (a small crate on the workstation) reads the files. Without `F2CAP`
//! nothing here runs: the call site tests `enabled()` once per subject and the default closure is
//! used; the capture is out of line (`#[cold]`, `#[inline(never)]`).

use super::x_seed_fused::x_scan_two_hit;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

struct Config {
    dir: PathBuf,
    max_subjects: usize,
}

static DONE: AtomicBool = AtomicBool::new(false);

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

/// True while the capture is open (once it has written its files the default closure runs).
#[inline]
pub(super) fn enabled() -> bool {
    config().is_some() && !DONE.load(Ordering::Relaxed)
}

#[derive(Default)]
struct State {
    started: bool,
    lookup_ptr: usize,
    thread: Option<std::thread::ThreadId>,
    // block
    query: Vec<u8>,
    frame_bases: Vec<i32>,
    x_drops: Vec<i32>,
    cutoffs: Vec<i32>,
    args_word_size: usize,
    word_size: i32,
    // per subject
    subjects: Vec<u8>,
    index: Vec<u8>,
    cands: Vec<u8>,
    nsub: u64,
    ncand: u64,
    diag_mask: i32,
    window: i32,
    table_len: usize,
}

static STATE: Mutex<Option<State>> = Mutex::new(None);

/// Called before `scan` for every subject of the capture: the block's constants (kept from the
/// first call of a block).
#[cold]
#[inline(never)]
pub(super) fn set_block(
    query: &[u8],
    frame_bases: &[i32],
    x_drops: &[i32],
    cutoffs: &[i32],
    args_word_size: usize,
    word_size: i32,
) {
    let mut guard = STATE.lock().unwrap();
    let st = guard.get_or_insert_with(State::default);
    if st.started {
        return;
    }
    st.query = query.to_vec();
    st.frame_bases = frame_bases.to_vec();
    st.x_drops = x_drops.to_vec();
    st.cutoffs = cutoffs.to_vec();
    st.args_word_size = args_word_size;
    st.word_size = word_size;
}

fn as_bytes<T: Copy>(v: &[T]) -> &[u8] {
    // SAFETY: plain-old-data slices (u32, i32), read as bytes.
    unsafe { std::slice::from_raw_parts(v.as_ptr().cast::<u8>(), std::mem::size_of_val(v)) }
}

fn write_file(dir: &std::path::Path, name: &str, bytes: &[u8]) {
    let mut f = std::fs::File::create(dir.join(name)).expect("F2CAP: create file");
    f.write_all(bytes).expect("F2CAP: write file");
}

#[cold]
#[inline(never)]
fn finish(st: &mut State, cfg: &Config) {
    std::fs::create_dir_all(&cfg.dir).expect("F2CAP: create dir");
    write_file(&cfg.dir, "query.bin", &st.query);
    let mut ctx: Vec<i32> = Vec::new();
    for i in 0..st.frame_bases.len() {
        ctx.extend_from_slice(&[st.frame_bases[i], st.x_drops[i], st.cutoffs[i]]);
    }
    write_file(&cfg.dir, "ctx.bin", as_bytes(&ctx));
    write_file(&cfg.dir, "subjects.bin", &st.subjects);
    write_file(&cfg.dir, "subjects.idx", &st.index);
    write_file(&cfg.dir, "cands.bin", &st.cands);
    let mut m = String::new();
    m += &format!("nsub={}\n", st.nsub);
    m += &format!("candidates={}\n", st.ncand);
    m += &format!("nctx={}\n", st.frame_bases.len());
    m += &format!("query_len={}\n", st.query.len());
    m += &format!("args_word_size={}\n", st.args_word_size);
    m += &format!("word_size={}\n", st.word_size);
    m += &format!("diag_mask={}\n", st.diag_mask);
    m += &format!("window={}\n", st.window);
    m += &format!("table_len={}\n", st.table_len);
    std::fs::write(cfg.dir.join("meta_end.txt"), m).expect("F2CAP: meta end");
    DONE.store(true, Ordering::Relaxed);
    eprintln!(
        "[F2CAP] wrote {} subjects, {} candidates, {} contexts, query {} to {}",
        st.nsub,
        st.ncand,
        st.frame_bases.len(),
        st.query.len(),
        cfg.dir.display()
    );
}

/// `x_scan_two_hit` with the capture; the same arguments, the same effects. `candidate` is the
/// default closure's body and returns [ctx_idx, flags, q_start, s_start, length, score, s_last, 0]
/// (flags: bits 0-1 = 0 context early return / 1 extension None / 2 extension Some, bit 2 =
/// right_extend, bit 3 = init HSP saved).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn scan<F>(
    lookup: &BlastAaLookupTable,
    aa_seq: &[u8],
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct) -> [i32; 8],
{
    let cfg = config().as_ref().expect("F2CAP config");
    let mut guard = STATE.lock().unwrap();
    let st = guard.get_or_insert_with(State::default);
    let ptr = lookup as *const BlastAaLookupTable as usize;
    if !st.started {
        st.started = true;
        st.lookup_ptr = ptr;
        st.thread = Some(std::thread::current().id());
        st.table_len = diag_array.len();
    }
    if st.lookup_ptr != ptr || st.thread != Some(std::thread::current().id()) {
        // The next query block (or another thread): the capture ends before this subject.
        finish(st, cfg);
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
            |q: i32, s: i32, lh: i32, d: i32, e: &mut DiagStruct| {
                candidate(q, s, lh, d, e);
            },
        );
        return;
    }
    st.diag_mask = diag_mask;
    st.window = window;
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
            let x = candidate(q, s, lh, d, e);
            let after = e.raw_bits();
            for v in [q, s, lh, d] {
                recs.extend_from_slice(&v.to_le_bytes());
            }
            recs.extend_from_slice(&before.to_le_bytes());
            recs.extend_from_slice(&after.to_le_bytes());
            for v in &x[..7] {
                recs.extend_from_slice(&v.to_le_bytes());
            }
            n += 1;
        },
    );
    // idx record (24 bytes): subject offset u64, aa_seq len u32, aa_len u32, diag_offset i32, n u32.
    st.index.extend_from_slice(&(st.subjects.len() as u64).to_le_bytes());
    st.index.extend_from_slice(&(aa_seq.len() as u32).to_le_bytes());
    st.index.extend_from_slice(&(subject_length as u32).to_le_bytes());
    st.index.extend_from_slice(&diag_offset.to_le_bytes());
    st.index.extend_from_slice(&n.to_le_bytes());
    st.subjects.extend_from_slice(aa_seq);
    st.cands.extend_from_slice(&recs);
    st.ncand += u64::from(n);
    st.nsub += 1;
    if st.nsub as usize >= cfg.max_subjects {
        finish(st, cfg);
    }
}
