//! With LOSAT_X_BPCLUT: the fused seed loop of blastp-fast
//! (`x_seed_fused_compressed.rs`) on the flat compressed lookup table
//! (`tblastx/lookup/x_clut_direct.rs`).
//!
//! The loop is `x_scan_two_hit_compressed` with the flat table's scan (`for_each_cell`) in place of
//! `x_scan_cells`: the same subject offsets and cells in the same order, and each cell's hits as its
//! slice, which holds them in the order of `x_cell_hits` (NCBI's scan order; the permutation is
//! proved in `x_clut_direct.rs`). The diagonal test, the prefetch of a later hit's entry and the
//! closure are those of `x_scan_two_hit_compressed`, so the argument of `x_seed_fused_compressed.rs`
//! carries over unchanged.
//!
//! EXPERIMENT (LOSAT_X_BPCLUTSCANBF): the scan lists the offsets without a branch on the presence
//! bit (`for_each_cell::<true, _>`, the form of `x_seed_fused::x_scan_cells_bf`, or
//! `x_list_pairs_lean` below); the list, and so the cells and their order, are the same. That
//! instance is out of line with a closure of its own type (`x_scan_flat_bf`), so the instance with
//! the branch, which `=0` takes, is compiled as without it. Every flat table takes it: until
//! 2026-11-02 a table with fewer than 500 used cells per 10,000 (`LOSAT_X_BPCLUTSCANBFMIN`) kept
//! the loop with the branch natively, measured when the list was `for_each_cell`'s (x1.07 the
//! cycles at 234 per 10,000, #437); with the loops below, single queries of 10 to 739 residues x
//! UniRef50 volume 00, whose tables kept it, took x0.60-0.68 the cycles of the whole search without
//! the branch (line SF S3, family bpf-seed-loop). A switch of the strict set (`utils/xswitch.rs`):
//! on unless it is `0`.
//!
//! The loops after the branch-free list, on every target (until 2026-11-02 the switches
//! LOSAT_X_BPFSEED, LOSAT_X_BPFSTAGE and LOSAT_X_BPFLEAN, folded here by line SF S3, family
//! bpf-seed-loop: each `=0` was slower on every input measured, #472, #543, #551; WebAssembly
//! since 2026-11-18, below). On a table that takes the branch-free
//! list, a subject where `XHitBf::fits` holds runs the loop of `x_scan_flat_fseed_loop`, which
//! differs from `x_scan_flat_bf` in two ways: each hit goes first through blastp's branch-free
//! two-hit test (`XHitBf::test`, the test of LOSAT_X_BPHITBF, which does not act here; its argument
//! in `x_seed_fused.rs` holds here, as the diagonal table, its offsets and the extension's stores
//! are those of blastp's loop), and a flag-1 entry or a candidate goes on to `x_test_hit`,
//! unchanged; and no hit reads ahead a later hit's diagonal entry (a cache hint). The cells, their
//! order, the hits and their order are those of `x_scan_flat_bf`, the table holds after each hit
//! what `x_test_hit` leaves there, and the extension is called for the same hits with the same
//! values and entries.
//!
//! Where `XFlatCompressedLookup::x_pairs_fit` holds as well (a subject of a word or more, with
//! blastp-fast's 15-letter alphabet and words of 2 to 6 letters), the subject runs
//! `x_scan_flat_lean`, which differs from `x_scan_flat_fseed_loop` in how it reads, not in what it
//! tests. The scan is `x_list_pairs_lean`, two offsets per step in a loop of its own, out of line
//! (the same list). Each cell's hits are staged in scan order, in runs of [`X_RUN`] copied whole
//! (so that no loop ends at a cell's last hit, a mispredict per cell), as `q - s` and
//! `s + diagonal offset`, by `x_lean_copy` (out of line, without a call), which stops when fewer
//! than [`X_RUN`] of the [`X_STAGE_CAP`] places are left or the list ends; one loop
//! (`x_stage_test`, out of line) then gives the staged hits in order to the test of
//! `x_scan_flat_fseed_loop`: the diagonal is `(q - s) & mask` and the query offset `(q - s) + s`.
//! A run's second half of 4 hits is read from `k + 4` only where the cell has more than 4 hits
//! left, else from `k` again (places past the run's length, written and left unread), so that a
//! cell of few hits does not read the line after them. So the hits, their order, the tests, the
//! entries and the calls are those of `x_scan_flat_fseed_loop` (which makes its whole list before
//! any test, too); the staged values do not depend on the table's entries. Each loop is out of line
//! with a closure of its own type, so the loops without it keep their code.
//!
//! WebAssembly takes the same loops since 2026-11-18 (line W2). Until then wasm32 kept
//! `x_scan_flat_bf` on every subject (`x_fseed_on`), as in Node the branch-free test had taken
//! x1.06-1.09 the time of 3 proteins x UniRef50, and the staged loop x1.13-1.18 on 3 proteins and
//! on 42 residues. Measured again on dev ad3fb1927 (the serial WASI command in Node 26.8.2, one
//! process, 7 rounds in rotating order, one core of ws node 1; blastp-fast, `-evalue 0.001`,
//! proteins of one species of an OrthoFinder set against the 6,600 proteins of another; every
//! report identical): 8,000, 1,000 and 100 query proteins took x0.878, x0.896 and x0.889 the time
//! of the whole search (7/7 rounds faster each), 10 and 1 proteins and one of 1,001 residues
//! x0.918, x0.957 and x0.946 (runs of 0.04-0.08 s). In callgrind (the compiled code alone, dev
//! 7d968922d) the instructions were x1.015 at 1,000 proteins and x0.91-0.98 on the smaller cases,
//! the mispredicted branches of its model x0.72 at 1,000 proteins (the seed loop's own 19.3 M ->
//! 1.9 M), x0.76-0.86 on 10 and 100 proteins and a long one, x1.03 of 0.31 M on one protein: the
//! time saved comes with fewer mispredicted branches, not fewer instructions.

use std::cell::Cell;

use super::x_seed_fused::{x_prefetch, x_test_hit, XHitBf, X_PREFETCH_DISTANCE};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_prefetch as x_prefetch_ptr, XFlatCompressedLookup, X_HITS_AHEAD, X_STARTS_AHEAD,
};

thread_local! {
    // No NCBI counterpart: the scan's list of (subject offset, backbone index), reused across
    // subjects.
    static X_CELLS: Cell<Option<Vec<(i32, u32)>>> = const { Cell::new(None) };
    // No NCBI counterpart: the staged hits of `x_scan_flat_lean`, reused across subjects.
    static X_STAGE: Cell<Option<Box<XStage>>> = const { Cell::new(None) };
}

// No NCBI counterpart: counters; they do not change any value NCBI computes.
/// The subjects each loop of [`x_scan_two_hit_compressed_flat`] scanned (counted with the `xstats`
/// feature, for the line at exit): [`x_scan_flat_lean`], [`x_scan_flat_fseed_loop`],
/// [`x_scan_flat_bf`] and the loop with the branch.
static X_FLAT_SUBJECTS: [std::sync::atomic::AtomicU64; 4] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 4];

// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
/// With LOSAT_X_STATS, the `[X_STATS] BPFLAT lean=N fseed=N bf=N branch=N` line at exit: the
/// subjects each loop on the flat compressed table scanned (only when one of them is not zero).
pub(crate) fn x_bp_flat_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_none() {
        return;
    }
    let [lean, fseed, bf, branch] = X_FLAT_SUBJECTS
        .each_ref()
        .map(|n| n.load(std::sync::atomic::Ordering::Relaxed));
    if lean + fseed + bf + branch > 0 {
        eprintln!("[X_STATS] BPFLAT lean={lean} fseed={fseed} bf={bf} branch={branch}");
    }
}

// No NCBI counterpart: reads the LOSAT_X_BPCLUTSCANBF switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_clut_scan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPCLUTSCANBF"))
}

/// `x_seed_fused_compressed::x_scan_two_hit_compressed` on the flat table (same arguments).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-518
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_compressed_flat<const CALLER: u8, F>(
    lookup: &XFlatCompressedLookup,
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
    // EXPERIMENT (LOSAT_X_BPCLUTSCANBF): the offsets listed without a branch on the presence bit.
    if x_bp_clut_scan_bf() {
        // The branch-free two-hit test and no read-ahead of the diagonal entries, where the test
        // gives the reference's outcomes (module header).
        if let Some(bf) = x_fseed_test(diag_offset, subject_length, window, word_size) {
            x_scan_flat_fseed::<CALLER, _>(
                lookup,
                subject,
                subject_length,
                diag_array,
                diag_mask,
                diag_offset,
                window,
                word_size,
                bf,
                &mut candidate,
            );
            return;
        }
        crate::utils::xstats::add(&X_FLAT_SUBJECTS[2], 1);
        x_scan_flat_bf(
            lookup,
            subject,
            subject_length,
            diag_array,
            diag_mask,
            diag_offset,
            window,
            word_size,
            &mut candidate,
        );
        return;
    }
    crate::utils::xstats::add(&X_FLAT_SUBJECTS[3], 1);
    x_scan_flat_in::<false, F>(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        &mut candidate,
    );
}

/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): [`x_scan_two_hit_compressed_flat`] with the offsets listed
/// without a branch on the presence bit. Out of line, and with a closure of its own type, so that
/// the instance with the branch is compiled as without this one (`x_seed_fused::x_scan_two_hit_bf`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_bf<F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut own = |query_offset, subject_offset, last_hit, diff, entry: &mut DiagStruct| {
        candidate(query_offset, subject_offset, last_hit, diff, entry)
    };
    x_scan_flat_in::<true, _>(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        &mut own,
    );
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's outcomes.
/// The branch-free test of a subject of a table with the branch-free list where [`XHitBf::fits`]
/// holds, which then takes [`x_scan_flat_fseed`].
fn x_fseed_test(
    diag_offset: i32,
    subject_length: usize,
    window: i32,
    word_size: i32,
) -> Option<XHitBf> {
    XHitBf::new(diag_offset, subject_length, window, word_size)
}

/// The loop of a subject with the test `bf`: [`x_scan_flat_lean`] where the subject fits its scan
/// (`x_pairs_fit`), else [`x_scan_flat_fseed_loop`].
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_fseed<const CALLER: u8, F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    if lookup.x_pairs_fit(subject, subject_length) {
        crate::utils::xstats::add(&X_FLAT_SUBJECTS[0], 1);
        // SAFETY: x_pairs_fit holds (tested above).
        unsafe {
            x_scan_flat_lean::<CALLER, _>(
                lookup,
                subject,
                subject_length,
                diag_array,
                diag_mask,
                diag_offset,
                window,
                word_size,
                bf,
                candidate,
            )
        };
        return;
    }
    crate::utils::xstats::add(&X_FLAT_SUBJECTS[1], 1);
    x_scan_flat_fseed_loop(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        bf,
        candidate,
    );
}

/// [`x_scan_flat_bf`] with the branch-free two-hit test `bf` (made
/// for the subject, [`XHitBf::new`]) first and no read-ahead of the diagonal entries (module
/// header).
/// Out of line, and with a closure of its own type, so that the loops without it are compiled as
/// without this one (see [`x_scan_flat_bf`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at
// `x_seed_fused::x_test_hit`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_fseed_loop<F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    let mut cells = X_CELLS.with(Cell::take).unwrap_or_default();
    lookup.for_each_cell::<true, _>(subject, subject_length, &mut cells, |s, cell_hits| {
        let here = bf.here(s);
        for &query_offset in cell_hits {
            let diag_coord = ((query_offset - s) & diag_mask) as usize;
            // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
            let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
            if bf.test(entry, here) {
                continue;
            }
            if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size) {
                candidate(query_offset, s, last_hit, diff, entry);
            }
        }
    });
    X_CELLS.with(|slot| slot.set(Some(cells)));
}

// No NCBI counterpart: sizes of LOSAT's staging; they do not change any value NCBI computes.
/// The hits copied per step of [`x_lean_copy`] (a cell of more hits takes several steps), and the
/// staged hits tested per call of [`x_stage_test`].
const X_RUN: usize = 8;
const X_STAGE_CAP: usize = 1024;
// The copy of a run after a flush writes X_RUN places of the stage.
const _: () = assert!(X_RUN <= X_STAGE_CAP);

// No NCBI counterpart: the staged hits of LOSAT's loop after the branch-free list.
/// Staged hit `i` of [`x_scan_flat_lean`] is `diffs[i] = q - s` and `heres[i] = s + diagonal
/// offset` (two's complement).
struct XStage {
    diffs: Vec<u32>,
    heres: Vec<u32>,
}

/// The test of [`x_scan_flat_fseed_loop`] on staged hits in
/// order (`diffs[i] = q - s`, `heres[i] = s + diag_offset`), out of line so that its loop keeps
/// its own registers.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at
// `x_seed_fused::x_test_hit`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_stage_test<F>(
    diffs: &[u32],
    heres: &[u32],
    diag_array: &mut [DiagStruct],
    mask: u32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_array.len() > mask as usize);
    for (&diff, &here) in diffs.iter().zip(heres) {
        // SAFETY: diff & mask <= mask < diag_array.len() (asserted above).
        let entry = unsafe { diag_array.get_unchecked_mut((diff & mask) as usize) };
        if bf.test(entry, here) {
            continue;
        }
        let s = here.wrapping_sub(diag_offset as u32) as i32;
        if let Some((last_hit, gap)) = x_test_hit(entry, s, diag_offset, window, word_size) {
            candidate((diff as i32).wrapping_add(s), s, last_hit, gap, entry);
        }
    }
}

/// [`x_scan_flat_fseed_loop`] with the list of `x_list_pairs_lean` and each cell's hits staged by
/// the copy [`x_lean_copy`] before the test (module header): each call of the copy fills the empty
/// stage from where the last one stopped, and [`x_stage_test`] tests what it staged, until the list
/// ends. The read-ahead of the slice bounds and of the hits is that of `for_each_cell`, without the
/// tests: the list holds [`X_STARTS_AHEAD`] entries of cell 0 past its end while the loop runs.
///
/// # Safety
///
/// `lookup.x_pairs_fit(subject, subject_length)` holds (the scan reads without bounds checks, and
/// the slice bounds of every listed index exist).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (cited at
// `x_scan_two_hit_compressed_flat`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn x_scan_flat_lean<const CALLER: u8, F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    let mut cells = X_CELLS.with(Cell::take).unwrap_or_default();
    // SAFETY: x_pairs_fit holds (this function's contract).
    unsafe { lookup.x_list_pairs_lean::<CALLER>(subject, subject_length, &mut cells) };
    let listed = cells.len();
    cells.extend_from_slice(&[(0, 0); X_STARTS_AHEAD]);
    // The stage made here, not by a helper (#551 measured this code).
    let mut stage = X_STAGE.with(Cell::take).unwrap_or_else(|| {
        Box::new(XStage {
            diffs: vec![0; X_STAGE_CAP],
            heres: vec![0; X_STAGE_CAP],
        })
    });
    let (diffs, heres) = (&mut stage.diffs[..], &mut stage.heres[..]);
    let starts = &lookup.starts[..];
    let hits = &lookup.hits[..];
    let mask = diag_mask as u32;
    let mut at = XLeanAt {
        n: 0,
        k: 0,
        end: 0,
        s: 0,
    };
    loop {
        // SAFETY: the list's indices are below 15^word length < starts.len() - 1 (x_pairs_fit),
        // and `cells` holds X_STARTS_AHEAD entries past them (above).
        let staged =
            unsafe { x_lean_copy(&cells, listed, starts, hits, diffs, heres, bf, &mut at) };
        x_stage_test(
            &diffs[..staged],
            &heres[..staged],
            diag_array,
            mask,
            diag_offset,
            window,
            word_size,
            bf,
            candidate,
        );
        if at.n >= listed && at.k >= at.end {
            break;
        }
    }
    cells.truncate(listed);
    X_STAGE.with(|slot| slot.set(Some(stage)));
    X_CELLS.with(|slot| slot.set(Some(cells)));
}

// No NCBI counterpart: where LOSAT's staging copy stopped.
/// Where [`x_lean_copy`] stopped: the next cell `n` of the list and,
/// of the cell before it, the hits `k..end` not yet staged and its subject offset `s`.
struct XLeanAt {
    n: usize,
    k: usize,
    end: usize,
    s: i32,
}

/// The copy of [`x_scan_flat_lean`], out of line and without the test's call so that its state
/// stays in registers: from `at`, stages runs ([`x_lean_run`]) at the start of the stage until
/// fewer than [`X_RUN`] places are left or the first `listed` cells of `cells` are done, leaves in
/// `at` where it stopped and returns the staged count. It reads ahead, [`X_STARTS_AHEAD`] cells on,
/// the slice bound and the one after it, and the first hit of the cell [`X_HITS_AHEAD`] on.
///
/// # Safety
///
/// Every index in `cells` is below `starts.len() - 1` (the slice bounds are read without bounds
/// checks), as in the list of `x_list_pairs_lean` where `x_pairs_fit` holds and the zero entries
/// past it.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (cited at
// `x_scan_two_hit_compressed_flat`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn x_lean_copy(
    cells: &[(i32, u32)],
    listed: usize,
    starts: &[u32],
    hits: &[i32],
    diffs: &mut [u32],
    heres: &mut [u32],
    bf: XHitBf,
    at: &mut XLeanAt,
) -> usize {
    assert!(cells.len() >= listed + X_STARTS_AHEAD);
    assert!(diffs.len() >= X_STAGE_CAP && heres.len() >= X_STAGE_CAP);
    let (mut n, mut k, mut end, mut s) = (at.n, at.k, at.end, at.s);
    let diffs = diffs.as_mut_ptr();
    let heres = heres.as_mut_ptr();
    let mut staged = 0usize;
    'cells: loop {
        let here = bf.here(s);
        while k < end {
            if staged + X_RUN > X_STAGE_CAP {
                break 'cells;
            }
            // SAFETY: staged + X_RUN <= X_STAGE_CAP, the length of both arrays (above).
            let take = unsafe { x_lean_run(hits, k, end, s, here, diffs, heres, staged) };
            staged += take;
            k += take;
        }
        if n >= listed {
            break;
        }
        // SAFETY: n + X_STARTS_AHEAD < cells.len() (asserted above; X_HITS_AHEAD <=
        // X_STARTS_AHEAD), and every index in `cells` is below starts.len() - 1 (this function's
        // contract).
        unsafe {
            let (cell_s, index) = *cells.get_unchecked(n);
            // Cache hints only, as in `for_each_cell`, and the bound after `ahead`'s.
            let ahead = cells.get_unchecked(n + X_STARTS_AHEAD).1 as usize;
            x_prefetch_ptr(starts.as_ptr().add(ahead));
            x_prefetch_ptr(starts.as_ptr().add(ahead + 1));
            let later = cells.get_unchecked(n + X_HITS_AHEAD).1 as usize;
            x_prefetch_ptr(
                hits.as_ptr()
                    .wrapping_add(*starts.get_unchecked(later) as usize),
            );
            s = cell_s;
            k = *starts.get_unchecked(index as usize) as usize;
            end = *starts.get_unchecked(index as usize + 1) as usize;
        }
        n += 1;
    }
    *at = XLeanAt { n, k, end, s };
    staged
}

// x_lean_run reads a run in two halves of 4.
const _: () = assert!(X_RUN == 8);

/// Stages the run of hits `k..min(end, k + X_RUN)` of the cell at subject offset `s`
/// (`here = bf.here(s)`) at place `staged` and returns its length. Where [`X_RUN`] hits from `k`
/// are in the array the run is copied whole in two halves of 4, the second read from `k + 4` where
/// the cell has more than 4 hits left and else from `k` again; the places past the run's length
/// are written and later overwritten or left unread. Else (the last hits of the array) it is
/// copied hit by hit.
///
/// # Safety
///
/// `diffs` and `heres` have `staged + X_RUN` places (`k < end` is the caller's, for the length).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (cited at
// `x_scan_two_hit_compressed_flat`)
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn x_lean_run(
    hits: &[i32],
    k: usize,
    end: usize,
    s: i32,
    here: u32,
    diffs: *mut u32,
    heres: *mut u32,
    staged: usize,
) -> usize {
    let take = (end - k).min(X_RUN);
    if k + X_RUN <= hits.len() {
        // SAFETY: both reads end at or before hi_at + 4 <= k + X_RUN <= hits.len() (tested),
        // and the places are this function's contract.
        unsafe {
            let p = hits.as_ptr();
            let lo = std::ptr::read_unaligned(p.add(k).cast::<[i32; 4]>());
            let hi_at = k + 4 * usize::from(end - k > 4);
            let hi = std::ptr::read_unaligned(p.add(hi_at).cast::<[i32; 4]>());
            std::ptr::write_unaligned(
                diffs.add(staged).cast::<[u32; 4]>(),
                lo.map(|q| q.wrapping_sub(s) as u32),
            );
            std::ptr::write_unaligned(
                diffs.add(staged + 4).cast::<[u32; 4]>(),
                hi.map(|q| q.wrapping_sub(s) as u32),
            );
            std::ptr::write_unaligned(heres.add(staged).cast::<[u32; X_RUN]>(), [here; X_RUN]);
        }
    } else {
        for (j, &q) in hits[k..k + take].iter().enumerate() {
            // SAFETY: j < take <= X_RUN (this function's contract).
            unsafe {
                *diffs.add(staged + j) = q.wrapping_sub(s) as u32;
                *heres.add(staged + j) = here;
            }
        }
    }
    take
}

/// The loop of [`x_scan_two_hit_compressed_flat`]; `BRANCH_FREE` is `for_each_cell`'s.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (cited at
// `x_scan_two_hit_compressed_flat`)
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_in<const BRANCH_FREE: bool, F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    let prefetch_distance = X_PREFETCH_DISTANCE;
    let mut cells = X_CELLS.with(Cell::take).unwrap_or_default();
    lookup.for_each_cell::<BRANCH_FREE, _>(subject, subject_length, &mut cells, |s, cell_hits| {
        for (i, &query_offset) in cell_hits.iter().enumerate() {
            // Read ahead the entry of a later hit of this cell (a cache hint; see `x_scan_two_hit`).
            if let Some(&ahead) = cell_hits.get(i + prefetch_distance) {
                x_prefetch(diag_array, ((ahead - s) & diag_mask) as usize);
            }
            let diag_coord = ((query_offset - s) & diag_mask) as usize;
            // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
            let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
            if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size) {
                candidate(query_offset, s, last_hit, diff, entry);
            }
        }
    });
    X_CELLS.with(|slot| slot.set(Some(cells)));
}

#[cfg(test)]
#[path = "x_seed_fused_clut_tests.rs"]
mod tests;
