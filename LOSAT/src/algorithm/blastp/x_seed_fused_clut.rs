//! EXPERIMENT (LOSAT_X_BPSEED with LOSAT_X_BPCLUT): the fused seed loop of blastp-fast
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
//! EXPERIMENT (LOSAT_X_BPCLUTSCANBF): on a table with at least LOSAT_X_BPCLUTSCANBFMIN used cells
//! per 10,000 (`XFlatCompressedLookup::x_scan_bf_pays`) the scan lists the offsets without a branch
//! on the presence bit (`for_each_cell::<true, _>`, the form of LOSAT_X_BPSCANBF's
//! `x_seed_fused::x_scan_cells_bf`); the list, and so the cells and their order, are the same. That
//! instance is out of line with a closure of its own type (`x_scan_flat_bf`), so the instance with
//! the branch, which the other tables take, is compiled as without it. A switch of the strict set
//! (`utils/xswitch.rs`): on unless it is `0`.
//!
//! EXPERIMENT (LOSAT_X_BPFSEED): on a table that takes the branch-free list, a subject where
//! `XHitBf::fits` holds runs the loop of `x_scan_flat_fseed_loop`, which differs in two ways: each hit
//! goes first through blastp's branch-free two-hit test (`XHitBf::test`, the test of
//! LOSAT_X_BPHITBF, which does not act here; its argument in `x_seed_fused.rs` holds here, as the
//! diagonal table, its offsets and the extension's stores are those of blastp's loop), and a
//! flag-1 entry or a candidate goes on to
//! `x_test_hit`, unchanged; and no hit reads ahead a later hit's diagonal entry (a cache hint).
//! The cells, their order, the hits and their order are those of `x_scan_flat_bf`, the table
//! holds after each hit what `x_test_hit` leaves there, and the extension is called for the same
//! hits with the same values and entries. The loop is out of line with a closure of its own type,
//! so the loops without it keep their code. A sparser table keeps the loop with the branch: there
//! few positions have hits, and the test gained nothing on the workstation. WebAssembly keeps its
//! loops (`x_fseed_on`): in Node the loop took x1.06-1.09 the time of 3 proteins x UniRef50. A
//! switch of the strict set: on unless it is `0`.
//!
//! EXPERIMENT (LOSAT_X_BPFSTAGE): a subject that takes the loop of LOSAT_X_BPFSEED, where
//! `XFlatCompressedLookup::x_pairs_fit` holds, runs `x_scan_flat_stage` instead, which differs in
//! three ways. The scan is out of line and two offsets per step (`x_list_pairs`), with the same
//! list. Each cell's hits are staged in scan order, in runs of [`X_RUN`] copied whole (so that no
//! loop ends at a cell's last hit, a mispredict per cell), as `q - s` and `s + diagonal offset`;
//! when fewer than [`X_RUN`] of the [`X_STAGE_CAP`] places are left, and at the end, one loop
//! (`x_stage_test`, out of line) gives the staged hits in order to the test of
//! `x_scan_flat_fseed_loop`: the diagonal is `(q - s) & mask` and the query offset `(q - s) + s`.
//! So the hits, their order, the tests, the entries and the calls are those of
//! `x_scan_flat_fseed_loop` (which makes its whole list before any test, too); the staged values do
//! not depend on the table's entries. Where `x_pairs_fit` does not hold (a subject shorter than a
//! word) the subject takes `x_scan_flat_fseed_loop`; the chooser is `x_scan_flat_fseed`. A switch
//! of the strict set: on unless it is `0`; it acts only where LOSAT_X_BPFSEED's loop runs.

use std::cell::Cell;

use super::x_seed_fused::{x_prefetch, x_prefetch_distance, x_test_hit, XHitBf};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_prefetch as x_prefetch_ptr, XFlatCompressedLookup, X_HITS_AHEAD, X_STARTS_AHEAD,
};

thread_local! {
    // No NCBI counterpart: the scan's list of (subject offset, backbone index), reused across
    // subjects.
    static X_CELLS: Cell<Option<Vec<(i32, u32)>>> = const { Cell::new(None) };
    // No NCBI counterpart: LOSAT_X_BPFSTAGE's staged hits, reused across subjects.
    static X_STAGE: Cell<Option<Box<XStage>>> = const { Cell::new(None) };
}

// No NCBI counterpart: reads the LOSAT_X_BPFSTAGE switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPFSTAGE): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_fstage() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPFSTAGE"))
}

// No NCBI counterpart: reads the LOSAT_X_BPFSEED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPFSEED): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_fseed() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPFSEED"))
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's outcomes.
/// EXPERIMENT (LOSAT_X_BPFSEED): whether tables with the branch-free list take
/// [`x_scan_flat_fseed`]: the switch is on and the target is not wasm32 (module header).
pub(crate) fn x_fseed_on() -> bool {
    !cfg!(target_arch = "wasm32") && x_bp_fseed()
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's outcomes.
/// EXPERIMENT (LOSAT_X_BPFSEED): whether the subjects of `lookup` may take [`x_scan_flat_fseed`]
/// (the rule of [`x_scan_two_hit_compressed_flat`] before [`XHitBf::fits`]).
pub(crate) fn x_fseed_table(lookup: &XFlatCompressedLookup) -> bool {
    x_clut_scan_bf(lookup) && x_fseed_on()
}

// No NCBI counterpart: reads the LOSAT_X_BPCLUTSCANBF switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_clut_scan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPCLUTSCANBF"))
}

// No NCBI counterpart: chooses LOSAT's list form; both list the same offsets.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): whether the fused loop lists `lookup`'s offsets without a
/// branch on the presence bit: the switch is on and the table has enough cells used.
pub(crate) fn x_clut_scan_bf(lookup: &XFlatCompressedLookup) -> bool {
    x_bp_clut_scan_bf() && lookup.x_scan_bf_pays()
}

/// `x_seed_fused_compressed::x_scan_two_hit_compressed` on the flat table (same arguments).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518
// ```c
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//
//             /* calculate the diagonal associated with this query-subject pair
//              */
//
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_compressed_flat<F>(
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
    // EXPERIMENT (LOSAT_X_BPCLUTSCANBF): on a table with enough cells used, the offsets listed
    // without a branch on the presence bit (`for_each_cell`).
    if x_clut_scan_bf(lookup) {
        // EXPERIMENT (LOSAT_X_BPFSEED): the branch-free two-hit test and no read-ahead of the
        // diagonal entries, where the test gives the reference's outcomes (module header).
        if let Some(bf) = x_fseed_on()
            .then(|| x_fseed_test(diag_offset, subject_length, window, word_size))
            .flatten()
        {
            x_scan_flat_fseed(
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
/// EXPERIMENT (LOSAT_X_BPFSEED): the branch-free test of a subject of a table with the branch-free
/// list where [`XHitBf::fits`] holds, which then takes [`x_scan_flat_fseed`]; the side is recorded
/// for `LOSAT_X_REACH`.
fn x_fseed_test(
    diag_offset: i32,
    subject_length: usize,
    window: i32,
    word_size: i32,
) -> Option<XHitBf> {
    let bf = XHitBf::new(diag_offset, subject_length, window, word_size);
    crate::utils::xreach::x_reach!("LOSAT_X_BPFSEED", bf.is_some());
    bf
}

/// EXPERIMENT (LOSAT_X_BPFSEED, LOSAT_X_BPFSTAGE): the loop of a subject with the test `bf`:
/// [`x_scan_flat_stage`] where LOSAT_X_BPFSTAGE is on and the subject fits its scan, else
/// [`x_scan_flat_fseed_loop`]; the side is recorded for `LOSAT_X_REACH`.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_fseed<F>(
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
    if x_bp_fstage() {
        let fits = lookup.x_pairs_fit(subject, subject_length);
        crate::utils::xreach::x_reach!("LOSAT_X_BPFSTAGE", fits);
        if fits {
            // SAFETY: x_pairs_fit holds (tested above).
            unsafe {
                x_scan_flat_stage(
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
    }
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

/// EXPERIMENT (LOSAT_X_BPFSEED): [`x_scan_flat_bf`] with the branch-free two-hit test `bf` (made
/// for the subject, [`XHitBf::new`]) first and no read-ahead of the diagonal entries (module
/// header).
/// Out of line, and with a closure of its own type, so that the loops without it are compiled as
/// without this one (see [`x_scan_flat_bf`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at
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
/// EXPERIMENT (LOSAT_X_BPFSTAGE): the hits copied per step of [`x_scan_flat_stage`] (a cell of
/// more hits takes several steps), and the staged hits tested per call of [`x_stage_test`].
const X_RUN: usize = 8;
const X_STAGE_CAP: usize = 1024;
// The copy of a run after a flush writes X_RUN places of the stage.
const _: () = assert!(X_RUN <= X_STAGE_CAP);

// No NCBI counterpart: LOSAT_X_BPFSTAGE's staged hits.
/// EXPERIMENT (LOSAT_X_BPFSTAGE): staged hit `i` of [`x_scan_flat_stage`] is `diffs[i] = q - s`
/// and `heres[i] = s + diagonal offset` (two's complement).
struct XStage {
    diffs: Vec<u32>,
    heres: Vec<u32>,
}

impl XStage {
    fn new() -> Box<Self> {
        Box::new(Self {
            diffs: vec![0; X_STAGE_CAP],
            heres: vec![0; X_STAGE_CAP],
        })
    }
}

/// EXPERIMENT (LOSAT_X_BPFSTAGE): [`x_scan_flat_fseed_loop`] with the list made by
/// `x_list_pairs` and each cell's hits staged before the test (module header). The read-ahead of
/// the slice bounds and of the hits is that of `for_each_cell`, without the tests: the list holds
/// [`X_STARTS_AHEAD`] entries of cell 0 past its end while the loop runs.
///
/// # Safety
///
/// `lookup.x_pairs_fit(subject, subject_length)` holds (the scan reads without bounds checks, and
/// the slice bounds of every listed index exist).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (quoted at
// `x_scan_two_hit_compressed_flat`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn x_scan_flat_stage<F>(
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
    unsafe { lookup.x_list_pairs(subject, subject_length, &mut cells) };
    let listed = cells.len();
    cells.extend_from_slice(&[(0, 0); X_STARTS_AHEAD]);
    let mut stage = X_STAGE.with(Cell::take).unwrap_or_else(XStage::new);
    let XStage { diffs, heres } = &mut *stage;
    let starts = &lookup.starts[..];
    let hits = &lookup.hits[..];
    let mask = diag_mask as u32;
    let mut staged = 0usize;
    for n in 0..listed {
        // SAFETY: n + X_STARTS_AHEAD < cells.len() (the entries above, and X_HITS_AHEAD <=
        // X_STARTS_AHEAD); every index of the list is below 15^word length < starts.len() - 1
        // (x_pairs_fit).
        let (s, start, end) = unsafe {
            let (s, index) = *cells.get_unchecked(n);
            // Cache hints only, as in `for_each_cell`.
            let ahead = cells.get_unchecked(n + X_STARTS_AHEAD).1 as usize;
            x_prefetch_ptr(starts.as_ptr().add(ahead));
            let later = cells.get_unchecked(n + X_HITS_AHEAD).1 as usize;
            x_prefetch_ptr(
                hits.as_ptr()
                    .wrapping_add(*starts.get_unchecked(later) as usize),
            );
            (
                s,
                *starts.get_unchecked(index as usize) as usize,
                *starts.get_unchecked(index as usize + 1) as usize,
            )
        };
        let here = bf.here(s);
        let mut k = start;
        loop {
            if staged + X_RUN > X_STAGE_CAP {
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
                staged = 0;
            }
            let take = (end - k).min(X_RUN);
            if k + X_RUN <= hits.len() {
                // SAFETY: k + X_RUN <= hits.len() (tested) and staged + X_RUN <= X_STAGE_CAP, the
                // length of both arrays (above). The places past `take` are written and later
                // overwritten or left unread.
                unsafe {
                    let run = std::ptr::read_unaligned(hits.as_ptr().add(k).cast::<[i32; X_RUN]>());
                    std::ptr::write_unaligned(
                        diffs.as_mut_ptr().add(staged).cast::<[u32; X_RUN]>(),
                        run.map(|q| q.wrapping_sub(s) as u32),
                    );
                    std::ptr::write_unaligned(
                        heres.as_mut_ptr().add(staged).cast::<[u32; X_RUN]>(),
                        [here; X_RUN],
                    );
                }
            } else {
                for (j, &q) in hits[k..k + take].iter().enumerate() {
                    diffs[staged + j] = q.wrapping_sub(s) as u32;
                    heres[staged + j] = here;
                }
            }
            staged += take;
            k += take;
            if k >= end {
                break;
            }
        }
    }
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
    cells.truncate(listed);
    X_STAGE.with(|slot| slot.set(Some(stage)));
    X_CELLS.with(|slot| slot.set(Some(cells)));
}

/// EXPERIMENT (LOSAT_X_BPFSTAGE): the test of [`x_scan_flat_fseed_loop`] on staged hits in
/// order (`diffs[i] = q - s`, `heres[i] = s + diag_offset`), out of line so that its loop keeps
/// its own registers.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at
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

/// The loop of [`x_scan_two_hit_compressed_flat`]; `BRANCH_FREE` is `for_each_cell`'s.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (quoted at
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
    let prefetch_distance = x_prefetch_distance();
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
mod tests {
    use super::*;
    use crate::algorithm::blastp::x_seed_fused_compressed::x_scan_two_hit_compressed;
    use crate::algorithm::tblastx::lookup::compressed::try_build_blosum62_compressed_lookup;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_try_build_flat_compressed_lookup, XFlatBuild,
    };

    // A small deterministic generator (xorshift64).
    struct XorShift(u64);

    impl XorShift {
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x % n as u64) as usize
        }
    }

    // NCBISTDAA letters inside the 15-letter compressed alphabet and outside it (NULLB, X, U, *, O).
    const LETTERS: [u8; 23] = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 27,
    ];
    const OUTSIDE: [u8; 5] = [0, 21, 24, 25, 26];

    // A stand-in for the extension: it logs the call and writes the entry as the reference closure
    // does (a flag and a later `last_hit` after a right extension, else this hit).
    fn stand_in(
        log: &mut Vec<(i32, i32, i32, i32)>,
        diag_offset: i32,
        query_offset: i32,
        subject_offset: i32,
        last_hit: i32,
        diff: i32,
        entry: &mut DiagStruct,
    ) {
        log.push((query_offset, subject_offset, last_hit, diff));
        match (query_offset ^ subject_offset) & 3 {
            0 => {
                entry.set_flag(1);
                entry.set_last_hit(subject_offset + 7 + diag_offset);
            }
            1 => {}
            _ => entry.set_last_hit(subject_offset + diag_offset),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    //             /* If the reset bit is set, an extension just happened. */
    //             if (diag_array[diag_coord].flag) {
    //             ...
    //                 if (diff >= window) {
    //             ...
    //                 if (diff < wordsize) {
    //                     continue;
    //                 }
    // ```
    // Queries of few letters and low thresholds (cells of many hits, several lookup segments),
    // subjects with letters outside the alphabet and copies of the query, small diagonal tables: the
    // loop on the flat table makes the candidates and leaves the diagonal table of the loop of
    // `x_seed_fused_compressed.rs` on the reference's table, subject after subject, with the list
    // of the branch and the branch-free one (LOSAT_X_BPCLUTSCANBF), and so does the loop of
    // LOSAT_X_BPFSEED.
    #[test]
    fn the_flat_loop_is_the_fused_loop_on_the_reference_table() {
        let mut rng = XorShift(0x510e_527f_ade6_82d1);
        for (threshold, letters) in [
            (1100, &[1u8, 11, 17, 10][..]),
            (900, &[1u8, 11, 17][..]),
            (1300, &LETTERS[..]),
        ] {
            let mut query = Vec::new();
            let mut locations = Vec::new();
            for length in [120usize, 75, 200, 33] {
                let left = query.len() as i32;
                query.extend((0..length).map(|_| letters[rng.below(letters.len())]));
                // Two segments (a masked run between them).
                let cut = left + length as i32 / 2;
                locations.push((left, cut - 3));
                locations.push((cut, left + length as i32 - 1));
                query.push(0);
            }
            let tables = [0..4, 6..8];
            let reference =
                try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                    .expect("word size 5")
                    .expect("not refused");
            let flat =
                match x_try_build_flat_compressed_lookup(5, threshold, &query, &locations, &tables)
                {
                    Some(XFlatBuild::Table(table)) => table,
                    _ => panic!("flat table"),
                };
            assert!(reference.longest_chain > 5, "cells with overflow lists");
            let mut subjects = Vec::new();
            for outside_one_in in [0usize, 50, 6, 2] {
                for _ in 0..6 {
                    let length = 1 + rng.below(300);
                    let mut subject: Vec<u8> = (0..length)
                        .map(|_| {
                            if outside_one_in != 0 && rng.below(outside_one_in) == 0 {
                                OUTSIDE[rng.below(OUTSIDE.len())]
                            } else {
                                LETTERS[rng.below(LETTERS.len())]
                            }
                        })
                        .collect();
                    if length > 60 && rng.below(2) == 0 {
                        let from = rng.below(query.len() - 60);
                        let to = rng.below(length - 50);
                        subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                    }
                    subject.push(0);
                    subjects.push(subject);
                }
            }
            for diag_len in [64usize, 1024, 4096] {
                let window = 40;
                let diag_mask = diag_len as i32 - 1;
                let mut flat_diag = vec![DiagStruct::default(); diag_len];
                let mut bf_diag = vec![DiagStruct::default(); diag_len];
                let mut fseed_diag = vec![DiagStruct::default(); diag_len];
                let mut stage_diag = vec![DiagStruct::default(); diag_len];
                let mut reference_diag = vec![DiagStruct::default(); diag_len];
                let mut diag_offset = window;
                let mut total = 0usize;
                for (n, subject) in subjects.iter().enumerate() {
                    let length = subject.len() - 1;
                    // Both instances, whichever the switch and the table's density choose.
                    let mut flat_log = Vec::new();
                    x_scan_flat_in::<false, _>(
                        &flat,
                        subject,
                        length,
                        &mut flat_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut flat_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    let mut bf_log = Vec::new();
                    x_scan_flat_bf(
                        &flat,
                        subject,
                        length,
                        &mut bf_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut bf_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    // LOSAT_X_BPFSEED (the subjects are in XHitBf's range).
                    let bf = XHitBf::new(diag_offset, length, window, 5).expect("in range");
                    let mut fseed_log = Vec::new();
                    x_scan_flat_fseed_loop(
                        &flat,
                        subject,
                        length,
                        &mut fseed_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        bf,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut fseed_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    // LOSAT_X_BPFSTAGE on the subjects of a word or more (they keep their
                    // sentinel); a shorter one takes the loop without the switch.
                    assert_eq!(flat.x_pairs_fit(subject, length), length >= 5);
                    let mut stage_log = Vec::new();
                    if length >= 5 {
                        // SAFETY: x_pairs_fit holds (asserted above).
                        unsafe {
                            x_scan_flat_stage(
                                &flat,
                                subject,
                                length,
                                &mut stage_diag,
                                diag_mask,
                                diag_offset,
                                window,
                                5,
                                bf,
                                &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                                    stand_in(
                                        &mut stage_log,
                                        diag_offset,
                                        q,
                                        s,
                                        last_hit,
                                        diff,
                                        entry,
                                    )
                                },
                            )
                        };
                    } else {
                        x_scan_flat_fseed_loop(
                            &flat,
                            subject,
                            length,
                            &mut stage_diag,
                            diag_mask,
                            diag_offset,
                            window,
                            5,
                            bf,
                            &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                                stand_in(&mut stage_log, diag_offset, q, s, last_hit, diff, entry)
                            },
                        );
                    }
                    let mut reference_log = Vec::new();
                    x_scan_two_hit_compressed(
                        &reference,
                        subject,
                        length,
                        &mut reference_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        |q, s, last_hit, diff, entry| {
                            stand_in(&mut reference_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    assert_eq!(flat_log, reference_log, "subject {n}: candidates");
                    assert_eq!(bf_log, reference_log, "subject {n}: branch-free candidates");
                    assert_eq!(fseed_log, reference_log, "subject {n}: BPFSEED candidates");
                    assert_eq!(stage_log, reference_log, "subject {n}: BPFSTAGE candidates");
                    for (diag, form) in [
                        (&flat_diag, "branch"),
                        (&bf_diag, "branch-free"),
                        (&fseed_diag, "BPFSEED"),
                        (&stage_diag, "BPFSTAGE"),
                    ] {
                        assert!(
                            diag.iter()
                                .zip(&reference_diag)
                                .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                            "subject {n}: diagonal table ({form})"
                        );
                    }
                    total += flat_log.len();
                    diag_offset += length as i32 + window;
                }
                assert!(total > 0, "the test reaches the extension");
            }
        }
    }

    // LOSAT_X_BPFSTAGE: the cell whose slice ends the hits array is copied hit by hit (no run of
    // X_RUN reads past the array), with the outcomes of the loop without the switch.
    #[test]
    fn the_staged_loop_copies_the_last_cell_of_the_hits() {
        let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
        let mut query = Vec::new();
        let mut locations = Vec::new();
        for length in [400usize, 300] {
            let left = query.len() as i32;
            query.extend((0..length).map(|_| LETTERS[rng.below(LETTERS.len())]));
            locations.push((left, left + length as i32 - 1));
            query.push(0);
        }
        let flat = match x_try_build_flat_compressed_lookup(
            5,
            1100,
            &query,
            &locations,
            &[0..locations.len()],
        ) {
            Some(XFlatBuild::Table(table)) => table,
            _ => panic!("flat table"),
        };
        // The last cell with hits, its slice at the end of `hits`, and a word that reaches it: the
        // letter of digit k (index = sum of digit k * 15^k) sits k letters into the word.
        let last = (0..flat.starts.len() - 1)
            .rev()
            .find(|&i| flat.starts[i + 1] > flat.starts[i])
            .expect("a used cell");
        assert_eq!(flat.starts[last + 1] as usize, flat.hits.len());
        let scale = 15i32.pow(4);
        let letter = |digit: i32| {
            (0u8..=255)
                .find(|&b| flat.scaled_compress_table_u8[usize::from(b)] == digit * scale)
                .expect("a letter for every digit")
        };
        let word: Vec<u8> = (0..5)
            .map(|k| letter(last as i32 / 15i32.pow(k) % 15))
            .collect();
        let mut subject = Vec::new();
        for _ in 0..6 {
            subject.extend((0..7).map(|_| LETTERS[rng.below(LETTERS.len())]));
            subject.extend_from_slice(&word);
        }
        let length = subject.len();
        subject.push(0);
        assert!(flat.x_pairs_fit(&subject, length));
        let mut cells = Vec::new();
        // SAFETY: x_pairs_fit holds (asserted above).
        unsafe { flat.x_list_pairs(&subject, length, &mut cells) };
        assert!(
            cells.iter().any(|&(_, index)| index as usize == last),
            "the subject reaches it"
        );
        let (window, diag_len, diag_offset) = (40, 4096usize, 40);
        let diag_mask = diag_len as i32 - 1;
        let bf = XHitBf::new(diag_offset, length, window, 5).expect("in range");
        let mut got_diag = vec![DiagStruct::default(); diag_len];
        let mut want_diag = got_diag.clone();
        let (mut got, mut want) = (Vec::new(), Vec::new());
        // SAFETY: x_pairs_fit holds (asserted above).
        unsafe {
            x_scan_flat_stage(
                &flat,
                &subject,
                length,
                &mut got_diag,
                diag_mask,
                diag_offset,
                window,
                5,
                bf,
                &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                    stand_in(&mut got, diag_offset, q, s, last_hit, diff, entry)
                },
            )
        };
        x_scan_flat_fseed_loop(
            &flat,
            &subject,
            length,
            &mut want_diag,
            diag_mask,
            diag_offset,
            window,
            5,
            bf,
            &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                stand_in(&mut want, diag_offset, q, s, last_hit, diff, entry)
            },
        );
        assert_eq!(got, want, "candidates");
        assert!(
            got_diag
                .iter()
                .zip(&want_diag)
                .all(|(a, b)| a.raw_bits() == b.raw_bits()),
            "diagonal table"
        );
    }

    // LOSAT_X_BPFSEED's choice: on a table with the branch-free list the fused loop
    // (`x_scan_two_hit_compressed_flat`) gives the candidates and the diagonal table of
    // `x_scan_flat_bf`, on subjects where `XHitBf::fits` holds (the loop of `x_scan_flat_fseed`) and
    // at the first diagonal offset where it does not (`x_scan_flat_bf` itself).
    #[test]
    fn the_fused_loop_takes_either_loop_with_the_same_outcomes() {
        let mut rng = XorShift(0x9b05_688c_2b3e_6c1f);
        let mut query = Vec::new();
        let mut locations = Vec::new();
        for length in [700usize, 500, 800] {
            let left = query.len() as i32;
            query.extend((0..length).map(|_| LETTERS[rng.below(LETTERS.len())]));
            locations.push((left, left + length as i32 - 1));
            query.push(0);
        }
        let flat = match x_try_build_flat_compressed_lookup(
            5,
            1100,
            &query,
            &locations,
            &[0..locations.len()],
        ) {
            Some(XFlatBuild::Table(table)) => table,
            _ => panic!("flat table"),
        };
        assert!(flat.x_scan_bf_pays(), "a table with the branch-free list");
        let window = 40;
        let diag_len = 4096usize;
        let diag_mask = diag_len as i32 - 1;
        let mut total = 0usize;
        // 3,000 letters stage several times X_STAGE_CAP hits (LOSAT_X_BPFSTAGE's flushes).
        for length in [50usize, 300, 3000] {
            let mut subject: Vec<u8> = (0..length)
                .map(|_| LETTERS[rng.below(LETTERS.len())])
                .collect();
            let from = rng.below(query.len() - 50);
            subject[..40].copy_from_slice(&query[from..from + 40]);
            subject.push(0);
            if length == 3000 {
                let mut staged = 0;
                let mut cells = Vec::new();
                flat.for_each_cell::<true, _>(&subject, length, &mut cells, |_, hits| {
                    staged += hits.len();
                });
                assert!(staged > 3 * X_STAGE_CAP, "{staged} hits: several flushes");
            }
            let edge = (1 << 30) - 2 * window - length as i32;
            for (diag_offset, fits) in [(window, true), (edge - 1, true), (edge, false)] {
                assert_eq!(XHitBf::fits(diag_offset, length, window, 5), fits);
                let mut got_diag = vec![DiagStruct::default(); diag_len];
                let mut want_diag = got_diag.clone();
                let mut got = Vec::new();
                x_scan_two_hit_compressed_flat(
                    &flat,
                    &subject,
                    length,
                    &mut got_diag,
                    diag_mask,
                    diag_offset,
                    window,
                    5,
                    |q, s, last_hit, diff, entry| {
                        stand_in(&mut got, diag_offset, q, s, last_hit, diff, entry)
                    },
                );
                let mut want = Vec::new();
                x_scan_flat_bf(
                    &flat,
                    &subject,
                    length,
                    &mut want_diag,
                    diag_mask,
                    diag_offset,
                    window,
                    5,
                    &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                        stand_in(&mut want, diag_offset, q, s, last_hit, diff, entry)
                    },
                );
                assert_eq!(
                    got, want,
                    "length {length}, offset {diag_offset}: candidates"
                );
                assert!(
                    got_diag
                        .iter()
                        .zip(&want_diag)
                        .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                    "length {length}, offset {diag_offset}: diagonal table"
                );
                total += got.len();
            }
        }
        assert!(total > 0, "the test reaches the extension");
    }
}
