//! EXPERIMENT (LOSAT_X_BXDIAG / LOSAT_X_BXDIAGSHADOW): the BLASTX two-hit seed stage with the
//! subject scan and the diagonal test in one loop, and the hits of one lookup cell tested in AVX2
//! batches of eight (by default the batches cross the cells' borders: see the LOSAT_X_BXCROSS
//! paragraph below, where they rest on a conflict test instead).
//!
//! The reference (`word_finder`, and `x_word_finder_direct` of `x_seed_scan.rs`) has the
//! scanner copy every word hit into the offset-pair buffer and then reads the buffer back for the
//! two-hit test (NCBI's `s_BlastAaScanSubject` and `s_BlastAaWordFinder_TwoHit`). With the switch,
//! blastp's fused scan (`blastp::x_seed_fused::x_scan_two_hit_with`: blastp's fused loop with
//! LOSAT_X_BPDIAG's cell test) reads the hits of each lookup cell in place, in the same order, and
//! applies the outcomes of the test that touch only the hit's diagonal entry; a hit that passes the
//! test goes to the closure of `x_word_finder_fused`, which runs the rest of the reference's
//! statements unchanged (the context check, the two-hit extension, the save, the entry's update).
//!
//! Why the outcome is the reference's: blastx's diagonal entries are NCBI's `DiagStruct`, as
//! blastp's are, the test is blastp's statement for statement (aa_ungapped.c:513-552), and the
//! argument of `x_seed_fused`'s module header holds as it stands: a hit reads and writes only the
//! entry of its diagonal `(q - s) & mask`, the closure writes only that entry and the hit list, and
//! the hits are visited in the reference's order (subject offset ascending, a cell's hits in their
//! stored order). The query offsets of one cell are distinct and the table is at least the six-frame
//! query's length plus the window (`Diagonals::new`), so the hits of one subject offset fall on
//! distinct entries, which the batches of eight within a cell rest on (the batches that cross cells
//! rest on the conflict test of LOSAT_X_BXCROSS, below). The values are those of the bit-field and
//! of the wrapping i32 arithmetic, in both loops. The closure takes a candidate's context from
//! blastp's `XContextTable` (made once per query batch and thread), which gives the index of the
//! reference's binary search.
//!
//! Where it acts: the two-hit stage (`-window_size` above 0) with the standard lookup table
//! (`-word_size 3`, the default), when nothing traces the seeds; elsewhere the reference runs. `LOSAT_X_BXDIAG` is a
//! switch of the strict set: on unless it is `0`. When it
//! is not set (or empty) it acts only where the cell test is LOSAT_X_BXCROSS's (an x86-64 CPU with
//! AVX2 and POPCNT), the cell test its rule was timed with; set to a value, it runs on any CPU. With
//! a value the cell test is LOSAT_X_BXCROSS's where that runs, otherwise the masked AVX2 batch on an
//! x86-64 CPU with AVX2 and the per-hit loop elsewhere, as LOSAT_X_BPDIAG chooses;
//! `LOSAT_X_BXDIAG=per-hit` and `LOSAT_X_BXDIAG=scalar` force the per-hit loop or the scalar batch
//! (diagnostics; the tests use them). `LOSAT_X_BXDIAGSHADOW` (any value) also runs the reference on
//! a copy of the diagonal table for every subject and stops (panics) if the hit lists or the tables
//! differ, and prints a summary at exit. With `LOSAT_TIMING=1` a run prints
//! `[TIMING] blastx_diag: mode=<avx2-cross-defer|avx2-cross|avx2-defer-masked|scalar|per-hit> subjects=N cells=N
//! batched_hits=N candidates=N` (the subjects the fused loop scanned, the cells it tested in
//! batches, their hits, and the hits passed to the closure; zeros but the subjects with `per-hit`).
//!
//! The masked AVX2 batch is `XDiagMode::Avx2DeferMasked`: blastp's LOSAT_X_BPDEFER (the calls of a
//! cell's candidates after its batches) with the lanes' entry coordinates from scalar loads and the
//! last batch of fewer than eight hits as a masked AVX2 batch; the argument that the calls, their
//! order and the table are the per-hit loop's is in `x_seed_fused`'s module header. It runs only with
//! a value of `LOSAT_X_BXDIAG` or `LOSAT_X_BXDIAGSHADOW` where LOSAT_X_BXCROSS's cell test does not
//! (`LOSAT_X_BXCROSS=0`, or a CPU without POPCNT). Until 2026-11-02 it had a switch and a rule of its own (LOSAT_X_BXDEFER,
//! LOSAT_X_BXDEFERMIN at 7 hits per used cell), and those places fused tables from 9 hits per used
//! cell; no default table reached the rule, as 9 is above 7, and every CPU with AVX2 has POPCNT.
//!
//! A query batch whose lookup table has fewer hits per used backbone cell than
//! `LOSAT_X_BXDIAGMIN` (hundredths, built-in 670; `0` fuses every batch) is scanned by the direct
//! loop instead: with few hits per cell the fused loop's per-cell work costs more than the copy it
//! saves. With LOSAT_X_BXCROSS's cell test they broke even near 6.7 on 2026-10-24 (about 6.46
//! against the fly proteome, 6.95 against Swiss-Prot; the cross test against the direct loop: 4.05
//! x1.13, 5.34 x1.06, 6.33 x1.01, 6.63 x0.99, 7.34 x0.95, 8.38 x0.90), hence 670. LOSAT_X_BXSCANBYTE
//! then made the direct loop faster, and on 2026-11-01 they break even near 6.9 and 7.5 (fly 6.33
//! x1.03, 6.77 x1.01, 7.34 x0.96; Swiss-Prot 6.33 x1.09, 6.77 x1.05, 7.34 x1.01); the value was not
//! moved. No one form wins on both sides: against the fly proteome the fused loop reads x1.77 at
//! 1.76 hits per used cell and the direct loop x1.28 at 13.32; the fused loop on presence bytes with
//! one register for a short cell, and the direct loop's pairs tested by the cross-cell batches,
//! still read x1.43 and x1.37 at 1.76 (line SF S3 s4, family bx-seed-loop). Both loops give the
//! reference's hits and table, so the choice changes only the speed; it is made once per thread and
//! table build
//! (`x_fused_pays`; a table built from the pooled vectors of the one before is a new build),
//! and once for all the threads that search one table's subjects in tasks
//! (`x_shared_table_ratio`). The timing line then also gives `min_hits_per_cell=` (the cutoff),
//! `tables_fused=` and `tables_other=` (the choices made) and
//! `hits_per_cell=LEAST-MOST` (hundredths).
//!
//! The hits of a table the fused loop scans are stored in 16 bits (`BxLookup::Narrow`,
//! `BackboneCell16`: 8-byte cells and a `u16` overflow array, half the bytes of the 16-byte cells
//! and `Int4` overflow): the builder chooses the width after its count pass, from the same
//! hits-per-used-cell test as the choice of the loop (`x_narrow_cutoff`) and the length of the
//! lookup query (at most 65,535 letters, so every query offset and every cell's count fits). The
//! fused loop is generic over the table's cell, and widens each hit to an `i32` before any
//! arithmetic (`(q - s) & mask`, `q - diff`), so the tables, calls and reports are those of the
//! table of `Int4` hits. Every other table (the direct loop's, `-window_size 0`, the shadow mode,
//! a longer query batch, the compressed table of `-task blastx-fast`) keeps `Int4` hits.
//!
//! EXPERIMENT (LOSAT_X_BXCROSS): the AVX2 batch is `XCellTest::Cross`:
//! the subject's hits staged and tested in windows of eight cut across the cells' borders, a
//! window with a repeated entry in parts cut before each hit that repeats an entry of its part, the
//! last hits of the subject one after the other; the argument that the calls, their order and the
//! table are the per-hit loop's is in `x_seed_fused`'s module header (it rests on the closures of
//! `x_word_finder_fused` and `x_word_finder_fused_sent` writing only their hit's entry). The
//! stage keeps no query offset: it is recovered as `(coord + s) & mask`, exact because the table
//! made by `Diagonals::new` for the query batch has at least its length plus the window entries,
//! so `mask` is at least every query offset of the batch's lookup
//! (`x_seed_fused::x_scan_cells_cross` asserts it per subject). A
//! switch of the strict set (since 2026-10-23): on unless it is `0`. On, it is the cell test
//! of every table the fused loop takes where the CPU has AVX2 and POPCNT; `LOSAT_X_BXCROSS=0`, or a
//! CPU without POPCNT, leaves the direct loop on every table (with a value of `LOSAT_X_BXDIAG` or
//! `LOSAT_X_BXDIAGSHADOW`, the masked AVX2 batch above). `LOSAT_X_BXDIAGMIN` still chooses the tables (above), and `LOSAT_X_BXDIAGSHADOW`
//! compares it with the reference like the other cell tests. The batches run at 1,024 staged
//! hits. The timing line names `mode=avx2-cross` (`mode=avx2-cross-defer` with
//! LOSAT_X_BXCROSSDEFER, below) and ends with `cross_splits=N` (the windows of
//! eight staged hits with a repeated entry, tested in parts; counted per subject scan like the
//! other counts).
//!
//! EXPERIMENT (LOSAT_X_BXCROSSDEFER): `XCellTest::Cross { defer: true }`: LOSAT_X_BXCROSS's drain
//! defers each batch's calls to the batch that next reads one of their entries, or to the end of
//! the drain (`x_seed_fused`'s module header has the form and why the calls, their order and the
//! table are unchanged). A switch of the strict set (since 2026-10-25): on unless it is `0`;
//! `LOSAT_X_BXCROSSDEFER=0` gives LOSAT_X_BXCROSS's drain with the calls at once. It acts only
//! where LOSAT_X_BXCROSS's cell test runs; the timing line then names `mode=avx2-cross-defer`.
//!
//! LOSAT_X_BXSCANBF (`x_seed_scan.rs`; graduated 2026-10-26, no variable, always on): the fused
//! loop finds the cells with hits by tblastn's branch-free candidate loop (`x_scan_two_hit_with`
//! with `branch_free: true`, blastp's `x_scan_cells_bf`): the same cells in the same order.
//!
//! EXPERIMENT (LOSAT_X_BXUNGAPSENT, `x_ungap_sent.rs`): for a subject on that path the closure's
//! two-hit extension is the sentinel walks, in an instance of its own (`x_word_finder_fused_sent`,
//! never inlined; the other statements unchanged), so the loop without it keeps its code.

use super::x_seed_scan::XKeptPresence;
use super::x_ungap_sent::XSent;
use super::*;
use crate::algorithm::blastp::x_seed_fused::{
    x_cell_test, x_scan_two_hit_with, XCellTest, XContextTable,
};
use crate::algorithm::tblastx::lookup::{x_hits_per_cell_hundredths, XCell};
use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::OnceLock;

thread_local! {
    /// EXPERIMENT (LOSAT_X_BXDIAGMIN): the last lookup table this thread chose
    /// a loop for (its build number, `XKeptPresence::build`) and its hits per used cell (see
    /// [`x_table_ratio`]).
    static X_CHOICE: std::cell::Cell<(u64, u64)> = const { std::cell::Cell::new((0, 0)) };
    /// EXPERIMENT (LOSAT_X_BXDIAG): the context table of the query batch this thread searches,
    /// with the frame bases and the span it was made for (see [`x_with_context_table`]).
    static X_CONTEXTS: RefCell<Option<(Vec<i32>, usize, XContextTable)>> =
        const { RefCell::new(None) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:211-232
// ```c
// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
// {
//     Int4 m=0, b=0, e=0, size=0;
//
//     size = A->last_context+1;
//     ...
//     while (b < e - 1) {
//         m = (b + e) / 2;
//         if (A->contexts[m].query_offset > n)
//             e = m;
//         else
//             b = m;
//     }
//     return b;
// }
// ```
/// EXPERIMENT (LOSAT_X_BXDIAG): runs `f` with blastp's context table (`XContextTable`: a coarse
/// table and a short walk over the frame bases) for `batch`, the query batch whose lookup table
/// `kept` goes with. It gives the index the reference's binary search gives, the last context
/// whose offset is at most the query offset (`XContextTable::get`). The table is made once per
/// lookup table and kept with it (LOSAT_X_BXCTXONCE), or with that switch `=0` once per query batch
/// and thread, reused only while the batch's frame bases and span are the ones it was made for
/// (compared on every call), so a new batch gets a new table even if its contexts lie at the same
/// address.
fn x_with_context_table<R>(
    batch: &PreparedQueryBatch,
    kept: &XKeptPresence,
    f: impl FnOnce(&XContextTable) -> R,
) -> R {
    let span = batch.contexts.last().map_or(0, |c| c.offset + c.length);
    // Dispatch point (LOSAT_X_BXCTXONCE, in the strict set): the table kept with the lookup table,
    // made by the first thread that scans with it. A lookup table is built from one query batch
    // (`Lookup::new`) and scanned only with that batch, and its companion is made with it, so the
    // frame bases and the span are those of `batch` on every call; a thread that takes tasks of
    // several tables in turn (LOSAT_X_BXUNITCAP, the subject flow) then makes none again.
    if x_ctx_once() {
        let table = kept.contexts.get_or_init(|| {
            let bases: Vec<i32> = batch.contexts.iter().map(|c| c.offset as i32).collect();
            XContextTable::new(&bases, span)
        });
        return f(table);
    }
    X_CONTEXTS.with(|cell| {
        let mut cell = cell.borrow_mut();
        let fits = cell.as_ref().is_some_and(|(bases, s, _)| {
            *s == span
                && bases.len() == batch.contexts.len()
                && bases
                    .iter()
                    .zip(&batch.contexts)
                    .all(|(&b, c)| b == c.offset as i32)
        });
        if !fits {
            let bases: Vec<i32> = batch.contexts.iter().map(|c| c.offset as i32).collect();
            let table = XContextTable::new(&bases, span);
            *cell = Some((bases, span, table));
        }
        f(&cell.as_ref().expect("made above").2)
    })
}

// No NCBI counterpart: reads the LOSAT_X_BXDIAG switches once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): 0 = off, 1 = the fused loop, 2 = shadow
/// (`LOSAT_X_BXDIAGSHADOW`: the fused loop and the reference on a copy of the table, compared).
/// A switch of the strict set: on unless it is `0`; on
/// by default (not set, or empty) only with LOSAT_X_BXCROSS's cell test ([`x_cell_mode`]).
pub(crate) fn x_bxdiag_mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        let named = std::env::var("LOSAT_X_BXDIAG").is_ok_and(|v| !v.trim().is_empty());
        if std::env::var_os("LOSAT_X_BXDIAGSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_BXDIAG") && (named || x_cross_batches()) {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: chooses LOSAT's code path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): whether the cell test is LOSAT_X_BXCROSS's, where the default acts.
fn x_cross_batches() -> bool {
    match x_cell_mode() {
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { .. } => true,
        _ => false,
    }
}

// No NCBI counterpart: reads the LOSAT_X_BXCROSS switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXCROSS): whether the AVX2 batch is `XCellTest::Cross` where the CPU
/// has POPCNT (module header); a switch of the strict set, on unless it is `0`.
fn x_bx_cross() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCROSS"))
}

// No NCBI counterpart: reads the LOSAT_X_BXCROSSDEFER switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): whether `XCellTest::Cross` defers its calls (module
/// header); a switch of the strict set, on unless it is `0`.
fn x_bx_cross_defer() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCROSSDEFER"))
}

// No NCBI counterpart: reads the LOSAT_X_BXCTXONCE switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXCTXONCE): whether a query batch's context table is made once per lookup
/// table and kept with it ([`x_with_context_table`]); a switch of the strict set, on unless it is
/// `0`.
fn x_ctx_once() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCTXONCE"))
}

// No NCBI counterpart: chooses LOSAT's code path once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): the cell test of the fused loop (the switch's value names one).
fn x_cell_mode() -> XCellTest {
    static MODE: OnceLock<XCellTest> = OnceLock::new();
    *MODE.get_or_init(|| {
        x_cell_test(
            std::env::var("LOSAT_X_BXDIAG").unwrap_or_default().as_str(),
            x_bx_cross(),
            x_bx_cross_defer(),
        )
    })
}

// No NCBI counterpart: reads LOSAT_X_BXDIAGMIN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the fewest hits per used lookup cell, in hundredths, for
/// which the fused loop scans a query batch's subjects; a batch whose table has fewer is scanned by
/// the direct loop. Built-in 670, near where the two loops take the same time (module header). `0`
/// fuses every batch; a value that is not a number takes the built-in one.
fn x_min_hits_per_cell() -> u64 {
    static MIN: OnceLock<u64> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BXDIAGMIN")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(670)
    })
}

// No NCBI counterpart: chooses the width of the standard table's hits (tblastx/lookup/backbone.rs
// `BackboneCell16`); it does not change any value NCBI computes.
/// The cutoff for the width of the standard table's hits: `Some(LOSAT_X_BXDIAGMIN)` where the fused
/// loop scans the tables of this search, `None` elsewhere. The build stores the hits in 16 bits
/// when the lookup query fits and the table's hits per used cell reach the cutoff, which is the
/// test [`x_fused_pays`] makes of the finished table: a `BxLookup::Narrow` table is one the fused
/// loop scans (the fused loop is the one for `-window_size` above 0 when `x_bxdiag_mode()` is 1,
/// [`x_word_finder_diag`]); the direct loop scans the tables with fewer hits per cell, the one-hit
/// stage, and `LOSAT_X_BXDIAGSHADOW`'s (mode 2) every table, all with `Int4` hits.
pub(crate) fn x_narrow_cutoff(options: &ResolvedOptions) -> Option<u64> {
    (x_bxdiag_mode() == 1 && options.window_size > 0).then(x_min_hits_per_cell)
}

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's loop; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the table's hits per used cell, in hundredths (0 for a table
/// without hits).
fn x_hits_per_cell<K: XCell>(table: &BlastAaLookupTable<K>) -> u64 {
    let (hits, used) = K::x_hits_and_used(&table.backbone);
    x_hits_per_cell_hundredths(hits, used)
}

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's loop; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): `table`'s hits per used cell
/// ([`x_hits_per_cell`]), worked out once per table and thread (`build` = the table's build number
/// identifies it), and whether this call worked it out.
fn x_table_ratio<K: XCell>(table: &BlastAaLookupTable<K>, build: u64) -> (u64, bool) {
    X_CHOICE.with(|choice| {
        let (last, ratio) = choice.get();
        if last == build {
            return (ratio, false);
        }
        let ratio = x_hits_per_cell(table);
        choice.set((build, ratio));
        (ratio, true)
    })
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's hits and table.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): whether the fused loop scans the subjects of `table`'s query
/// batch: its hits per used cell are at least [`x_min_hits_per_cell`]. Counted for the timing line
/// once per table and thread.
fn x_fused_pays<K: XCell>(table: &BlastAaLookupTable<K>, build: u64) -> bool {
    let (ratio, first) = x_table_ratio(table, build);
    if first {
        x_count_choice(ratio);
    }
    ratio >= x_min_hits_per_cell()
}

// No NCBI counterpart: counters of the timing line.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): counts one choice of the loop for a table of `ratio` hits per
/// used cell (see [`x_fused_pays`]).
fn x_count_choice(ratio: u64) {
    let fused = ratio >= x_min_hits_per_cell();
    X_TABLES[usize::from(!fused)].fetch_add(1, Relaxed);
    X_RATIO_MIN.fetch_min(ratio, Relaxed);
    X_RATIO_MAX.fetch_max(ratio, Relaxed);
}

// No NCBI counterpart: chooses LOSAT's loop once for the threads that share a table's subjects
// (LOSAT_X_BXUNITCAP's tasks and the subject flow, blastx/x_bx_mem.rs); it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): `lookup`'s hits per used cell where the fused loop may scan its
/// subjects, worked out (and counted) on this thread, to be given to the threads that search some
/// of them ([`x_take_table_ratio`]); `None` where [`x_word_finder_diag`] does not ask for it. A
/// thread that takes tasks of several tables in turn would otherwise work it out again at each
/// change of table (its cache holds one).
#[cfg(feature = "parallel")]
pub(crate) fn x_shared_table_ratio(lookup: &BxLookup, options: &ResolvedOptions) -> Option<u64> {
    let (ratio, first) = match lookup {
        BxLookup::Wide(Lookup::Standard(table, kept)) => {
            x_table_ratio_if_fused(table, kept, options)?
        }
        BxLookup::Narrow(NarrowLookup { table, kept }) => {
            x_table_ratio_if_fused(table, kept, options)?
        }
        BxLookup::Wide(Lookup::Compressed(_)) => return None,
    };
    if first {
        x_count_choice(ratio);
    }
    Some(ratio)
}

// No NCBI counterpart: see `x_shared_table_ratio` (speed only).
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): makes `ratio` (from [`x_shared_table_ratio`] for `lookup`) this
/// thread's choice for `lookup`'s table, so that the thread does not work it out again.
#[cfg(feature = "parallel")]
pub(crate) fn x_take_table_ratio(lookup: &BxLookup, ratio: Option<u64>) {
    let kept = match lookup {
        BxLookup::Wide(Lookup::Standard(_, kept)) => kept,
        BxLookup::Narrow(NarrowLookup { kept, .. }) => kept,
        BxLookup::Wide(Lookup::Compressed(_)) => return,
    };
    if let Some(ratio) = ratio {
        X_CHOICE.with(|choice| choice.set((kept.build, ratio)));
    }
}

// No NCBI counterpart: see `x_shared_table_ratio` (speed only).
/// [`x_table_ratio`] of a standard table where the fused loop may scan its subjects, else `None`.
#[cfg(feature = "parallel")]
fn x_table_ratio_if_fused<K: XCell>(
    table: &BlastAaLookupTable<K>,
    kept: &XKeptPresence,
    options: &ResolvedOptions,
) -> Option<(u64, bool)> {
    if x_bxdiag_mode() != 1 || options.window_size <= 0 {
        return None;
    }
    Some(x_table_ratio(table, kept.build))
}

// No NCBI counterpart: counters of the timing line and of the shadow summary.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the (thread, table) choices for the fused loop (`[0]`) and
/// for the other loop (`[1]`), and the least and greatest hits per used cell (hundredths) seen.
static X_TABLES: [AtomicU64; 2] = [AtomicU64::new(0), AtomicU64::new(0)];
static X_RATIO_MIN: AtomicU64 = AtomicU64::new(u64::MAX);
static X_RATIO_MAX: AtomicU64 = AtomicU64::new(0);
static X_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static X_CELLS: AtomicU64 = AtomicU64::new(0);
static X_HITS: AtomicU64 = AtomicU64::new(0);
static X_CANDIDATES: AtomicU64 = AtomicU64::new(0);
/// EXPERIMENT (LOSAT_X_BXCROSS): the windows of eight staged hits with a repeated entry, tested in
/// parts, added with the other counts once per subject scan ([`x_add_counts`]).
static X_SPLITS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_SHORT: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: a LOSAT timing line.
/// EXPERIMENT (LOSAT_X_BXDIAG): with the switch on, the `[TIMING] blastx_diag` line of the run.
pub(crate) fn x_write_timing(out: &mut dyn std::io::Write) -> std::io::Result<()> {
    if x_bxdiag_mode() != 0 {
        write!(
            out,
            "[TIMING] blastx_diag: mode={} subjects={} cells={} batched_hits={} candidates={} min_hits_per_cell={} tables_fused={} tables_other={} hits_per_cell={}-{}",
            x_cell_mode().name(),
            X_SUBJECTS.load(Relaxed),
            X_CELLS.load(Relaxed),
            X_HITS.load(Relaxed),
            X_CANDIDATES.load(Relaxed),
            x_min_hits_per_cell(),
            X_TABLES[0].load(Relaxed),
            X_TABLES[1].load(Relaxed),
            X_RATIO_MIN.load(Relaxed).min(X_RATIO_MAX.load(Relaxed)),
            X_RATIO_MAX.load(Relaxed),
        )?;
        // EXPERIMENT (LOSAT_X_BXCROSS): the windows tested in parts.
        #[cfg(target_arch = "x86_64")]
        if matches!(x_cell_mode(), XCellTest::Cross { .. }) {
            write!(out, " cross_splits={}", X_SPLITS.load(Relaxed))?;
        }
        writeln!(out)?;
    }
    Ok(())
}

// No NCBI counterpart: prints the shadow-mode counters; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGSHADOW): one line at exit with what the shadow mode compared.
pub(crate) fn x_print_shadow_summary() {
    if x_bxdiag_mode() == 2 {
        eprintln!(
            "[X_SHADOW] BXDIAG subjects_compared={} (hit lists identical; diagonal tables identical, not compared for {} subjects shorter than a word less one letter)",
            X_SHADOW_SUBJECTS.load(Relaxed),
            X_SHADOW_SHORT.load(Relaxed)
        );
    }
}

/// EXPERIMENT (LOSAT_X_BXDIAG): the fused loop for one subject where the switch acts (module
/// header), `None` elsewhere (the one-hit stage, the compressed lookup table, a table with fewer
/// hits per used cell than LOSAT_X_BXDIAGMIN outside the shadow mode). With the shadow
/// switch the reference runs on a copy of the table and both results are compared.
///
/// A subject shorter than a word less one letter (one residue with word size 3) has one scan
/// position whose word runs past the sentinel: the reference reads that byte unchecked and the
/// fused scan reads it as a sentinel (`x_seed_fused::x_scan_cells`). With one scan position no
/// hit is extended, and the entry such a hit leaves lies a window and the subject's length behind
/// the next subject's offsets, where every later hit reads it as no previous hit; so the hit lists
/// are compared and the tables, which may differ in that entry, are not (as blastp's
/// fused loop does for the tables of subjects shorter than a word).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_word_finder_diag<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<Option<Vec<InitHsp>>> {
    let Lookup::Standard(table, kept) = lookup else {
        if x_bxdiag_mode() != 2 {
            crate::utils::xreach::x_reach!("LOSAT_X_BXDIAG", false);
        }
        return Ok(None);
    };
    if options.window_size <= 0 {
        if x_bxdiag_mode() != 2 {
            crate::utils::xreach::x_reach!("LOSAT_X_BXDIAG", false);
        }
        return Ok(None);
    }
    let mode = x_cell_mode();
    if x_bxdiag_mode() != 2 {
        // EXPERIMENT (LOSAT_X_BXDIAGMIN): a batch with few hits per cell goes to the other loop.
        let pays = x_fused_pays(table, kept.build);
        crate::utils::xreach::x_reach!("LOSAT_X_BXDIAG", pays);
        if !pays {
            return Ok(None);
        }
        let (hits, counts) = x_word_finder_fused_or_sent(
            scoring, mode, batch, parameters, options, subject, table, kept, diagonals, sent,
        )?;
        x_add_counts(counts);
        return Ok(Some(hits));
    }
    let mut copy = Diagonals {
        entries: diagonals.entries.clone(),
        offset: diagonals.offset,
        window: diagonals.window,
    };
    let reference = word_finder_any(
        scoring,
        batch,
        parameters,
        options,
        subject,
        lookup,
        &mut copy,
        |_| {},
    )?;
    let (fused, counts) = x_word_finder_fused_or_sent(
        scoring, mode, batch, parameters, options, subject, table, kept, diagonals, sent,
    )?;
    x_add_counts(counts);
    if fused != reference {
        let at = fused
            .iter()
            .zip(&reference)
            .position(|(a, b)| a != b)
            .unwrap_or(fused.len().min(reference.len()));
        panic!(
            "LOSAT_X_BXDIAGSHADOW: hit lists differ at {at} (fused {} / reference {}): fused {:?} reference {:?}",
            fused.len(),
            reference.len(),
            fused.get(at),
            reference.get(at)
        );
    }
    assert_eq!(
        copy.offset, diagonals.offset,
        "LOSAT_X_BXDIAGSHADOW: diagonal offset differs"
    );
    if subject.len() < table.word_length as usize {
        X_SHADOW_SHORT.fetch_add(1, Relaxed);
    } else if let Some(at) = copy
        .entries
        .iter()
        .zip(&diagonals.entries)
        .position(|(a, b)| a != b)
    {
        panic!(
            "LOSAT_X_BXDIAGSHADOW: diagonal entry {at} differs: fused {:#x} reference {:#x}",
            diagonals.entries[at].raw_bits(),
            copy.entries[at].raw_bits()
        );
    }
    X_SHADOW_SUBJECTS.fetch_add(1, Relaxed);
    Ok(Some(fused))
}

/// EXPERIMENT (LOSAT_X_BXDIAG): the fused loop for one subject on the table of 16-bit hits
/// (`x_word_finder_diag` for a [`NarrowLookup`]). The table was built to be fused (the build's
/// test is [`x_fused_pays`]'s on the finished table, and [`x_narrow_cutoff`] its options'), so it
/// is never handed to the direct loop; where the fused loop does not act (`-window_size 0` or
/// another `LOSAT_X_BXDIAG` mode, which no table of 16-bit hits is built for) the reference word
/// finder scans it.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606 (as in
// `x_word_finder_diag`)
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_word_finder_narrow<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &NarrowLookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<Vec<InitHsp>> {
    let NarrowLookup { table, kept } = lookup;
    if options.window_size <= 0 || x_bxdiag_mode() != 1 {
        return word_finder_narrow(
            scoring,
            batch,
            parameters,
            options,
            subject,
            lookup,
            diagonals,
            |_| {},
        );
    }
    // The choice is counted for the timing line (`tables_fused=`) as the direct loop's tables are.
    let pays = x_fused_pays(table, kept.build);
    debug_assert!(pays, "a table of 16-bit hits is built to be fused");
    crate::utils::xreach::x_reach!("LOSAT_X_BXDIAG", true);
    let (hits, counts) = x_word_finder_fused_or_sent(
        scoring,
        x_cell_mode(),
        batch,
        parameters,
        options,
        subject,
        table,
        kept,
        diagonals,
        sent,
    )?;
    x_add_counts(counts);
    Ok(hits)
}

// No NCBI counterpart: adds one subject's counts to the timing line's totals.
/// EXPERIMENT (LOSAT_X_BXDIAG): one subject's counts of [`x_word_finder_fused`] (the cells, their
/// hits, the candidates, and LOSAT_X_BXCROSS's windows tested in parts) added to the totals.
fn x_add_counts(counts: [u64; 4]) {
    X_SUBJECTS.fetch_add(1, Relaxed);
    if counts[0] != 0 {
        X_CELLS.fetch_add(counts[0], Relaxed);
        X_HITS.fetch_add(counts[1], Relaxed);
        X_CANDIDATES.fetch_add(counts[2], Relaxed);
        if counts[3] != 0 {
            X_SPLITS.fetch_add(counts[3], Relaxed);
        }
    }
}

// No NCBI counterpart: chooses the instance of the fused loop; both give the reference's hits and
// table.
/// EXPERIMENT (LOSAT_X_BXUNGAPSENT): [`x_word_finder_fused`], or with `sent` its instance with the
/// sentinel walks, [`x_word_finder_fused_sent`] (with the keys of the matrix of `sent`).
#[allow(clippy::too_many_arguments)]
fn x_word_finder_fused_or_sent<K: XCell, S: BxScoring>(
    scoring: S,
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable<K>,
    kept: &XKeptPresence,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<(Vec<InitHsp>, [u64; 4])> {
    match sent {
        Some(sent) => x_word_finder_fused_sent(
            mode, batch, parameters, options, subject, table, kept, diagonals, sent,
        ),
        None => x_word_finder_fused(
            scoring, mode, batch, parameters, options, subject, table, kept, diagonals,
        ),
    }
}

/// EXPERIMENT (LOSAT_X_BXDIAG): `x_word_finder_direct`'s two-hit stage for one subject on the
/// standard lookup table `table`, with the scan and the diagonal test fused (cell test `mode`;
/// the branch-free scan, BXSCANBF, graduated); also returns the cells tested in
/// batches, their hits, the candidates and LOSAT_X_BXCROSS's windows tested in parts.
/// `subject` starts at the first residue and keeps the sentinel after the last one.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
// ```c
//         wordsize = lookup->word_length;
//     }
//
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
// ```
#[allow(clippy::too_many_arguments)]
fn x_word_finder_fused<K: XCell, S: BxScoring>(
    scoring: S,
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable<K>,
    kept: &XKeptPresence,
    diagonals: &mut Diagonals,
) -> Result<(Vec<InitHsp>, [u64; 4])> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let offset = diagonals.offset;
    let mask = diagonals.entries.len() as i32 - 1;
    // Lookup::scan's set-up, with length = subject.len() - 1 as word_finder passes it; the
    // fused scan runs from 0 to max(length - word length, 0), the reference's range with a
    // two-hit window.
    let length = subject.len().saturating_sub(1);
    let mut hits = Vec::new();
    let counts = x_with_context_table(batch, kept, |contexts| {
        x_scan_two_hit_with(
            mode,
            true,
            table,
            subject,
            length,
            &mut diagonals.entries,
            mask,
            offset,
            options.window_size,
            word,
            |q, s, previous, diff, entry| {
                // From here to the end of the closure: the statements of x_word_finder_direct after
                // the diagonal tests (`previous` = last_hit - diag_offset, `diff` = s - previous).
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560
                // ```c
                //                 curr_context = BSearchContextInfo(query_offset, query_info);
                // ```
                // The table's index is the binary search's (x_with_context_table); an offset outside
                // the table (none from the lookup) takes the binary search.
                let context = contexts.get(q).unwrap_or_else(|| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                });
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
                // ```c
                //                 if (query_offset - diff <
                //                     query_info->contexts[curr_context].query_offset) {
                //
                //                     /* there was no last hit for this diagnol; start a new hit */
                //                     diag_array[diag_coord].last_hit =
                //                         subject_offset + diag_offset;
                //                     continue;
                //                 }
                // ```
                if q - diff < batch.contexts[context].offset as i32 {
                    entry.set_last_hit(s.wrapping_add(offset));
                    return;
                }
                let p = &parameters[context];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                // ```c
                //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
                //                                               last_hit + wordsize,
                //                                               subject_offset, query_offset,
                //                                               cutoffs->x_dropoff,
                //                                               &hsp_q, &hsp_s,
                //                                               &hsp_len, use_pssm,
                //                                               wordsize, &right_extend,
                // ```
                // The matrix's walks: BLOSUM62's `extend_two_hit_blosum62`, or another matrix's
                // (`bx_two_hit`).
                let Some(result) = bx_two_hit(
                    scoring,
                    query,
                    subject,
                    (previous + word) as usize,
                    s as usize,
                    q as usize,
                    p.word_xdrop,
                    word as usize,
                ) else {
                    return;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
                // ```c
                //                 if (score >= cutoffs->cutoff_score)
                //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                                      query_offset, subject_offset, hsp_len,
                //                                      score);
                //
                //                 /* If an extension to the right happened, reset the last hit
                //                    so that future hits to this diagonal must start over. */
                //
                //                 if (right_extend) {
                //                     diag_array[diag_coord].flag = 1;
                //                     diag_array[diag_coord].last_hit =
                //                         s_last_off - (wordsize - 1) + diag_offset;
                //                 }
                //                 /* Otherwise, make the present hit into the previous hit for
                //                    this diagonal */
                //                 else {
                //                     diag_array[diag_coord].last_hit =
                //                         subject_offset + diag_offset;
                //                 }
                // ```
                let u = result.ungapped_data;
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                }
                if result.right_extend {
                    entry.set_flag(1);
                    entry.set_last_hit((result.s_last_off - (word - 1)).wrapping_add(offset));
                } else {
                    entry.set_last_hit(s.wrapping_add(offset));
                }
            },
        )
    });
    // The rest of word_finder, unchanged.
    diagonals.finish_subject(length)?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok((hits, counts))
}

/// EXPERIMENT (LOSAT_X_BXUNGAPSENT): [`x_word_finder_fused`] with the two-hit extension of `sent`
/// (`x_ungap_sent.rs`), in a function of its own so that the loop without it keeps its code. The
/// closure's other statements are those of `x_word_finder_fused`'s, in the same order.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-606 (as in
// `x_word_finder_fused`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_fused_sent<K: XCell, S: BxScoring>(
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable<K>,
    kept: &XKeptPresence,
    diagonals: &mut Diagonals,
    sent: XSent<'_, S>,
) -> Result<(Vec<InitHsp>, [u64; 4])> {
    let word = options.word_size;
    let offset = diagonals.offset;
    let mask = diagonals.entries.len() as i32 - 1;
    let length = subject.len().saturating_sub(1);
    let mut hits = Vec::new();
    let counts = x_with_context_table(batch, kept, |contexts| {
        // The candidate closure with the extension `two_hit::<$shadow>`.
        macro_rules! scan {
            ($shadow:literal) => {
                x_scan_two_hit_with(
                    mode,
                    true,
                    table,
                    subject,
                    length,
                    &mut diagonals.entries,
                    mask,
                    offset,
                    options.window_size,
                    word,
                    |q, s, previous, diff, entry| {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560-573
                        // (as in `x_word_finder_fused`)
                        let context = contexts.get(q).unwrap_or_else(|| {
                            batch
                                .contexts
                                .partition_point(|c| c.offset as i32 <= q)
                                .saturating_sub(1)
                        });
                        if q - diff < batch.contexts[context].offset as i32 {
                            entry.set_last_hit(s.wrapping_add(offset));
                            return;
                        }
                        let p = &parameters[context];
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                        // (`x_ungap_sent.rs`)
                        let Some(result) = sent.two_hit::<$shadow>(
                            (previous + word) as usize,
                            s as usize,
                            q as usize,
                            p.word_xdrop,
                        ) else {
                            return;
                        };
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
                        // (as in `x_word_finder_fused`)
                        let u = result.ungapped_data;
                        if u.score >= p.word_cutoff {
                            hits.push(InitHsp {
                                q_seed: q,
                                s_seed: s,
                                q_start: u.q_start,
                                s_start: u.s_start,
                                length: u.length,
                                score: u.score,
                            });
                        }
                        if result.right_extend {
                            entry.set_flag(1);
                            entry.set_last_hit(
                                (result.s_last_off - (word - 1)).wrapping_add(offset),
                            );
                        } else {
                            entry.set_last_hit(s.wrapping_add(offset));
                        }
                    },
                )
            };
        }
        if sent.shadow() {
            scan!(true)
        } else {
            scan!(false)
        }
    });
    // The rest of word_finder, unchanged.
    diagonals.finish_subject(length)?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok((hits, counts))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::x_seed_fused::XDiagMode;
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    /// The cell tests this machine can run (the masked AVX2 batch, and with LOSAT_X_BXCROSS, its
    /// calls deferred with LOSAT_X_BXCROSSDEFER and not), each with the branch-free scan (BXSCANBF,
    /// graduated).
    fn modes() -> Vec<XCellTest> {
        let mut modes = vec![
            XCellTest::Diag(XDiagMode::PerHit),
            XCellTest::Diag(XDiagMode::Scalar),
        ];
        for (cross, defer) in [(false, false), (true, false), (true, true)] {
            let auto = x_cell_test("", cross, defer);
            if auto != XCellTest::Diag(XDiagMode::PerHit) {
                modes.push(auto);
            }
        }
        modes
    }

    /// A search of `query` (nucleotide FASTA) with `-window_size window` (default options
    /// otherwise), its statistics set for subjects of `subject_length` residues; its standard table
    /// has `Int4` hits, or with `narrow` 16-bit hits.
    fn search(
        query: &[u8],
        window: i32,
        subject_length: usize,
        narrow: bool,
    ) -> (
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        BxLookup,
    ) {
        search_with(query, window, subject_length, narrow, "BLOSUM62", &[])
    }

    /// [`search`] with `-matrix matrix` (its default gap costs and threshold, its score block) and
    /// the options `extra`; another matrix than BLOSUM62 takes its build (`Lookup::x_new_of`, item
    /// SC6b).
    fn search_with(
        query: &[u8],
        window: i32,
        subject_length: usize,
        narrow: bool,
        matrix: &str,
        extra: &[&str],
    ) -> (
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        BxLookup,
    ) {
        let q = parse_fasta(query, false, false).unwrap();
        let window = window.to_string();
        let mut args = vec![
            "losat",
            "blastx",
            "-query",
            "query",
            "-subject",
            "subject",
            "-window_size",
            &window,
            "-matrix",
            matrix,
        ];
        args.extend_from_slice(extra);
        let crate::cli::Commands::Blastx(args) =
            crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                .unwrap()
                .command
        else {
            unreachable!()
        };
        let options = args.resolve().unwrap();
        let batch = prepare_queries(&q, &options).unwrap();
        let mut params = parameters::score_block_for(&batch, options.scoring);
        parameters::effective_lengths(
            &batch,
            &mut params,
            &options,
            subject_length as i64,
            1,
            None,
        )
        .unwrap();
        let lookup = Lookup::x_new_for(&batch, &params, &options, 1, narrow.then_some(0)).unwrap();
        assert_eq!(matches!(lookup, BxLookup::Narrow(_)), narrow);
        parameters::subject_parameters(
            &batch,
            &mut params,
            &options,
            subject_length as i64,
            1,
            subject_length,
            None,
        )
        .unwrap();
        (batch, params, options, lookup)
    }

    /// The subject records of a protein FASTA, encoded with the sentinel after the last residue.
    fn subjects(fasta: &[u8]) -> Vec<Vec<u8>> {
        parse_fasta(fasta, true, false)
            .unwrap()
            .iter()
            .map(|record| {
                let mut subject = record
                    .sequence
                    .iter()
                    .copied()
                    .map(crate::utils::matrix::aa_char_to_ncbistdaa)
                    .collect::<Vec<_>>();
                subject.push(0);
                subject
            })
            .collect()
    }

    /// Runs every subject through the reference and through the fused loop (cell test `mode`;
    /// the branch-free scan), each on its own table starting from
    /// `start(table)`, and compares the hit lists and the tables after each subject; returns the
    /// hits and the counts of the fused loop.
    fn compare(
        search: &(
            PreparedQueryBatch,
            Vec<ContextParameters>,
            ResolvedOptions,
            BxLookup,
        ),
        subjects: &[Vec<u8>],
        mode: XCellTest,
        start: impl Fn(&mut Diagonals),
    ) -> (usize, [u64; 4]) {
        match &search.3 {
            BxLookup::Wide(Lookup::Standard(table, kept)) => {
                compare_on(search, table, kept, subjects, mode, start)
            }
            BxLookup::Narrow(NarrowLookup { table, kept }) => {
                compare_on(search, table, kept, subjects, mode, start)
            }
            BxLookup::Wide(Lookup::Compressed(_)) => panic!("the standard lookup table"),
        }
    }

    /// [`compare`] for the standard table `table` of the search's lookup, of either width.
    fn compare_on<K: XCell>(
        (batch, params, options, lookup): &(
            PreparedQueryBatch,
            Vec<ContextParameters>,
            ResolvedOptions,
            BxLookup,
        ),
        table: &BlastAaLookupTable<K>,
        kept: &XKeptPresence,
        subjects: &[Vec<u8>],
        mode: XCellTest,
        start: impl Fn(&mut Diagonals),
    ) -> (usize, [u64; 4]) {
        let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        let mut reference = Diagonals::new(length, options.window_size).unwrap();
        let mut fused = Diagonals::new(length, options.window_size).unwrap();
        // EXPERIMENT (LOSAT_X_BXUNGAPSENT): the instances with the sentinel walks, without and
        // with the shadow's comparison, each on a table of its own.
        let mut sent = [
            Diagonals::new(length, options.window_size).unwrap(),
            Diagonals::new(length, options.window_size).unwrap(),
        ];
        start(&mut reference);
        start(&mut fused);
        sent.iter_mut().for_each(&start);
        let (mut hits, mut counts) = (0, [0u64; 4]);
        let matrix = options.scoring;
        for (i, subject) in subjects.iter().enumerate() {
            // BLOSUM62's reference and fused loop, or the instances of another matrix (item SC6b).
            let (expected, (actual, c)) = if matrix == ScoringMatrix::Blosum62 {
                (
                    word_finder_bx(
                        batch,
                        params,
                        options,
                        subject,
                        lookup,
                        &mut reference,
                        |_| {},
                    )
                    .unwrap(),
                    x_word_finder_fused(
                        Blosum62Bx, mode, batch, params, options, subject, table, kept, &mut fused,
                    )
                    .unwrap(),
                )
            } else {
                let scoring = MatrixBx::of(matrix);
                (
                    word_finder_bx_of(
                        scoring,
                        batch,
                        params,
                        options,
                        subject,
                        lookup,
                        &mut reference,
                        |_| {},
                    )
                    .unwrap(),
                    x_word_finder_fused(
                        scoring, mode, batch, params, options, subject, table, kept, &mut fused,
                    )
                    .unwrap(),
                )
            };
            let run = format!("{matrix:?} {mode:?} subject {i}");
            assert_eq!(actual, expected, "{run}");
            assert_eq!(fused.offset, reference.offset, "{run}");
            assert!(fused.entries == reference.entries, "{run}");
            let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
            for (table_sent, shadow) in sent.iter_mut().zip([false, true]) {
                let (got, c2) = if matrix == ScoringMatrix::Blosum62 {
                    let x = XSent::for_test(Blosum62Bx, &batch.sequence_start, &lead, shadow);
                    x_word_finder_fused_sent(
                        mode, batch, params, options, subject, table, kept, table_sent, x,
                    )
                } else {
                    let x =
                        XSent::for_test(MatrixBx::of(matrix), &batch.sequence_start, &lead, shadow);
                    x_word_finder_fused_sent(
                        mode, batch, params, options, subject, table, kept, table_sent, x,
                    )
                }
                .unwrap();
                assert_eq!(got, expected, "{run} sent shadow={shadow}");
                assert_eq!(c2, c, "{run} sent shadow={shadow}");
                assert_eq!(table_sent.offset, reference.offset, "{run} sent");
                assert!(table_sent.entries == reference.entries, "{run} sent");
            }
            hits += actual.len();
            for (total, n) in counts.iter_mut().zip(c) {
                *total += n;
            }
        }
        (hits, counts)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // The session F fixture (the frozen NCBI INIT rows of tests/blastx_f_diagonal) at its two-hit
    // windows, from a fresh table, a table near INT4_MAX / 4, a reset one and one whose offset is
    // just below 2^30 (the 31-bit `last_hit` wraps): the fused loop gives the reference's hits and
    // table with every cell test.
    #[test]
    fn the_fused_loop_matches_the_reference_on_the_session_f_fixture() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let query = std::fs::read(dir.join("query.fna")).unwrap();
        let subject = subjects(&std::fs::read(dir.join("subject.faa")).unwrap());
        let length = subject[0].len() - 1;
        for (window, narrow) in [1, 39, 40, 41]
            .into_iter()
            .flat_map(|window| [(window, false), (window, true)])
        {
            let s = search(&query, window, length, narrow);
            let flagged = |last_hit: i32| {
                let mut entry = DiagStruct::default();
                entry.set_last_hit(last_hit);
                entry.set_flag(1);
                entry
            };
            for mode in modes() {
                let repeated = [subject[0].clone(), subject[0].clone(), subject[0].clone()];
                compare(&s, &repeated, mode, |_| {});
                compare(&s, &repeated, mode, |d| {
                    d.offset = 536_870_910;
                    d.entries.fill(flagged(123));
                });
                compare(&s, &repeated, mode, |d| {
                    d.offset = 536_870_911;
                    d.entries.fill(flagged(536_870_900));
                    d.finish_subject(83).unwrap();
                });
                compare(&s, &repeated, mode, |d| {
                    d.offset = (1 << 30) - 64;
                    d.entries.fill(DiagStruct::clear(window));
                });
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:211-232
    // ```c
    // Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
    // {
    //     Int4 m=0, b=0, e=0, size=0;
    //
    //     size = A->last_context+1;
    //     ...
    //     while (b < e - 1) {
    //         m = (b + e) / 2;
    //         if (A->contexts[m].query_offset > n)
    //             e = m;
    //         else
    //             b = m;
    //     }
    //     return b;
    // }
    // ```
    // The context table of a batch gives the binary search's index at every query offset, and a
    // second batch (other queries) gets its own table.
    #[test]
    fn the_context_table_gives_the_binary_search_index_for_every_offset() {
        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let genome = parse_fasta(
            &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
            false,
            false,
        )
        .unwrap();
        let letters = &genome[0].sequence;
        for pieces in [[1_500, 2_200, 900, 3_100], [4_000, 333, 2_500, 1_200]] {
            let mut query = Vec::new();
            for (i, &n) in pieces.iter().enumerate() {
                query.extend_from_slice(format!(">p{i}\n").as_bytes());
                query.extend_from_slice(&letters[i * 50_000..i * 50_000 + n]);
                query.push(b'\n');
            }
            let (batch, _, _, lookup) = search(&query, 40, 500, false);
            let BxLookup::Wide(Lookup::Standard(_, kept)) = &lookup else {
                panic!("a standard table");
            };
            let span = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            assert!(batch.contexts.len() >= 24);
            x_with_context_table(&batch, kept, |table| {
                for q in 0..span as i32 {
                    let expected = batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1);
                    assert_eq!(table.get(q), Some(expected), "q={q}");
                }
                assert_eq!(table.get(span as i32), None);
            });
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    //             /* If the reset bit is set, an extension just happened. */
    //             if (diag_array[diag_coord].flag) {
    //                 /* If we've already extended past this hit, skip it. */
    //                 if ((Int4) (subject_offset + diag_offset) <
    //                     diag_array[diag_coord].last_hit) {
    //                     continue;
    //                 }
    // ```
    // 6 kb of the AP027131 genome against the first 60 of its proteins, one table for all of
    // them: many lookup cells with eight hits or more (so the batches run, not only the scalar
    // rest of a cell), extensions with and without a right extension, and the entries a subject
    // leaves for the next. Windows 40 (the default) and 400, and a run starting just below 2^30.
    #[test]
    fn the_fused_loop_matches_the_reference_on_many_cells_and_subjects() {
        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let genome = parse_fasta(
            &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
            false,
            false,
        )
        .unwrap();
        let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
        let query = format!(">piece\n{piece}\n").into_bytes();
        let all = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
        let subjects = &all[..60];
        let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
        for (window, narrow) in [(40, false), (40, true), (400, false), (400, true)] {
            let s = search(&query, window, mean, narrow);
            // The lookup has cells of a full batch (eight hits) and more.
            let cells: Vec<usize> = match &s.3 {
                BxLookup::Wide(Lookup::Standard(table, _)) => {
                    table.backbone.iter().map(|c| c.count()).collect()
                }
                BxLookup::Narrow(NarrowLookup { table, .. }) => {
                    table.backbone.iter().map(|c| c.count()).collect()
                }
                BxLookup::Wide(Lookup::Compressed(_)) => panic!("the standard lookup table"),
            };
            assert!(cells.iter().any(|&n| n >= 16));
            for mode in modes() {
                let (hits, counts) = compare(&s, subjects, mode, |_| {});
                assert!(hits > 0, "window={window} {mode:?}: no hit");
                if mode == XCellTest::Diag(XDiagMode::PerHit) {
                    assert_eq!(counts, [0, 0, 0, 0]);
                } else {
                    assert!(counts[0] > 0 && counts[2] > 0, "{mode:?} {counts:?}");
                    assert!(counts[1] > counts[0], "{mode:?} {counts:?}");
                }
                compare(&s, subjects, mode, |d| {
                    d.offset = (1 << 30) - 3_000;
                    d.entries.fill(DiagStruct::clear(window));
                });
            }
        }
    }

    /// The 6 kb of AP027131 as a query and its first 40 proteins as subjects (with their mean
    /// length) of the tests of the other matrices.
    fn other_matrix_inputs() -> (Vec<u8>, Vec<Vec<u8>>, usize) {
        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let genome = parse_fasta(
            &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
            false,
            false,
        )
        .unwrap();
        let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
        let query = format!(">piece\n{piece}\n").into_bytes();
        let mut subjects = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
        subjects.truncate(40);
        let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
        (query, subjects, mean)
    }

    /// BLASTX's other matrices, gapped, and IDENTITY, which an ungapped search takes with word size
    /// 3 (its byte table has keys for the sentinel walks).
    const OTHER_MATRICES: [(&str, &[&str]); 8] = [
        ("BLOSUM45", &[]),
        ("BLOSUM50", &[]),
        ("BLOSUM80", &[]),
        ("BLOSUM90", &[]),
        ("PAM30", &[]),
        ("PAM70", &[]),
        ("PAM250", &[]),
        ("IDENTITY", &["-ungapped", "-comp_based_stats", "0"]),
    ];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-606 (as above), with
    // `matrix->data` of the search's matrix (`s_BlastAaExtendTwoHit`).
    // The test above with every other matrix of BLASTX's table and IDENTITY (item SC6b): their
    // tables of either width, the fused loop with their walks (`bx_two_hit`) and the sentinel walks
    // with their keys give the reference's hits and diagonal table.
    #[test]
    fn the_fused_loop_matches_the_reference_with_other_matrices() {
        let (query, subjects, mean) = other_matrix_inputs();
        for (matrix, extra) in OTHER_MATRICES {
            for narrow in [false, true] {
                let s = search_with(&query, 40, mean, narrow, matrix, extra);
                assert_ne!(s.2.scoring, ScoringMatrix::Blosum62);
                for mode in modes() {
                    let (hits, counts) = compare(&s, &subjects, mode, |_| {});
                    assert!(hits > 0, "{matrix} {mode:?}: no hit");
                    if mode != XCellTest::Diag(XDiagMode::PerHit) {
                        assert!(counts[2] > 0, "{matrix} {mode:?} {counts:?}");
                    }
                }
            }
        }
    }

    // No NCBI counterpart (speed only): the table of 16-bit hits is the table of 32-bit hits of the
    // same query in every property that chooses the loop (its hits per used cell are the build's
    // choice of the width, `x_narrow_cutoff`), and the entry point of the seed stage scans it, with
    // the fused loop or (where that does not act) the reference scan, to the reference's hits and
    // diagonal table.
    #[test]
    fn a_table_of_16_bit_hits_is_scanned_to_the_reference_hits_by_the_seed_stage() {
        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let genome = parse_fasta(
            &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
            false,
            false,
        )
        .unwrap();
        let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
        let query = format!(">piece\n{piece}\n").into_bytes();
        let all = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
        let subjects = &all[..40];
        let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
        let wide = search(&query, 40, mean, false);
        let narrow = search(&query, 40, mean, true);
        let (
            BxLookup::Wide(Lookup::Standard(wide_table, _)),
            BxLookup::Narrow(NarrowLookup {
                table: narrow_table,
                ..
            }),
        ) = (&wide.3, &narrow.3)
        else {
            panic!("the tables of the two widths");
        };
        let ratio = x_hits_per_cell(wide_table);
        assert!(ratio > 0);
        assert_eq!(x_hits_per_cell(narrow_table), ratio);
        // The build's choice at the table's own ratio and one hundredth above it.
        for (cutoff, narrower) in [(ratio, true), (ratio + 1, false)] {
            let (batch, params, options, _) = &narrow;
            let query_batch = Lookup::x_new_as(batch, params, options, 1, Some(cutoff)).unwrap();
            assert_eq!(matches!(query_batch, BxLookup::Narrow(_)), narrower);
        }
        let (batch, params, options, lookup) = &narrow;
        let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        let mut reference = Diagonals::new(length, options.window_size).unwrap();
        let mut scanned = Diagonals::new(length, options.window_size).unwrap();
        let mut hits = 0;
        for subject in subjects {
            let expected = word_finder_bx(
                batch,
                params,
                options,
                subject,
                lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            let actual = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                Blosum62Bx,
                batch,
                params,
                options,
                subject,
                lookup,
                &mut scanned,
                None,
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert!(scanned.offset == reference.offset && scanned.entries == reference.entries);
            hits += actual.len();
        }
        assert!(hits > 0);
    }

    // No NCBI counterpart (speed only): the entry point of another matrix's seed stage
    // (`x_word_finder_scan` with `MatrixBx`, item SC6b) scans its tables as the search builds them
    // (16-bit hits where `x_narrow_cutoff` lets the fused loop take them), with and without the
    // sentinel walks, and with `-window_size 0` (the one-hit stage, the direct loop), to the
    // reference's hits and diagonal table.
    #[test]
    fn another_matrix_s_seed_stage_gives_the_reference_hits() {
        let (query, subjects, mean) = other_matrix_inputs();
        let mut widths = Vec::new();
        for (matrix, extra) in OTHER_MATRICES {
            for window in [40, 0] {
                let (batch, params, options, _) =
                    search_with(&query, window, mean, false, matrix, extra);
                let lookup =
                    Lookup::x_new_of(&batch, &params, &options, 1, x_narrow_cutoff(&options))
                        .unwrap();
                if window > 0 {
                    widths.push((matrix, matches!(lookup, BxLookup::Narrow(_))));
                } else {
                    assert!(
                        matches!(lookup, BxLookup::Wide(_)),
                        "{matrix}: one-hit table"
                    );
                }
                let scoring = MatrixBx::of(options.scoring);
                let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
                // The sentinel walks only with a two-hit window (`query_fits_of`).
                for with_sent in [false, window > 0] {
                    let mut reference = Diagonals::new(length, options.window_size).unwrap();
                    let mut scanned = Diagonals::new(length, options.window_size).unwrap();
                    let mut hits = 0;
                    for subject in &subjects {
                        let expected = word_finder_bx_of(
                            scoring,
                            &batch,
                            &params,
                            &options,
                            subject,
                            &lookup,
                            &mut reference,
                            |_| {},
                        )
                        .unwrap();
                        let lead: Vec<u8> =
                            std::iter::once(0).chain(subject.iter().copied()).collect();
                        let sent = with_sent
                            .then(|| XSent::for_test(scoring, &batch.sequence_start, &lead, false));
                        let actual =
                            crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                                scoring,
                                &batch,
                                &params,
                                &options,
                                subject,
                                &lookup,
                                &mut scanned,
                                sent,
                            )
                            .unwrap();
                        let run = format!("{matrix} window={window} sent={with_sent}");
                        assert_eq!(actual, expected, "{run}");
                        assert!(
                            scanned.offset == reference.offset
                                && scanned.entries == reference.entries,
                            "{run}"
                        );
                        hits += actual.len();
                    }
                    assert!(hits > 0, "{matrix} window={window}");
                }
            }
        }
        // Where the fused loop acts (`x_bxdiag_mode`: AVX2 and POPCNT, or a value of
        // LOSAT_X_BXDIAG), some of these tables take 16-bit hits and the others the direct loop.
        if x_bxdiag_mode() == 1 {
            assert!(widths.iter().any(|w| w.1), "{widths:?}");
            assert!(widths.iter().any(|w| !w.1), "{widths:?}");
        }
        eprintln!("16-bit tables: {widths:?}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-606 (as above), with the
    // word length of the search's lookup table (`wordsize = lookup->word_length`).
    // No NCBI counterpart (speed only): the seed stage of the word sizes other than 3 and 5 (item SC9)
    // scans the tables the search builds (`Lookup::x_new_for` with `x_narrow_cutoff`: for 2 and 4
    // the lookup table, of 16-bit hits where the fused loop takes it, and for 6 and 7 the
    // compressed table), with BLOSUM62 and another matrix, two-hit and one-hit, to the reference's
    // hits and diagonal table.
    #[test]
    fn every_word_size_s_seed_stage_gives_the_reference_hits() {
        fn seed_stage<S: BxScoring>(
            scoring: S,
            search: &(PreparedQueryBatch, Vec<ContextParameters>, ResolvedOptions),
            lookup: &BxLookup,
            subjects: &[Vec<u8>],
            run: &str,
        ) -> usize {
            let (batch, params, options) = search;
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let mut reference = Diagonals::new(length, options.window_size).unwrap();
            let mut scanned = Diagonals::new(length, options.window_size).unwrap();
            let mut hits = 0;
            for subject in subjects {
                let expected = word_finder_bx_of(
                    scoring,
                    batch,
                    params,
                    options,
                    subject,
                    lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap();
                let actual = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                    scoring,
                    batch,
                    params,
                    options,
                    subject,
                    lookup,
                    &mut scanned,
                    None,
                )
                .unwrap();
                assert_eq!(actual, expected, "{run}");
                assert!(
                    scanned.offset == reference.offset && scanned.entries == reference.entries,
                    "{run}"
                );
                hits += actual.len();
            }
            hits
        }
        let (query, subjects, mean) = other_matrix_inputs();
        let mut tables = Vec::new();
        // Each word size with BLOSUM62 and PAM30, two-hit or one-hit (a table of 4 letters
        // takes seconds to build in a debug build).
        for (matrix, word, window) in [
            ("BLOSUM62", "2", 40),
            ("BLOSUM62", "2", 0),
            ("PAM30", "2", 40),
            ("BLOSUM62", "4", 40),
            ("PAM30", "4", 0),
            ("BLOSUM62", "6", 40),
            ("PAM30", "6", 0),
            ("BLOSUM62", "7", 0),
            ("PAM30", "7", 40),
        ] {
            let (batch, params, options, _) =
                search_with(&query, window, mean, false, matrix, &["-word_size", word]);
            let lookup =
                Lookup::x_new_for(&batch, &params, &options, 1, x_narrow_cutoff(&options)).unwrap();
            let kind = match &lookup {
                BxLookup::Narrow(l) => ("16-bit", l.table.word_length),
                BxLookup::Wide(Lookup::Standard(l, _)) => ("standard", l.word_length),
                BxLookup::Wide(Lookup::Compressed(l)) => {
                    ("compressed", x_compressed_word_length(l))
                }
            };
            assert_eq!(kind.1.to_string(), word, "{matrix} {word}");
            assert_eq!(kind.0 == "compressed", options.compressed_lookup());
            let run = format!("{matrix} -word_size {word} window={window} {}", kind.0);
            let search = (batch, params, options);
            let hits = if matrix == "BLOSUM62" {
                seed_stage(Blosum62Bx, &search, &lookup, &subjects, &run)
            } else {
                seed_stage(
                    MatrixBx::of(search.2.scoring),
                    &search,
                    &lookup,
                    &subjects,
                    &run,
                )
            };
            assert!(hits > 0, "{run}");
            tables.push(run);
        }
        // Where the fused loop acts, the two-hit tables of 2 letters take 16-bit hits.
        if x_bxdiag_mode() == 1 {
            assert!(tables.iter().any(|t| t.contains("16-bit")), "{tables:?}");
        }
        eprintln!("tables: {tables:?}");
    }

    // No NCBI counterpart (memory only): a search of two threads or more builds its tables of 2 and 3
    // letters pooled and gives their vectors back when they are dropped; a table of 4 letters (the
    // joint builder, which takes no kept vectors, item SC9) is not marked pooled, so that the pool
    // does not keep every dropped table of the search.
    #[test]
    fn a_table_of_four_letters_gives_no_vectors_back() {
        let (query, _, mean) = other_matrix_inputs();
        for (matrix, word) in [
            ("BLOSUM62", "2"),
            ("BLOSUM62", "4"),
            ("PAM30", "3"),
            ("PAM30", "4"),
        ] {
            let extra = ["-word_size", word, "-num_threads", "8"];
            let (batch, params, options, _) = search_with(&query, 40, mean, false, matrix, &extra);
            let lookup = Lookup::x_new_for(&batch, &params, &options, 1, None).unwrap();
            let BxLookup::Wide(Lookup::Standard(table, kept)) = &lookup else {
                panic!("{matrix} {word}: the standard table of Int4 hits");
            };
            assert_eq!(table.word_length.to_string(), word);
            let pooled = crate::algorithm::tblastx::lookup::x_bx_pooled_for(8);
            assert_eq!(kept.pooled, pooled && word != "4", "{matrix} {word}");
        }
    }

    // No NCBI counterpart (speed only): the choice of the loop names a table by its build number,
    // so a table at the address and with the vectors of the one before it is worked out again.
    #[test]
    fn the_loop_choice_is_worked_out_again_for_each_table_build() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let query = std::fs::read(dir.join("query.fna")).unwrap();
        let (_, _, _, lookup) = search(&query, 40, 500, false);
        let BxLookup::Wide(Lookup::Standard(table, kept)) = &lookup else {
            panic!("a standard table");
        };
        let next = crate::algorithm::blastx::seed::x_seed_scan::XKeptPresence::new();
        assert!(kept.build > 0 && next.build > kept.build);
        let (ratio, first) = x_table_ratio(table, kept.build);
        assert!(first);
        assert_eq!(x_table_ratio(table, kept.build), (ratio, false));
        // The same table under the next build number, as a pooled rebuild would give it.
        assert_eq!(x_table_ratio(table, next.build), (ratio, true));
        assert_eq!(x_table_ratio(table, next.build), (ratio, false));
    }
}
