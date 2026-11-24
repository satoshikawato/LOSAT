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
//! (`-word_size` 3, the default, and 2 and 4); elsewhere the reference runs. `LOSAT_X_BXDIAG` is a
//! switch of the strict set: on unless it is `0`. When it
//! is not set (or empty) it acts only where the cell test is LOSAT_X_BXCROSS's (an x86-64 CPU with
//! AVX2 and POPCNT), the cell test its rule was timed with, or LOSAT_X_BXCROSS128's (the same
//! batches in 128 bits: the WebAssembly build, and x86-64 with SSSE3, SSE4.1 and POPCNT but no AVX2);
//! set to a value, it runs on any CPU. With a value the cell test is one of those where it runs,
//! otherwise the masked AVX2 batch on an x86-64 CPU with AVX2 and the per-hit loop elsewhere, as
//! LOSAT_X_BPDIAG chooses;
//! `LOSAT_X_BXDIAG=per-hit` and `LOSAT_X_BXDIAG=scalar` force the per-hit loop or the scalar batch
//! (diagnostics; the tests use them). `LOSAT_X_BXDIAGSHADOW` (any value) also runs the reference on
//! a copy of the diagonal table for every subject and stops (panics) if the hit lists or the tables
//! differ, and prints a summary at exit. With `LOSAT_TIMING=1` a run prints
//! `[TIMING] blastx_diag: mode=<avx2-cross-defer|avx2-cross|w128-cross|avx2-defer-masked|scalar|per-hit> subjects=N cells=N
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
//! a longer query batch, the compressed table of `-word_size` 5 to 7 and `-task blastx-fast`, and
//! the table of 4-letter words, which the joint builder makes with `Int4` hits and the fused loop
//! scans as it is) keeps `Int4` hits.
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
//! of every table the fused loop takes where the CPU has AVX2 and POPCNT, and in 128 bits
//! (LOSAT_X_BXCROSS128, below) where it does not; `LOSAT_X_BXCROSS=0`, or a CPU without either,
//! leaves the direct loop on every table (with a value of `LOSAT_X_BXDIAG` or
//! `LOSAT_X_BXDIAGSHADOW`, the masked AVX2 batch above, or the per-hit loop without AVX2).
//! `LOSAT_X_BXDIAGMIN` still chooses the tables (above), and `LOSAT_X_BXDIAGSHADOW` compares it
//! with the reference like the other cell tests. The batches run at 1,024 staged hits. The timing
//! line names `mode=avx2-cross` (`mode=avx2-cross-defer` with
//! LOSAT_X_BXCROSSDEFER, below) and ends with `cross_splits=N` (the windows of
//! eight staged hits with a repeated entry, which these batches test in parts and
//! LOSAT_X_BXCROSS128's hit by hit; counted per subject scan like the other counts).
//!
//! EXPERIMENT (LOSAT_X_BXCROSSDEFER): `XCellTest::Cross { defer: true }`: LOSAT_X_BXCROSS's drain
//! defers each batch's calls to the batch that next reads one of their entries, or to the end of
//! the drain (`x_seed_fused`'s module header has the form and why the calls, their order and the
//! table are unchanged). A switch of the strict set (since 2026-10-25): on unless it is `0`;
//! `LOSAT_X_BXCROSSDEFER=0` gives LOSAT_X_BXCROSS's drain with the calls at once. It acts only
//! where LOSAT_X_BXCROSS's cell test runs; the timing line then names `mode=avx2-cross-defer`.
//!
//! EXPERIMENT (LOSAT_X_BXCROSS128): `XCellTest::Cross128`: LOSAT_X_BXCROSS's batches in two 128-bit
//! vectors of four lanes, each batch's calls at once and a window with a repeated entry tested hit
//! by hit, where the AVX2 batches do not run (the
//! WebAssembly build; x86-64 with SSSE3, SSE4.1 and POPCNT but no AVX2); `x_seed_fused`'s module
//! header has the form and why the calls, their order and the table are those of the AVX2 batches.
//! A switch of the strict set (class platform): on unless it is `0`; `2` runs it in place of the
//! AVX2 batches on any x86-64 CPU with SSSE3, SSE4.1 and POPCNT (diagnostic; the tests). With it the
//! fused loop runs by default in the WebAssembly build too (`x_bxdiag_mode`), on the tables
//! `LOSAT_X_BXDIAGMIN` gives it, and they take 16-bit hits there as on x86-64. The timing line then
//! names `mode=w128-cross` and ends with `cross_splits=N`.
//!
//! LOSAT_X_BXSCANBF (`x_seed_scan.rs`; graduated 2026-10-26, no variable, always on): the fused
//! loop finds the cells with hits by tblastn's branch-free candidate loop (`x_scan_two_hit_with`
//! with `branch_free: true`, blastp's `x_scan_cells_bf`): the same cells in the same order.
//!
//! EXPERIMENT (LOSAT_X_BXUNGAPSENT, `x_ungap_sent.rs`): for a subject on that path the closure's
//! two-hit extension is the sentinel walks, in an instance of its own (`x_word_finder_fused_sent`,
//! never inlined; the other statements unchanged), so the loop without it keeps its code.

use super::x_seed_scan::XKeptPresence;
use super::x_ungap_sent::{by_word_size, XSent};
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:211-242
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
/// by default (not set, or empty) only with a cell test of batches across the cells' borders,
/// LOSAT_X_BXCROSS's or LOSAT_X_BXCROSS128's ([`x_cross_batches`]).
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
/// EXPERIMENT (LOSAT_X_BXDIAG): whether the cell test is LOSAT_X_BXCROSS's or LOSAT_X_BXCROSS128's
/// (the batches across the cells' borders), where the default acts.
fn x_cross_batches() -> bool {
    match x_cell_mode() {
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { .. } => true,
        // EXPERIMENT (LOSAT_X_BXCROSS128): the same batches in 128 bits.
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XCellTest::Cross128 => true,
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

// No NCBI counterpart: reads the LOSAT_X_BXCROSS128 switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXCROSS128): 0 = off, 1 = LOSAT_X_BXCROSS's batches in 128-bit vectors
/// (`XCellTest::Cross128`) where the AVX2 batches do not run (wasm32 simd128; x86-64 with SSSE3,
/// SSE4.1 and POPCNT but not AVX2), 2 = those in place of the AVX2 batches too (diagnostic: the
/// tests run them end to end). A switch of the strict set, on unless it is `0`.
fn x_bx_cross128() -> usize {
    static ON: OnceLock<usize> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BXCROSS128", 1..=2, 1).unwrap_or(0))
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
            x_bx_cross128(),
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
/// EXPERIMENT (LOSAT_X_BXCROSS): the windows of eight staged hits with a repeated entry (tested in
/// parts by LOSAT_X_BXCROSS's AVX2 batches, hit by hit by LOSAT_X_BXCROSS128's), added with the
/// other counts once per subject scan ([`x_add_counts`]).
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
        // EXPERIMENT (LOSAT_X_BXCROSS, LOSAT_X_BXCROSS128): the windows with a repeated entry.
        if x_cross_batches() {
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
        x_db_masks: None,
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
/// hits, the candidates, and the cross-cell batches' windows with a repeated entry) added to the
/// totals.
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
        // The instance of the word size, chosen once for the subject (item SC9b-1).
        Some(sent) => by_word_size!(
            options.word_size,
            W => x_word_finder_fused_sent::<K, S, W>(
                mode, batch, parameters, options, subject, table, kept, diagonals, sent,
            )
        ),
        None => x_word_finder_fused(
            scoring, mode, batch, parameters, options, subject, table, kept, diagonals,
        ),
    }
}

/// EXPERIMENT (LOSAT_X_BXDIAG): `x_word_finder_direct`'s two-hit stage for one subject on the
/// standard lookup table `table`, with the scan and the diagonal test fused (cell test `mode`;
/// the branch-free scan, BXSCANBF, graduated); also returns the cells tested in
/// batches, their hits, the candidates and the cross-cell batches' windows with a repeated entry.
/// `subject` starts at the first residue and keeps the sentinel after the last one.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
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
                // The table's index is the binary search's (x_with_context_table); an offset outside
                // the table (none from the lookup) takes the binary search.
                let context = contexts.get(q).unwrap_or_else(|| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                });
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
                // (unsigned, as in `seed::word_finder_over`)
                if (q.wrapping_sub(diff) as u32) < batch.contexts[context].offset as u32 {
                    entry.set_last_hit(s.wrapping_add(offset));
                    return;
                }
                let p = &parameters[context];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                // The matrix's walks: BLOSUM62's `extend_two_hit_blosum62`, or another matrix's
                // (`bx_two_hit`).
                let Some(result) = bx_two_hit(
                    scoring,
                    query,
                    subject,
                    previous + word,
                    s as usize,
                    q as usize,
                    p.word_xdrop,
                    word as usize,
                ) else {
                    return;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
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
/// closure's other statements are those of `x_word_finder_fused`'s, in the same order. `W` is
/// `options.word_size` as a const (2 to 4), the walks' instance: the caller chooses it once per
/// subject (`by_word_size!`), and the closure still reads the word size from `options`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-606 (as in
// `x_word_finder_fused`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_fused_sent<K: XCell, S: BxScoring, const W: usize>(
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
                        if (q.wrapping_sub(diff) as u32) < batch.contexts[context].offset as u32 {
                            entry.set_last_hit(s.wrapping_add(offset));
                            return;
                        }
                        let p = &parameters[context];
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                        // (`x_ungap_sent.rs`)
                        let Some(result) = sent.two_hit::<$shadow, W>(
                            previous + word,
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
#[path = "x_seed_diag_tests.rs"]
mod tests;
