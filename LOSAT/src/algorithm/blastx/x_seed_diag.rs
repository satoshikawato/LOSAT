//! EXPERIMENT (LOSAT_X_BXDIAG / LOSAT_X_BXDIAGSHADOW): the BLASTX two-hit seed stage with the
//! subject scan and the diagonal test in one loop, and the hits of one lookup cell tested in AVX2
//! batches of eight (by default the batches cross the cells' borders: see the LOSAT_X_BXCROSS
//! paragraph below, where they rest on a conflict test instead).
//!
//! The reference (`word_finder`, and `x_word_finder_direct` of `x_seed_scan.rs`) has the
//! scanner copy every word hit into the offset-pair buffer and then reads the buffer back for the
//! two-hit test (NCBI's `s_BlastAaScanSubject` and `s_BlastAaWordFinder_TwoHit`). With the switch,
//! blastp's fused scan (`blastp::x_seed_fused::x_scan_two_hit_with`: LOSAT_X_BPSEED's loop with
//! LOSAT_X_BPDIAG's cell test) reads the hits of each lookup cell in place, in the same order, and
//! applies the outcomes of the test that touch only the hit's diagonal entry; a hit that passes the
//! test goes to the closure of [`x_word_finder_fused`], which runs the rest of the reference's
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
//! is not set (or empty) it acts only where the cell test is the AVX2 batch (an x86-64 CPU with
//! AVX2), the cell test it was timed with; set to a value, it runs on any CPU. The cell test is the
//! AVX2 batch on an x86-64 CPU with AVX2 and the per-hit loop elsewhere, as LOSAT_X_BPDIAG chooses;
//! `LOSAT_X_BXDIAG=per-hit` and `LOSAT_X_BXDIAG=scalar` force the per-hit loop or the scalar batch
//! (diagnostics; the tests use them). `LOSAT_X_BXDIAGSHADOW` (any value) also runs the reference on
//! a copy of the diagonal table for every subject and stops (panics) if the hit lists or the tables
//! differ, and prints a summary at exit. With `LOSAT_TIMING=1` a run prints
//! `[TIMING] blastx_diag: mode=<avx2-defer-masked|avx2-cross-defer|avx2-cross|avx2|scalar|per-hit> subjects=N cells=N
//! batched_hits=N candidates=N` (the subjects the fused loop scanned, the cells it tested in
//! batches, their hits, and the hits passed to the closure; zeros but the subjects with `per-hit`).
//!
//! EXPERIMENT (LOSAT_X_BXDEFER): the AVX2 batch is `XDiagMode::Avx2DeferMasked`: blastp's
//! LOSAT_X_BPDEFER (the calls of a cell's candidates after its batches) with the lanes' entry
//! coordinates from scalar loads and the last batch of fewer than eight hits as a masked AVX2 batch;
//! the argument that the calls, their order and the table are the AVX2 batch's is in
//! `x_seed_fused`'s module header. A switch of the strict set: on unless it is `0`;
//! `LOSAT_X_BXDEFER=0` gives the AVX2 batch (`mode=avx2`).
//! It acts where LOSAT_X_BXCROSS's cell test (below) does not run (`LOSAT_X_BXCROSS=0`, or a CPU
//! without POPCNT), on the tables with at least `LOSAT_X_BXDEFERMIN` hits per used backbone cell
//! (hundredths, built-in 700; `0` every table); a fused table with fewer takes the AVX2 batch, as
//! there the masked batch and the buffer cost more than they save (the fused loop forced on the
//! first 1,500 nt of a 2.1-kb fly cDNA x the fly proteome, 4.91 hits per used cell: x1.033
//! cycles; the whole cDNA, 6.77: x0.996). The choice is made with LOSAT_X_BXDIAGMIN's, once per
//! thread and table; as LOSAT_X_BXDIAGMIN's built-in value with this cell test (900) is above
//! 700, no fused table falls under the rule by default (only with LOSAT_X_BXDIAGMIN set below
//! 700). The timing line then also gives `defer_min_hits_per_cell=` and `tables_avx2=` (the fused
//! tables that took the AVX2 batch).
//! `LOSAT_X_BXDIAGSHADOW` compares it with the reference like the other cell tests, on the tables
//! the rule gives it.
//!
//! A query batch whose lookup table has fewer hits per used backbone cell than
//! `LOSAT_X_BXDIAGMIN` (hundredths, built-in 670 with LOSAT_X_BXCROSS's cell test and 900 with
//! the others; `0` fuses every batch, a number set is the cutoff of every cell test) is scanned
//! by the direct loop instead: with few hits per cell the fused loop's per-cell work costs more
//! than the copy it saves. Since LOSAT_X_BXHITBF made the direct loop faster, the fused loop with
//! LOSAT_X_BXDEFER's cell test and the direct loop break even
//! near 9 hits per used cell (about 8.75 against the fly proteome, 10.7 against Swiss-Prot; one
//! core, user-space cycles of the fused loop with the cell test it takes against the direct loop,
//! fly cDNAs x the fly proteome: 4.91 hits per used cell x1.20, 6.77 x1.13, 8.38 x1.03, 9.37
//! x0.98, 13.32 x0.90); the built-in value was 500 before (measured against the direct loop
//! without LOSAT_X_BXHITBF: 1.76 x1.20, 4.91 x1.00, 6.77 x0.94, 800 cDNAs x yeast x0.72), and
//! `LOSAT_X_BXDIAGMIN=500` gives that rule. With LOSAT_X_BXCROSS's cell test they break even near
//! 6.7 (about 6.46 against the fly proteome, 6.95 against Swiss-Prot; the cross test against the
//! direct loop: 4.05 x1.13, 5.34 x1.06, 6.33 x1.01, 6.63 x0.99, 7.34 x0.95, 8.38 x0.90), hence 670
//! for it (`LOSAT_X_BXDIAGMIN=900` gives the rule before). Both loops give the reference's hits
//! and table, so the choice changes only the speed; it is made once per thread and table build
//! ([`x_fused_pays`]; a table built from the pooled vectors of the one before is a new build),
//! and once for all the threads that search one table's subjects in tasks
//! ([`x_shared_table_ratio`]). The timing line then also gives `min_hits_per_cell=` (the cutoff of
//! the cell test that ran), `tables_fused=` and `tables_other=` (the choices made) and
//! `hits_per_cell=LEAST-MOST` (hundredths).
//!
//! EXPERIMENT (LOSAT_X_BXCROSS): the AVX2 batch is `XCellTest::Cross`:
//! the subject's hits staged and tested in windows of eight cut across the cells' borders, a
//! window with a repeated entry in parts cut before each hit that repeats an entry of its part, the
//! last hits of the subject one after the other; the argument that the calls, their order and the
//! table are the per-hit loop's is in `x_seed_fused`'s module header (it rests on the closures of
//! [`x_word_finder_fused`] and [`x_word_finder_fused_sent`] writing only their hit's entry). The
//! stage keeps no query offset: it is recovered as `(coord + s) & mask`, exact because the table
//! made by `Diagonals::new` for the query batch has at least its length plus the window entries,
//! so `mask` is at least every query offset of the batch's lookup
//! (`x_seed_fused::x_scan_cells_cross` asserts it per subject). A
//! switch of the strict set (since 2026-10-23): on unless it is `0`; `LOSAT_X_BXCROSS=0` gives
//! LOSAT_X_BXDEFER's choice above. On, it is the cell test
//! of every table the fused loop takes where the CPU has AVX2 and POPCNT, in place of
//! LOSAT_X_BXDEFER's and of LOSAT_X_BXDEFERMIN's rule (`tables_avx2=0`); without POPCNT the fused
//! loop takes LOSAT_X_BXDEFER's choice. `LOSAT_X_BXDIAGMIN` still chooses the tables (its cutoff
//! against the direct loop is a rule of its own; with this cell test the two loops break even
//! lower, and its built-in value is 670 instead of 900: above), and `LOSAT_X_BXDIAGSHADOW`
//! compares it with the reference like the other cell tests. The batches run at 1,024 staged
//! hits. The timing line names `mode=avx2-cross` and ends with `cross_splits=N` (the windows of
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
    x_cell_test, x_scan_two_hit_with, XCellTest, XContextTable, XDiagMode,
};
use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::OnceLock;

thread_local! {
    /// EXPERIMENT (LOSAT_X_BXDIAGMIN, LOSAT_X_BXDEFERMIN): the last lookup table this thread chose
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
/// by default (not set, or empty) only with the AVX2 batch ([`x_cell_mode`]).
pub(crate) fn x_bxdiag_mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        let named = std::env::var("LOSAT_X_BXDIAG").is_ok_and(|v| !v.trim().is_empty());
        if std::env::var_os("LOSAT_X_BXDIAGSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_BXDIAG") && (named || x_avx2_batches()) {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: chooses LOSAT's code path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): whether the cell test is the AVX2 batch, where the default acts.
fn x_avx2_batches() -> bool {
    match x_cell_mode() {
        #[cfg(target_arch = "x86_64")]
        XCellTest::Diag(XDiagMode::Avx2 | XDiagMode::Avx2DeferMasked) | XCellTest::Cross { .. } => {
            true
        }
        _ => false,
    }
}

// No NCBI counterpart: reads the LOSAT_X_BXDEFER switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXDEFER): whether the AVX2 batch is `XDiagMode::Avx2DeferMasked` (module
/// header); a switch of the strict set, on unless it is `0`.
fn x_bx_defer() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXDEFER"))
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
            x_bx_defer(),
            x_bx_cross(),
            x_bx_cross_defer(),
        )
    })
}

// No NCBI counterpart: reads LOSAT_X_BXDIAGMIN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the fewest hits per used lookup cell, in hundredths, for
/// which the fused loop with the cell test `mode` scans a query batch's subjects; a batch whose
/// table has fewer is scanned by the direct loop. Built-in, near where the two loops take the
/// same time: 670 with LOSAT_X_BXCROSS's cell test, 900 with the others (since LOSAT_X_BXHITBF;
/// 500 before it). A number set is the cutoff of every cell test, `0` fuses every batch; a value
/// that is not a number takes the built-in one.
fn x_min_hits_per_cell(mode: XCellTest) -> u64 {
    static MIN: OnceLock<Option<u64>> = OnceLock::new();
    let set = *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BXDIAGMIN")
            .ok()
            .and_then(|v| v.trim().parse().ok())
    });
    set.unwrap_or(match mode {
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { .. } => 670,
        XCellTest::Diag(_) => 900,
    })
}

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's loop; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the table's hits per used cell, in hundredths (0 for a table
/// without hits).
fn x_hits_per_cell(table: &BlastAaLookupTable) -> u64 {
    let (mut hits, mut used) = (0u64, 0u64);
    for cell in &table.backbone {
        if cell.num_used > 0 {
            hits += cell.num_used as u64;
            used += 1;
        }
    }
    (hits * 100).checked_div(used).unwrap_or(0)
}

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's loop; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN, LOSAT_X_BXDEFERMIN): `table`'s hits per used cell
/// ([`x_hits_per_cell`]), worked out once per table and thread (`build` = the table's build number
/// identifies it), and whether this call worked it out.
fn x_table_ratio(table: &BlastAaLookupTable, build: u64) -> (u64, bool) {
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
/// batch: its hits per used cell are at least [`x_min_hits_per_cell`] of `mode`. Counted for the
/// timing line once per table and thread, with the cell test `mode` would take for it
/// (LOSAT_X_BXDEFERMIN, [`x_table_mode`]).
fn x_fused_pays(table: &BlastAaLookupTable, build: u64, mode: XCellTest) -> bool {
    let (ratio, first) = x_table_ratio(table, build);
    if first {
        x_count_choice(ratio, mode);
    }
    ratio >= x_min_hits_per_cell(mode)
}

// No NCBI counterpart: counters of the timing line.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): counts one choice of the loop for a table of `ratio` hits per
/// used cell (see [`x_fused_pays`]).
fn x_count_choice(ratio: u64, mode: XCellTest) {
    let fused = ratio >= x_min_hits_per_cell(mode);
    X_TABLES[usize::from(!fused)].fetch_add(1, Relaxed);
    X_RATIO_MIN.fetch_min(ratio, Relaxed);
    X_RATIO_MAX.fetch_max(ratio, Relaxed);
    if fused && x_table_mode(mode, ratio) != mode {
        X_TABLES_AVX2.fetch_add(1, Relaxed);
    }
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
pub(crate) fn x_shared_table_ratio(lookup: &Lookup, options: &ResolvedOptions) -> Option<u64> {
    let Lookup::Standard(table, kept) = lookup else {
        return None;
    };
    if x_bxdiag_mode() == 0 || options.window_size <= 0 {
        return None;
    }
    let (ratio, first) = x_table_ratio(table, kept.build);
    if first {
        x_count_choice(ratio, x_cell_mode());
    }
    Some(ratio)
}

// No NCBI counterpart: see `x_shared_table_ratio` (speed only).
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): makes `ratio` (from [`x_shared_table_ratio`] for `lookup`) this
/// thread's choice for `lookup`'s table, so that the thread does not work it out again.
#[cfg(feature = "parallel")]
pub(crate) fn x_take_table_ratio(lookup: &Lookup, ratio: Option<u64>) {
    if let (Lookup::Standard(_, kept), Some(ratio)) = (lookup, ratio) {
        X_CHOICE.with(|choice| choice.set((kept.build, ratio)));
    }
}

// No NCBI counterpart: reads LOSAT_X_BXDEFERMIN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDEFERMIN): the fewest hits per used lookup cell, in hundredths, for
/// which the fused loop takes LOSAT_X_BXDEFER's cell test (built-in 700); a table with fewer takes
/// the AVX2 batch. With few hits per cell the masked last batch and the buffer cost more than they
/// save (the fused loop on 4.91 hits per used cell x1.033 cycles, 6.77 x0.996). `0` takes it for
/// every table; a value that is not a number takes the built-in one.
fn x_min_hits_per_cell_defer() -> u64 {
    static MIN: OnceLock<u64> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BXDEFERMIN")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(700)
    })
}

// No NCBI counterpart: chooses LOSAT's cell test; every cell test gives the reference's hits and
// table.
/// EXPERIMENT (LOSAT_X_BXDEFERMIN): the cell test for a table of `ratio` hits per used cell:
/// `mode`, but the AVX2 batch instead of LOSAT_X_BXDEFER's below [`x_min_hits_per_cell_defer`].
/// LOSAT_X_BXCROSS's cell test is kept on every table (module header).
fn x_table_mode(mode: XCellTest, ratio: u64) -> XCellTest {
    #[cfg(target_arch = "x86_64")]
    if mode == XCellTest::Diag(XDiagMode::Avx2DeferMasked) && ratio < x_min_hits_per_cell_defer() {
        return XCellTest::Diag(XDiagMode::Avx2);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = ratio;
    mode
}

// No NCBI counterpart: counters of the timing line and of the shadow summary.
/// EXPERIMENT (LOSAT_X_BXDIAGMIN): the (thread, table) choices for the fused loop (`[0]`) and
/// for the other loop (`[1]`), and the least and greatest hits per used cell (hundredths) seen.
static X_TABLES: [AtomicU64; 2] = [AtomicU64::new(0), AtomicU64::new(0)];
/// EXPERIMENT (LOSAT_X_BXDEFERMIN): the (thread, table) choices for the fused loop that took the
/// AVX2 batch instead of LOSAT_X_BXDEFER's cell test.
static X_TABLES_AVX2: AtomicU64 = AtomicU64::new(0);
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
            "[TIMING] blastx_diag: mode={} subjects={} cells={} batched_hits={} candidates={} min_hits_per_cell={} tables_fused={} tables_other={} hits_per_cell={}-{} defer_min_hits_per_cell={} tables_avx2={}",
            x_cell_mode().name(),
            X_SUBJECTS.load(Relaxed),
            X_CELLS.load(Relaxed),
            X_HITS.load(Relaxed),
            X_CANDIDATES.load(Relaxed),
            x_min_hits_per_cell(x_cell_mode()),
            X_TABLES[0].load(Relaxed),
            X_TABLES[1].load(Relaxed),
            X_RATIO_MIN.load(Relaxed).min(X_RATIO_MAX.load(Relaxed)),
            X_RATIO_MAX.load(Relaxed),
            x_min_hits_per_cell_defer(),
            X_TABLES_AVX2.load(Relaxed),
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
/// are compared and the tables, which may differ in that entry, are not (as LOSAT_X_BPSEEDSHADOW
/// does for blastp).
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
pub(crate) fn x_word_finder_diag(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_>>,
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
        let pays = x_fused_pays(table, kept.build, mode);
        crate::utils::xreach::x_reach!("LOSAT_X_BXDIAG", pays);
        if !pays {
            return Ok(None);
        }
        // EXPERIMENT (LOSAT_X_BXDEFERMIN): BXDEFER's cell test only with enough hits per cell.
        let table_mode = x_table_mode(mode, x_table_ratio(table, kept.build).0);
        #[cfg(target_arch = "x86_64")]
        if mode == XCellTest::Diag(XDiagMode::Avx2DeferMasked) {
            crate::utils::xreach::x_reach!("LOSAT_X_BXDEFER", table_mode == mode);
        }
        let (hits, counts) = x_word_finder_fused_or_sent(
            table_mode, batch, parameters, options, subject, table, kept, diagonals, sent,
        )?;
        x_add_counts(counts);
        return Ok(Some(hits));
    }
    let mut copy = Diagonals {
        entries: diagonals.entries.clone(),
        offset: diagonals.offset,
        window: diagonals.window,
    };
    let reference = word_finder(
        batch,
        parameters,
        options,
        subject,
        lookup,
        &mut copy,
        |_| {},
    )?;
    let (fused, counts) = x_word_finder_fused_or_sent(
        x_table_mode(mode, x_table_ratio(table, kept.build).0),
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
/// sentinel walks, [`x_word_finder_fused_sent`].
#[allow(clippy::too_many_arguments)]
fn x_word_finder_fused_or_sent(
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable,
    kept: &XKeptPresence,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_>>,
) -> Result<(Vec<InitHsp>, [u64; 4])> {
    match sent {
        Some(sent) => x_word_finder_fused_sent(
            mode, batch, parameters, options, subject, table, kept, diagonals, sent,
        ),
        None => x_word_finder_fused(
            mode, batch, parameters, options, subject, table, kept, diagonals,
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
fn x_word_finder_fused(
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable,
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
                let Some(result) = extend_two_hit_blosum62(
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
fn x_word_finder_fused_sent(
    mode: XCellTest,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    table: &BlastAaLookupTable,
    kept: &XKeptPresence,
    diagonals: &mut Diagonals,
    sent: XSent<'_>,
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
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    /// The cell tests this machine can run (the AVX2 batch with LOSAT_X_BXDEFER and without, and
    /// with LOSAT_X_BXCROSS, its calls deferred with LOSAT_X_BXCROSSDEFER and not), each with the
    /// branch-free scan (BXSCANBF, graduated).
    fn modes() -> Vec<XCellTest> {
        let mut modes = vec![
            XCellTest::Diag(XDiagMode::PerHit),
            XCellTest::Diag(XDiagMode::Scalar),
        ];
        for (masked, cross, defer) in [
            (true, false, false),
            (false, false, false),
            (false, true, false),
            (false, true, true),
        ] {
            let auto = x_cell_test("", masked, cross, defer);
            if auto != XCellTest::Diag(XDiagMode::PerHit) {
                modes.push(auto);
            }
        }
        modes
    }

    /// A search of `query` (nucleotide FASTA) with `-window_size window` (default options
    /// otherwise), its statistics set for subjects of `subject_length` residues.
    fn search(
        query: &[u8],
        window: i32,
        subject_length: usize,
    ) -> (
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        Lookup,
    ) {
        let q = parse_fasta(query, false, false).unwrap();
        let window = window.to_string();
        let args = [
            "losat",
            "blastx",
            "-query",
            "query",
            "-subject",
            "subject",
            "-window_size",
            &window,
        ];
        let crate::cli::Commands::Blastx(args) =
            crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                .unwrap()
                .command
        else {
            unreachable!()
        };
        let options = args.resolve().unwrap();
        let batch = prepare_queries(&q, &options).unwrap();
        let mut params = parameters::score_block(&batch);
        parameters::effective_lengths(
            &batch,
            &mut params,
            &options,
            subject_length as i64,
            1,
            None,
        )
        .unwrap();
        let lookup = Lookup::new(&batch, &params, &options).unwrap();
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
        (batch, params, options, lookup): &(
            PreparedQueryBatch,
            Vec<ContextParameters>,
            ResolvedOptions,
            Lookup,
        ),
        subjects: &[Vec<u8>],
        mode: XCellTest,
        start: impl Fn(&mut Diagonals),
    ) -> (usize, [u64; 4]) {
        let Lookup::Standard(table, kept) = lookup else {
            panic!("the standard lookup table")
        };
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
        for (i, subject) in subjects.iter().enumerate() {
            let expected = word_finder(
                batch,
                params,
                options,
                subject,
                lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            let (actual, c) = x_word_finder_fused(
                mode, batch, params, options, subject, table, kept, &mut fused,
            )
            .unwrap();
            let run = format!("{mode:?} subject {i}");
            assert_eq!(actual, expected, "{run}");
            assert_eq!(fused.offset, reference.offset, "{run}");
            assert!(fused.entries == reference.entries, "{run}");
            let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
            for (table_sent, shadow) in sent.iter_mut().zip([false, true]) {
                let x = XSent::for_test(&batch.sequence_start, &lead, shadow);
                let (got, c2) = x_word_finder_fused_sent(
                    mode, batch, params, options, subject, table, kept, table_sent, x,
                )
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
        for window in [1, 39, 40, 41] {
            let s = search(&query, window, length);
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
            let (batch, _, _, lookup) = search(&query, 40, 500);
            let Lookup::Standard(_, kept) = &lookup else {
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
        for window in [40, 400] {
            let s = search(&query, window, mean);
            // The lookup has cells of a full batch (eight hits) and more.
            let Lookup::Standard(table, _) = &s.3 else {
                panic!("the standard lookup table")
            };
            assert!(table.backbone.iter().any(|cell| cell.num_used >= 16));
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited above);
    // no NCBI counterpart for the table (buffer reuse, LOSAT_X_BXDIAGHOLD): a table of two and eight
    // times the entries (`Diagonals::x_holds`), first used for another, longer query batch's
    // subjects, gives the reference's hits (and offsets moved on alike), subject after subject, in the reference
    // loop and the fused loop of every cell test: 6 kb of the AP027131 genome against 60 of its
    // proteins (a 9 kb piece before), windows 40 and 400.
    #[test]
    fn a_longer_table_gives_the_same_hits() {
        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let genome = parse_fasta(
            &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
            false,
            false,
        )
        .unwrap();
        let piece = |from: usize, n: usize| {
            let letters = String::from_utf8(genome[0].sequence[from..from + n].to_vec()).unwrap();
            format!(">piece\n{letters}\n").into_bytes()
        };
        let all = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
        let subjects = &all[..60];
        let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
        for window in [40, 400] {
            let s = search(&piece(200_000, 6_000), window, mean);
            let before = search(&piece(400_000, 9_000), window, mean);
            let length = s.0.contexts.last().map(|c| c.offset + c.length).unwrap();
            let longer_length = before
                .0
                .contexts
                .last()
                .map(|c| c.offset + c.length)
                .unwrap();
            for factor in [2, 8] {
                let mut reference = Diagonals::new(length, window).unwrap();
                let mut tables: Vec<Diagonals> = (0..=modes().len())
                    .map(|_| Diagonals::new(longer_length * factor, window).unwrap())
                    .collect();
                for table in tables.iter_mut() {
                    assert!(table.x_holds(length, window) && table.x_holds(longer_length, window));
                    assert!(table.entries.len() > reference.entries.len());
                    for subject in &subjects[..20] {
                        word_finder(
                            &before.0,
                            &before.1,
                            &before.2,
                            subject,
                            &before.3,
                            table,
                            |_| {},
                        )
                        .unwrap();
                    }
                }
                let (batch, params, options, lookup) = &s;
                let Lookup::Standard(table, kept) = lookup else {
                    panic!("the standard lookup table")
                };
                // The earlier subjects moved the longer tables' offsets on; each subject then moves
                // every table's offset by the same amount.
                let ahead: Vec<i32> = tables.iter().map(|t| t.offset - reference.offset).collect();
                for (i, subject) in subjects.iter().enumerate() {
                    let expected = word_finder(
                        batch,
                        params,
                        options,
                        subject,
                        lookup,
                        &mut reference,
                        |_| {},
                    )
                    .unwrap();
                    let run = format!("window {window} x{factor} subject {i}");
                    let got = word_finder(
                        batch,
                        params,
                        options,
                        subject,
                        lookup,
                        &mut tables[0],
                        |_| {},
                    )
                    .unwrap();
                    assert_eq!(got, expected, "{run}");
                    assert_eq!(tables[0].offset - reference.offset, ahead[0], "{run}");
                    for ((mode, longer), &ahead) in
                        modes().into_iter().zip(&mut tables[1..]).zip(&ahead[1..])
                    {
                        let (got, _) = x_word_finder_fused(
                            mode, batch, params, options, subject, table, kept, longer,
                        )
                        .unwrap();
                        assert_eq!(got, expected, "{run} {mode:?}");
                        assert_eq!(longer.offset - reference.offset, ahead, "{run} {mode:?}");
                    }
                }
            }
        }
    }

    // No NCBI counterpart (speed only): the choice of the loop names a table by its build number,
    // so a table at the address and with the vectors of the one before it is worked out again.
    #[test]
    fn the_loop_choice_is_worked_out_again_for_each_table_build() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let query = std::fs::read(dir.join("query.fna")).unwrap();
        let (_, _, _, lookup) = search(&query, 40, 500);
        let Lookup::Standard(table, kept) = &lookup else {
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
