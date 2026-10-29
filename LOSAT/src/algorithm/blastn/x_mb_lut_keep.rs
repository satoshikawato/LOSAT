//! EXPERIMENT (LOSAT_X_MBLUTKEEP; in the strict set since 2026-10-19: on unless it is `0`)
//! : the array of one element per cell of a query batch's lookup table kept from one query batch
//! to the next of the batch loop, instead of a fresh allocation per batch that is freed at the end
//! of its batch: the hashtable of `build_mb_lookup` (the megablast table of the megablast and
//! blastn tasks, zeroed for the next build), the `offsets` of a table laid out as CSR
//! (LOSAT_X_BNCSR, `x_bn_csr.rs`; LOSAT_X_DCCSR, `x_dc_csr.rs`) and the chain heads of a
//! discontiguous table of one template (`build_disc_mb_lookup_kept`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262,1376
//! ```c
//!    mb_lt->hashtable = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
//! ...
//!    sfree(mb_lt->hashtable);
//! ```
//!
//! What it changes. Each query batch builds its table (`BlastMBLookupTableNew`) and frees it with
//! the batch (`LookupTableWrapFree`, `BlastMBLookupTableDestruct`). A table of 4^12 cells is 64
//! MiB: with a fresh allocation every batch maps it anew (glibc maps blocks this large and unmaps
//! them when they are freed), its pages fault while the fill writes them (or in the prefault round
//! of LOSAT_X_MBLUTPREFAULT), and freeing it unmaps the pages and flushes their entries from the
//! TLBs of every CPU the search ran on. Here the batch loop holds a `Keep`; a table built for the
//! loop goes to it when its batch drops the table (`MbLookupTable`'s `Drop`), and the next build
//! takes it when it is at least as long as the build's table: it writes 0 to the cells of its table
//! that earlier builds wrote (on the search pool when the build runs on a pool of two threads or
//! more; the other cells are 0 already) and the build uses its first cells as its table. The
//! longest table is kept: NCBI's `CBatchSizeMixer` alternates full batches with small probe
//! batches (megablast of 31,168 D. melanogaster cDNAs: 18 tables of 4^12 cells and 19 of 4^8 or
//! 4^11), and a probe batch's table is the start of the kept one. The keep records the longest
//! start of the table that builds wrote since it was last zeroed (`dirty`), so a probe batch's
//! build zeroes its own cells only and leaves the rest of the full batch's cells to the next build
//! that uses them (a search that ends with probe batches never zeroes them; megablast of 1,775
//! cDNAs builds tables of 4^11, 4^8 and 4^12 cells). A build whose table is longer than the kept
//! one frees the kept one first and allocates a fresh table, which the keep takes back instead.
//! The loop frees the `Keep` and its table when it ends.
//!
//! A table laid out as CSR has no hashtable; its `offsets` (a cell's first position, one element
//! per cell and one or two more) is built in the kept table instead (`take_buffer`: the build
//! writes every element it uses, so the table is not zeroed for it) and goes back to the keep when
//! the batch drops the CSR (`give_back_buffer`); a discontiguous table of one template takes its
//! chain heads from the keep as `build_mb_lookup` takes its hashtable. These arrays do not depend
//! on the batch's letters, only on its table's cells (4^11 or 4^12), so one kept table serves the
//! search's batches whichever layout each one takes: `CBatchSizeMixer`'s full batches are laid out
//! as CSR and its probe batches as chains (dc-megablast of 31,168 cDNAs: 18 CSR builds and 19
//! chain builds of 4^11 cells). A fresh table has room for two elements more than its cells, so a
//! later CSR can take it.
//!
//! Why it matters for memory. glibc maps a block of at least its mmap threshold (128 KiB at the
//! start) and unmaps it when it is freed, but the `free` of a mapped block also raises the
//! threshold to the block's size (up to 32 MiB), and the heap then serves every later allocation
//! below it. A 16 MiB array freed after each batch put the threshold up after the first batch,
//! the heap took the batches' other buffers, and the free space in its arena stayed resident: the
//! dc-megablast above at one thread peaked at an anonymous RSS of 84-92 MiB depending on the
//! build (where glibc's arenas fall), against 80.8 MiB for the memory floor of plan v4 line MF.
//! With the CSR and the chains on the kept table the threshold no longer moves at a batch's end:
//! 79.5 MiB, with fewer page faults than before (183 k against 191 k).
//!
//! Why the table is NCBI's. The fill reads a cell before it writes it (`hashtable[ecode] == 0`, the
//! cell's head) and relies on `calloc`'s zeros; on a kept table every cell of the build's table is
//! 0 before the fill starts (its cells that earlier builds wrote are set to 0, and no build writes
//! past its own cells), so the fill sees what it sees in a fresh table and writes the same cells
//! (the other arrays, `next_pos`, `pv_array` and the counters of `longest_chain`, are allocated as
//! before). Every reader takes cells by a word's index, below `1 << 2 * lut_word_length` (the
//! table's cells, `hashtable[..table_size]` wherever the build passes the table on), so the
//! elements past the table are never read. A discontiguous table's chain heads are the kept
//! table's first cells, cut to the table's length. A CSR build writes 0 to the elements of its
//! `offsets` before it counts the positions into them and reads no element it has not written
//! (`x_dc_csr::build`), so a kept buffer gives the CSR of a fresh one. No floating-point operation
//! is involved.
//!
//! The prefault round. On a kept table the round of LOSAT_X_MBLUTPREFAULT does not run: it exists
//! to map the pages of a fresh table by writing them before the fill reads them (its module
//! documentation, `x_mb_lookup_once.rs`), and a kept table's pages are mapped already (written by
//! an earlier fill or by the zeroing, which writes every cell of the build: the cells past those
//! that builds wrote too, as a packed CSR of LOSAT_X_BNCSRPACK writes only the start of a table).
//! Its writes are zeros over zeros, so the arrays are the same either way.
//!
//! Where. Only the batch loop of `run_in_pool` keeps a table, for the batches it searches (and the
//! query chunks of a split batch). Until 2026-10-28 it kept none while the query-batch slots of
//! LOSAT_X_QBATCHPAR (removed then) could run. The tables of PAIRPAR's sub-lookups, the discontiguous tables of two templates and the tables
//! of the tests are not kept.
//!
//! Memory. Between two batches the loop holds the longest table built so far, and no build
//! allocates a second table beside it (a longer one replaces it). A probe batch after a full batch
//! keeps the full batch's table resident where the table was freed before; its peak (the kept
//! table, the probe batch's `next_pos` and the rest of a small batch) stays below the full
//! batch's (the same table, a full batch's `next_pos` and the rest), so the search's peak is the
//! full batch's as before. With LOSAT_X_BNLUTDROP a batch drops its tables before its traceback,
//! and the kept table stays resident through it (16 MiB of 4^11 cells, 64 MiB of 4^12); the
//! peaks above are measured with it.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_mb_lut_keep: kept= fresh= csr= zero=` gives the loop's
//! builds of a hashtable or chain heads on the kept table, its builds on a fresh one (of either
//! layout), its CSR builds whose `offsets` is the kept table, and the time spent zeroing (also
//! inside `x_mb_lookup_split`'s `alloc=`).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// No NCBI counterpart: LOSAT_X_MBLUTKEEP (read once).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTKEEP"))
}

/// No NCBI counterpart (scheduling only): the elements one job of the pool's zeroing writes at
/// least (256 KiB).
const ZERO_CHUNK: usize = 1 << 16;

static KEPT: AtomicU64 = AtomicU64::new(0);
static FRESH: AtomicU64 = AtomicU64::new(0);
static CSR: AtomicU64 = AtomicU64::new(0);
static ZERO_NANOS: AtomicU64 = AtomicU64::new(0);

/// The kept table: as long as the largest table built so far; its first `dirty` elements may
/// hold values of earlier fills (the longest start that builds wrote since it was last zeroed),
/// the others are 0.
struct Kept {
    table: Vec<u32>,
    dirty: usize,
}

/// The table one batch loop keeps between its batches (LOSAT_X_MBLUTKEEP).
#[derive(Default)]
pub(crate) struct Keep {
    kept: Mutex<Option<Kept>>,
}

impl Keep {
    /// A keep for one batch loop, when LOSAT_X_MBLUTKEEP is on.
    pub(crate) fn for_loop() -> Option<Arc<Keep>> {
        on().then(|| Arc::new(Keep::default()))
    }

    /// A table for a build of `cells` cells, its first `cells` elements 0, and the number of its
    /// first elements that may hold values once the build has filled it: the kept one when it has
    /// at least that many elements (`true`; the build uses its first `cells`), with its first
    /// `cells` elements written with 0 (the written elements past them stay as they are and stay
    /// counted); else a fresh `vec![0; cells]` (`false`), after the kept one is freed. `on_pool`:
    /// the build runs on a pool of two threads or more, which writes the zeros.
    pub(super) fn take(&self, cells: usize, on_pool: bool) -> (Vec<u32>, bool, usize) {
        let kept = self.kept.lock().unwrap_or_else(|e| e.into_inner()).take();
        // The reference side is a kept table too small for this build (the first build of a loop
        // has none on either path).
        let too_small = kept.as_ref().is_some_and(|kept| kept.table.len() < cells);
        let Some(Kept { mut table, dirty }) = kept.filter(|kept| kept.table.len() >= cells) else {
            if too_small {
                crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", false);
            }
            FRESH.fetch_add(1, Ordering::Relaxed);
            // Room for the `offsets` of a later CSR build (two elements more, `x_bn_csr`), which
            // takes the table from the keep (`take_buffer`).
            let mut table = vec![0u32; cells + 2];
            table.truncate(cells);
            return (table, false, cells);
        };
        let start = super::x_mb_lookup_par::timing_on().then(std::time::Instant::now);
        // Every element of the build's cells, also those past `dirty`, which no build wrote: a
        // packed CSR (LOSAT_X_BNCSRPACK) writes only the start of its `offsets`, so the rest of a
        // table it gave back may be pages never mapped. The fill reads a cell before it writes it,
        // and on such a page the read maps the zero page and the write copies it, a TLB flush of
        // every CPU the search runs on; writing first maps the page as the prefault round does on
        // a fresh table. Elsewhere `dirty` covers the cells already, and this writes the same.
        zero(&mut table[..cells], on_pool);
        if let Some(start) = start {
            ZERO_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", true);
        KEPT.fetch_add(1, Ordering::Relaxed);
        (table, true, dirty.max(cells))
    }

    /// The kept table as the buffer of a CSR's `offsets` (`x_dc_csr::x_offsets_buffer`; the
    /// build writes every element it uses): the table, its length and the number of its first
    /// elements that earlier builds wrote, when it has room for `len` elements (a discontiguous
    /// table's heads are the first cells of a CSR's `offsets`, one shorter); else an empty vector,
    /// after the kept one is freed (as `take`). The CSR gives its `offsets` back when it is
    /// dropped (`give_back_buffer`), so a loop whose batches alternate between CSR and chains
    /// keeps one table for both.
    pub(super) fn take_buffer(&self, len: usize) -> (Vec<u32>, usize, usize) {
        let kept = self.kept.lock().unwrap_or_else(|e| e.into_inner()).take();
        match kept {
            Some(Kept { table, dirty }) if table.capacity() >= len => {
                crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", true);
                CSR.fetch_add(1, Ordering::Relaxed);
                let length = table.len();
                (table, length, dirty)
            }
            kept => {
                if kept.is_some() {
                    crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", false);
                }
                FRESH.fetch_add(1, Ordering::Relaxed);
                drop(kept);
                (Vec::new(), 0, 0)
            }
        }
    }

    /// Takes back a CSR's `offsets` (`take_buffer`'s table of `length` elements, or a vector of its
    /// own when `length` is 0), with every element the build wrote counted as written.
    pub(super) fn give_back_buffer(&self, mut table: Vec<u32>, length: usize, dirty: usize) {
        let used = table.len();
        if length > used {
            // SAFETY: the elements past `used` are those of the kept table of `length` elements
            // (initialized, unchanged: the build wrote its first `used` only, in its capacity).
            unsafe { table.set_len(length) };
        }
        self.give_back(table, dirty.max(used));
    }

    /// Takes back the table of a dropped `MbLookupTable`, whose first `dirty` elements may hold
    /// values (`take`'s third value): kept when nothing is kept or the kept one is shorter, else
    /// freed.
    pub(super) fn give_back(&self, table: Vec<u32>, dirty: usize) {
        let mut slot = self.kept.lock().unwrap_or_else(|e| e.into_inner());
        if slot
            .as_ref()
            .is_none_or(|kept| kept.table.len() < table.len())
        {
            let dirty = dirty.min(table.len());
            let old = slot.replace(Kept { table, dirty });
            drop(slot);
            drop(old);
        }
    }
}

/// Writes 0 to every element, with the pool's threads when `on_pool`.
fn zero(table: &mut [u32], on_pool: bool) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if on_pool {
        use rayon::prelude::*;
        table
            .par_chunks_mut(ZERO_CHUNK)
            .for_each(|chunk| chunk.fill(0));
        return;
    }
    let _ = (on_pool, ZERO_CHUNK);
    table.fill(0);
}

/// No NCBI counterpart: printed by `main` at exit (through `x_mb_lookup_print_summary`): with
/// LOSAT_TIMING the `[TIMING] x_mb_lut_keep` line when a loop kept tables. It prints counters and
/// changes no value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if super::x_mb_lookup_par::timing_on() && load(&KEPT) + load(&FRESH) + load(&CSR) > 0 {
        eprintln!(
            "[TIMING] x_mb_lut_keep: kept={} fresh={} csr={} zero={:.4}s",
            load(&KEPT),
            load(&FRESH),
            load(&CSR),
            load(&ZERO_NANOS) as f64 * 1e-9
        );
    }
}

#[cfg(test)]
mod tests {
    use super::super::x_mb_lookup_par::tests::{random_case, Lcg};
    use super::super::{build_mb_lookup, x_mb_narrow, MbLookupTable};
    use super::*;

    /// The kept table's length and the number of its first elements that may be dirty.
    fn kept(keep: &Keep) -> Option<(usize, usize)> {
        let kept = keep.kept.lock().unwrap();
        kept.as_ref().map(|kept| (kept.table.len(), kept.dirty))
    }

    /// The build's table (its first cells) and the other arrays and fields equal the fresh
    /// build's (both of 3-byte elements, LOSAT_X_MB24, or both not).
    fn assert_same(
        built: &MbLookupTable,
        fresh: &MbLookupTable,
        cells: usize,
        next_pos_len: usize,
        what: &str,
    ) {
        assert_eq!(built.x_narrow.is_some(), fresh.x_narrow.is_some(), "{what}");
        let (built_heads, built_next) = built.x_wide(cells, next_pos_len);
        let (fresh_heads, fresh_next) = fresh.x_wide(cells, next_pos_len);
        assert!(built_heads == fresh_heads, "{what}: hashtable");
        assert!(built_next == fresh_next, "{what}: next_pos");
        assert!(built.pv_array == fresh.pv_array, "{what}: pv_array");
        assert_eq!(built.pv_array_bts, fresh.pv_array_bts, "{what}");
        assert_eq!(built.longest_chain, fresh.longest_chain, "{what}");
    }

    // A keep takes back the longest table with the start that builds wrote, gives it out to
    // builds of no more cells with their cells set to 0 (the written elements past them stay and
    // stay counted), and frees it for a longer build.
    #[test]
    fn the_keep_holds_the_longest_table_and_gives_it_out_zeroed() {
        let keep = Keep::default();
        let (table, from_keep, dirty) = keep.take(8, false);
        assert_eq!((from_keep, dirty), (false, 8));
        assert_eq!(table, vec![0; 8]);
        keep.give_back(vec![7; 8], 8);
        keep.give_back(vec![5; 4], 4);
        assert_eq!(kept(&keep), Some((8, 8)), "a shorter table is freed");
        let (mut table, from_keep, dirty) = keep.take(4, false);
        assert!(from_keep, "a shorter build takes the kept table");
        assert_eq!(
            dirty, 8,
            "the written elements past the build's cells stay counted"
        );
        assert_eq!(
            table,
            [0, 0, 0, 0, 7, 7, 7, 7],
            "only the build's cells are set to 0"
        );
        assert_eq!(kept(&keep), None);
        table[..2].fill(3);
        keep.give_back(table, dirty);
        assert_eq!(kept(&keep), Some((8, 8)));
        let (table, from_keep, dirty) = keep.take(8, false);
        assert_eq!((from_keep, dirty), (true, 8));
        assert_eq!(table, vec![0; 8]);
        keep.give_back(table, 4);
        assert_eq!(kept(&keep), Some((8, 4)));
        let (table, from_keep, dirty) = keep.take(2, false);
        assert_eq!((from_keep, dirty), (true, 4));
        keep.give_back(table, dirty);
        let (table, from_keep, _) = keep.take(16, false);
        assert!(!from_keep, "a longer build allocates");
        assert_eq!(table, vec![0; 16]);
        assert_eq!(kept(&keep), None, "the shorter kept table is freed first");
        keep.give_back(vec![9; 1 << 20], 1 << 20);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(4)
                .build()
                .unwrap();
            let (table, from_keep, dirty) = pool.install(|| keep.take(1 << 12, true));
            assert_eq!((from_keep, dirty), (true, 1 << 20));
            assert!(
                table[..1 << 12].iter().all(|&cell| cell == 0),
                "zeroed on the pool"
            );
            assert!(
                table[1 << 12..].iter().all(|&cell| cell == 9),
                "past the build's cells"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262
    // ```c
    //    mb_lt->hashtable = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
    // ```
    // Random queries (ambiguity codes, masks, both cell orders) of tables of 4^4 to 4^12 cells in
    // a row through one keep, one after another and on pools of 2 and 5 threads (the pool fill
    // with the `resident` flag; these queries are too short for the prefault round, which the
    // whole-search tests cover): every build equals the same build on a fresh table, every build
    // of no more cells than the kept table takes it (filled by earlier builds), and the elements
    // past the start that builds wrote stay 0. With LOSAT_X_MB24 the newest-first tables of 4^12
    // cells take the words of their 3-byte heads, fewer than 4^12. With LOSAT_X_BNCSRPACK a
    // newest-first table of 4^8 to 4^11 cells is CSR at any size, so those tables are built in
    // reading order here, which keeps their chains (the switch off builds them newest first too).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn builds_on_a_kept_table_are_the_fresh_builds() {
        const LUTS: [usize; 8] = [11, 8, 11, 11, 4, 12, 9, 12];
        let mut rng = Lcg(0x6b65_6570);
        for threads in [1usize, 2, 5] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let keep = Arc::new(Keep::default());
                let (mut reused, mut reused_dirty) = (0, 0);
                for case in 0..LUTS.len() * 2 {
                    let (queries, offsets, masks, word_length, _) = random_case(&mut rng);
                    let lut = LUTS[case % LUTS.len()];
                    let word_length = word_length.max(lut);
                    let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                    let ascending = case % 3 == 1
                        || (super::super::x_bn_csr::x_will_pack(false) && (8..=11).contains(&lut));
                    let fresh = build_mb_lookup(
                        &refs,
                        &offsets,
                        word_length,
                        lut,
                        &masks,
                        None,
                        0,
                        1_000_000,
                        ascending,
                        None,
                    );
                    let next_pos_len = queries
                        .iter()
                        .zip(&offsets)
                        .map(|(q, &o)| o.max(0) as usize + q.len())
                        .max()
                        .unwrap_or(0)
                        + 1;
                    let cells = if fresh.x_narrow.is_some() {
                        x_mb_narrow::words_for(1 << (2 * lut))
                    } else {
                        1 << (2 * lut)
                    };
                    let before = keep.kept.lock().unwrap().as_ref().map(|kept| {
                        assert!(
                            kept.table[kept.dirty..].iter().all(|&cell| cell == 0),
                            "case {case}: an element past the dirty cells is not 0"
                        );
                        (
                            kept.table.len(),
                            kept.dirty,
                            kept.table[..kept.dirty.min(cells)].iter().any(|&c| c != 0),
                        )
                    });
                    let built = build_mb_lookup(
                        &refs,
                        &offsets,
                        word_length,
                        lut,
                        &masks,
                        None,
                        0,
                        1_000_000,
                        ascending,
                        Some(&keep),
                    );
                    let dirty = match before {
                        Some((len, dirty, written)) if len >= cells => {
                            assert_eq!(kept(&keep), None, "case {case}: the kept table taken");
                            assert_eq!(built.x_table().len(), len, "case {case}");
                            reused += 1;
                            reused_dirty += usize::from(written);
                            dirty.max(cells)
                        }
                        _ => {
                            assert_eq!(built.x_table().len(), cells, "case {case}: fresh");
                            cells
                        }
                    };
                    assert_same(
                        &built,
                        &fresh,
                        1 << (2 * lut),
                        next_pos_len,
                        &format!("case {case} lut {lut} threads {threads}"),
                    );
                    let len = built.x_table().len();
                    drop(built);
                    assert_eq!(kept(&keep), Some((len, dirty)), "case {case}: given back");
                }
                // Every build but the first (4^11) and the first 4^12, which is longer than the
                // kept table: 14 of 16. With LOSAT_X_MB24 the first 4^12 (case 5) is a 3-byte
                // table, shorter than a 4-byte one, so the next 4^12 (case 7, ascending: 4-byte)
                // is longer than the kept table too: 13.
                let expected = if x_mb_narrow::on() { 13 } else { 14 };
                assert_eq!(
                    reused, expected,
                    "threads {threads}: builds on a kept table"
                );
                assert!(
                    reused_dirty >= 7,
                    "threads {threads}: {reused_dirty} filled tables kept"
                );
            });
        }
    }
}
