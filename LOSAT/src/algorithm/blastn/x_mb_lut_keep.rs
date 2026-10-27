//! EXPERIMENT (LOSAT_X_MBLUTKEEP; in the strict set since 2026-10-19: on unless it is `0`)
//! : the hashtable of `build_mb_lookup` (the megablast table
//! of the megablast and blastn tasks) kept from one query batch to the next of the batch loop and
//! zeroed for the next build, instead of a fresh zeroed allocation per batch that is freed at the
//! end of its batch.
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
//! Why the table is NCBI's. The fill reads a cell before it writes it (`hashtable[ecode] == 0`, the
//! cell's head) and relies on `calloc`'s zeros; on a kept table every cell of the build's table is
//! 0 before the fill starts (its cells that earlier builds wrote are set to 0, and no build writes
//! past its own cells), so the fill sees what it sees in a fresh table and writes the same cells
//! (the other arrays, `next_pos`, `pv_array` and the counters of `longest_chain`, are allocated as
//! before). Every reader takes cells by a word's index, below `1 << 2 * lut_word_length` (the
//! table's cells, `hashtable[..table_size]` wherever the build passes the table on), so the
//! elements past the table are never read. No floating-point operation is involved.
//!
//! The prefault round. On a kept table the round of LOSAT_X_MBLUTPREFAULT does not run: it exists
//! to map the pages of a fresh table by writing them before the fill reads them (its module
//! documentation, `x_mb_lookup_once.rs`), and a kept table's pages are mapped already (written by
//! the zeroing or by an earlier fill). Its writes are zeros over zeros, so the arrays are the same
//! either way.
//!
//! Where. Only the batch loop of `run_in_pool` keeps a table, for the batches it searches itself
//! (and the query chunks of a split batch), and only while no LOSAT_X_QBATCHPAR slot can run:
//! after the driver turned off, or before the first batch is committed where no slot can start
//! (`x_qbatch_par::Driver::may_speculate`: no `LOSAT_X_QBATCHPARSLOTS`, and a budget of 0 or one
//! thread). A slot holds a whole batch while others run, and QBATCHPAR's budget counts one table
//! per batch in flight, so a table kept between the loop's batches would add to the slots' peak.
//! The tables of PAIRPAR's sub-lookups, of the discontiguous tables
//! (`build_disc_mb_lookup`) and of the tests are not kept.
//!
//! Memory. Between two batches the loop holds the longest table built so far, and no build
//! allocates a second table beside it (a longer one replaces it). A probe batch after a full batch
//! keeps the full batch's table resident where the table was freed before; its peak (the kept
//! table, the probe batch's `next_pos` and the rest of a small batch) stays below the full
//! batch's (the same table, a full batch's `next_pos` and the rest), so the search's peak is the
//! full batch's as before.
//!
//! The chains (EXPERIMENT LOSAT_X_MBCHAINKEEP; in the strict set from the pull request that added
//! it, 2026-10-27: on unless it is `0`). Beside the heads a table has its chains (NCBI's
//! `next_pos`, one element a query position: 4 bytes, or 3 in a table of LOSAT_X_MB24,
//! `x_mb_narrow.rs`), allocated zeroed per build (`calloc`, `blast_nalookup.c:979`) and freed with
//! the table after the preliminary stage (LOSAT_X_BNLUTDROP). Megablast of 31,168 D. melanogaster
//! cDNAs against a database of 3 Mb builds 18 tables of 4^12 cells with chains of up to 28.6 MiB:
//! once glibc's dynamic mmap threshold has risen past them (a freed mapped block raises it), each
//! batch's chains come from the heap of the main arena, other blocks fill the hole they leave, and
//! the next batch's chains grow the heap again, which glibc keeps (peak `RssAnon` 159.8 MiB at one
//! thread and about 170 at 16, against 131 and 146 with the threshold fixed at 1 MiB, the floor of
//! that retention). With the switch the keep is the chains' workspace too: a dropped table's chains
//! go back to it while a later table may follow in the loop (a later batch, `x_batches_follow`,
//! which the loop sets before each batch, or a later query chunk of a split batch,
//! `x_chunks_follow`, which the chunk loop sets before each chunk), and the next build takes them
//! as its `next_pos` without zeroing them, extended with zeros when it needs more elements (a
//! longer one is used as it is: the readers index it by offset). A search holds one block of
//! chains from its first build to its last: no hole, no regrowth. The last table frees its chains
//! as before, and when no table follows, the kept chains are freed after the preliminary stage of
//! the batch (`x_after_prelim`), so that the last traceback holds none. Between the builds the
//! workspace stays resident: through the tracebacks of the batches before the last, which held no
//! chains before (megablast of E. coli Sakai and 2 Mb of EDL933 against MG1655, two batches: peak
//! `RssAnon` +5.0 MiB at one thread, +8.1 at 16).
//!
//! Why the chains are NCBI's although they are not zeroed. A chain is read from a head (0 or the
//! 1-based offset of a word of this build) and then only at the offsets of the words it visits
//! (`next_pos[q_off]`), and every fill writes the element of every word's offset (the link to the
//! cell's previous head, or the re-threaded link of `ascending_cells`) before any reader runs, so
//! every element a reader reads holds the value of NCBI's `calloc`ed array; the others are never
//! read and may hold anything. The fills write an element whole: 4-byte elements by plain stores;
//! 3-byte elements on a pool as their exact 3 bytes (`store_exact`), serially as one 4-byte store
//! whose fourth byte, the first byte of the next element, comes before that element's own store
//! whenever that element is a word (the words come in increasing order of offset). The heads are
//! zeroed by the keep as before. No floating-point operation is involved.
//!
//! Where. The builds of the loop's keep (above: no LOSAT_X_QBATCHPAR slot can run), except a
//! table laid out as CSR (LOSAT_X_BNCSR: it has no chains, and the keep frees the kept chains with
//! the kept table) and a build with `BLEMIR_DEBUG` (its counters print the length of `next_pos`).
//! The discontiguous tables are not kept.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_mb_lut_keep: kept= fresh= zero= chains_kept=` gives the
//! loop's builds on the kept table, its builds on a fresh one, the time spent zeroing (also inside
//! `x_mb_lookup_split`'s `alloc=`) and the builds on kept chains.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// No NCBI counterpart: LOSAT_X_MBLUTKEEP (read once).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTKEEP"))
}

/// No NCBI counterpart: LOSAT_X_MBCHAINKEEP (read once; module documentation, "The chains").
pub(crate) fn chains_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBCHAINKEEP"))
}

/// No NCBI counterpart (scheduling only): the elements one job of the pool's zeroing writes at
/// least (256 KiB).
const ZERO_CHUNK: usize = 1 << 16;

static KEPT: AtomicU64 = AtomicU64::new(0);
static FRESH: AtomicU64 = AtomicU64::new(0);
static ZERO_NANOS: AtomicU64 = AtomicU64::new(0);
static CHAINS_KEPT: AtomicU64 = AtomicU64::new(0);

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
    /// LOSAT_X_MBCHAINKEEP: the words of the chains a dropped table gave back, their elements as
    /// the fills left them.
    chains: Mutex<Option<Vec<u32>>>,
    /// LOSAT_X_MBCHAINKEEP: whether a later batch of the loop follows the current one
    /// (`x_batches_follow`).
    batches_follow: AtomicBool,
    /// LOSAT_X_MBCHAINKEEP: whether a later query chunk of the current batch follows the current
    /// one (`x_chunks_follow`).
    chunks_follow: AtomicBool,
}

impl Keep {
    /// A keep for one batch loop, when LOSAT_X_MBLUTKEEP is on.
    pub(crate) fn for_loop() -> Option<Arc<Keep>> {
        on().then(|| Arc::new(Keep::default()))
    }

    /// A table for a build of `cells` cells, its first `cells` elements 0, and the number of its
    /// first elements that may hold values once the build has filled it: the kept one when it has
    /// at least that many elements (`true`; the build uses its first `cells`), with the elements of
    /// its first `cells` that earlier fills wrote set to 0 (the written elements past them stay as
    /// they are and stay counted); else a fresh `vec![0; cells]` (`false`), after the kept one is
    /// freed. `on_pool`: the build runs on a pool of two threads or more, which writes the zeros.
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
            return (vec![0u32; cells], false, cells);
        };
        let start = super::x_mb_lookup_par::timing_on().then(std::time::Instant::now);
        zero(&mut table[..dirty.min(cells)], on_pool);
        if let Some(start) = start {
            ZERO_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTKEEP", true);
        KEPT.fetch_add(1, Ordering::Relaxed);
        (table, true, dirty.max(cells))
    }

    /// Frees the kept table and chains, if any: a build laid out as CSR (LOSAT_X_BNCSR,
    /// `x_bn_csr.rs`) takes nothing from the keep, and the loop should not hold a table beside it;
    /// the next build with chains allocates a fresh table, which the keep takes back as before.
    pub(super) fn release(&self) {
        let old = self.kept.lock().unwrap_or_else(|e| e.into_inner()).take();
        drop(old);
        self.release_chains();
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

    /// LOSAT_X_MBCHAINKEEP (module documentation, "The chains"): whether a later batch of the
    /// loop follows the one it is about to search.
    pub(crate) fn x_batches_follow(&self, follow: bool) {
        self.batches_follow.store(follow, Ordering::Relaxed);
    }

    /// LOSAT_X_MBCHAINKEEP: whether a later query chunk of the batch follows the one about to be
    /// searched (`false` after the last).
    pub(crate) fn x_chunks_follow(&self, follow: bool) {
        self.chunks_follow.store(follow, Ordering::Relaxed);
    }

    /// LOSAT_X_MBCHAINKEEP: the end of a batch's preliminary stage: the kept chains are freed
    /// when no later table may take them.
    pub(crate) fn x_after_prelim(&self) {
        if !self.follow() {
            self.release_chains();
        }
    }

    /// LOSAT_X_MBCHAINKEEP: whether a later batch or chunk follows.
    fn follow(&self) -> bool {
        self.batches_follow.load(Ordering::Relaxed) || self.chunks_follow.load(Ordering::Relaxed)
    }

    fn release_chains(&self) {
        let old = self.chains.lock().unwrap_or_else(|e| e.into_inner()).take();
        drop(old);
    }

    /// LOSAT_X_MBCHAINKEEP: the words for a build's chains, `len` of them or more: the kept ones,
    /// extended to `len` with zeros when shorter (their other words as earlier fills left them),
    /// else `vec![0; len]` as before.
    pub(super) fn chains(&self, len: usize) -> Vec<u32> {
        let kept = self.chains.lock().unwrap_or_else(|e| e.into_inner()).take();
        let Some(mut chains) = kept else {
            return vec![0u32; len];
        };
        if chains.len() < len {
            chains.reserve_exact(len - chains.len());
            chains.resize(len, 0);
        }
        CHAINS_KEPT.fetch_add(1, Ordering::Relaxed);
        chains
    }

    /// LOSAT_X_MBCHAINKEEP: takes the chains of a dropped table when a later batch or chunk
    /// follows; otherwise leaves them to be freed with the table, as before.
    pub(super) fn give_back_chains(&self, chains: &mut Vec<u32>) {
        if !chains_on() {
            return;
        }
        let follow = self.follow();
        crate::utils::xreach::x_reach!("LOSAT_X_MBCHAINKEEP", follow);
        if follow {
            let old = self
                .chains
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .replace(std::mem::take(chains));
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
    if super::x_mb_lookup_par::timing_on() && load(&KEPT) + load(&FRESH) > 0 {
        eprintln!(
            "[TIMING] x_mb_lut_keep: kept={} fresh={} zero={:.4}s chains_kept={}",
            load(&KEPT),
            load(&FRESH),
            load(&ZERO_NANOS) as f64 * 1e-9,
            load(&CHAINS_KEPT)
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
    // cells take the words of their 3-byte heads, fewer than 4^12.
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
                    let ascending = case % 3 == 1;
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

    /// The words of the kept chains, if any (LOSAT_X_MBCHAINKEEP).
    fn kept_chains(keep: &Keep) -> Option<usize> {
        keep.chains.lock().unwrap().as_ref().map(Vec::len)
    }

    /// The walk of every cell of the first `cells` whose head is not 0, as the readers walk it
    /// (`for_each_hit` without the PV test; `x_mb_narrow::Narrow::each` for 3-byte elements); a
    /// walk of more than `limit` links fails.
    fn walks(table: &MbLookupTable, cells: usize, limit: usize) -> Vec<(usize, Vec<u32>)> {
        let heads: Vec<usize> = match table.x_narrow.as_ref() {
            Some(narrow) => {
                let words = &narrow.heads()[..x_mb_narrow::words_for(cells)];
                let mut heads = std::collections::BTreeSet::new();
                for (w, _) in words.iter().enumerate().filter(|(_, &word)| word != 0) {
                    heads.extend((4 * w / 3..=(4 * w + 3) / 3).filter(|&cell| cell < cells));
                }
                heads.into_iter().collect()
            }
            None => (0..cells)
                .filter(|&cell| table.hashtable[cell] != 0)
                .collect(),
        };
        heads
            .into_iter()
            .filter_map(|cell| {
                let mut walk = Vec::new();
                let mut push = |q_off: u32| {
                    walk.push(q_off);
                    assert!(walk.len() <= limit, "cell {cell}: more than {limit} links");
                };
                match table.x_narrow.as_ref() {
                    Some(narrow) => narrow.each(cell, push),
                    None => {
                        let mut q_off = table.hashtable[cell];
                        while q_off != 0 {
                            push(q_off);
                            q_off = table.next_pos[q_off as usize];
                        }
                    }
                }
                (!walk.is_empty()).then_some((cell, walk))
            })
            .collect()
    }

    // LOSAT_X_MBCHAINKEEP (module documentation, "The chains"): random queries (ambiguity codes,
    // masks) of tables of 4^9, 4^11 and 4^12 cells (3-byte elements, LOSAT_X_MB24, and 4-byte
    // ones: 4^9, 4^11 and 4^12 in `ascending_cells` order) built in a row through one keep while a
    // later batch follows, at one thread (the serial fills: 3-byte links stored as 4 bytes) and on
    // pools of 2 and 5 threads (the pool fills: exact 3-byte stores). The first build takes kept
    // chains of every byte 0xFF (each element past every chain) and fewer words than it needs
    // (extended with zeros), the others the chains the build before filled, of the other width or
    // longer: every walk, every head and the PV array equal the fresh build's, and the chains go
    // back. A build with no later batch or chunk frees its chains with the table, a later chunk
    // keeps them, and the end of a batch's preliminary stage frees the kept ones when nothing
    // follows. With the switch off no chains are kept.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn builds_on_kept_chains_walk_the_fresh_chains() {
        const LUTS: [(usize, bool); 6] = [
            (12, false),
            (11, false),
            (12, true),
            (9, false),
            (12, false),
            (11, true),
        ];
        let mut rng = Lcg(0x6368_6e73);
        let cases: Vec<_> = LUTS
            .iter()
            .map(|&(lut, ascending)| (lut, ascending, random_case(&mut rng)))
            .collect();
        for threads in [1usize, 2, 5] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let keep = Arc::new(Keep::default());
                keep.x_batches_follow(true);
                if chains_on() {
                    *keep.chains.lock().unwrap() = Some(vec![u32::MAX; 1_000]);
                }
                let build = |case: usize, keep: Option<&Arc<Keep>>| {
                    let (lut, ascending, (queries, offsets, masks, word_length, _)) = &cases[case];
                    let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                    let next_pos_len = queries
                        .iter()
                        .zip(offsets)
                        .map(|(q, &o)| o.max(0) as usize + q.len())
                        .max()
                        .unwrap_or(0)
                        + 1;
                    let table = build_mb_lookup(
                        &refs,
                        offsets,
                        (*word_length).max(*lut),
                        *lut,
                        masks,
                        None,
                        0,
                        1_000_000,
                        *ascending,
                        keep,
                    );
                    (table, 1usize << (2 * lut), next_pos_len)
                };
                let mut taken = 0;
                for case in 0..cases.len() {
                    let what = format!("case {case} threads {threads}");
                    let (fresh, cells, next_pos_len) = build(case, None);
                    assert!(fresh.x_csr.is_none(), "{what}: a CSR table");
                    let before = kept_chains(&keep);
                    let (built, _, _) = build(case, Some(&keep));
                    let narrow = built.x_narrow.as_ref();
                    assert_eq!(narrow.is_some(), fresh.x_narrow.is_some(), "{what}");
                    let (chain, need, head_words) = match narrow {
                        Some(narrow) => (
                            narrow.chain(),
                            x_mb_narrow::words_for(next_pos_len),
                            x_mb_narrow::words_for(cells),
                        ),
                        None => (&built.next_pos[..], next_pos_len, cells),
                    };
                    match before {
                        Some(len) => {
                            assert_eq!(chain.len(), len.max(need), "{what}: the kept chains");
                            taken += 1;
                        }
                        None => assert_eq!(chain.len(), need, "{what}: fresh chains"),
                    }
                    assert!(
                        built.x_table()[..head_words] == fresh.x_table()[..head_words],
                        "{what}: heads"
                    );
                    assert!(built.pv_array == fresh.pv_array, "{what}: pv_array");
                    assert!(
                        walks(&built, cells, next_pos_len) == walks(&fresh, cells, next_pos_len),
                        "{what}: walks"
                    );
                    let len = chain.len();
                    drop(built);
                    assert_eq!(
                        kept_chains(&keep),
                        chains_on().then_some(len),
                        "{what}: given back"
                    );
                }
                if chains_on() {
                    assert_eq!(taken, cases.len(), "threads {threads}");
                }
                // The last batch: its chains go with its table.
                keep.x_batches_follow(false);
                drop(build(0, Some(&keep)));
                assert_eq!(kept_chains(&keep), None, "threads {threads}: last batch");
                // A later chunk of the last batch: they go back; after the last chunk and the
                // batch's preliminary stage they go.
                keep.x_chunks_follow(true);
                drop(build(1, Some(&keep)));
                assert_eq!(kept_chains(&keep).is_some(), chains_on(), "a chunk follows");
                keep.x_after_prelim();
                assert_eq!(kept_chains(&keep).is_some(), chains_on(), "a chunk follows");
                keep.x_chunks_follow(false);
                keep.x_after_prelim();
                assert_eq!(
                    kept_chains(&keep),
                    None,
                    "threads {threads}: nothing follows"
                );
            });
        }
    }
}
