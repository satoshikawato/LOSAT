//! EXPERIMENT (LOSAT_X_MBLUTKEEP; in the strict set since 2026-10-19: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set): the hashtable of `build_mb_lookup` (the megablast table
//! of the megablast and blastn tasks) kept from one query batch to the next of the batch loop, and
//! zeroed for the next table of the same size, instead of a fresh zeroed allocation per batch that
//! is freed at the end of its batch.
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
//! loop goes to it when its batch drops the table (`MbLookupTable`'s `Drop`), and the next build of
//! a table of the same number of cells takes it and writes 0 to every cell before the fill (on the
//! search pool when the build runs on a pool of two threads or more). The largest table is kept:
//! NCBI's `CBatchSizeMixer` alternates full batches with small probe batches (megablast of 31,168
//! D. melanogaster cDNAs: 18 tables of 4^12 cells and 19 of 4^8 or 4^11), and a smaller table is
//! allocated and freed as before while the large one waits. At the end of the loop the `Keep` and
//! its table are freed.
//!
//! Why the table is NCBI's. The fill reads a cell before it writes it (`hashtable[ecode] == 0`, the
//! cell's head) and relies on `calloc`'s zeros; a kept table has every cell set to 0 before the
//! fill starts, so the fill sees what it sees in a fresh table and writes the same `hashtable`
//! (the other arrays, `next_pos`, `pv_array` and the counters of `longest_chain`, are allocated as
//! before). The number of cells is the table's (`1 << 2 * lut_word_length`): a kept table is taken
//! only when its length is that number. No floating-point operation is involved.
//!
//! The prefault round. On a kept table the round of LOSAT_X_MBLUTPREFAULT does not run: it exists
//! to map the pages of a fresh table by writing them before the fill reads them (its module
//! documentation, `x_mb_lookup_once.rs`), and a kept table's pages are mapped already (written by
//! the zeroing just before). Its writes are zeros over zeros, so the arrays are the same either
//! way.
//!
//! Where. Only the batch loop of `run_in_pool` keeps a table, for the batches it searches itself
//! (and the query chunks of a split batch), and only while no LOSAT_X_QBATCHPAR slot can run
//! (the driver's phase is off): a slot holds a whole batch while others run, and QBATCHPAR's
//! budget counts one table per batch in flight, so a table kept between the loop's batches would
//! add to the slots' peak. The tables of LOSAT_X_PAIRPAR's sub-lookups, of the discontiguous
//! tables (`build_disc_mb_lookup`) and of the tests are not kept.
//!
//! Memory. Between two batches the loop holds the largest table built so far; during a batch
//! that builds a smaller table both are allocated, the kept one not touched. A full batch's peak
//! (its table, `next_pos` and the rest of the batch) is above a probe batch's peak with the kept
//! table beside it, so the search's peak is the full batch's as before; on the workstation the
//! max RSS of megablast of 31,168 cDNAs fell (248 -> 235 MiB at 16 threads), the table no longer
//! unmapped and mapped again between batches.
//!
//! Shadow and timing: with `LOSAT_X_MBLUTKEEPSHADOW` every build on a kept table is built again
//! on a fresh table (`vec![0; cells]`, the same fill) and the arrays and `longest_chain` are
//! compared (`[X_SHADOW] MBLUTKEEP tables= cells= diffs=0 identical` at exit; a difference
//! panics). With LOSAT_TIMING, `[TIMING] x_mb_lut_keep: kept= fresh= zero=` gives the builds on a
//! kept table, the builds of the loop on a fresh one, and the time spent zeroing kept tables
//! (also inside `x_mb_lookup_split`'s `alloc=`).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// No NCBI counterpart: LOSAT_X_MBLUTKEEP (read once).
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTKEEP"))
}

/// No NCBI counterpart: LOSAT_X_MBLUTKEEPSHADOW (read once).
pub(super) fn shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_MBLUTKEEPSHADOW").is_some())
}

/// No NCBI counterpart (scheduling only): the elements one job of the pool's zeroing writes at
/// least (256 KiB).
const ZERO_CHUNK: usize = 1 << 16;

static KEPT: AtomicU64 = AtomicU64::new(0);
static FRESH: AtomicU64 = AtomicU64::new(0);
static ZERO_NANOS: AtomicU64 = AtomicU64::new(0);
static SHADOW_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);

/// The table one batch loop keeps between its batches (LOSAT_X_MBLUTKEEP): the largest hashtable
/// its batches dropped, holding the values of its last fill.
#[derive(Default)]
pub(crate) struct Keep {
    table: Mutex<Option<Vec<u32>>>,
}

impl Keep {
    /// A keep for one batch loop, when LOSAT_X_MBLUTKEEP is on.
    pub(crate) fn for_loop() -> Option<Arc<Keep>> {
        on().then(|| Arc::new(Keep::default()))
    }

    /// A zeroed table of `cells` cells: the kept one when it has that many (`true`), else a fresh
    /// allocation (`false`), which leaves the kept one where it is. `on_pool`: the build runs on a
    /// pool of two threads or more, which writes the zeros.
    pub(super) fn take(&self, cells: usize, on_pool: bool) -> (Vec<u32>, bool) {
        let kept = {
            let mut slot = self.table.lock().unwrap_or_else(|e| e.into_inner());
            match slot.as_ref() {
                Some(table) if table.len() == cells => slot.take(),
                _ => None,
            }
        };
        let Some(mut table) = kept else {
            FRESH.fetch_add(1, Ordering::Relaxed);
            return (vec![0u32; cells], false);
        };
        let start = super::x_mb_lookup_par::timing_on().then(std::time::Instant::now);
        zero(&mut table, on_pool);
        if let Some(start) = start {
            ZERO_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        KEPT.fetch_add(1, Ordering::Relaxed);
        (table, true)
    }

    /// Takes back the table of a dropped `MbLookupTable`: kept when nothing is kept or the kept
    /// one has fewer cells, else freed.
    pub(super) fn give_back(&self, table: Vec<u32>) {
        let mut slot = self.table.lock().unwrap_or_else(|e| e.into_inner());
        if slot.as_ref().is_none_or(|kept| kept.len() < table.len()) {
            let old = slot.replace(table);
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

/// The shadow: the arrays (`hashtable`, `next_pos`, `pv_array`) and the longest chain (which the
/// counters of `helper_array` give) of a build on a kept table (`kept`) equal those of the same
/// build on a fresh table (`fresh`); a difference panics.
pub(super) fn shadow_check(kept: [&[u32]; 3], fresh: [&[u32]; 3], chains: (usize, usize)) {
    const NAMES: [&str; 3] = ["hashtable", "next_pos", "pv_array"];
    let mut cells = 0u64;
    for ((name, a), b) in NAMES.iter().zip(kept).zip(fresh) {
        assert_eq!(
            a.len(),
            b.len(),
            "LOSAT_X_MBLUTKEEPSHADOW: {name} has {} elements on the kept table, {} on a fresh one",
            a.len(),
            b.len()
        );
        if let Some(i) = (0..a.len()).find(|&i| a[i] != b[i]) {
            panic!(
                "LOSAT_X_MBLUTKEEPSHADOW: {name}[{i}] is {} on the kept table, {} on a fresh one",
                a[i], b[i]
            );
        }
        cells += a.len() as u64;
    }
    assert_eq!(
        chains.0, chains.1,
        "LOSAT_X_MBLUTKEEPSHADOW: longest_chain {} on the kept table, {} on a fresh one",
        chains.0, chains.1
    );
    SHADOW_TABLES.fetch_add(1, Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(cells, Ordering::Relaxed);
}

/// No NCBI counterpart: printed by `main` at exit (through `x_mb_lookup_print_summary`): the
/// shadow line with LOSAT_X_MBLUTKEEPSHADOW, and with LOSAT_TIMING the `[TIMING] x_mb_lut_keep`
/// line when a loop kept tables. They print counters and change no value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if shadow() {
        eprintln!(
            "[X_SHADOW] MBLUTKEEP tables={} cells={} diffs=0 identical",
            load(&SHADOW_TABLES),
            load(&SHADOW_CELLS)
        );
    }
    if super::x_mb_lookup_par::timing_on() && load(&KEPT) + load(&FRESH) > 0 {
        eprintln!(
            "[TIMING] x_mb_lut_keep: kept={} fresh={} zero={:.4}s",
            load(&KEPT),
            load(&FRESH),
            load(&ZERO_NANOS) as f64 * 1e-9
        );
    }
}

#[cfg(test)]
mod tests {
    use super::super::x_mb_lookup_par::tests::{random_case, Lcg};
    use super::super::{build_mb_lookup, MbLookupTable};
    use super::*;

    /// The arrays and fields of one build.
    fn arrays(t: &MbLookupTable) -> (Vec<u32>, Vec<u32>, Vec<u32>, usize, usize) {
        (
            t.hashtable.clone(),
            t.next_pos.clone(),
            t.pv_array.clone(),
            t.pv_array_bts,
            t.longest_chain,
        )
    }

    fn kept_len(keep: &Keep) -> Option<usize> {
        keep.table.lock().unwrap().as_ref().map(Vec::len)
    }

    // A keep takes back the largest table, gives it out zeroed only for a table of its size, and
    // leaves it where it is for any other size.
    #[test]
    fn the_keep_holds_the_largest_table_and_gives_it_out_zeroed() {
        let keep = Keep::default();
        let (table, kept) = keep.take(8, false);
        assert!(!kept);
        assert_eq!(table, vec![0; 8]);
        keep.give_back(vec![7; 8]);
        keep.give_back(vec![5; 4]);
        assert_eq!(kept_len(&keep), Some(8), "a smaller table is freed");
        let (table, kept) = keep.take(4, false);
        assert!(!kept);
        assert_eq!(table, vec![0; 4]);
        assert_eq!(
            kept_len(&keep),
            Some(8),
            "another size leaves the kept table"
        );
        keep.give_back(vec![9; 16]);
        assert_eq!(kept_len(&keep), Some(16), "a larger table replaces it");
        let (table, kept) = keep.take(16, false);
        assert!(kept);
        assert_eq!(table, vec![0; 16]);
        assert_eq!(kept_len(&keep), None);
        let (_, kept) = keep.take(16, false);
        assert!(!kept, "the kept table is out");
        keep.give_back(vec![3; 1 << 20]);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let (table, kept) = pool.install(|| keep.take(1 << 20, true));
        assert!(kept);
        assert!(table.iter().all(|&cell| cell == 0), "zeroed on the pool");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262
    // ```c
    //    mb_lt->hashtable = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
    // ```
    // Random queries (ambiguity codes, masks, both cell orders) of tables of 4^4 to 4^12 cells in
    // a row through one keep, one after another and on pools of 2 and 5 threads (the pool fill and
    // its prefault round, skipped on a kept table): every build equals the same build on a fresh
    // table, and the builds whose size the keep holds take the kept table, filled by the build
    // before.
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
                let mut reused = 0;
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
                    let cells = 1usize << (2 * lut);
                    let dirty = keep
                        .table
                        .lock()
                        .unwrap()
                        .as_ref()
                        .filter(|table| table.len() == cells)
                        .map(|table| table.iter().any(|&cell| cell != 0));
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
                    if let Some(dirty) = dirty {
                        assert_eq!(kept_len(&keep), None, "case {case}: the kept table taken");
                        reused += 1;
                        assert!(
                            dirty || case < LUTS.len(),
                            "case {case}: a filled table kept"
                        );
                    }
                    assert!(
                        arrays(&built) == arrays(&fresh),
                        "case {case} lut {lut} threads {threads}: the build on the keep differs"
                    );
                    drop(built);
                    assert!(kept_len(&keep) >= Some(cells), "case {case}: given back");
                }
                // 4^11 twice, then 4^12 (kept from then on) once in the first round and twice in
                // the second.
                assert_eq!(reused, 5, "threads {threads}: builds on a kept table");
            });
        }
    }
}
