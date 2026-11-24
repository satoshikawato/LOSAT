//! EXPERIMENT (LOSAT_X_MBLUTKEEP; in the strict set since 2026-10-19: on unless it is `0`)
//! : the array of one element per cell of a query batch's lookup table kept from one query batch
//! to the next of the batch loop, instead of a fresh allocation per batch that is freed at the end
//! of its batch: the hashtable of `build_mb_lookup` (the megablast table of the megablast and
//! blastn tasks, zeroed for the next build), the `offsets` of a table laid out as CSR
//! (LOSAT_X_BNCSR, `x_bn_csr.rs`; LOSAT_X_DCCSR, `x_dc_csr.rs`) and the chain heads of a
//! discontiguous table of one template (`build_disc_mb_lookup_kept`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262,1376
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
//! No build follows the preliminary stage of the loop's last batch, so the keep frees its table
//! there (`Keep::release`, called by `x_lut_drop::after_prelim` when the batch's tables drop); the
//! loop frees the `Keep` when it ends.
//!
//! A table laid out as CSR has no hashtable; its `offsets` (a cell's first position, one element
//! per cell and one or two more) is built in the kept table instead (`take_buffer`: the build
//! writes every element it uses, so the table is not zeroed for it) and goes back to the keep when
//! the batch drops the CSR (`give_back_buffer`); a discontiguous table of one template takes its
//! chain heads from the keep as `build_mb_lookup` takes its hashtable. These arrays do not depend
//! on the batch's letters, only on its table's cells (4^11 or 4^12), so one kept table serves the
//! search's batches whichever layout each one takes: `CBatchSizeMixer`'s full batches are laid out
//! as CSR and its probe batches as chains (dc-megablast of 31,168 cDNAs: 18 CSR builds and 19
//! chain builds of 4^11 cells; at 4^12 cells LOSAT_X_MBCSR's size rule can take the small batches
//! and leave the full ones as chains). A fresh table has room for two elements more than its cells, so a
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
//! and the kept table stays resident through it (16 MiB of 4^11 cells, 64 MiB of 4^12; what a
//! packed CSR wrote, less), except in the loop's last batch, whose traceback runs without it (no
//! build is left to take it: the traceback of NCBI's last batch holds no table either); the
//! peaks above are measured with the kept table resident. A search of one batch (a `-subject`
//! search of one genome) has the table only in its preliminary stage, and a peak in its traceback
//! falls by the table's resident pages (2026-11-17, line EB s29, E. coli Sakai x EDL933
//! `-subject`, VmHWM: blastn at one thread 95.8 -> 90.3 MiB, a packed CSR's start; dc-megablast at
//! 16 threads 147.9 -> 129.0 MiB, the chain heads of 4^11 cells). Megablast's peak is in its scan,
//! where the table is read, and does not move.
//!
//! EXPERIMENT (LOSAT_X_CSRQSKEEP, in the strict set from 2026-11-23: on unless it is `0`): the
//! keep also holds the `qs` of a discontiguous CSR whose `offsets` goes back to it (the CSR's
//! second array, one element per word of the batch's table: 37 MB for a full batch of the
//! dc-megablast of 31,168 cDNAs, above glibc's largest mmap threshold of 32 MiB, so each build
//! mapped it, faulted its pages in during the place pass and unmapped it when the batch dropped
//! the CSR; at one thread 91 % of that search's page faults were in the dc CSR build). NCBI
//! allocates the chains that `qs` stands for (`next_pos`) per batch and frees them with the table,
//! as it does the hashtable.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:685,979,1377
//!
//! A discontiguous CSR build that takes its `offsets` through the keep (`take_buffer`,
//! `build_disc_mb_lookup_kept`) also takes the kept `qs` (`take_qs`, an empty vector when none is
//! kept) and uses it where it allocated `qs` before (`qs_buffer`, the switch's decision point:
//! `x_dc_csr::build`, `x_disc_lookup_par::csr_on_pool` and `csr_on_pool_steps`) when its capacity
//! holds the elements the build allocates; else it frees the kept one first and allocates a fresh
//! one, as `take_buffer` does with a table too short. The CSR gives its `qs` back when it is
//! dropped (`x_dc_csr::Csr::x_kept_with_qs`, `give_back_qs`, after its `offsets`): kept when no
//! `qs` is kept or the kept one has less room, else freed with the CSR. No build follows the
//! loop's last tables, so before they drop (`x_lut_drop::after_prelim` of its last batch, or the
//! last query chunk of a split last batch) the keep frees a kept `qs` and takes none back from
//! then on (`close_qs`; `release` closes it too): the last CSR frees its `qs` with its table, where
//! and when it did before the switch, and a search of one batch runs as before it. The switch has
//! no rule of its own: it acts wherever a CSR of LOSAT_X_DCCSR takes its `offsets` through the
//! keep of LOSAT_X_MBLUTKEEP.
//!
//! The contiguous CSR builds (LOSAT_X_BNCSR, `x_bn_csr`, `x_mb_lookup_once`) keep their own `qs`
//! (`Csr::x_kept`: their `offsets` alone goes back to the keep). Their CSR tables are those of
//! 4^11 cells or fewer (at 4^12 cells only LOSAT_X_MBCSR's side, beside a large subject), while
//! megablast's full batches build tables of 4^12 cells: in the megablast of the 31,168 cDNAs the
//! CSR batches are `CBatchSizeMixer`'s small ones (packed heads of LOSAT_X_BNCSRPACK) beside full
//! batches of chains, so a kept `qs` (4 MiB for a packed head of 4^11 cells) idled through the
//! full batches' scans, where the search's peak is, and a `qs` that small is served again from
//! glibc's heap without faults anyway: with it (2026-11-23, line FD, workstation node 1, dev
//! 5541322b2 against its merge 141e00ebc, the 31,168 cDNAs against 3 Mb of the D. melanogaster
//! genome, `-subject`) megablast at 16 threads took 0.757 s against 0.748 s with VmHWM +2.90 %
//! (137,004 against 133,146 KiB; `=0` +0.10 %), at one thread 3.935 s against 3.916 s with max RSS
//! +2.5 %, and blastn at 16 threads 4.718 s against 4.651 s (VmHWM -1.6 %): no gain for more
//! memory. The discontiguous CSR batches are the full
//! ones (LOSAT_X_DCCSR's size rule) and their probe batches chains, so the kept `qs` goes from a
//! full batch to the next: dc-megablast at 16 threads x0.957 (VmHWM +0.11 %), 8 x0.967 (+0.28
//! %), one x0.990; 5,000 of the cDNAs at 16 threads x0.965 (+0.28 %).
//!
//! Why the CSR is the same. Every discontiguous build writes each element of `qs` that it uses
//! before it reads one: its place pass writes each of the places `0 .. words` once (the runs
//! partition them; on the pool each owner writes its own places); the elements past `words` are
//! never read (`assemble` sets the length to `words`). So a kept `qs` is not zeroed, and the CSR
//! is the CSR of a fresh one. No floating-point operation is involved.
//!
//! Memory. No build allocates a `qs` beside the kept one (a kept `qs` too small is freed before
//! the fresh allocation). Between the drop of a CSR and the next CSR build the loop holds the
//! `qs` that was freed at the drop before: through the batch's traceback (LOSAT_X_BNLUTDROP), as
//! the kept table, except in the loop's last batch (`close_qs`), and through the chain builds of
//! `CBatchSizeMixer`'s probe batches; a smaller batch's CSR uses the largest kept `qs` (the
//! resident pages that earlier builds wrote) where it allocated its own before. So the search's
//! peak stays the one of the batch that built the largest `qs`, which holds it itself, unless the
//! memory of another stage beside the kept `qs` (a traceback, a probe batch) passes it.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_mb_lut_keep: kept= fresh= csr= zero= released= qs_kept=
//! qs_fresh=` gives the loop's builds of a hashtable or chain heads on the kept table, its builds
//! on a fresh one (of either layout), its CSR builds whose `offsets` is the kept table, the time
//! spent zeroing (also inside `x_mb_lookup_split`'s `alloc=`), the kept tables freed after a
//! loop's last preliminary stage (`Keep::release`), and (LOSAT_X_CSRQSKEEP) its discontiguous CSR
//! builds whose `qs` is the kept one and those that allocated their `qs` (none kept, or the kept
//! one too small; both 0 with the switch off).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// No NCBI counterpart: LOSAT_X_MBLUTKEEP (read once).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTKEEP"))
}

/// No NCBI counterpart: LOSAT_X_CSRQSKEEP (read once).
fn qs_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_CSRQSKEEP"))
}

/// No NCBI counterpart (scheduling only): the elements one job of the pool's zeroing writes at
/// least (256 KiB).
const ZERO_CHUNK: usize = 1 << 16;

static KEPT: AtomicU64 = AtomicU64::new(0);
static FRESH: AtomicU64 = AtomicU64::new(0);
static CSR: AtomicU64 = AtomicU64::new(0);
static ZERO_NANOS: AtomicU64 = AtomicU64::new(0);
static RELEASED: AtomicU64 = AtomicU64::new(0);
static QS_KEPT: AtomicU64 = AtomicU64::new(0);
static QS_FRESH: AtomicU64 = AtomicU64::new(0);

/// The kept table: as long as the largest table built so far; its first `dirty` elements may
/// hold values of earlier fills (the longest start that builds wrote since it was last zeroed),
/// the others are 0.
struct Kept {
    table: Vec<u32>,
    dirty: usize,
}

/// EXPERIMENT (LOSAT_X_CSRQSKEEP): the kept `qs` (the one with the most room of the dropped CSRs
/// that gave theirs back, or an empty vector) and whether the keep takes no more (`closed`: the
/// loop's last tables are being dropped, `Keep::close_qs`).
#[derive(Default)]
struct KeptQs {
    qs: Vec<u32>,
    closed: bool,
}

/// The table one batch loop keeps between its batches (LOSAT_X_MBLUTKEEP).
#[derive(Default)]
pub(crate) struct Keep {
    kept: Mutex<Option<Kept>>,
    /// EXPERIMENT (LOSAT_X_CSRQSKEEP): the `qs` a discontiguous CSR build takes (`take_qs`).
    qs: Mutex<KeptQs>,
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

    /// EXPERIMENT (LOSAT_X_CSRQSKEEP): for a discontiguous CSR build that takes its `offsets`
    /// through the keep (`take_buffer`), the kept `qs` (an empty vector when none is kept), which
    /// the build passes to `qs_buffer`; `None` with the switch off (the build allocates its `qs` as
    /// before).
    pub(super) fn take_qs(&self) -> Option<Vec<u32>> {
        self.take_qs_when(qs_on())
    }

    /// `take_qs` with the switch's state given (`on`).
    fn take_qs_when(&self, on: bool) -> Option<Vec<u32>> {
        on.then(|| std::mem::take(&mut self.qs.lock().unwrap_or_else(|e| e.into_inner()).qs))
    }

    /// EXPERIMENT (LOSAT_X_CSRQSKEEP): takes back the `qs` of a CSR dropped after its `offsets` went
    /// back (`x_dc_csr::Csr`'s `Drop`, for a CSR made with `x_kept_with_qs`): kept when no `qs` is
    /// kept or the kept one has less room (which is freed), else left in `qs`, which the CSR frees;
    /// always left there with the switch off, as before the switch, and once the keep is closed
    /// (`close_qs`, `release`), so the loop's last CSR frees its `qs` itself, where it did before.
    pub(super) fn give_back_qs(&self, qs: &mut Vec<u32>) {
        self.give_back_qs_when(qs, qs_on());
    }

    /// `give_back_qs` with the switch's state given (`on`).
    fn give_back_qs_when(&self, qs: &mut Vec<u32>, on: bool) {
        if !on {
            return;
        }
        let mut slot = self.qs.lock().unwrap_or_else(|e| e.into_inner());
        if !slot.closed && slot.qs.capacity() < qs.capacity() {
            let old = std::mem::replace(&mut slot.qs, std::mem::take(qs));
            drop(slot);
            drop(old);
        }
    }

    /// EXPERIMENT (LOSAT_X_CSRQSKEEP): no CSR build follows: frees the kept `qs` and keeps no `qs`
    /// from here on (`give_back_qs`). Called before the loop's last tables drop (the last batch's,
    /// `x_lut_drop::after_prelim`, or the last query chunk's of a split last batch), so the last
    /// CSR frees its `qs` with its table, where and when it did before the switch. No effect with
    /// the switch off.
    pub(crate) fn close_qs(&self) {
        self.close_qs_when(qs_on());
    }

    /// `close_qs` with the switch's state given (`on`).
    fn close_qs_when(&self, on: bool) {
        if !on {
            return;
        }
        let mut slot = self.qs.lock().unwrap_or_else(|e| e.into_inner());
        slot.closed = true;
        let old = std::mem::take(&mut slot.qs);
        drop(slot);
        drop(old);
    }

    /// Frees the kept table, when the loop has no batch left to build one on it: after the
    /// preliminary stage of its last batch (`x_lut_drop::after_prelim`), where NCBI's traceback
    /// stage frees the batch's table, instead of when the loop ends; and closes the keep to `qs`
    /// (`close_qs`, LOSAT_X_CSRQSKEEP; NCBI frees the batch's chains with its table), which holds
    /// none after it. Whether a table was kept.
    pub(crate) fn release(&self) -> bool {
        self.close_qs();
        let kept = self.kept.lock().unwrap_or_else(|e| e.into_inner()).take();
        let released = kept.is_some();
        if released {
            RELEASED.fetch_add(1, Ordering::Relaxed);
        }
        drop(kept);
        released
    }
}

/// EXPERIMENT (LOSAT_X_CSRQSKEEP): the decision point of the switch, where a CSR build allocated its
/// `qs` with room for `len` elements: `kept` (`Keep::take_qs`) emptied when it is a kept `qs`
/// (not the empty vector of a keep that holds none) whose capacity holds them; else, after `kept`
/// is freed, a fresh `Vec::with_capacity(len)`; `None` (no keep, or the switch off) is that fresh
/// vector alone, as before the switch. Each build writes every element of it that it uses before
/// it reads one (module comment), so `kept` is not zeroed.
pub(super) fn qs_buffer(kept: Option<Vec<u32>>, len: usize) -> Vec<u32> {
    let Some(mut qs) = kept else {
        return Vec::with_capacity(len);
    };
    if qs.capacity() > 0 && qs.capacity() >= len {
        crate::utils::xreach::x_reach!("LOSAT_X_CSRQSKEEP", true);
        QS_KEPT.fetch_add(1, Ordering::Relaxed);
        qs.clear();
        return qs;
    }
    // The reference side: the loop keeps no `qs` (its first CSR build) or one too small for this
    // build.
    crate::utils::xreach::x_reach!("LOSAT_X_CSRQSKEEP", false);
    QS_FRESH.fetch_add(1, Ordering::Relaxed);
    drop(qs);
    Vec::with_capacity(len)
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
            "[TIMING] x_mb_lut_keep: kept={} fresh={} csr={} zero={:.4}s released={} qs_kept={} \
             qs_fresh={}",
            load(&KEPT),
            load(&FRESH),
            load(&CSR),
            load(&ZERO_NANOS) as f64 * 1e-9,
            load(&RELEASED),
            load(&QS_KEPT),
            load(&QS_FRESH)
        );
    }
}

#[cfg(test)]
#[path = "x_mb_lut_keep_tests.rs"]
mod tests;
