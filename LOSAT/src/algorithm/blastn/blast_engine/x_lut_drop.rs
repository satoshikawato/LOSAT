//! EXPERIMENT (LOSAT_X_BNLUTDROP; in the strict set since 2026-10-24: on unless it is `0`): a
//! query batch's lookup tables dropped when its preliminary stage ends, before its traceback,
//! instead of at the end of the batch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:201-211
//! ```c++
//!     if (is_phi) {
//!         _ASSERT(m_InternalData->m_LookupTable);
//!     ...
//!     }
//!     else
//!     {
//!         m_InternalData->m_LookupTable.Reset(NULL);
//!     }
//! ```
//!
//! What it changes. `search_query_batch` builds the batch's tables (`build_lookup_tables`: the
//! megablast table, the direct or the small table of the blastn tasks, the discontiguous tables),
//! scans the subjects with them (the preliminary stage) and then runs the traceback of the HSPs
//! that the preliminary stage kept. Without the switch the tables live until the batch returns, so
//! the traceback's own memory (the subjects' letters, the gapped alignments, the HSP lists) comes
//! on top of them; for megablast of a 5.5 Mb query (E. coli Sakai x MG1655, `-subject` and `-db`)
//! the batch's peak was in the traceback: the 64 MiB hashtable, the chains (`next_pos`, 4 bytes a
//! position of both strands, 42 MiB) and the traceback's 21 MiB (plan v4 line MF, `REPORT.md`
//! §2.3). With the switch the batch drops them where NCBI's traceback stage drops its table: the
//! `-subject` and `-db` traceback reads the kept HSPs (`process_subject` with the preliminary
//! hits, which goes past the scan), and the lookup tables are passed to `process_subject` by the
//! preliminary stage's calls only (`None` in the traceback's), so the compiler proves that no
//! traceback reads them. A LOSAT_X_MBLUTKEEP table goes back to the batch loop's keep when the
//! tables drop (`x_mb_lut_keep.rs`), as at the end of the batch before: the keep is zeroed and
//! taken by the next build as before; what the traceback no longer holds is the rest (the chains,
//! the PV array, and every table that is not kept). A batch split into query chunks has no table
//! (its chunks' calls stop after their preliminary stage and drop theirs when they return).
//!
//! Why the reports are the same. The tables are read only by the scan of the preliminary stage
//! (`collect_prelim_hits_for_chunk`); the batch's values derived from them before the scan (the
//! offset array size, the shape of the subject search) are plain numbers kept in their own
//! variables. Dropping memory changes no value that the search computes. No floating-point
//! operation is involved.
//!
//! Memory: the peak of a batch whose traceback held the tables on top of its own memory moves to
//! the scan (megablast of Sakai x MG1655, one thread: 135.7 -> 131.3 MiB of `RssAnon`); the slots
//! of LOSAT_X_QBATCHPAR, which run whole batches side by side, hold their tables for a shorter time.
//! Where glibc's dynamic mmap threshold rises when the chains are freed (a mapped block of up to 32
//! MiB), it now rises one traceback earlier, whose medium buffers then come from the heap: megablast
//! of 31,168 cDNAs at one thread reads 2.2 MiB more (none with the threshold fixed).
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_bn_lut_drop: batches=` gives the number of batches whose
//! tables were dropped before their traceback.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::algorithm::blastn::coordination::LookupTables;

/// No NCBI counterpart: LOSAT_X_BNLUTDROP (read once).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNLUTDROP"))
}

// No NCBI counterpart: the batches whose tables were dropped before their traceback.
static DROPPED: AtomicU64 = AtomicU64::new(0);

/// The end of a batch's preliminary stage: with the switch the batch's tables are dropped here and
/// `None` is returned; without it they are returned, for the caller to hold until the batch ends
/// (the previous lifetime).
pub(super) fn after_prelim(tables: LookupTables) -> Option<LookupTables> {
    if !on() {
        return Some(tables);
    }
    if tables.two_stage_lookup.is_some()
        || tables.pv_direct_lookup.is_some()
        || tables.na_lookup.is_some()
    {
        DROPPED.fetch_add(1, Ordering::Relaxed);
    }
    drop(tables);
    None
}

/// No NCBI counterpart: with LOSAT_TIMING, the `[TIMING] x_bn_lut_drop` line `main` prints at
/// exit (when a batch dropped its tables).
pub(super) fn print_summary() {
    let dropped = DROPPED.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && dropped > 0 {
        eprintln!("[TIMING] x_bn_lut_drop: batches={dropped}");
    }
}
