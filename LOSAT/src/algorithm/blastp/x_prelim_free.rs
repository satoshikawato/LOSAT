//! EXPERIMENT (LOSAT_X_BPPFREE): blastp's preliminary hit lists hold no more memory than their
//! HSPs need, and each query's list is freed once its composition-based redo is done.
//!
//! The preliminary stage keeps, for every query of a search, the HSP lists of up to the
//! preliminary hit list size of subjects (`BlastpHitList`), and the redo of the traceback stage
//! reads them query by query. Two things keep more memory than NCBI does:
//! - a subject's HSP list for a query is collected from all its preliminary HSPs and then reaped
//!   by the preliminary E-value (`reap_hsplist_by_evalue`, `s_Blast_HSPListReapByPrelimEvalue`);
//!   the vector keeps the room of the reaped HSPs for as long as the list is kept;
//! - the redo reads each query's list in place, and every list is kept until the search ends,
//!   while NCBI's redo (`Blast_RedoAlignmentCore_MT`, below) frees each HSP list it has read
//!   from the stream once it is redone. The peak of a block of queries with many matches is the
//!   end of the redo, where the preliminary lists of all queries and the redone alignments of
//!   all queries are held together.
//!
//! With the switch on, the vector of a reaped list gives back its spare room
//! (`Vec::shrink_to_fit`), and the serial redo frees a query's preliminary list after its last
//! match has been redone. Neither changes a value: the reaped list keeps the same HSPs in the
//! same order, and nothing reads a query's preliminary list after its redo (the redo of a query
//! reads only its own list; the early termination test reads the redone heaps of the batch,
//! `redone_matches`, not the preliminary lists). The parallel redo paths are unchanged (they read
//! the lists again when they replay the kept matches).
//!
//! Off unless `LOSAT_X_BPPFREE` is set. Not in the strict set. With `LOSAT_TIMING=1` a search
//! with the switch prints `[TIMING] x_bppfree: freed_lists=<n>` (the lists the serial redo freed).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3398-3399
//! ```c
//!         while (BlastHSPStreamRead(hsp_stream, &localMatch)
//!                 != kBlastHSPStream_Eof) {
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3734-3738
//! ```c
//! match_loop_cleanup:
//!                 if (seqSrc) {
//!                     localMatch = Blast_HSPListFree(localMatch);
//!                 } else {
//! ```

// No NCBI counterpart: reads the LOSAT_X_BPPFREE switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPPFREE): off unless set to any value.
pub(crate) fn x_bp_prelim_free() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BPPFREE").is_some())
}
