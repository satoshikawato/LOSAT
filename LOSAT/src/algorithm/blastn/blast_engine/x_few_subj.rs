//! EXPERIMENT (LOSAT_X_BNFEWSUBJ, a strict-set switch): a `-subject` search with fewer 5 Mb work
//! units than threads searches its subjects one at a time, and a subject's batch of fewer chunks
//! than threads its chunks one at a time, each on the whole pool with the subject scan in strips
//! (`x_pair_seed.rs`, LOSAT_X_PAIRPAR and LOSAT_X_PAIRPARDISC).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//!        Int4 stat_length;
//!        if (seq_arg.oid == BLAST_SEQSRC_ERROR) {
//!            status = BLASTERR_SEQSRC;
//!            break;
//!        }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-484
//! ```c
//!     while (TRUE) {
//!         status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
//!                                        dbseq_chunk_overlap);
//!         if (status == SUBJECT_SPLIT_DONE) break;
//!         if (status == SUBJECT_SPLIT_NO_RANGE) continue;
//! ```
//!
//! Without the switch, a search of two subjects or more at two threads or more searches each
//! subject as one task of the pool (`use_parallel`, `search_subjects`), and a subject of two
//! chunks or more (`MAX_DBSEQ_LEN`) searches the chunks of each batch as tasks of the pool
//! (`reuse_prelim_hits` false). The strips of the seed stage run only for a chunk that one subject's
//! owner searches by itself (`one_subject_owner`), so neither shape has them: the Sakai genome
//! against MG1655 and its 92 kb plasmid used two of 16 threads for MG1655's scan (dc-megablast
//! 2.22 s against 0.53 s with MG1655 alone), and against EDL933 (5.53 Mb, two chunks) two threads too
//! (line G `PROBE_1.md` #4).
//!
//! With the switch on, a search whose strips would run (`x_pair_seed::takes_owner_chunks`, on a
//! thread of the search pool) takes the shape of a one-subject search when its work is in fewer
//! units than threads (`serial_subjects`): the subjects in OID order on the calling thread
//! (`use_parallel` false), each searched as a single `-subject` is (strips, LOSAT_X_QGROUPPAR,
//! the speculative traceback). And in that shape, or with one subject, a batch of fewer chunks than
//! threads (`serial_chunks`) searches its chunks in order with `reuse_prelim_hits`, the statements
//! of the one-chunk batch and of the one-thread loop.
//!
//! Why the reports do not change: both shapes are ones LOSAT already runs with the same bytes as one
//! thread. A search of one subject at two threads or more is the first (the subject loop of
//! `search_subjects` is the same code for any number of subjects, and `sink` receives the subjects
//! in OID order in both shapes: the parallel one collects them and then hands them over in order).
//! A subject's result is a function of the subject and of the search's fixed inputs: the scratch
//! that the loop carries from one subject to the next is reset per chunk, as the one-thread loop
//! relies on. The chunks of a batch searched in order are the one-thread loop's statements (search
//! with `reuse_prelim_hits`, merge, `-subject_besthit`, buffer back), and the parallel batch merges
//! its chunks' lists in the same order. The parts that come on with the shape (strips, per-query
//! gapped stages, speculative traceback) have their own shadows and gates.
//!
//! The rule (scheduling only; no NCBI counterpart): `units` counts the 5 Mb chunks of all subjects,
//! as LOSAT_X_QBATCHPAR does (`x_qbatch_par_off_of`). With `units >= threads` the subject tasks
//! fill the pool already, and the shape stays.

use std::sync::OnceLock;

use super::MAX_DBSEQ_LEN;

// No NCBI counterpart: reads the LOSAT_X_BNFEWSUBJ switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BNFEWSUBJ): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNFEWSUBJ"))
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// The 5 Mb work units (`MAX_DBSEQ_LEN` chunks, without the overlap) of subjects of these lengths.
pub(super) fn units(lengths: impl Iterator<Item = usize>) -> usize {
    lengths.map(|length| length.div_ceil(MAX_DBSEQ_LEN)).sum()
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Whether a `-subject` search of `subjects` subjects of `units` work units on `threads` threads
/// searches its subjects one at a time on the whole pool: the switch on, `strips` (the seed stage's
/// strips would take an unmasked chunk of the search on this thread), two subjects or more, and
/// fewer units than threads.
pub(super) fn serial_subjects(subjects: usize, units: usize, threads: usize, strips: bool) -> bool {
    on() && strips && subjects >= 2 && units < threads
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Whether a batch of `chunks` chunks of one subject, which the subject's owner searches on a pool
/// of `threads` threads, searches its chunks one at a time: the switch on, `strips` (as for
/// `serial_subjects`, and the subject has no lowercase mask), two chunks or more, and fewer chunks
/// than threads.
pub(super) fn serial_chunks(chunks: usize, threads: usize, strips: bool) -> bool {
    on() && strips && chunks >= 2 && chunks < threads
}

#[cfg(test)]
mod tests {
    use super::*;

    // Units are the 5 Mb chunks of each subject, without the overlap (as LOSAT_X_QBATCHPAR counts).
    #[test]
    fn units_count_the_5_mb_chunks_of_each_subject() {
        assert_eq!(units([].into_iter()), 0);
        assert_eq!(units([1].into_iter()), 1);
        assert_eq!(units([MAX_DBSEQ_LEN].into_iter()), 1);
        assert_eq!(units([MAX_DBSEQ_LEN + 1].into_iter()), 2);
        // MG1655 and pO157; EDL933; both genomes.
        assert_eq!(units([4_641_652, 92_721].into_iter()), 2);
        assert_eq!(units([5_528_445].into_iter()), 2);
        assert_eq!(units([4_641_652, 5_528_445].into_iter()), 3);
    }

    // The boundaries of the rule (the switch is on by default in a test run, as in the strict set,
    // unless the test runs with LOSAT_X_DEFAULT=0 or LOSAT_X_BNFEWSUBJ=0).
    #[test]
    fn the_rule_needs_strips_two_items_and_fewer_units_than_threads() {
        if !on() {
            assert!(!serial_subjects(2, 2, 16, true));
            assert!(!serial_chunks(2, 16, true));
            return;
        }
        assert!(serial_subjects(2, 2, 3, true));
        assert!(!serial_subjects(2, 2, 2, true));
        assert!(!serial_subjects(2, 2, 16, false));
        assert!(!serial_subjects(1, 1, 16, true));
        assert!(!serial_subjects(16, 16, 16, true));
        assert!(serial_subjects(15, 15, 16, true));
        assert!(serial_chunks(2, 3, true));
        assert!(!serial_chunks(2, 2, true));
        assert!(!serial_chunks(1, 16, true));
        assert!(!serial_chunks(2, 16, false));
    }
}
