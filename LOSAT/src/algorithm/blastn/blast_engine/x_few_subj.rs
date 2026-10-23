//! EXPERIMENT (LOSAT_X_BNFEWSUBJ, a strict-set switch): a `-subject` search whose longest 5 Mb
//! work unit is a large part of its letters searches its subjects one at a time, and a subject's
//! batch of chunks the same way its chunks one at a time, each on the whole pool with the subject
//! scan in strips (`x_pair_seed.rs`, PAIRPAR, LOSAT_X_PAIRPARDISC and
//! LOSAT_X_PAIRPARDISCARR).
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
//! (`reuse_prelim_hits` false). The strips of the seed stage run only for a chunk that one
//! subject's owner searches by itself (`one_subject_owner`), so neither shape has them: the Sakai
//! genome against MG1655 and its 92 kb plasmid used two of 16 threads for MG1655's scan
//! (dc-megablast 2.22 s against 0.53 s with MG1655 alone), and against EDL933 (5.53 Mb, two
//! chunks) two threads too (line G `PROBE_1.md` #4).
//!
//! With the switch on, a search whose strips would run (`x_pair_seed::takes_owner_chunks`, on a
//! thread of the search pool) takes the shape of a one-subject search when the rule below says so
//! (`serial_subjects`): the subjects in OID order on the calling thread (`use_parallel` false),
//! each searched as a single `-subject` is (strips, LOSAT_X_QGROUPPAR, the speculative traceback).
//! And in that shape, or with one subject, a batch of chunks for which the rule says so
//! (`serial_chunks`) searches its chunks in order with `reuse_prelim_hits`, the statements of the
//! one-chunk batch and of the one-thread loop.
//!
//! Why the reports do not change: both shapes are ones LOSAT already runs with the same bytes as
//! one thread. The subjects one at a time on the calling thread are the one-thread loop: the
//! subject loop of `search_subjects` is the same code for any number of subjects and threads, it
//! starts from one scratch (`new_state`) and searches the subjects in OID order, so every subject
//! meets the scratch history it meets at one thread (the diagonal hash, kept across chunks with an
//! advancing offset, included); only the stages inside a chunk use the pool, each with its own
//! exactness argument and shadow. `sink` receives the subjects in OID order in both shapes (the
//! parallel one hands them over in order: collected, or streamed with LOSAT_X_BNSTREAM). The
//! chunks of a batch searched in order are the one-thread loop's statements (search with
//! `reuse_prelim_hits`, merge, `-subject_besthit`, buffer back), and the parallel batch merges its
//! chunks' lists in the same order. The parts that come on with the shape (strips, per-query
//! gapped stages, speculative traceback) have their own shadows and gates.
//!
//! The rule (scheduling only; no NCBI counterpart): the subject tasks take at least the time of
//! their longest work unit (a 5 Mb chunk, `MAX_DBSEQ_LEN`, or a shorter subject: the chunks of a
//! subject's batch are tasks of the pool too), and the one-subject shape about the time of all the
//! letters over the strips' speed-up S(p) at p threads. The shape changes when the longest unit's
//! letters times Ŝ(p) are more than all the letters (`serial`), and a batch of chunks the same way
//! with its chunks. Ŝ(p) (`SPEEDUP_TENTHS`) is the lowest S(p) = one-thread wall / wall of the
//! shape that line G's probe measured for the kind of table, at the largest of 2, 4, 8 and 16
//! threads not above p, floored to tenths (workstation, node 1, 2026-10-20: Sakai against MG1655 +
//! pO157, EDL933, MG1655 + EDL933, and MG1655 in 3 and in 10 parts; dc-megablast 1.75 / 2.56 /
//! 3.24 / 3.40 with EDL933, blastn 1.74 / 2.94 / 4.18 / 4.75, megablast 1.51 / 2.06 / 2.55 / 2.74).
//! Every cell of the probe that the rule changes was faster in the new shape; MG1655 in 10 parts,
//! where the subject tasks fill the pool (x1.01-1.48 slower in the new shape), and 2,000 16S genes
//! (x5.2 slower) keep theirs.

use std::sync::OnceLock;

use super::MAX_DBSEQ_LEN;

// No NCBI counterpart: reads the LOSAT_X_BNFEWSUBJ switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BNFEWSUBJ): a strict-set switch, on unless it is `0`.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNFEWSUBJ"))
}

/// No NCBI counterpart: the kind of seed stage whose strips' speed-up the rule uses (a row of
/// `SPEEDUP_TENTHS`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Kind {
    /// A discontiguous table (dc-megablast).
    Disc,
    /// A contiguous table of `-task blastn`.
    Blastn,
    /// Any other contiguous table (megablast; the lowest row).
    Other,
}

impl Kind {
    // No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
    pub(super) fn of(disc: bool, task: &str) -> Kind {
        if disc {
            Kind::Disc
        } else if task == "blastn" {
            Kind::Blastn
        } else {
            Kind::Other
        }
    }
}

/// No NCBI counterpart: Ŝ(p) in tenths per `Kind` at 2, 4, 8 and 16 threads (see the module
/// comment).
const SPEEDUP_TENTHS: [[u64; 4]; 3] = [[17, 25, 32, 34], [17, 29, 41, 47], [15, 20, 25, 27]];

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Ŝ(p) in tenths: the row of `kind` at the largest of 2, 4, 8 and 16 threads not above
/// `threads`, and 10 (no speed-up) below two threads.
fn speedup_tenths(kind: Kind, threads: usize) -> u64 {
    let column = match threads {
        0 | 1 => return 10,
        2..=3 => 0,
        4..=7 => 1,
        8..=15 => 2,
        _ => 3,
    };
    SPEEDUP_TENTHS[kind as usize][column]
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Whether work of `letters` letters whose longest unit has `longest` takes the one-subject shape
/// on `threads` threads: `longest` x Ŝ(p) > `letters` (the switch on).
fn serial(longest: usize, letters: usize, threads: usize, kind: Kind) -> bool {
    on() && (longest as u64).saturating_mul(speedup_tenths(kind, threads))
        > (letters as u64).saturating_mul(10)
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Whether a `-subject` search of subjects of these lengths on `threads` threads searches them one
/// at a time on the whole pool (the caller has checked that the strips would run): two subjects or
/// more, and `serial` with the longest unit = the longest subject up to `MAX_DBSEQ_LEN`.
pub(super) fn serial_subjects(
    lengths: impl Iterator<Item = usize>,
    threads: usize,
    kind: Kind,
) -> bool {
    let (mut subjects, mut longest, mut letters) = (0usize, 0usize, 0usize);
    for length in lengths {
        subjects += 1;
        longest = longest.max(length.min(MAX_DBSEQ_LEN));
        letters = letters.saturating_add(length);
    }
    subjects >= 2 && serial(longest, letters, threads, kind)
}

// No NCBI counterpart: the shape of the search (scheduling); it changes no value NCBI computes.
/// Whether a batch of chunks of these lengths (with their overlaps), which one subject's owner
/// searches on a pool of `threads` threads, searches its chunks one at a time (the caller has
/// checked that the strips would run and that the subject has no lowercase mask): two chunks or
/// more, and `serial` with the longest chunk.
pub(super) fn serial_chunks(
    lengths: impl Iterator<Item = usize>,
    threads: usize,
    kind: Kind,
) -> bool {
    let (mut chunks, mut longest, mut letters) = (0usize, 0usize, 0usize);
    for length in lengths {
        chunks += 1;
        longest = longest.max(length);
        letters = letters.saturating_add(length);
    }
    chunks >= 2 && serial(longest, letters, threads, kind)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MG1655: usize = 4_641_652;
    const PO157: usize = 92_721;
    const EDL933: usize = 5_528_445;

    // Ŝ(p) is the row of the kind at the largest probed thread count not above p.
    #[test]
    fn speedup_takes_the_largest_probed_thread_count_not_above_p() {
        assert_eq!(speedup_tenths(Kind::Disc, 1), 10);
        assert_eq!(speedup_tenths(Kind::Disc, 2), 17);
        assert_eq!(speedup_tenths(Kind::Disc, 3), 17);
        assert_eq!(speedup_tenths(Kind::Blastn, 4), 29);
        assert_eq!(speedup_tenths(Kind::Blastn, 15), 41);
        assert_eq!(speedup_tenths(Kind::Other, 16), 27);
        assert_eq!(speedup_tenths(Kind::Other, 128), 27);
        assert_eq!(Kind::of(true, "dc-megablast"), Kind::Disc);
        assert_eq!(Kind::of(false, "blastn"), Kind::Blastn);
        assert_eq!(Kind::of(false, "megablast"), Kind::Other);
        assert_eq!(Kind::of(false, "blastn-short"), Kind::Other);
    }

    // The cells of the probe (the switch is on by default in a test run, as in the strict set,
    // unless the test runs with LOSAT_X_BNFEWSUBJ=0).
    #[test]
    fn the_rule_follows_the_probe() {
        let subjects = |lengths: &[usize], threads, kind| {
            serial_subjects(lengths.iter().copied(), threads, kind)
        };
        let chunks = |lengths: &[usize], threads, kind| {
            serial_chunks(lengths.iter().copied(), threads, kind)
        };
        if !on() {
            assert!(!subjects(&[MG1655, PO157], 16, Kind::Disc));
            assert!(!chunks(
                &[MAX_DBSEQ_LEN, EDL933 - MAX_DBSEQ_LEN],
                16,
                Kind::Disc
            ));
            return;
        }
        // One subject, or one thread: never.
        assert!(!subjects(&[MG1655], 16, Kind::Disc));
        assert!(!chunks(&[MAX_DBSEQ_LEN], 16, Kind::Disc));
        assert!(!subjects(&[MG1655, PO157], 1, Kind::Disc));
        // MG1655 + pO157 and EDL933's two chunks: every kind from two threads.
        for kind in [Kind::Disc, Kind::Blastn, Kind::Other] {
            assert!(subjects(&[MG1655, PO157], 2, kind));
            assert!(chunks(&[MAX_DBSEQ_LEN, EDL933 - MAX_DBSEQ_LEN], 2, kind));
        }
        // MG1655 + EDL933 (longest unit a 5 Mb chunk): from 4 threads, megablast from 8.
        assert!(!subjects(&[MG1655, EDL933], 3, Kind::Disc));
        assert!(subjects(&[MG1655, EDL933], 4, Kind::Disc));
        assert!(subjects(&[MG1655, EDL933], 4, Kind::Blastn));
        assert!(!subjects(&[MG1655, EDL933], 7, Kind::Other));
        assert!(subjects(&[MG1655, EDL933], 8, Kind::Other));
        // MG1655 in 3 parts: dc-megablast and blastn from 8 threads, megablast never.
        let third = [MG1655 / 3 + 1, MG1655 / 3, MG1655 / 3];
        assert!(!subjects(&third, 7, Kind::Disc));
        assert!(subjects(&third, 8, Kind::Disc));
        assert!(subjects(&third, 8, Kind::Blastn));
        assert!(!subjects(&third, 64, Kind::Other));
        // MG1655 in 10 parts, a 10 Mb subject's two equal chunks at two threads: never.
        assert!(!subjects(&[MG1655 / 10; 10], 64, Kind::Blastn));
        assert!(!chunks(&[MAX_DBSEQ_LEN, MAX_DBSEQ_LEN], 2, Kind::Disc));
        assert!(chunks(&[MAX_DBSEQ_LEN, MAX_DBSEQ_LEN], 4, Kind::Disc));
        // Letters past u64 do not wrap.
        assert!(!subjects(&[usize::MAX, usize::MAX], 16, Kind::Disc));
    }
}
