//! A `-db` search of TBLASTN reads each database sequence where it searches it (graduated
//! 2026-10-27; the switch was `LOSAT_X_TBNDBSTREAM`), instead of a copy of the whole database made
//! before the first query batch. This module holds TBLASTN's stages over the stream; the stream
//! itself, with the check, the reads and the exactness argument, is `common/x_db_stream.rs`.
//!
//! The copy (`read_stores`, deleted) packed every sequence of the database into anonymous memory
//! before the first query batch, and every stage read the copy until the reports: 23.9 MiB for
//! the 100 Mb of the C. elegans genome, 32.2 MiB for a zebrafish cDNA set, alive at the peak of
//! every thread count (v4-D12 s8 probe, workstation), where NCBI maps the sequence file. The
//! reports then cloned the packed sequences they show (`database::reported_subjects`).
//!
//! TBLASTN's stages over the stream:
//!
//! - The preliminary stage of a subject reads its sequence when it makes the subject's plan
//!   (`search_seed::x_scan_plan`, `XPrelimSubject::Stream`) and drops it once the plan is made:
//!   the frames are translated from it and the stage reads nothing else of the sequence. A plan
//!   that the search keeps between query batches (`XSubjectSide`) is not made again, and its
//!   sequence is not read again. The plans of short subjects are made from the bytes of their
//!   block (`XDbStream::with_block`), without a store: the preliminary stage asks for the subjects
//!   of a range in ascending order on each thread, so most short subjects come without a read or a
//!   copy.
//! - The composition redo and the traceback of `-comp_based_stats 0` read a subject of their
//!   lists when they first align it and hold it until the stage has done its last match
//!   (`XStageSubjects`, LOSAT_X_TBNSTAGEDROP; with `=0` until the stage ends); the reports read
//!   the sequences they show, one read each, when a report reads their letters
//!   (`database::reported_subjects`; the other reports read the lengths of the stream).
//! - A sequence whose plan has no residues (LOSAT_X_TBNUNITXLATE: a first frame of 1 MiB or more,
//!   in a search that keeps no plan of it) is read before its plan and held while its units
//!   translate it, as a packed subject of the copy is (`search_seed::x_lean_subject`,
//!   `search_gapped::x_preliminary_protein_hsps_by_unit`).
//! - While a stream is open (`x_frames_pooled`), the frames of 1 MiB or more of a dropped subject
//!   plan go to `X_BYTES` and the next frames and long reads take them with their pages resident
//!   (`x_subject_translate::x_frame_buffer`, `x_reuse_frame`), as BLASTN's database reads do. The
//!   frames go idle when their subject is searched, after the peak of the preliminary search, so
//!   the composition redo and the reports read their sequences into them, and at one thread each
//!   chromosome's frames are those of the one before (C. elegans genome: 6 chromosomes of 13.8 to
//!   20.9 Mb, frames of 4.6 to 7 MB). A frame takes an idle buffer that touched at most twice its
//!   bytes, since the pool does not count the pages of a lent buffer. No store is given to the
//!   pool: the stores read for the plans of 8 and 16 threads, made at once, would stay idle across
//!   the peak.
//!
//! A read that fails during the search is an `XReadError` (`common/x_db_stream.rs`): the search of
//! the batches is then made again one batch at a time, so that the batch that reads it fails after
//! the reports of the batches before, as with the check's error (`tblastn_app.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1380-1386 (the redo fetches
//! each match's subject)
//! ```c
//!         if( program_number == eBlastTypeTblastn ) {
//!             seq_info->seq_arg.encoding = eBlastEncodingNcbi4na;
//!         } else {
//!             seq_info->seq_arg.encoding = eBlastEncodingProtein;
//!         }
//!         if (BlastSeqSrcGetSequence(seqSrc, &seq_info->seq_arg) >= 0) {
//! ```
//!
//! It acts for every TBLASTN search with `-db`. The shadow mode (`LOSAT_X_TBNDBSTREAMSHADOW=1`,
//! which compared the check and every read with the copy) went with the copy.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Result;

use crate::algorithm::common::x_db_stream::XDbStream;
use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};

// No NCBI counterpart: reads the LOSAT_X_TBNSTAGEDROP variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSTAGEDROP): a switch of the strict set (`utils/xswitch.rs`): a stage's
/// subject is let go after its last match (`XStageSubjects::done`) instead of at the stage's end.
pub(super) fn x_tbn_stage_drop() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSTAGEDROP"))
}

/// The subjects that a stage after the preliminary search reads
/// (the composition redo, the traceback of `-comp_based_stats 0`): with a `-db` stream, each
/// subject of the stage's lists read when the stage first aligns it (`XDbStream::read_shared`)
/// and held until the stage has passed its last match (LOSAT_X_TBNSTAGEDROP, `done`) or, with
/// `LOSAT_X_TBNSTAGEDROP=0`, until the stage ends; otherwise the search's subjects.
pub(super) enum XStageSubjects<'a> {
    All(super::x_subject_store::SubjectLetters<'a>),
    Read {
        stream: &'a XDbStream<'a>,
        /// The subjects of the lists, ascending.
        oids: Vec<usize>,
        /// The store of each, once read, until its last match is done.
        stores: Vec<Mutex<Option<Arc<PackedNucleotide>>>>,
        /// The matches of each that are not done (LOSAT_X_TBNSTAGEDROP); empty with the switch
        /// off, where the stores stay until the stage ends.
        left: Vec<AtomicUsize>,
    },
}

/// A subject of a stage, held while a match reads it (`XStageSubjects::held`).
pub(super) enum XHeld<'a> {
    Nt(SubjectNt<'a>),
    Store(Arc<PackedNucleotide>),
}

impl XHeld<'_> {
    /// The subject for the readers of windows.
    pub(super) fn nt(&self) -> SubjectNt<'_> {
        match self {
            Self::Nt(nt) => *nt,
            Self::Store(store) => SubjectNt::Packed(store),
        }
    }
}

impl<'a> XStageSubjects<'a> {
    /// The subjects of a stage whose lists name the subjects `oids` (one item per match).
    pub(super) fn new(
        subjects: super::x_subject_store::SubjectLetters<'a>,
        oids: impl Iterator<Item = usize>,
    ) -> Self {
        let super::x_subject_store::SubjectLetters::Stream(stream) = subjects else {
            return Self::All(subjects);
        };
        let mut all: Vec<usize> = oids.collect();
        all.sort_unstable();
        let mut oids: Vec<usize> = Vec::new();
        let mut counts: Vec<usize> = Vec::new();
        for oid in all {
            if oids.last() == Some(&oid) {
                *counts.last_mut().expect("a count per subject") += 1;
            } else {
                oids.push(oid);
                counts.push(1);
            }
        }
        let stores = oids.iter().map(|_| Mutex::new(None)).collect();
        let left = if x_tbn_stage_drop() {
            counts.into_iter().map(AtomicUsize::new).collect()
        } else {
            Vec::new()
        };
        Self::Read {
            stream,
            oids,
            stores,
            left,
        }
    }

    /// Subject `oid` for the readers of windows (`SubjectLetters::nt`), read the first time.
    pub(super) fn held(&self, oid: usize) -> Result<XHeld<'_>> {
        match self {
            Self::All(subjects) => Ok(XHeld::Nt(subjects.nt(oid))),
            Self::Read {
                stream,
                oids,
                stores,
                ..
            } => {
                let at = oids
                    .binary_search(&oid)
                    .expect("every subject of the stage's lists is listed");
                let mut slot = stores[at].lock().unwrap_or_else(|e| e.into_inner());
                if let Some(store) = slot.as_ref() {
                    return Ok(XHeld::Store(Arc::clone(store)));
                }
                // Two threads may ask at once: the stream reads it once (`read_shared`).
                drop(slot);
                let store = stream.read_shared(oid)?;
                slot = stores[at].lock().unwrap_or_else(|e| e.into_inner());
                Ok(XHeld::Store(Arc::clone(slot.get_or_insert(store))))
            }
        }
    }

    /// A guard of one match of subject `oid` that calls `done` when it is dropped. Declared before
    /// the match's `held`, it runs after the match let the subject go.
    pub(super) fn x_match(&self, oid: usize) -> XMatchDone<'_, 'a> {
        XMatchDone { stage: self, oid }
    }

    /// One match of subject `oid` is done (redone, or passed by the early termination): after its
    /// last one the stage lets the subject go (LOSAT_X_TBNSTAGEDROP); a match that reads it later
    /// would read it again. No NCBI counterpart: NCBI's redo fetches each match's subject and
    /// releases it after the match (`blast_kappa.c` `s_MatchingSequenceRelease`).
    pub(super) fn done(&self, oid: usize) {
        let Self::Read {
            oids, stores, left, ..
        } = self
        else {
            return;
        };
        if left.is_empty() {
            return;
        }
        let at = oids
            .binary_search(&oid)
            .expect("every subject of the stage's lists is listed");
        let before = left[at].fetch_sub(1, Ordering::AcqRel);
        debug_assert!(
            before > 0,
            "more matches done than listed for subject {oid}"
        );
        if before == 1 {
            let store = stores[at].lock().unwrap_or_else(|e| e.into_inner()).take();
            drop(store);
        }
    }
}

/// Calls `XStageSubjects::done` for one match when it is dropped (`XStageSubjects::x_match`).
pub(super) struct XMatchDone<'s, 'a> {
    stage: &'s XStageSubjects<'a>,
    oid: usize,
}

impl Drop for XMatchDone<'_, '_> {
    fn drop(&mut self) {
        self.stage.done(self.oid);
    }
}
