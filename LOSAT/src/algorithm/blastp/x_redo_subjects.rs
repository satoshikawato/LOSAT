//! The composition redo of a `-db` search reads the subject of each match from the database when
//! it redoes the match, as NCBI does, instead of a copy of every listed subject made before the
//! redo (graduated 2026-10-27; the switch was `LOSAT_X_BPREDODBSUBJ`). A copy of the residues of
//! every subject of the preliminary hit lists (the listed subjects) would, with the UniRef50
//! database, list 70,701 subjects of 50.8 MB of residues for the first query batch (91 queries) of
//! 100 at 32 threads, all in the main heap until the batch's reports end, where the search reaches
//! its peak of anonymous memory (v4-D12 s3 probe, workstation); the final hit lists hold 38,230 of
//! those subjects.
//!
//! `ListedSubjects::x_list` keeps the ordinal ids and the lengths of the listed
//! subjects (read from the index: the estimate of `LOSAT_X_BPREDOLPT`, `x_job_costs`, reads the
//! lengths before a window's matches are redone). Every redo path (the one-thread loop, the two
//! query-parallel paths and `x_redo_par.rs`) and the first Newton problems of the Newton lanes
//! read a match's subject with `SeqDb::read_sequence` (`XRedoSubjects::read`): the stored residues
//! that the copy holds, so the redo computes the same values. Each read is dropped after its match:
//! the first Newton problems read a match ahead of its redo, and the redo reads it again (keeping
//! the reads ahead for the redo raised the maximum resident set of the UniRef50 `uniref50.04`
//! search at 8 threads from 53.6 to 60.3 MiB, medians of 3, against 53.9 MiB for the copy: v4-D12
//! s5m3, workstation). The one-thread loop reads after the early termination test, as NCBI does,
//! and after the first Newton problems; the parallel paths read the matches that they compute ahead
//! of the test's replay. A read that finds a pyrrolysine (O) sets the subject's flag
//! (`ListedSubjects::x_has_o`).
//! The block that the read returns holds the sentinel bytes of the sequence file on both sides of
//! the residues, where the copy holds `NCBISTDAA_GAP`; the redo reads the residues only, and SEG
//! makes its own copy of a range (`build_subject_range_data`). After the redo,
//! `ListedSubjects::x_read_shown` reads the subjects of the final hit lists again for the reports
//! (NCBI reads them again when it formats the Seq-aligns): their lengths are kept for every listed
//! subject, and their residues where a report reads them: where every report is tabular rows of the
//! default fields or Seq-aligns, which read the letters only to recount the HSPs whose query or
//! subject has an O (`recount_pyrrolysine_identities`), a subject with a pyrrolysine (O) and every
//! subject of an HSP whose query has an O (the others are not read again: the redo of a match reads
//! its subject before the subject can reach a final hit list, so the flags of the subjects with an
//! O are set); every subject of the final hit lists when a report renders the alignments (outfmt 0,
//! XML, tabular fields other than the default ones, the hits of the library API); every listed
//! subject when the reports make records (`ListedSubjects::read_titles`: `LOSAT_X_BPLISTEDLEAN=0`,
//! ABI v1's titles). No report reads the letters of a subject left without residues.
//!
//! A reader error of the redo's reads is the error that `ListedSubjects::read` gives for the same
//! read (`database_error`); it is reached in the redo, where NCBI reaches its own reads, and only
//! for a subject that the redo reads. Only a database changed after its preliminary read can fail
//! there.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1381-1385, quoted above
//! `ListedSubjects::read`: NCBI's redo fetches each match's sequence from the database in the
//! encoding of the search.
//!
//! It acts for every task of BLASTP with `-db`.

use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::Result;

use super::encoding::EncodedProtein;
use crate::seqdb::{SeqDb, SequenceBlock};

/// The stored residue of pyrrolysine (O) (`blast_engine.rs` `STORED_O`).
const STORED_O: u8 = 26;

// No NCBI counterpart: how LOSAT's redo reaches its subjects.
/// The subjects that the redo reads, by the index that its hit lists name them with.
#[derive(Clone, Copy)]
pub(crate) enum XRedoSubjects<'a> {
    /// The encoded subject set (`-subject`).
    Encoded(&'a [EncodedProtein]),
    /// A database's listed subjects, by rank, read when a match is redone.
    Database(XDbSubjects<'a>),
}

/// A database's listed subjects as the redo reads them.
#[derive(Clone, Copy)]
pub(crate) struct XDbSubjects<'a> {
    pub(crate) db: &'a SeqDb,
    /// The listed ordinal ids, by rank (`ListedSubjects::ordinals`).
    pub(crate) ordinals: &'a [u32],
    /// Their numbers of residues, by rank (`ListedSubjects::x_lengths`).
    pub(crate) lengths: &'a [u32],
    /// The program that the refusals name.
    pub(crate) program: &'static str,
    /// Whether a read found a pyrrolysine (O), by rank (`ListedSubjects::x_has_o`). The redo of a
    /// match reads its subject before it can reach a final hit list.
    pub(crate) has_o: &'a [AtomicBool],
}

/// A subject as the redo reads it: its residues, without the sentinels.
pub(crate) enum XRedoSubject<'a> {
    Encoded(&'a [u8]),
    /// The sequence read from the database.
    Read(SequenceBlock),
}

impl XRedoSubject<'_> {
    /// The residues, without the sentinels.
    #[inline]
    pub(crate) fn residues(&self) -> &[u8] {
        match self {
            Self::Encoded(residues) => residues,
            Self::Read(block) => block.protein(0),
        }
    }
}

impl<'a> XRedoSubjects<'a> {
    /// The number of subjects.
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Encoded(subjects) => subjects.len(),
            Self::Database(subjects) => subjects.ordinals.len(),
        }
    }

    /// The number of residues of the subject `s_idx`, if there is one.
    pub(crate) fn aa_len(&self, s_idx: usize) -> Option<usize> {
        match self {
            Self::Encoded(subjects) => subjects.get(s_idx).map(|subject| subject.aa_len),
            Self::Database(subjects) => subjects.lengths.get(s_idx).map(|&length| length as usize),
        }
    }

    /// The subject `s_idx` (which the subjects hold).
    #[inline]
    pub(crate) fn read(&self, s_idx: u32) -> Result<XRedoSubject<'a>> {
        match self {
            Self::Encoded(subjects) => {
                let subject = &subjects[s_idx as usize];
                Ok(XRedoSubject::Encoded(
                    &subject.aa_seq[1..subject.aa_seq.len() - 1],
                ))
            }
            Self::Database(subjects) => subjects.read(s_idx).map(XRedoSubject::Read),
        }
    }
}

impl XDbSubjects<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1381-1385, quoted above
    // `ListedSubjects::read`: the redo fetches each match's sequence from the database.
    /// Reads the subject of rank `rank`.
    #[inline(never)]
    fn read(&self, rank: u32) -> Result<SequenceBlock> {
        let ordinal = self.ordinals[rank as usize];
        let block = self.db.read_sequence(ordinal).map_err(|error| {
            crate::algorithm::common::database::database_error(error, self.program)
        })?;
        if x_has_stored_o(block.protein(0)) {
            self.has_o[rank as usize].store(true, Ordering::Relaxed);
        }
        Ok(block)
    }
}

/// Whether stored residues hold a pyrrolysine (O).
pub(crate) fn x_has_stored_o(residues: &[u8]) -> bool {
    residues.contains(&STORED_O)
}
