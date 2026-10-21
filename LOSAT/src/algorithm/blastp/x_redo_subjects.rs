//! EXPERIMENT (LOSAT_X_BPREDODBSUBJ / LOSAT_X_BPREDODBSUBJSHADOW): the composition redo of a `-db`
//! search reads the subject of each match from the database when it redoes the match, as NCBI
//! does, instead of a copy of every listed subject made before the redo.
//!
//! Without the switch, `database::ListedSubjects::read` copies the residues of every subject of the
//! preliminary hit lists (the listed subjects) into an `EncodedProtein` before the redo, and the
//! redo and the reports read the copies. With the UniRef50 database, the first query batch (91
//! queries) of 100 at 32 threads lists 70,701 subjects of 50.8 MB of residues, all in the main heap
//! until the batch's reports end, where the search reaches its peak of anonymous memory (v4-D12 s3
//! probe, workstation); the final hit lists hold 38,230 of those subjects.
//!
//! With the switch on, `ListedSubjects::x_list` keeps the ordinal ids and the lengths of the listed
//! subjects (read from the index: the estimate of `LOSAT_X_BPREDOLPT`, `x_job_costs`, reads the
//! lengths before a window's matches are redone). Every redo path (the one-thread loop, the two
//! query-parallel paths and `x_redo_par.rs`) and the first Newton problems of `LOSAT_X_NEWTONLANES`
//! read a match's subject with `SeqDb::read_sequence` (`XRedoSubjects::read`): the stored residues
//! that the copy holds, so the redo computes the same values. The redo of a match takes the subject
//! that the first Newton problems read ahead for it (`XPrefetchedSubjects`, which drops what no
//! redo took at their next call), so a match is read once. The one-thread loop reads after the
//! early termination test, as NCBI does, or takes the subject that the first Newton problems read
//! ahead of the test; the parallel paths read the matches that they compute ahead of the test's
//! replay. A read that finds a pyrrolysine (O) sets the subject's flag (`ListedSubjects::x_has_o`).
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
//! its shadow mode, ABI v1's titles). No report reads the letters of a subject left without
//! residues.
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
//! `LOSAT_X_BPREDODBSUBJ` is a switch of the strict set: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. It acts for every task of BLASTP with `-db`.
//!
//! `LOSAT_X_BPREDODBSUBJSHADOW=1` (shadow mode, implies the switch): `ListedSubjects::x_list` makes
//! the copy of `read` too, first, and compares the lengths with it, and every read of the redo and
//! of `x_read_shown` is compared with the copy (the length and every residue), and a subject of the
//! final hit lists left without residues without a read has no O in the copy. The first difference
//! stops the run (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] BPREDODBSUBJ searches= lengths= reads= residues= identical` (`searches`: listings;
//! `lengths`: lengths compared; `reads`: reads compared; `residues`: residues compared).

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use super::encoding::EncodedProtein;
use crate::seqdb::{SeqDb, SequenceBlock};

/// The stored residue of pyrrolysine (O) (`blast_engine.rs` `STORED_O`).
const STORED_O: u8 = 26;

/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): where the redo of a `-db` search reads its subjects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum XRedoDbSubj {
    /// The copy of every listed subject (`ListedSubjects::read`).
    Off,
    /// The database, per match (`ListedSubjects::x_list`).
    On,
    /// `On`, with the copy of `Off` made too and compared.
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPREDODBSUBJ / LOSAT_X_BPREDODBSUBJSHADOW switches once;
// it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set; `LOSAT_X_BPREDODBSUBJSHADOW` (any value) selects the
/// shadow mode.
pub(crate) fn x_bp_redo_db_subj() -> XRedoDbSubj {
    static MODE: OnceLock<XRedoDbSubj> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPREDODBSUBJSHADOW").is_some() {
            XRedoDbSubj::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BPREDODBSUBJ") {
            XRedoDbSubj::On
        } else {
            XRedoDbSubj::Off
        }
    })
}

// No NCBI counterpart: how LOSAT's redo reaches its subjects.
/// The subjects that the redo reads, by the index that its hit lists name them with.
#[derive(Clone, Copy)]
pub(crate) enum XRedoSubjects<'a> {
    /// The encoded subject set (`-subject`), or the copy of a database's listed subjects
    /// (`database::ListedSubjects::encoded`).
    Encoded(&'a [EncodedProtein]),
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): a database's listed subjects, by rank, read when a match
    /// is redone.
    Database(XDbSubjects<'a>),
}

/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): a database's listed subjects as the redo reads them.
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
    /// The shadow mode's copy (`ListedSubjects::read`), by rank.
    pub(crate) shadow: Option<&'a [EncodedProtein]>,
}

/// A subject as the redo reads it: its residues, without the sentinels.
pub(crate) enum XRedoSubject<'a> {
    Encoded(&'a [u8]),
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the sequence read from the database.
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

// No NCBI counterpart: keeps reads for their second reader; no value is computed differently.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the subjects of a database that the first Newton problems of
/// `LOSAT_X_NEWTONLANES` read for the next matches (`x_newton_lanes_prefetch_matches`), which the
/// redo of those matches takes instead of reading them again.
#[derive(Default)]
pub(crate) struct XPrefetchedSubjects<'a> {
    read: Vec<(u32, XRedoSubject<'a>)>,
}

impl<'a> XPrefetchedSubjects<'a> {
    /// Drops the subjects read ahead that no redo took.
    pub(crate) fn clear(&mut self) {
        self.read.clear();
    }

    /// Reads the subject `s_idx` of `subjects` (a database's) ahead, for its redo.
    pub(crate) fn read_ahead(
        &mut self,
        subjects: &XRedoSubjects<'a>,
        s_idx: u32,
    ) -> Result<&XRedoSubject<'a>> {
        self.read.push((s_idx, subjects.read(s_idx)?));
        Ok(&self.read.last().expect("just pushed").1)
    }

    /// The subject `s_idx` for its redo: the one read ahead, or a read of `subjects`.
    #[inline]
    pub(crate) fn take_or_read(
        &mut self,
        subjects: &XRedoSubjects<'a>,
        s_idx: u32,
    ) -> Result<XRedoSubject<'a>> {
        match self.read.iter().position(|&(read, _)| read == s_idx) {
            Some(at) => Ok(self.read.swap_remove(at).1),
            None => subjects.read(s_idx),
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
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): reads the subject of rank `rank` (shadow mode: compared
    /// with the copy).
    #[inline(never)]
    fn read(&self, rank: u32) -> Result<SequenceBlock> {
        let ordinal = self.ordinals[rank as usize];
        let block = self.db.read_sequence(ordinal).map_err(|error| {
            crate::algorithm::common::database::database_error(error, self.program)
        })?;
        if x_has_stored_o(block.protein(0)) {
            self.has_o[rank as usize].store(true, Ordering::Relaxed);
        }
        if let Some(copy) = self.shadow {
            x_shadow_compare_residues(
                rank as usize,
                ordinal,
                &copy[rank as usize],
                block.protein(0),
            );
        }
        Ok(block)
    }
}

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_LENGTHS: AtomicU64 = AtomicU64::new(0);
static X_READS: AtomicU64 = AtomicU64::new(0);
static X_RESIDUES: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BPREDODBSUBJ searches={} lengths={} reads={} residues={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_LENGTHS.load(Ordering::Relaxed),
        X_READS.load(Ordering::Relaxed),
        X_RESIDUES.load(Ordering::Relaxed),
        if identical { " identical" } else { "" }
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJSHADOW): the shadow mode's line at exit.
pub(crate) fn x_bp_redo_db_subj_print_summary() {
    if x_bp_redo_db_subj() == XRedoDbSubj::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BPREDODBSUBJ {message}");
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJSHADOW): the residues read for the subject of rank `rank`
/// against its copy: the length and every residue.
pub(crate) fn x_shadow_compare_residues(
    rank: usize,
    ordinal: u32,
    copy: &EncodedProtein,
    residues: &[u8],
) {
    if copy.aa_len != residues.len() || copy.aa_seq[1..=copy.aa_len] != *residues {
        x_shadow_difference(format!(
            "rank {rank} (ordinal {ordinal}): {} residues read, {} in the copy, or another residue",
            residues.len(),
            copy.aa_len
        ));
    }
    X_READS.fetch_add(1, Ordering::Relaxed);
    X_RESIDUES.fetch_add(residues.len() as u64, Ordering::Relaxed);
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJSHADOW): the lengths of `ListedSubjects::x_list` against the
/// copy, rank by rank.
pub(crate) fn x_shadow_compare_lengths(ordinals: &[u32], lengths: &[u32], copy: &[EncodedProtein]) {
    X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    if lengths.len() != copy.len() || ordinals.len() != copy.len() {
        x_shadow_difference(format!(
            "ranks: lengths {} / copy {} / listed {}",
            lengths.len(),
            copy.len(),
            ordinals.len()
        ));
    }
    for (rank, (&length, subject)) in lengths.iter().zip(copy).enumerate() {
        if length as usize != subject.aa_len {
            x_shadow_difference(format!(
                "rank {rank} (ordinal {}): length {length}, {} in the copy",
                ordinals[rank], subject.aa_len
            ));
        }
    }
    X_LENGTHS.fetch_add(lengths.len() as u64, Ordering::Relaxed);
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDODBSUBJSHADOW): a subject of the final hit lists that the reports keep
/// without its residues because no read found an O has no O in its copy.
pub(crate) fn x_shadow_check_no_o(rank: usize, ordinal: u32, copy: &EncodedProtein) {
    if x_has_stored_o(&copy.aa_seq[1..=copy.aa_len]) {
        x_shadow_difference(format!(
            "rank {rank} (ordinal {ordinal}): no read found an O, the copy holds one"
        ));
    }
}

/// Whether stored residues hold a pyrrolysine (O).
pub(crate) fn x_has_stored_o(residues: &[u8]) -> bool {
    residues.contains(&STORED_O)
}
