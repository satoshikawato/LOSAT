//! The Bioseq that NCBI's places give for a subject's Seq-align id where that is another database
//! sequence, the holder (line DBID hr-2): every record writer reads it for the subject.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:380-420
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:350-420
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:120-180
//! The writers take the subject's Bioseq handle from the scope by the Seq-align id
//! (`GetBioseqHandle`): where the scope's places ([`crate::report::db_subject::ScopeWalk`]) give the
//! holder's blob, the names, titles, taxids, length and letters that they write are the
//! holder's; the Seq-align's coordinates, scores and gaps stay the subject's. A program writes it
//! by putting this record, length and letters in the subject's place of its reported subjects
//! (`algorithm::common::database::ShownHolders`), so every writer of the program reads them
//! through the one table it already reads.

use std::sync::Arc;

use crate::report::db_subject::{HitLoad, ParsedSubject};
use crate::seqdb::DbSequence;
use crate::seqid::ParsedSeqId;

/// The holder's Bioseq as the writers read it for one subject: one per (holder, B-line, Seq-align
/// id), made where the program reads its reported subjects; the program reads the holder's
/// length and letters by its ordinal id ([`ShownBioseq::ordinal`]).
#[derive(Clone, Debug)]
pub struct ShownBioseq {
    /// The holder's record viewed by its loaded line (the B-line) with the subject's Seq-align id
    /// ([`ParsedSubject::shown_by`]).
    pub record: Arc<ParsedSubject>,
    /// The holder, whose length (B-len) and letters the writers read.
    pub holder: DbSequence,
}

impl ShownBioseq {
    /// The Bioseq of the holder whose record is `holder_record` (its own record as the reports
    /// read it, of sequence `holder`), loaded by its line `line`, for the Seq-align id
    /// `alignment_id`. `Err` as [`ParsedSubject::shown_by`]'s.
    pub fn new(
        holder_record: &ParsedSubject,
        holder: DbSequence,
        line: usize,
        alignment_id: &ParsedSeqId,
    ) -> Result<Self, String> {
        let record = ParsedSubject::shown_by(holder_record, line, alignment_id)?;
        Ok(Self {
            record: Arc::new(record),
            holder,
        })
    }

    /// The ordinal id of the holder, by which the program reads its letters and length.
    pub(crate) fn ordinal(&self) -> u32 {
        self.holder.ordinal()
    }
}

/// What the record writers of a program read for one reported subject with parsed Seq-ids, the
/// one carrier of every program: the subject's own record, its loads' answer and its own letters
/// and length, or where its places give a holder ([`ShownBioseq`], made by
/// `algorithm::common::database::ShownHolders::shown`) the holder's record, loaded, with the
/// holder's letters and length (read by [`Written::ordinal`] from the program's letter table;
/// BLASTN's and BLASTP's rows then pass the holder's end through
/// [`crate::report::row_letters::shown_subject_row`]). Each program puts it in the subject's place
/// of its reported subjects and checks it as the subject's own
/// (`algorithm::common::database::check_reported_parsed`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Written<'s> {
    /// The record that the writers write.
    pub(crate) record: &'s Arc<ParsedSubject>,
    /// The ordinal id of the sequence whose letters and length the writers read.
    pub(crate) ordinal: u32,
    /// What the reports do at the subject's `Hit`: the holder's Bioseq is loaded.
    pub(crate) load: HitLoad,
    /// Whether the record is a holder's.
    pub(crate) held: bool,
}

impl<'s> Written<'s> {
    /// What the writers read for the subject `ordinal` with the record `subject`, whose holder
    /// they write where `shown` is given; `load`, the subject's own loads' answer.
    pub(crate) fn of(
        shown: Option<&'s ShownBioseq>,
        (ordinal, subject): (u32, &'s Arc<ParsedSubject>),
        load: impl FnOnce() -> HitLoad,
    ) -> Self {
        match shown {
            Some(bioseq) => Self {
                record: &bioseq.record,
                ordinal: bioseq.ordinal(),
                load: HitLoad::Loaded,
                held: true,
            },
            None => Self {
                record: subject,
                ordinal,
                load: load(),
                held: false,
            },
        }
    }

    /// Whether a query's scope finds the written sequence's Seq-align id, of `seen` (whether it
    /// finds the subject's): never for a holder, whose record `ShownHolders::shown` refuses where
    /// a query's scope finds the subject's or the holder's id.
    pub(crate) fn seen_by_scope(&self, seen: bool) -> bool {
        seen && !self.held
    }
}
