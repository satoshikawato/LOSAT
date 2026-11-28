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

use crate::report::db_subject::ParsedSubject;
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

    /// This Bioseq whose description line (outfmt 0 to 4) writes the title of the holder's line
    /// `line` ([`ParsedSubject::described_by`]). `Err` as that function's.
    pub fn described_by(self, line: usize) -> Result<Self, String> {
        let record = Arc::unwrap_or_clone(self.record).described_by(line)?;
        Ok(Self {
            record: Arc::new(record),
            holder: self.holder,
        })
    }

    /// The ordinal id of the holder, by which the program reads its letters and length.
    pub(crate) fn ordinal(&self) -> u32 {
        self.holder.ordinal()
    }
}
