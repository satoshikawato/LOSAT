//! The subjects of a BLASTX search, by oid: the records of `-subject`, or the stored
//! sequences of a protein database (`-db`, `algorithm::common::database::StoredProteins`,
//! oid = ordinal id). The search, the runtime and the composition redo read them through
//! this one type; the reports read records (`FastaRecord`) in both cases.

use std::borrow::Cow;

use anyhow::{anyhow, Result};

use super::input::FastaRecord;
use super::x_db_stream::{XBlock, XCursor, XPass, XStream};
use crate::algorithm::common::database::StoredProteins;
use crate::utils::matrix::aa_char_to_ncbistdaa;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
// ```c
// #define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */
// ```
/// The length that a database answers for its shortest sequence
/// (`BlastSeqSrcGetMinSeqLen`, seqdbimpl.cpp:128-136: "Do not bother scanning the db").
pub(crate) const BLAST_SEQSRC_MINLENGTH: usize = 10;

/// The subjects of a search (see the module documentation).
#[derive(Clone, Copy)]
pub(crate) enum Subjects<'a> {
    /// `-subject`: the records of the subject file.
    Fasta(&'a [FastaRecord]),
    /// `-db`: the stored sequences of the database.
    Database(&'a StoredProteins),
    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): `-db` read in blocks (`x_db_stream.rs`), with the block of
    /// the current pass or task.
    Stream(XStream<'a>),
    /// `-db` whose sequences could not be read (`database::BlastxDatabase::read_stores`):
    /// a query batch without a valid context is reported as NCBI reports it (not
    /// searched); the first batch with one fails with `UnreadableSubjects`, which the
    /// caller replaces by the reader's error (decision 112).
    Unreadable,
}

/// The error of a query batch with a valid context whose subjects are `Unreadable`.
#[derive(Debug)]
pub(crate) struct UnreadableSubjects;

impl std::fmt::Display for UnreadableSubjects {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the database's sequences could not be read")
    }
}

impl std::error::Error for UnreadableSubjects {}

impl<'a> From<&'a [FastaRecord]> for Subjects<'a> {
    fn from(records: &'a [FastaRecord]) -> Self {
        Subjects::Fasta(records)
    }
}

impl<'a> From<&'a Vec<FastaRecord>> for Subjects<'a> {
    fn from(records: &'a Vec<FastaRecord>) -> Self {
        Subjects::Fasta(records)
    }
}

impl<'a> From<&'a StoredProteins> for Subjects<'a> {
    fn from(stores: &'a StoredProteins) -> Self {
        Subjects::Database(stores)
    }
}

impl<'a> Subjects<'a> {
    /// The number of subjects.
    pub(crate) fn len(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.len(),
            Subjects::Database(stores) => stores.len(),
            Subjects::Stream(stream) => stream.db.len(),
            Subjects::Unreadable => 0,
        }
    }

    /// Whether the subjects are `Unreadable`.
    pub(crate) fn is_unreadable(&self) -> bool {
        matches!(self, Subjects::Unreadable)
    }

    /// Whether there is no subject.
    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The length of subject `oid`.
    pub(crate) fn length(&self, oid: usize) -> usize {
        match self {
            Subjects::Fasta(records) => records[oid].sequence.len(),
            Subjects::Database(stores) => stores.length(oid),
            Subjects::Stream(stream) => stream.length(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` followed by the byte 0 (what the
    /// preliminary stage scans): a record's letters encoded (`aa_char_to_ncbistdaa`), a
    /// database sequence's stored bytes as they are (NCBI searches, redoes and shows the
    /// stored bytes of a database, decision 108).
    pub(crate) fn encoded_with_sentinel(&self, oid: usize) -> Cow<'a, [u8]> {
        match *self {
            Subjects::Fasta(records) => {
                let mut encoded = encode(&records[oid]);
                encoded.push(0);
                Cow::Owned(encoded)
            }
            Subjects::Database(stores) => Cow::Borrowed(stores.residues_with_sentinel(oid)),
            Subjects::Stream(stream) => stream.encoded_with_sentinel(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `encoded_with_sentinel` with the byte 0 before the
    /// residues too (a database's bytes as they are stored: the sentinel of the sequence before).
    pub(crate) fn x_encoded_with_sentinels(&self, oid: usize) -> Cow<'a, [u8]> {
        match *self {
            Subjects::Fasta(records) => {
                let record = &records[oid];
                let mut encoded = Vec::with_capacity(record.sequence.len() + 2);
                encoded.push(0);
                encoded.extend(record.sequence.iter().copied().map(aa_char_to_ncbistdaa));
                encoded.push(0);
                Cow::Owned(encoded)
            }
            Subjects::Database(stores) => Cow::Borrowed(stores.x_residues_with_sentinels(oid)),
            Subjects::Stream(stream) => stream.x_encoded_with_sentinels(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` (what the linking, the redo and the
    /// traceback read), as `encoded_with_sentinel` without the byte 0.
    pub(crate) fn encoded(&self, oid: usize) -> Cow<'a, [u8]> {
        match *self {
            Subjects::Fasta(records) => Cow::Owned(encode(&records[oid])),
            Subjects::Database(stores) => {
                let with_sentinel = stores.residues_with_sentinel(oid);
                Cow::Borrowed(&with_sentinel[..with_sentinel.len() - 1])
            }
            Subjects::Stream(stream) => stream.encoded(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The total length (Int8): the sum of the records' lengths, or the database's
    /// total of its index files (`BlastSeqSrcGetTotLen`).
    pub(crate) fn total_length(&self) -> Result<i64> {
        match self {
            Subjects::Fasta(records) => records
                .iter()
                .try_fold(0i64, |sum, s| {
                    sum.checked_add(i64::try_from(s.sequence.len()).ok()?)
                })
                .ok_or_else(|| anyhow!("BLASTX subject-set length exceeds Int8")),
            Subjects::Database(stores) => i64::try_from(stores.total_length())
                .map_err(|_| anyhow!("BLASTX database length exceeds Int8")),
            Subjects::Stream(stream) => i64::try_from(stream.db.total_length())
                .map_err(|_| anyhow!("BLASTX database length exceeds Int8")),
            Subjects::Unreadable => Ok(0),
        }
    }

    /// The shortest subject length that the hit-saving cutoff reads before it is raised
    /// to 10 (`BlastSeqSrcGetMinSeqLen`): the shortest record, or 10 for a database
    /// (D1; BLASTX translates the query, so the value is not divided by 3,
    /// blast_setup.c:968-975).
    pub(crate) fn shortest_length(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.iter().map(|s| s.sequence.len()).min().unwrap_or(0),
            Subjects::Database(_) | Subjects::Stream(_) | Subjects::Unreadable => {
                BLAST_SEQSRC_MINLENGTH
            }
        }
    }

    /// The length of the longest subject.
    pub(crate) fn longest_length(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.iter().map(|s| s.sequence.len()).max().unwrap_or(0),
            Subjects::Database(stores) => stores.max_length(),
            Subjects::Stream(stream) => stream.db.longest(),
            Subjects::Unreadable => 0,
        }
    }

    /// EXPERIMENT (LOSAT_X_BXTASKRES): the lengths of the subjects from `oid` on, at least one and
    /// at most `X_LENGTHS` of them (with `-db` read in blocks, to the end of `oid`'s volume, from
    /// the index: `XStreamDb::x_lengths_from`); an error of the index read is the reader's.
    pub(crate) fn x_lengths_from(&self, oid: usize) -> Result<Vec<u32>> {
        const X_LENGTHS: usize = 16_384;
        let end = self.len().min(oid.saturating_add(X_LENGTHS));
        match self {
            Subjects::Fasta(records) => Ok(records[oid..end]
                .iter()
                .map(|record| u32::try_from(record.sequence.len()).unwrap_or(u32::MAX))
                .collect()),
            Subjects::Database(stores) => Ok((oid..end)
                .map(|i| u32::try_from(stores.length(i)).unwrap_or(u32::MAX))
                .collect()),
            Subjects::Stream(stream) => stream.db.x_lengths_from(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` for the stages after the preliminary search (the
    /// composition redo, the traceback, the reports' records): `encoded`, or with
    /// `LOSAT_X_BXDBSTREAM` a read of the database whose error is the reader's error.
    pub(crate) fn read_encoded(&self, oid: usize) -> Result<Cow<'a, [u8]>> {
        match self {
            Subjects::Stream(stream) => stream.read_encoded(oid),
            _ => Ok(self.encoded(oid)),
        }
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the reader of one pass over the subjects in ascending
    /// order (`XSubjectPass::view`).
    pub(crate) fn x_pass(&self) -> XSubjectPass<'a> {
        match *self {
            Subjects::Stream(stream) => XSubjectPass::Stream(XPass::new(stream.db)),
            subjects => XSubjectPass::Whole(subjects),
        }
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the subjects `first..end` (non-empty) read at once, for
    /// `x_with_block`; `None` unless the subjects are read in blocks.
    pub(crate) fn x_read_range(&self, first: usize, end: usize) -> Result<Option<XBlock>> {
        match self {
            Subjects::Stream(stream) => stream.db.read_range(first, end).map(Some),
            _ => Ok(None),
        }
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the subjects with the block of `x_read_range`.
    pub(crate) fn x_with_block<'b>(&self, block: Option<&'b XBlock>) -> Subjects<'b>
    where
        'a: 'b,
    {
        match (*self, block) {
            (Subjects::Stream(stream), Some(block)) => Subjects::Stream(XStream {
                db: stream.db,
                block: Some(block),
            }),
            (subjects, _) => subjects,
        }
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): a pass of the `-ungapped` sink over the residues of the
    /// subjects with HSPs, when the subjects are read in blocks.
    pub(crate) fn x_residue_pass(&self) -> Option<XPass<'a>> {
        match self {
            Subjects::Stream(stream) => Some(XPass::new(stream.db)),
            _ => None,
        }
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the lengths that the sink reads of the subjects in its
    /// order, when they are read in blocks.
    pub(crate) fn x_cursor(&self) -> Option<XCursor<'a>> {
        match self {
            Subjects::Stream(stream) => Some(XCursor::new(stream.db)),
            _ => None,
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): the subjects of one pass over them in ascending order: the
/// subjects themselves, or the blocks of a database read in blocks.
pub(crate) enum XSubjectPass<'a> {
    Whole(Subjects<'a>),
    Stream(XPass<'a>),
}

impl XSubjectPass<'_> {
    /// The subjects to read subject `oid` from (`oid` at least that of the call before).
    #[inline]
    pub(crate) fn view(&mut self, oid: usize) -> Result<Subjects<'_>> {
        match self {
            XSubjectPass::Whole(subjects) => Ok(*subjects),
            XSubjectPass::Stream(pass) => pass.view(oid).map(Subjects::Stream),
        }
    }
}

/// The NCBIstdaa residues of the record `record`.
fn encode(record: &FastaRecord) -> Vec<u8> {
    record
        .sequence
        .iter()
        .copied()
        .map(aa_char_to_ncbistdaa)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastx::input::parse_fasta;

    // No NCBI counterpart: LOSAT's encoding of `-subject` records (`encoded_with_sentinel` is the
    // reference).
    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): a FASTA subject read with both sentinels is the record's
    /// encoding with the byte 0 added before it, for records of 0, 1 and more letters, lower case
    /// and letters outside the alphabet included.
    #[test]
    fn fasta_subjects_with_both_sentinels_add_the_byte_before() {
        let records =
            parse_fasta(b">a\nMKV*LLAx\n>b\nW\n>c\n\n>d\nbjouzBJOUZ-\n", true, false).unwrap();
        let subjects = Subjects::Fasta(&records);
        for oid in 0..records.len() {
            let both = subjects.x_encoded_with_sentinels(oid);
            assert_eq!(
                *both,
                [&[0u8][..], &subjects.encoded_with_sentinel(oid)].concat(),
                "{oid}"
            );
        }
    }
}
