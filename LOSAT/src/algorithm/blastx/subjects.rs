//! The subjects of a BLASTX search, by oid: the records of `-subject`, or the stored
//! sequences of a protein database (`-db`, `algorithm::common::database::StoredProteins`,
//! oid = ordinal id). The search, the runtime and the composition redo read them through
//! this one type; the reports read records (`FastaRecord`) in both cases.

use std::borrow::Cow;

use anyhow::{anyhow, Result};

use super::input::FastaRecord;
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
            Subjects::Database(_) | Subjects::Unreadable => BLAST_SEQSRC_MINLENGTH,
        }
    }

    /// The length of the longest subject.
    pub(crate) fn longest_length(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.iter().map(|s| s.sequence.len()).max().unwrap_or(0),
            Subjects::Database(stores) => stores.max_length(),
            Subjects::Unreadable => 0,
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
