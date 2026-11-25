//! The subjects of a BLASTX search, by oid: the records of `-subject`, or the sequences of a
//! protein database (`-db`, read in blocks where they are searched, `x_db_stream.rs`,
//! oid = ordinal id). The search, the runtime and the composition redo read them through
//! this one type; the reports read records (`FastaRecord`) in both cases.

use std::borrow::Cow;

use anyhow::{anyhow, Result};

use super::input::FastaRecord;
use super::x_db_stream::{XBlock, XCursor, XPass, XPastHead, XStream};
use crate::utils::matrix::aa_char_to_ncbistdaa;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
/// The length that a database answers for its shortest sequence
/// (`BlastSeqSrcGetMinSeqLen`, seqdbimpl.cpp:128-136: "Do not bother scanning the db").
pub(crate) const BLAST_SEQSRC_MINLENGTH: usize = 10;

/// The subjects of a search (see the module documentation).
#[derive(Clone, Copy)]
pub(crate) enum Subjects<'a> {
    /// `-subject`: the records of the subject file.
    Fasta(&'a [FastaRecord]),
    /// `-db` read in blocks (`x_db_stream.rs`), with the block of the current pass or task.
    Stream(XStream<'a>),
    /// `-db` refused before its sequences are read (`database::BlastxDatabase::x_open_stream`: a
    /// database without sequences, or a longest sequence above `MAX_DBSEQ_LEN` with the check's
    /// error; the search's reads raise the other errors): a query batch without a valid context
    /// is reported as NCBI reports it (not
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

impl<'a> Subjects<'a> {
    /// The database of a `-db` search whose sequences are read (`None` for `-subject`).
    pub(crate) fn seqdb(&self) -> Option<&'a crate::seqdb::SeqDb> {
        match self {
            Subjects::Stream(stream) => Some(stream.db.seqdb()),
            Subjects::Fasta(_) | Subjects::Unreadable => None,
        }
    }

    /// The number of subjects.
    pub(crate) fn len(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.len(),
            Subjects::Stream(stream) => stream.db.len(),
            Subjects::Unreadable => 0,
        }
    }

    /// The numbers that the parameters of the search read (`SubjectNumbers`): a database's
    /// statistics (`SeqDb::statistics`: `STATS_TOTLEN` or the total, `STATS_NSEQ` or
    /// `GetNumSeqs`, the length of OID 0 and one sequence for a list that kept none, and the
    /// Gumbel block's length before that fallback),
    /// `GetMinLength` and `BlastSeqSrcGetAvgSeqLen`; a FASTA set's total, records, shortest
    /// record (raised to `BLAST_SEQSRC_MINLENGTH`) and average.
    pub(crate) fn numbers(&self) -> Result<super::parameters::SubjectNumbers> {
        let shortest = self.shortest_length();
        Ok(match self {
            Subjects::Stream(stream) => {
                let statistics = stream.db.statistics();
                super::parameters::SubjectNumbers {
                    total: statistics.setup_total,
                    count: statistics.setup_count,
                    gumbel_total: statistics.gumbel_db_length,
                    shortest,
                    average: statistics.average,
                    request_count: statistics.request_count,
                }
            }
            Subjects::Fasta(_) | Subjects::Unreadable => {
                let total = self.total_length()?;
                let count = self.len() as i64;
                super::parameters::SubjectNumbers {
                    total,
                    count: i32::try_from(count).unwrap_or(i32::MAX),
                    gumbel_total: total,
                    shortest,
                    average: crate::stats::setup::fasta_average_length(total, self.len()),
                    request_count: i32::try_from(count).unwrap_or(i32::MAX),
                }
            }
        })
    }

    /// Whether the subjects are `Unreadable`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407, 1434 (cited at
    // `seqdb::DbStatistics`)
    /// Whether NCBI's preliminary search computes the statistics of each subject on its own: a
    /// database whose `BlastSeqSrcGetTotLen` is 0 (`DbStatistics::per_subject`). A FASTA set's
    /// total is the sum of its records (NCBI's "dbscan" mode).
    pub(crate) fn per_subject(&self) -> bool {
        matches!(self, Subjects::Stream(stream) if stream.db.statistics().per_subject())
    }

    pub(crate) fn is_unreadable(&self) -> bool {
        matches!(self, Subjects::Unreadable)
    }

    /// The error of the preliminary search's first read of the database after a batch's set-up
    /// checks, where a sequence file cannot be opened
    /// (`common::database::sequence_file_read_error`): NCBI's engine error, or LOSAT's refusal of
    /// an error that is not NCBI's; `None` for `-subject`.
    pub(crate) fn sequence_file_read_error(&self) -> Option<anyhow::Error> {
        use crate::algorithm::common::database::{sequence_file_read_error, PreliminaryReadError};
        let Subjects::Stream(stream) = self else {
            return None;
        };
        let error = sequence_file_read_error(stream.db.seqdb(), "BLASTX")?;
        Some(match error.downcast::<PreliminaryReadError>() {
            Ok(PreliminaryReadError(text)) => super::native::engine(text),
            Err(refusal) => refusal,
        })
    }

    /// `error`, raised by a search after the point where `Unreadable` subjects fail it
    /// (`search::search_core`): an `XPastHead` where the subjects are a stream whose reads have
    /// not checked every sequence, else `error` itself.
    pub(crate) fn x_past_head(&self, error: anyhow::Error) -> anyhow::Error {
        match self {
            Subjects::Stream(stream) if !stream.db.x_checked() => XPastHead(error).into(),
            _ => error,
        }
    }

    /// Whether the subjects are a stream whose reads have not checked every sequence and whose
    /// database fails the whole check (`XStreamDb::x_damaged`): such a database was refused before
    /// the search (`Unreadable`) when it was read before the first batch.
    pub(crate) fn x_damaged(&self) -> bool {
        match self {
            Subjects::Stream(stream) => !stream.db.x_checked() && stream.db.x_damaged(),
            _ => false,
        }
    }

    /// Whether there is no subject.
    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The length of subject `oid`.
    pub(crate) fn length(&self, oid: usize) -> usize {
        match self {
            Subjects::Fasta(records) => records[oid].sequence.len(),
            Subjects::Stream(stream) => stream.length(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` followed by the byte 0 (what the
    /// preliminary stage scans): a record's letters encoded (`search_residue`: a pyrrolysine
    /// (O) as X), a database sequence's stored bytes as they are (NCBI searches, redoes and
    /// shows the stored bytes of a database, decision 108).
    pub(crate) fn encoded_with_sentinel(&self, oid: usize) -> Cow<'a, [u8]> {
        match *self {
            Subjects::Fasta(records) => {
                let mut encoded = encode(&records[oid]);
                encoded.push(0);
                Cow::Owned(encoded)
            }
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
                encoded.extend(record.sequence.iter().copied().map(search_residue));
                encoded.push(0);
                Cow::Owned(encoded)
            }
            Subjects::Stream(stream) => stream.x_encoded_with_sentinels(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` (what the linking, the redo and the
    /// traceback read), as `encoded_with_sentinel` without the byte 0.
    pub(crate) fn encoded(&self, oid: usize) -> Cow<'a, [u8]> {
        match *self {
            Subjects::Fasta(records) => Cow::Owned(encode(&records[oid])),
            Subjects::Stream(stream) => stream.encoded(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    // No NCBI counterpart: the schedules (`x_subject_flow.rs`); it changes no value NCBI computes.
    /// The length of the searched sequences: the sum of the records' lengths, or the length of
    /// the database's included sequences (`SeqDb::included_length`); `None` when it cannot be
    /// read.
    #[cfg(feature = "parallel")]
    pub(crate) fn included_length(&self) -> Option<u64> {
        match self {
            Subjects::Fasta(records) => Some(
                records
                    .iter()
                    .map(|s| s.sequence.len() as u64)
                    .fold(0u64, u64::saturating_add),
            ),
            Subjects::Stream(stream) => stream.db.included_length().ok(),
            Subjects::Unreadable => None,
        }
    }

    /// The total length (Int8, `BlastSeqSrcGetTotLen`): the sum of the records' lengths, or the
    /// database's `GetTotalLength` (`SeqDb::total_length`: the walk of `LENGTH`, not
    /// `STATS_TOTLEN`), which the linking of an ungapped search reads (blast_engine.c:1407).
    pub(crate) fn total_length(&self) -> Result<i64> {
        match self {
            Subjects::Fasta(records) => records
                .iter()
                .try_fold(0i64, |sum, s| {
                    sum.checked_add(i64::try_from(s.sequence.len()).ok()?)
                })
                .ok_or_else(|| anyhow!("BLASTX subject-set length exceeds Int8")),
            // NCBI's `Int8` of the `Uint8` (c++/src/algo/blast/api/seqsrc_seqdb.cpp:183-189): negative for an alias
            // file's `LENGTH` of 2^63 or more, which takes the subject's own branch of the link
            // cutoffs (`linking::link_cutoffs`; line RG D1a-3).
            Subjects::Stream(stream) => Ok(stream.db.total_length() as i64),
            Subjects::Unreadable => Ok(0),
        }
    }

    /// The shortest subject length that the hit-saving cutoff reads (`BlastSeqSrcGetMinSeqLen`):
    /// the shortest record raised to `BLAST_SEQSRC_MINLENGTH` (seqsrc_multiseq.cpp:233-239), or
    /// a database's `CSeqDB::GetMinLength` (`SeqDb::min_length`: the alias files'
    /// `MIN_SEQ_LENGTH`, 10 where it is not positive, without a list, the shortest sequence NCBI
    /// scans with one, DI-8b DECISIONS 21; BLASTX translates the query, so the value is not
    /// divided by 3, blast_setup.c:968-975).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-975
    pub(crate) fn shortest_length(&self) -> i32 {
        match self {
            Subjects::Fasta(records) => i32::try_from(
                records
                    .iter()
                    .map(|s| s.sequence.len())
                    .min()
                    .unwrap_or(0)
                    .max(BLAST_SEQSRC_MINLENGTH),
            )
            .unwrap_or(i32::MAX),
            Subjects::Stream(stream) => stream.db.min_length(),
            Subjects::Unreadable => BLAST_SEQSRC_MINLENGTH as i32,
        }
    }

    /// The length of the longest subject.
    pub(crate) fn longest_length(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.iter().map(|s| s.sequence.len()).max().unwrap_or(0),
            Subjects::Stream(stream) => stream.db.longest(),
            Subjects::Unreadable => 0,
        }
    }

    /// EXPERIMENT (LOSAT_X_BXTASKRES): the lengths of the subjects from `oid` (below `len`) on, at
    /// least one and at most `X_LENGTHS` of them (with `-db` read in blocks, to the end of `oid`'s
    /// volume, from the index: `XStreamDb::x_lengths_from`); an error of the index read is the
    /// reader's.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_lengths_from(&self, oid: usize) -> Result<Vec<u32>> {
        const X_LENGTHS: usize = 16_384;
        let end = self.len().min(oid.saturating_add(X_LENGTHS));
        match self {
            Subjects::Fasta(records) => Ok(records[oid..end]
                .iter()
                .map(|record| u32::try_from(record.sequence.len()).unwrap_or(u32::MAX))
                .collect()),
            Subjects::Stream(stream) => stream.db.x_lengths_from(oid),
            Subjects::Unreadable => unreachable!("no subject is read from an unreadable database"),
        }
    }

    /// The NCBIstdaa residues of subject `oid` for the stages after the preliminary search (the
    /// composition redo, the traceback, the reports' records): `encoded`, or for a `-db` a read
    /// of the database whose error is the reader's error.
    pub(crate) fn read_encoded(&self, oid: usize) -> Result<Cow<'a, [u8]>> {
        match self {
            Subjects::Stream(stream) => stream.read_encoded(oid),
            _ => Ok(self.encoded(oid)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:328-353 and
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1631-1636
    // With `-db_hard_mask`, each fetch after the preliminary stage (the composition-based
    // statistics, the traceback) is a copy with `X` over the sequence's masks; the formatter
    // reads the stored letters.
    /// The residues of subject `oid` for the stages after the preliminary search that score it
    /// (the composition redo, the traceback): `read_encoded`, with `X` over the masks of
    /// `-db_hard_mask` (`options.db_mask`). The reports' records read `read_encoded`.
    pub(crate) fn read_searched(
        &self,
        oid: usize,
        options: &super::args::ResolvedOptions,
    ) -> Result<Cow<'a, [u8]>> {
        let encoded = self.read_encoded(oid)?;
        match options.db_mask.as_ref() {
            Some(mask) if mask.kind == crate::blastinput::db_mask::MaskKind::Hard => {
                self.x_hard_masked(oid, mask.id, encoded)
            }
            _ => Ok(encoded),
        }
    }

    /// `encoded` (the residues of database sequence `oid`) with `X` over its stored masks of the
    /// algorithm `algo_id` (`db_mask_search::mask_protein_letters`).
    #[cold]
    #[inline(never)]
    fn x_hard_masked(
        &self,
        oid: usize,
        algo_id: i32,
        encoded: Cow<'a, [u8]>,
    ) -> Result<Cow<'a, [u8]>> {
        let Some(db) = self.seqdb() else {
            return Ok(encoded);
        };
        let mut reader = crate::algorithm::common::db_mask_search::SubjectMasks::with_window(1);
        let masks = reader.of(
            db,
            oid as u32,
            algo_id,
            crate::algorithm::common::database::database_error_of,
            "BLASTX",
        )?;
        if masks.is_empty() {
            return Ok(encoded);
        }
        let mut residues = encoded.into_owned();
        crate::algorithm::common::db_mask_search::mask_protein_letters(&mut residues, masks);
        Ok(Cow::Owned(residues))
    }

    /// The reader of one pass over the subjects in ascending
    /// order (`XSubjectPass::view`).
    pub(crate) fn x_pass(&self) -> XSubjectPass<'a> {
        match *self {
            Subjects::Stream(stream) => XSubjectPass::Stream(XPass::new(stream.db)),
            subjects => XSubjectPass::Whole(subjects),
        }
    }

    /// The subjects `first..end` (non-empty) read at once, for
    /// `x_with_block`; `None` unless the subjects are read in blocks.
    pub(crate) fn x_read_range(&self, first: usize, end: usize) -> Result<Option<XBlock>> {
        match self {
            Subjects::Stream(stream) => stream.db.read_range(first, end).map(Some),
            _ => Ok(None),
        }
    }

    /// The subjects with the block of `x_read_range`.
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

    /// A pass of the `-ungapped` sink over the residues of the
    /// subjects with HSPs, when the subjects are read in blocks.
    pub(crate) fn x_residue_pass(&self) -> Option<XPass<'a>> {
        match self {
            Subjects::Stream(stream) => Some(XPass::new(stream.db)),
            _ => None,
        }
    }

    /// The lengths that the sink reads of the subjects in its
    /// order, when they are read in blocks.
    pub(crate) fn x_cursor(&self) -> Option<XCursor<'a>> {
        match self {
            Subjects::Stream(stream) => Some(XCursor::new(stream.db)),
            _ => None,
        }
    }
}

/// The subjects of one pass over them in ascending order: the
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
        .map(search_residue)
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:894-899
/// The NCBIstdaa residue that the search reads for the letter `letter` of a `-subject` record:
/// a pyrrolysine (O) is read as X, so it is identical to an X of the query in the Seq-aligns'
/// `num_ident` (`statistics::identity_and_positive`). The reports show the record's letters.
#[inline(always)]
fn search_residue(letter: u8) -> u8 {
    // `aa_char_to_ncbistdaa` of each ASCII byte with O as X: one load per letter, as before.
    const TABLE: [u8; 128] = {
        let mut table = [0; 128];
        let mut byte = 0;
        while byte < 128 {
            table[byte] = aa_char_to_ncbistdaa(byte as u8);
            byte += 1;
        }
        table[b'O' as usize] = aa_char_to_ncbistdaa(b'X');
        table
    };
    if letter < 128 {
        TABLE[letter as usize]
    } else {
        aa_char_to_ncbistdaa(letter)
    }
}

#[cfg(test)]
#[path = "subjects_tests.rs"]
mod tests;
