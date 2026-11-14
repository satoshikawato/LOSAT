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
    /// The number of subjects.
    pub(crate) fn len(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records.len(),
            Subjects::Stream(stream) => stream.db.len(),
            Subjects::Unreadable => 0,
        }
    }

    /// The numbers that the parameters of the search read (`SubjectNumbers`): a database's
    /// statistics (`SeqDb::statistics_numbers`: `STATS_TOTLEN` or the total, `STATS_NSEQ` or
    /// `GetNumSeqs`, the length of OID 0 and one sequence for a list that kept none),
    /// `GetMinLength` and `BlastSeqSrcGetAvgSeqLen`; a FASTA set's total, records, shortest
    /// record (raised to `BLAST_SEQSRC_MINLENGTH`) and average.
    pub(crate) fn numbers(&self) -> Result<super::parameters::SubjectNumbers> {
        let shortest = self.shortest_length();
        Ok(match self {
            Subjects::Stream(stream) => {
                let (total, count) = stream.db.statistics_numbers();
                super::parameters::SubjectNumbers {
                    total: i64::try_from(total)
                        .map_err(|_| anyhow!("BLASTX database length exceeds Int8"))?,
                    count: i64::from(count),
                    shortest,
                    average: stream.db.average_length(),
                }
            }
            Subjects::Fasta(_) | Subjects::Unreadable => {
                let total = self.total_length()?;
                let count = self.len() as i64;
                super::parameters::SubjectNumbers {
                    total,
                    count,
                    shortest,
                    average: (total / count.max(1)) as i32,
                }
            }
        })
    }

    /// Whether the subjects are `Unreadable`.
    pub(crate) fn is_unreadable(&self) -> bool {
        matches!(self, Subjects::Unreadable)
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:885-893,919,1236-1239
    // and c++/src/algo/blast/api/seqsrc_seqdb.cpp:375
    // ```c++
    //     buflen = CalculateSeqBufferLength(sv.size(), eBlastEncodingProtein);
    //     _ASSERT(buflen != 0);
    //     buf = buf_var = (Uint1*) malloc(sizeof(Uint1)*buflen);
    // ...
    //     *buf_var++ = GetSentinelByte(eBlastEncodingProtein);
    // ...
    //     *buf_var++ = GetSentinelByte(eBlastEncodingProtein);
    // ...
    //     case eBlastEncodingProtein:
    //         _ASSERT(sentinel == eSentinels);
    //         _ASSERT(strand == eNa_strand_unknown);
    //         retval = sequence_length + 2;
    // ...
    //         len = seqdb.GetSequence(oid, &buf);
    // ```
    /// The bytes that NCBI's buffer of subject `oid` holds after the sentinel byte that follows its
    /// residues: none for a record of `-subject` (`GetSequenceProtein` allocates the residues and a
    /// sentinel byte on each side, `CalculateSeqBufferLength`), and for a database sequence the
    /// rest of its volume's sequence file and of the file's last memory page, which NCBI maps
    /// whole (`XStreamDb::bytes_after`). The second item names the buffer in a refusal.
    pub(crate) fn buffer_after(&self, oid: usize) -> Result<(usize, &'static str)> {
        match self {
            Subjects::Fasta(_) => Ok((0, "subject")),
            Subjects::Stream(stream) => Ok((
                stream.db.bytes_after(oid)?,
                "database volume's last memory page",
            )),
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
            Subjects::Stream(stream) => i64::try_from(stream.db.total_length())
                .map_err(|_| anyhow!("BLASTX database length exceeds Int8")),
            Subjects::Unreadable => Ok(0),
        }
    }

    /// The shortest subject length that the hit-saving cutoff reads (`BlastSeqSrcGetMinSeqLen`):
    /// the shortest record raised to `BLAST_SEQSRC_MINLENGTH` (seqsrc_multiseq.cpp:233-239), or
    /// a database's `CSeqDB::GetMinLength` (`SeqDb::min_length`: the alias files'
    /// `MIN_SEQ_LENGTH`, 10 where it is not positive, without a list, the shortest sequence NCBI
    /// scans with one, DI-8b DECISIONS 21; BLASTX translates the query, so the value is not
    /// divided by 3, blast_setup.c:968-975).
    pub(crate) fn shortest_length(&self) -> usize {
        match self {
            Subjects::Fasta(records) => records
                .iter()
                .map(|s| s.sequence.len())
                .min()
                .unwrap_or(0)
                .max(BLAST_SEQSRC_MINLENGTH),
            Subjects::Stream(stream) => usize::try_from(stream.db.min_length()).unwrap_or(0),
            Subjects::Unreadable => BLAST_SEQSRC_MINLENGTH,
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
// ```c++
//     for (i = 0; i < sv.size(); i++) {
//         // Change unsupported residues to X
//         if (sv[i] == AMINOACID_TO_NCBISTDAA[(int)'O']) {
//             replaced_residues.push_back(i);
//             *buf_var++ = AMINOACID_TO_NCBISTDAA[(int)'X'];
// ```
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

    // NCBI reference: `search_residue` (blast_setup_cxx.cpp:894-899).
    /// A `-subject` record's pyrrolysine (O or o) is searched as X, through every encoding; the
    /// record, whose letters the reports show, keeps the O.
    #[test]
    fn a_fasta_subject_pyrrolysine_is_searched_as_x() {
        let records = parse_fasta(b">a\nSOSoXU\n", true, false).unwrap();
        let subjects = Subjects::Fasta(&records);
        let [s, x, u] = [b'S', b'X', b'U'].map(aa_char_to_ncbistdaa);
        assert_eq!(*subjects.encoded(0), [s, x, s, x, x, u]);
        assert_eq!(*subjects.encoded_with_sentinel(0), [s, x, s, x, x, u, 0]);
        assert_eq!(
            *subjects.x_encoded_with_sentinels(0),
            [0, s, x, s, x, x, u, 0]
        );
        assert_eq!(records[0].sequence, b"SOSOXU");
    }
}
