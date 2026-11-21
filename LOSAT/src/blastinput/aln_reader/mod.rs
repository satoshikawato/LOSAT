//! NCBI's alignment reader (`CAlnReader` of objtools/readers) as psiblast's `-in_msa` uses it:
//! `CAlnReader(in).SetClustal(eAlpha_Protein).Read(false, true)`, then `GetSeqs()` (the rows) and
//! `GetSeqEntry()` (the Bioseq of each row).
//!
//! - `lines`: the peek-ahead line input (`CPeekAheadStream`).
//! - `format_guess`: the guess of the format (`CAlnFormatGuesser`), all six tests.
//! - `clustal`: the Clustal scanner and the base scanner's verification. A scanner of another
//!   format is one more file that returns the ids and the rows (item PS6b of line PS); until then
//!   an input that NCBI reads as one of them is refused by name ([`AlnError::Unsupported`]).
//!
//! NCBI's reader reports every problem (error or warning) through an error reporter without a
//! listener, which throws `CObjReaderParseException`; `CPsiBlastInputClustalW::x_ReadAsciiMsa`
//! swallows that and `GetSeqEntry()` then throws "Seq_entry is not available until after Read()"
//! ([`AlnError::NotAvailable`]). So the read needs the exact decision of success or failure, not
//! the messages of the scanners.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_error_reporter.cpp:130-136
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:122-138
//! (No report of the readers has that text, so the branch never runs.)

mod clustal;
mod format_guess;
mod lines;

#[cfg(test)]
mod tests;

pub(crate) use format_guess::{guess_format, AlignFormat};
pub(crate) use lines::PeekAheadStream;

use crate::seqid::accguide::{identify_seq_id, parse_flags, F_NUC, F_PROT};
use crate::seqid::model::ParsedSeqId;
use crate::seqid::parse::parse_ids;

/// `CSeq_inst::EMol` of a Bioseq that the alignment reader makes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mol {
    /// `eMol_aa`: the letters are `iupacaa`.
    Aa,
    /// `eMol_na`, `eMol_dna` and `eMol_rna`: the letters are `iupacna`, packed by
    /// `CSeqportUtil::Pack`, which NCBI's PSSM code cannot convert to NCBIstdaa.
    Nucleotide,
}

/// One row of the alignment as `GetSeqEntry` makes its Bioseq.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AlnBioseq {
    /// `CSeq_id::ParseIDs` of the row's id string.
    pub(crate) ids: Vec<ParsedSeqId>,
    pub(crate) mol: Mol,
    /// The letters of the row without gaps, upper case (`m_SeqVec`).
    pub(crate) letters: Vec<u8>,
}

/// What `CAlnReader` gives `CPsiBlastInputClustalW`: `GetSeqs()` (the rows as in the file, with
/// their `-` gaps and their case) and the Bioseqs of `GetSeqEntry()`, in the order of the rows.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Alignment {
    pub(crate) rows: Vec<Vec<u8>>,
    pub(crate) bioseqs: Vec<AlnBioseq>,
}

/// Why `read_alignment` gives no alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum AlnError {
    /// NCBI's reader threw (a scanner or validation report) and `x_ReadAsciiMsa` swallowed it:
    /// `GetSeqEntry()` throws `CAlnReader::GetSeqEntry(): Seq_entry is not available until after
    /// Read()` (`CObjReaderParseException`, "BLAST query error").
    NotAvailable,
    /// `x_GetSequenceMolType`'s `sReportError` of a row with both `T` and `U`, thrown from
    /// `GetSeqEntry` (`CObjReaderParseException`, "BLAST query error").
    InvalidMolType,
    /// An input that NCBI reads as `format` (`AlignFormat`'s name), a format that LOSAT's reader
    /// does not read yet.
    Unsupported(&'static str),
    /// A state of the reader that no input gives: one of its two preconditions on the Seq-id
    /// parser does not hold ([`NO_SEQ_ID_WARNING`], [`AN_ID_OF_EVERY_ROW`]; named as LOSAT's
    /// internal error names it). Not a refusal: no alignment reaches it.
    Internal(&'static str),
}

/// The flags of `x_ParseAndValidateSeqIds` (`fParse_AnyLocal | fParse_RawText`, cited at
/// [`row_ids`]).
const ID_PARSE_FLAGS: u32 = parse_flags::ANY_LOCAL | parse_flags::RAW_TEXT;

/// Precondition 1 of [`row_ids`]: NCBI's Seq-id parser posts a warning only where
/// `fParse_PartialOK` is set (`CSeq_id::ParseIDs` and `CSeq_id::Set`, `seqid::parse`), and
/// [`ID_PARSE_FLAGS`] do not set it, so the parser of a row's id posts none (a post would be a
/// line of NCBI's standard error that this reader does not write).
pub(crate) const NO_SEQ_ID_WARNING: &str =
    "the Seq-id parser posted a warning without the flag that allows one";

/// Precondition 2 of [`row_ids`]: `CSeq_id::ParseIDs` gives no id and no exception only for a
/// string of white space, and a row's id is the first token of a line whose ends carry no white
/// space (`clustal::scan`), so every row has an id (`GetSeqEntry` reads `ids.front()`).
pub(crate) const AN_ID_OF_EVERY_ROW: &str = "the Seq-id parser read a row's id as no id";

/// The text of `x_GetSequenceMolType`'s error (the message that `CObjReaderParseException`
/// carries).
pub(crate) const INVALID_MOL_TYPE: &str = "Invalid Mol Type: U and T cannot appear in the same nucleotide sequence. Reinterpreting as protein.";

/// The text of the `CObjReaderParseException` that `GetSeqEntry()` throws for a read that did
/// not succeed.
pub(crate) const SEQ_ENTRY_NOT_AVAILABLE: &str =
    "CAlnReader::GetSeqEntry(): Seq_entry is not available until after Read()";

/// NCBI's name of the format in LOSAT's refusal.
fn format_name(format: AlignFormat) -> &'static str {
    match format {
        AlignFormat::Nexus => "NEXUS",
        AlignFormat::FastaGap => "FASTA",
        AlignFormat::Phylip => "PHYLIP",
        AlignFormat::Sequin => "Sequin",
        AlignFormat::MultAlign => "MultAlin",
        AlignFormat::Clustal => "Clustal",
        AlignFormat::Unknown => "unknown",
    }
}

/// `CAlnReader::Read(false, true)` and `GetSeqs()` and `GetSeqEntry()` of `bytes`, the whole
/// input.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:319-339
/// A format other than Clustal and the unknown one is refused by name (the scanners of PS6b);
/// the unknown format's scanner throws "Input file format not recognized.".
pub(crate) fn read_alignment(bytes: &[u8]) -> Result<Alignment, AlnError> {
    let mut input = PeekAheadStream::new(bytes);
    let format = guess_format(&mut input);
    match format {
        AlignFormat::Clustal => {}
        AlignFormat::Unknown => return Err(AlnError::NotAvailable),
        other => return Err(AlnError::Unsupported(format_name(other))),
    }
    let scanned = clustal::scan(&mut input).map_err(|_| AlnError::NotAvailable)?;
    let ids = verify_alignment_info(&scanned)?;
    entries(scanned.rows, ids)
}

/// `CAlnReader::x_VerifyAlignmentInfo`: at least two sequences, and the ids of the rows parsed
/// and checked ([`row_ids`]).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:377-403
fn verify_alignment_info(scanned: &clustal::Scanned) -> Result<Vec<Vec<ParsedSeqId>>, AlnError> {
    if scanned.rows.len() < 2 {
        return Err(AlnError::NotAvailable);
    }
    scanned.ids.iter().map(|id| row_ids(id)).collect()
}

/// `CAlnReader::x_ParseAndValidateSeqIds` of the id string of one row: its Seq-ids, at least one.
/// A parse exception and a report of the id validator are read failures (the error reporter
/// throws: no listener).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:349-368
/// (`flags` is 0 here, so `flags^fGenerateLocalIDs` is not 0 and the raw text flag is set.)
///
/// Two states cannot come out of the parser for a row's id, and each is LOSAT's internal error
/// ([`AlnError::Internal`]), not a refusal, because no alignment reaches it (review S3 of line PS
/// item PS6: 185 probed ids of every Seq-id type, and the parser's code): a warning
/// ([`NO_SEQ_ID_WARNING`]) and an empty list of ids ([`AN_ID_OF_EVERY_ROW`]). The tests hold both
/// preconditions (`the_seq_id_parser_posts_no_warning_under_the_readers_flags`,
/// `every_scanned_id_parses_to_an_id_or_fails`).
pub(crate) fn row_ids(id: &[u8]) -> Result<Vec<ParsedSeqId>, AlnError> {
    let mut warned = false;
    let parsed = parse_ids(id, ID_PARSE_FLAGS, &mut |_| warned = true)
        .map_err(|_| AlnError::NotAvailable)?;
    if warned {
        return Err(AlnError::Internal(NO_SEQ_ID_WARNING));
    }
    let ids: Vec<ParsedSeqId> = parsed.into_iter().flatten().collect();
    if ids.is_empty() {
        return Err(AlnError::Internal(AN_ID_OF_EVERY_ROW));
    }
    if crate::writedb::fasta_source::defline_ids::alignment_ids_are_reported(&ids) {
        return Err(AlnError::NotAvailable);
    }
    Ok(ids)
}

/// `CAlnReader::GetSeqAlign`'s `m_SeqVec` and `GetSeqEntry`: the upper-case letters of each row
/// without gaps, its Bioseq with the molecule of its first id (`IdentifyAccession`) or, for an id
/// that names none, of its letters (`x_GetSequenceMolType`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:749-769
/// `x_GetSequenceMolType` (the alphabet has 53 letters, so an undefined guess is protein):
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:670-698
/// (`sReportError` throws `CObjReaderParseException`: no listener. The letters are upper case with
/// no gap: `m_SeqVec`; the missing characters of the sequence info are none.) An id that names a
/// molecule decides before the letters are looked at.
fn entries(rows: Vec<Vec<u8>>, ids: Vec<Vec<ParsedSeqId>>) -> Result<Alignment, AlnError> {
    let mut bioseqs = Vec::with_capacity(rows.len());
    for (row, ids) in rows.iter().zip(ids) {
        let letters: Vec<u8> = row
            .iter()
            .filter(|&&c| c != b'-')
            .map(u8::to_ascii_uppercase)
            .collect();
        let info = identify_seq_id(
            &ids[0],
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |_| {},
        );
        let mol = if info & F_NUC != 0 {
            Mol::Nucleotide
        } else if info & F_PROT != 0 {
            Mol::Aa
        } else {
            match crate::writedb::fasta_source::reader::sequence_type(&letters) {
                Some(true) | None => Mol::Aa,
                Some(false) => {
                    if letters.contains(&b'T') && letters.contains(&b'U') {
                        return Err(AlnError::InvalidMolType);
                    }
                    Mol::Nucleotide
                }
            }
        };
        bioseqs.push(AlnBioseq { ids, mol, letters });
    }
    Ok(Alignment { rows, bioseqs })
}
