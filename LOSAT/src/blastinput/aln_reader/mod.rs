//! NCBI's alignment reader (`CAlnReader` of objtools/readers) as psiblast's `-in_msa` uses it:
//! a `CAlnReader` on the input, `SetClustal` with `eAlpha_Protein`, and `Read` with `false` and
//! `true`, then `GetSeqs()` (the rows) and `GetSeqEntry()` (the Bioseq of each row).
//!
//! - `lines`: the peek-ahead line input (`CPeekAheadStream`).
//! - `format_guess`: the guess of the format (`CAlnFormatGuesser`), all six tests.
//! - `base`: what every scanner shares (`CAlnScanner`'s checks, the sequence info, the id and
//!   line helpers).
//! - `clustal`, `fasta_gap`, `phylip`, `nexus`, `sequin`, `multalign`: one scanner per format; each
//!   returns the ids, the rows and the deflines of the file (item PS6b of line PS).
//! - `mods`: the defline of a row, its `[name=value]` modifiers and its title (`x_AddMods`).
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

mod base;
mod clustal;
pub(crate) mod descr;
mod fasta_gap;
mod format_guess;
mod lines;
pub(crate) mod mods;
mod multalign;
pub(crate) mod nexus;
pub(crate) mod phylip;
mod sequin;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_formats;

pub(crate) use format_guess::{guess_format, AlignFormat};
pub(crate) use lines::PeekAheadStream;
pub(crate) use mods::SetMol;

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
    /// The title descriptor that the row's defline gives (`x_AddMods`, `x_AddTitle`).
    pub(crate) title: Option<Vec<u8>>,
    /// What the defline's modifiers set on the Bioseq besides the title: descriptors, the protein
    /// feature and Seq-inst members (`in_msa::pssm_from_msa` writes them with the PSSM's query).
    pub(crate) content: descr::Content,
    /// The molecule that the defline's modifiers set on the Bioseq of a protein row where it is
    /// not protein (`molecule`, `mol-type`): NCBI carries it into the PSSM's query and the search.
    pub(crate) set_mol: Option<SetMol>,
}

/// What `CAlnReader` gives `CPsiBlastInputClustalW`: `GetSeqs()` (the rows as in the file, with
/// their `-` gaps and their case) and the Bioseqs of `GetSeqEntry()`, in the order of the rows.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Alignment {
    pub(crate) rows: Vec<Vec<u8>>,
    pub(crate) bioseqs: Vec<AlnBioseq>,
    /// The `ERR_POST` warnings of the modifiers of every row, in NCBI's order (messages only).
    pub(crate) warnings: Vec<Vec<u8>>,
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
    /// A modifier of a defline is an error (`CObjReaderParseException` of `CDefaultModErrorReporter`,
    /// "BLAST query error"), after the warnings that NCBI posted before it.
    Modifier(mods::ModsFailure),
    /// An exception that nothing in the reader catches, thrown while a defline's modifiers are
    /// applied, after the warnings that NCBI posted before it.
    ModifierException(Vec<Vec<u8>>, Box<crate::ncbi_exception::NcbiException>),
    /// A `std::exception` of the reader that nothing catches (`what()`): `Error: ` and its text,
    /// exit status 255.
    StdException(String),
    /// NCBI's reader does something here that LOSAT does not reproduce (it crashes, reads memory
    /// that the input does not give, or ends with an exception that nothing catches and LOSAT
    /// does not print): the cause of the refusal.
    Undefined(&'static str),
    /// An exception of NCBI's toolkit that nothing in the reader catches (`CStringException` of
    /// `NStr::StringToInt`): it ends the run as any uncaught exception, with its report.
    Exception(crate::ncbi_exception::NcbiException),
    /// A state of the reader that no input gives: one of its two preconditions on the Seq-id
    /// parser does not hold ([`NO_SEQ_ID_WARNING`], [`AN_ID_OF_EVERY_ROW`]; named as LOSAT's
    /// internal error names it). Not a refusal: no alignment reaches it.
    Internal(&'static str),
}

/// The flags of `x_ParseAndValidateSeqIds` (`fParse_AnyLocal` and `fParse_RawText`, cited at
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

/// `CAlnReader::Read` with `false` and `true`, and `GetSeqs()` and `GetSeqEntry()` of `bytes`, the whole
/// input.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:319-339
/// The unknown format's scanner throws "Input file format not recognized.".
pub(crate) fn read_alignment(bytes: &[u8]) -> Result<Alignment, AlnError> {
    let mut input = PeekAheadStream::new(bytes);
    let format = guess_format(&mut input);
    let scanned = match format {
        AlignFormat::Clustal => clustal::scan(&mut input),
        AlignFormat::FastaGap => fasta_gap::scan(&mut input),
        AlignFormat::Phylip => phylip::scan(&mut input),
        AlignFormat::Nexus => nexus::scan(&mut input),
        AlignFormat::Sequin => sequin::scan(&mut input),
        AlignFormat::MultAlign => multalign::scan(&mut input),
        AlignFormat::Unknown => return Err(AlnError::NotAvailable),
    };
    let scanned = scanned.map_err(|failure| match failure {
        base::ScanFailure::Failed => AlnError::NotAvailable,
        base::ScanFailure::Undefined(cause) => AlnError::Undefined(cause),
        base::ScanFailure::Exception(exception) => AlnError::Exception(exception),
        base::ScanFailure::StdException(what) => AlnError::StdException(what),
    })?;
    let ids = verify_alignment_info(&scanned)?;
    entries(scanned, ids)
}

/// `CAlnReader::x_VerifyAlignmentInfo`: at least two sequences, and the ids of the rows parsed
/// and checked ([`row_ids`]).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:377-403
///
/// The defline count: a scanner that read deflines must have read one for each id, else the
/// reporter's error throws (NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:405-427).
fn verify_alignment_info(scanned: &base::Scanned) -> Result<Vec<Vec<ParsedSeqId>>, AlnError> {
    if scanned.rows.len() < 2 {
        return Err(AlnError::NotAvailable);
    }
    let ids: Vec<Vec<ParsedSeqId>> = scanned
        .ids
        .iter()
        .map(|id| row_ids(id))
        .collect::<Result<_, _>>()?;
    if !scanned.deflines.is_empty() && scanned.deflines.len() != ids.len() {
        return Err(AlnError::NotAvailable);
    }
    Ok(ids)
}

/// `CAlnReader::x_ParseAndValidateSeqIds` of the id string of one row: its Seq-ids, at least one.
/// A parse exception and a report of the id validator are read failures (the error reporter
/// throws: no listener).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:349-368
/// (`flags` is 0 here, so `flags` XOR `fGenerateLocalIDs` is not 0 and the raw text flag is set.)
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
fn entries(scanned: base::Scanned, ids: Vec<Vec<ParsedSeqId>>) -> Result<Alignment, AlnError> {
    let base::Scanned {
        rows,
        deflines,
        info,
        ..
    } = scanned;
    let mut bioseqs = Vec::with_capacity(rows.len());
    for (row, ids) in rows.iter().zip(ids) {
        let letters: Vec<u8> = row
            .iter()
            .map(u8::to_ascii_uppercase)
            .filter(|c| !info.gap.contains(c))
            .collect();
        let identified = identify_seq_id(
            &ids[0],
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |_| {},
        );
        let mol = if identified & F_NUC != 0 {
            Mol::Nucleotide
        } else if identified & F_PROT != 0 {
            Mol::Aa
        } else {
            let guessed: Vec<u8> = letters
                .iter()
                .copied()
                .filter(|c| !info.missing.contains(c))
                .collect();
            match crate::writedb::fasta_source::reader::sequence_type(&guessed) {
                Some(true) | None => Mol::Aa,
                Some(false) => {
                    if guessed.iter().any(|c| matches!(c, b'T' | b't'))
                        && guessed.iter().any(|c| matches!(c, b'U' | b'u'))
                    {
                        return Err(AlnError::InvalidMolType);
                    }
                    Mol::Nucleotide
                }
            }
        };
        bioseqs.push(AlnBioseq {
            ids,
            mol,
            letters,
            title: None,
            content: descr::Content::default(),
            set_mol: None,
        });
    }
    // `x_AddMods` of each row's defline (`TruncateSpaces` of what the scanner kept): after every
    // row's molecule is known.
    let mut warnings: Vec<Vec<u8>> = Vec::new();
    if !deflines.is_empty() {
        for (bioseq, defline) in bioseqs.iter_mut().zip(&deflines) {
            let defline = crate::ncbi_str::truncate_spaces(defline);
            let applied =
                mods::apply_defline(defline, bioseq.mol == Mol::Aa).map_err(|stop| match stop {
                    mods::ModsStop::Failure(mut failure) => {
                        let mut all = std::mem::take(&mut warnings);
                        all.append(&mut failure.warnings);
                        failure.warnings = all;
                        AlnError::Modifier(failure)
                    }
                    mods::ModsStop::Exception(exception, mut posted) => {
                        let mut all = std::mem::take(&mut warnings);
                        all.append(&mut posted);
                        AlnError::ModifierException(all, exception)
                    }
                    mods::ModsStop::Undefined(cause) => AlnError::Undefined(cause),
                })?;
            warnings.extend(applied.warnings);
            bioseq.title = applied.title;
            bioseq.content = applied.effects.content;
            // A nucleotide row stays what its letters are (its Bioseq ends the run alike, whatever
            // the modifiers say: `f6b.m2.dna_mol_aa`).
            bioseq.set_mol = applied
                .effects
                .mol
                .filter(|set| bioseq.mol == Mol::Aa && *set != SetMol::Aa);
        }
    }
    Ok(Alignment {
        rows,
        bioseqs,
        warnings,
    })
}
