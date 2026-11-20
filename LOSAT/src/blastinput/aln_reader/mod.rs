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
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_error_reporter.cpp:111-135
//! ```c++
//!     string message(descr);
//!     if (!seqId.empty()) {
//!         message = "At ID \'" + seqId + "\': " + descr;
//!     }
//!     if (!mpEl) {
//!         NCBI_THROW2(CObjReaderParseException, eFormat, message, 0);
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/msa_pssm_input.cpp:112-134
//! ```c++
//!     CAlnReader reader(input_file);
//!     reader.SetClustal(CAlnReader::eAlpha_Protein);
//!     try {
//!         reader.Read(false, true);
//!     } catch (const CObjReaderParseException& e) {
//!         // Workaround to provide a more useful error message when repeated
//!         // Seq-IDs are encountered
//!         if ((e.GetErrCode() == CObjReaderParseException::eFormat) &&
//!             (NStr::Find(e.GetMsg(), "Not all sequences have same length") != NPOS)) {
//!             ...
//!         }
//!     }
//!     m_AsciiMsa = reader.GetSeqs();
//!     m_SeqEntry = reader.GetSeqEntry();
//! ```
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
    /// What NCBI reads and LOSAT's reader does not (named as the refusal's text).
    Refused(&'static str),
}

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
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:319-340
/// ```c++
/// void CAlnReader::Read(
///     bool /*guess*/,
///     bool generate_local_ids,
///     ncbi::objects::ILineErrorListener* /*pErrorListener*/)
/// {
///     // read the alignment stream
///     SAlignmentFile alignmentInfo;
///     try {
///         ReadAlignmentFile(
///             m_IS, generate_local_ids, m_UseNexusInfo, mSequenceInfo, alignmentInfo);
///         TReadFlags flags = 0;
///         x_VerifyAlignmentInfo(alignmentInfo, flags);
///     }
///     catch (const SShowStopper& showStopper) {
///         theErrorReporter->Fatal(showStopper);
///         return;
///     }
///     m_Dim = m_IdStrings.size();
///     m_ReadDone = true;
///     m_ReadSucceeded = true;
/// }
/// ```
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
/// (`x_ParseAndValidateSeqIds`) and checked.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:239-286
/// ```c++
///     const auto num_sequences = alignmentInfo.NumSequences();
///
///     if (num_sequences == 0) {
///         throw SShowStopper(...);
///     }
///
///     if (num_sequences == 1) {
///         throw SShowStopper(...);
///     }
///     ...
///     for (auto seqIdInfo : alignmentInfo.mIds) {
///         m_IdStrings.push_back(seqIdInfo.mData); // m_IdStrings is redundant and should be removed
///         TIdList ids;
///         x_ParseAndValidateSeqIds(seqIdInfo, flags, ids);
///         m_Ids.push_back(ids);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:216-237
/// ```c++
///     CSeq_id::TParseFlags parseFlags = CSeq_id::fParse_AnyLocal;
///     if (flags^fGenerateLocalIDs) {
///         parseFlags |= CSeq_id::fParse_RawText;
///     }
///
///     try {
///         CSeq_id::ParseIDs(ids, idString, parseFlags);
///     }
///     catch (...) {  // report an error and interpret the id string as a local ID
///         theErrorReporter->Error(...);   // throws
///         ids.push_back(Ref(new CSeq_id(CSeq_id::e_Local, idString)));
///     }
///
///     if (m_fValidateIds) {
///         m_fValidateIds(ids, seqIdInfo.mNumLine, theErrorReporter.get());
///     }
/// ```
/// (`flags` is 0 here, so `flags^fGenerateLocalIDs` is not 0 and the raw text flag is set.)
fn verify_alignment_info(scanned: &clustal::Scanned) -> Result<Vec<Vec<ParsedSeqId>>, AlnError> {
    if scanned.rows.len() < 2 {
        return Err(AlnError::NotAvailable);
    }
    let mut all = Vec::with_capacity(scanned.ids.len());
    for id in &scanned.ids {
        let mut warned = false;
        let parsed = parse_ids(
            id,
            parse_flags::ANY_LOCAL | parse_flags::RAW_TEXT,
            &mut |_| warned = true,
        )
        .map_err(|_| AlnError::NotAvailable)?;
        // A Seq-id warning is a post that LOSAT would have to write as NCBI does; none comes
        // from the ids of the corpus of line PS PS6.
        if warned {
            return Err(AlnError::Refused(
                "a sequence id for which NCBI's Seq-id parser posts a warning",
            ));
        }
        let ids: Vec<ParsedSeqId> = parsed.into_iter().flatten().collect();
        if crate::writedb::fasta_source::defline_ids::alignment_ids_are_reported(&ids) {
            return Err(AlnError::NotAvailable);
        }
        if ids.is_empty() {
            // `ids.front()` of an empty list in `GetSeqEntry`: undefined behaviour in NCBI.
            return Err(AlnError::Refused(
                "a sequence id that NCBI's Seq-id parser reads as no id",
            ));
        }
        all.push(ids);
    }
    Ok(all)
}

/// `CAlnReader::GetSeqAlign`'s `m_SeqVec` and `GetSeqEntry`: the upper-case letters of each row
/// without gaps, its Bioseq with the molecule of its first id (`IdentifyAccession`) or, for an id
/// that names none, of its letters (`x_GetSequenceMolType`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:640-680
/// ```c++
///     for (TNumrow row_i = 0; row_i < m_Dim; row_i++) {
///         const string& seq_str     = m_SeqVec[row_i];
///         auto pSubEntry = Ref(new CSeq_entry());
///
///         // seq-id(s)
///         auto& ids = pSubEntry->SetSeq().SetId();
///         ids = m_Ids[row_i];
///
///         // mol
///         CSeq_inst::EMol mol   = CSeq_inst::eMol_not_set;
///         CSeq_id::EAccessionInfo ai = ids.front()->IdentifyAccession();
///         if (ai & CSeq_id::fAcc_nuc) {
///             mol = CSeq_inst::eMol_na;
///         } else if (ai & CSeq_id::fAcc_prot) {
///             mol = CSeq_inst::eMol_aa;
///         } else {
///             const string seqId = ids.front()->AsFastaString();
///             mol = x_GetSequenceMolType(GetAlphabet(), seq_str, seqId, pErrorListener);
///         }
///         // seq-inst
///         auto pSeqInst = x_GetSeqInst(mol, seq_str);
/// ```
/// `x_GetSequenceMolType` (the alphabet has 53 letters, so an undefined guess is protein):
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:553-600
/// ```c++
///     auto formatGuess = CFormatGuess::SequenceType(seqChars.data(), seqChars.length());
///     if (formatGuess == CFormatGuess::eProtein) {
///         return CSeq_inst::eMol_aa;
///     }
///
///     //if alphabet contains full complement (26) of protein chars then
///     // it's definitely protein. It may also contain stop-codon characters:
///     if (formatGuess == CFormatGuess::eUndefined &&
///         alphabet.size() >= 2*26) {
///         return CSeq_inst::eMol_aa;
///     }
///
///     auto posFirstT = seqChars.find_first_of("Tt");
///     auto posFirstU = seqChars.find_first_of("Uu");
///     if (posFirstT != string::npos  &&  posFirstU != string::npos) {
///         string msg = "Invalid Mol Type: "
///                      "U and T cannot appear in the same nucleotide sequence. "
///                      "Reinterpreting as protein.";
///         sReportError(pErrorListener, ...);   // no listener: throws CObjReaderParseException
///         ...
///     }
///     return (posFirstU == string::npos  ?  CSeq_inst::eMol_dna : CSeq_inst::eMol_rna);
/// ```
/// (The letters are upper case with no gap: `m_SeqVec`; the missing characters of the sequence
/// info are none.) An id that names a molecule decides before the letters are looked at.
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
