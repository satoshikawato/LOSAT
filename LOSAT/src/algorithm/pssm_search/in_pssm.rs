//! The PSSM file of `-in_pssm` (psiblast, series PS item PS4a; tblastn's PSI-TBLASTN, item PS7): the file read
//! as NCBI's `CPsiBlastArgs` reads it for both programs, NCBI's errors for files it cannot read, NCBI's checks of
//! the PSSM (`CPsiBlastValidate::Pssm`), the query record made from the PSSM's Bioseq, the score rows
//! (`CScorematPssmConverter::GetScores`) and the statistical values of `PsiBlastSetupScoreBlock`, the frequency
//! ratios (item PS4b), the scores of a PSSM of frequency ratios only (the PSSM engine's score block), and LOSAT's
//! refusals of what the PSSM search does not read (DECISIONS 11 of line PS). The reader belongs to the PSSM search
//! (`PsiSearch`), which psiblast and tblastn both run; their own options stay with them
//! (`algorithm::psiblast::in_pssm`, `algorithm::tblastn::in_pssm`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1363
//! ```c++
//!     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
//!         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
//!         _ASSERT(m_Pssm.Empty());
//!         m_Pssm.Reset(new CPssmWithParameters);
//!         try {
//!             switch (CFormatGuess().Format(in)) {
//!             case CFormatGuess::eBinaryASN:
//!                 in >> MSerial_AsnBinary >> *m_Pssm;
//!                 break;
//!             case CFormatGuess::eTextASN:
//!                 in >> MSerial_AsnText >> *m_Pssm;
//!                 break;
//!             case CFormatGuess::eXml:
//!                 in >> MSerial_Xml >> *m_Pssm;
//!                 break;
//!             default:
//!                 NCBI_THROW(CInputException, eInvalidInput,
//!                            "Unsupported format for PSSM");
//!             }
//!         } catch (const CSerialException&) {
//!             string msg("Unrecognized format for PSSM in ");
//!             msg += args[kArgPSIInputChkPntFile].AsString() + " (must be ";
//!             msg += "PssmWithParameters)";
//!             NCBI_THROW(CInputException, eInvalidInput, msg);
//!         }
//!         _ASSERT(m_Pssm.NotEmpty());
//!     }
//! ```

use std::path::Path;

use anyhow::Result;

use super::rows::QueryPssm;
use super::{PsiSearch, PssmStats};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ScoringMatrix;
use crate::core::blast_util::{BLASTAA_SIZE, BLAST_SCORE_MIN};
use crate::core::psi::{KarlinBlk, PsiScoreBlk};
use crate::report::query_warnings::invalid_residues_error;
use crate::report::serial::pssm::{ObjectId, PssmFinalData, PssmWithParameters, SeqData};
use crate::seqid::ParsedSeqId;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-118
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'};
// ```
/// The letters of NCBIstdaa, the one table of the PSSM reader and of psiblast's PSSM files (the ASCII PSSM, the query
/// letters of a PSSM file; `algorithm::psiblast::in_pssm_files`).
pub(crate) const NCBISTDAA_TO_AMINOACID: &[u8; BLASTAA_SIZE] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
/// `BLAST_SCORE_MAX` (`INT2_MAX`, blast_stat.h:122): the largest score of one letter pair.
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

/// What `CFormatGuess::Format` makes of the file, as far as `-in_pssm` reads it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PssmFileFormat {
    /// `eTextASN`.
    TextAsn,
    /// A format that the `switch` sends to "Unsupported format for PSSM".
    Unsupported,
    /// Binary ASN.1 (`eBinaryASN`), XML (`eXml`) or a text that LOSAT does not classify as NCBI
    /// does: refused (binary and XML PSSMs are item PS9's).
    NotText,
}

/// `CFormatGuess::GuessFormat` for the formats a PSSM file can have. The compressed formats
/// come first in NCBI's order (`sm_CheckOrder`: eBam, eZip, eZstd, eGZip, eBZip2, eLzo), each
/// tested by its magic bytes; a file starting with `>` is never `eTextASN` (and NCBI's earlier
/// `eFasta` test takes it); `eTextASN` is `TestFormatTextAsn` on the test buffer. Of the other
/// files, those with only printing bytes and white space are taken as NCBI's unsupported formats
/// (they cannot be `eBinaryASN`, whose test is a non-printing byte, nor `eXml`, whose test is a
/// leading `<?xml`, `<!DOCTYPE` or `<Blast4-request>`); the rest are refused.
///
/// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1162-1183
/// ```c++
///     //
///     //  Criterion: Presence of any non-printing characters
///     //
///     EConfidence conf = eNo;
///     for (int i = 0;  i < m_iTestDataSize;  ++i) {
///         if ( !isgraph((unsigned char) m_pTestBuffer[i])  &&
///              !isspace((unsigned char) m_pTestBuffer[i]) )
///         {
///             if (m_pTestBuffer[i] == '\1') {
///                 conf = eMaybe;
///             } else {
///                 return true;
///             }
///         }
///     }
///     return (conf == eYes);
/// ```
fn guess_format(bytes: &[u8]) -> PssmFileFormat {
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:640-652
    // ```c++
    //     const streamsize k_TestBufferGranularity = 8096;
    //     ...
    //         m_iTestBufferSize = Multiplier * k_TestBufferGranularity;
    //         m_pTestBuffer = new char[ m_iTestBufferSize ];
    //         m_Stream.read( m_pTestBuffer, m_iTestBufferSize );
    //         m_iTestDataSize = m_Stream.gcount();
    //         if (m_iTestDataSize == 0) {
    //             ...
    //             return false; //empty file
    // ```
    if bytes.is_empty() {
        return PssmFileFormat::Unsupported;
    }
    // gzip (and BAM, gzip too), zip, zstd, bzip2: NCBI's magic-byte tests.
    const MAGICS: [&[u8]; 4] = [b"\x1f\x8b", b"PK\x03\x04", b"\x28\xb5\x2f\xfd", b"BZh"];
    if MAGICS.iter().any(|magic| bytes.starts_with(magic)) {
        return PssmFileFormat::Unsupported;
    }
    let buffer = &bytes[..bytes.len().min(8096)];
    if text_asn(buffer) {
        return PssmFileFormat::TextAsn;
    }
    let printing = buffer
        .iter()
        .all(|&c| c.is_ascii_graphic() || matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'));
    let trimmed = buffer.trim_ascii_start();
    let xml = [b"<?xml".as_slice(), b"<!doctype"].iter().any(|prefix| {
        trimmed.len() >= prefix.len() && trimmed[..prefix.len()].eq_ignore_ascii_case(prefix)
    }) || trimmed.starts_with(b"<Blast4-request>");
    if printing && !xml {
        PssmFileFormat::Unsupported
    } else {
        PssmFileFormat::NotText
    }
}

/// `TestFormatTextAsn`: at least 80 % of the test buffer letters, digits, white space or braces
/// (counted line by line, each non-empty line with one line end, as `EnsureStats` counts), and
/// the first line that is not blank or a comment of the form `Name ::= ...` with a letter first.
///
/// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1679-1726
/// ```c++
///     // reject obvious misfits:
///     if ( m_iTestDataSize == 0 || m_pTestBuffer[0] == '>' ) {
///         return false;
///     }
///
///     // criteria:
///     // at least 80% text-ish,
///     // 1st field of the first non-blank not comment line must start with letter.
///     // "::=" as the 2nd field of the first non-blank non comment line.
///     //
///     double dAlNumFraction = (double)(m_iStatsCountAlNumChars+m_iStatsCountBraces) /
///                             (double)m_iTestDataSize;
///     if ( dAlNumFraction < 0.80 ) {
///         return false;
///     }
///     ...
///     while ( ! TestBuffer.fail() ) {
///         vector<string> Fields;
///         NcbiGetline(TestBuffer, strLine, "\n\r");
///         NStr::Split(strLine, " \t", Fields, NStr::fSplit_Tokenize);
///         if ( IsAsnComment( Fields  ) ) {
///             continue;
///         }
///         return ( Fields.size() >= 2 && Fields[1] == "::=" && isalpha(Fields[0][0]));
///     }
///     return false;
/// ```
/// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:714-733
/// ```c++
///     while ( ! TestBuffer.fail() ) {
///         NcbiGetline( TestBuffer, strLine, "\r\n" );
/// // code in CFormatGuess::Format counts line ends
/// // so, we will count them here as well
///         if (!strLine.empty()) {
///             strLine += '\n';
///         }
///         size_t size = strLine.size();
///         ...
///             if ( type & (fAlpha | fDigit | fSpace) ) {
///                 ++m_iStatsCountAlNumChars;
///             }
///             else if (c == '{'  ||  c == '}') {
///                 ++m_iStatsCountBraces;
///             }
/// ```
fn text_asn(buffer: &[u8]) -> bool {
    if buffer.is_empty() || buffer[0] == b'>' {
        return false;
    }
    let mut counted = 0usize;
    for line in buffer.split(|&c| c == b'\n' || c == b'\r') {
        if line.is_empty() {
            continue;
        }
        counted += line
            .iter()
            .filter(|&&c| {
                c.is_ascii_alphanumeric() || c.is_ascii_whitespace() || c == b'{' || c == b'}'
            })
            .count()
            + 1;
    }
    if (counted as f64) / (buffer.len() as f64) < 0.80 {
        return false;
    }
    for line in buffer.split(|&c| c == b'\n' || c == b'\r') {
        let fields: Vec<&[u8]> = line
            .split(|&c| c == b' ' || c == b'\t')
            .filter(|field| !field.is_empty())
            .collect();
        if fields.is_empty() || fields[0].starts_with(b"--") {
            continue;
        }
        return fields.len() >= 2 && fields[1] == b"::=" && fields[0][0].is_ascii_alphabetic();
    }
    false
}

/// Who reads `-in_pssm`: the program that LOSAT's refusals name, and the options' matrix, with
/// which NCBI's PSSM engine computes the scores of a PSSM of frequency ratios only
/// (`scores_from_freq_ratios`). psiblast reads with `PSIBLAST` (BLOSUM62: its `-in_pssm` refuses
/// the other matrices, DECISIONS 11 of line PS); tblastn's PSI-TBLASTN with its own options
/// (`algorithm::tblastn::in_pssm`, DECISIONS 35).
#[derive(Clone, Copy, Debug)]
pub(crate) struct PssmReader<'a> {
    /// The program of LOSAT's refusals (`PSIBLAST`, `TBLASTN`).
    pub(crate) program: &'static str,
    /// The program's name in NCBI's messages (`[psiblast]`, `[tblastn]`).
    pub(crate) ncbi_name: &'static str,
    /// The options' matrix name (`opts.GetMatrixName()`).
    pub(crate) matrix_name: &'a str,
}

impl PssmReader<'static> {
    /// psiblast's reader.
    pub(crate) const PSIBLAST: Self = Self {
        program: "PSIBLAST",
        ncbi_name: "psiblast",
        matrix_name: "BLOSUM62",
    };
}

impl PssmReader<'_> {
    /// The refusal of the PSSM search's limits ("... is not supported by LOSAT's PSIBLAST").
    fn refuse(&self, what: impl std::fmt::Display) -> anyhow::Error {
        anyhow::anyhow!("{what} is not supported by LOSAT's {}", self.program)
    }
}

/// Reads `-in_pssm path`: NCBI's errors for a file it cannot read or guess, LOSAT's refusals of
/// what it does not read, then the PSSM search's query, rows and values, and the PSSM as NCBI's
/// object holds it when the search starts (`InPssmFile`, which psiblast's rounds
/// and PSSM files write back).
pub(crate) fn read_in_pssm(
    path: &Path,
    reader: PssmReader<'_>,
    diagnostics: &mut dyn std::io::Write,
) -> Result<(PsiSearch, InPssmFile)> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
    // ```c++
    // void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
    // {
    //     if ( !m_Ios ) {
    //         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
    //             "File is not accessible",AsString()));
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
    // ```c++
    //     if (AsString() == "-") {
    // #if defined(NCBI_OS_MSWIN)
    //         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
    // #endif
    //         m_Ios  = &cin;
    // ```
    // `-in_pssm -` is the standard input (`kArgPSIInputChkPntFile` is an `eInputFile`), as
    // `-query -` is; NCBI's messages name the file `-`.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:726-731
    // ```c++
    //             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
    //             if ( !fstrm->is_open() ) {
    //                 delete fstrm;
    //                 fstrm = NULL;
    //             } else {
    //                 m_DeleteFlag = true;
    // ```
    // Only a file that cannot be opened is "not accessible": one that opens and cannot be read (a
    // directory) gives the reader no bytes, which NCBI's format guess takes as an unsupported
    // format (probe `pf.dir` of line PS PS7).
    let bytes = if path.as_os_str() == "-" {
        let mut bytes = Vec::new();
        std::io::Read::read_to_end(&mut std::io::stdin().lock(), &mut bytes).map(|_| bytes)
    } else {
        std::fs::File::open(path).map(|mut file| {
            // The bytes the stream gave before it failed (none from a directory).
            let mut bytes = Vec::new();
            let _ = std::io::Read::read_to_end(&mut file, &mut bytes);
            bytes
        })
    }
    .map_err(|_| crate::cli::inaccessible("in_pssm", path))?;
    match guess_format(&bytes) {
        PssmFileFormat::TextAsn => {}
        PssmFileFormat::Unsupported => {
            return Err(crate::blastinput::app::options_error(
                "Unsupported format for PSSM",
            ))
        }
        PssmFileFormat::NotText => {
            return Err(reader.refuse(format!(
                "a PSSM file that is not ASN.1 text (binary ASN.1, XML or another format; {})",
                path.display()
            )))
        }
    }
    let mut bad_chars = Vec::new();
    let read =
        crate::blastinput::pssm_text::read_pssm_with_parameters_and_skips(&bytes, &mut bad_chars);
    // NCBI's reader stores `#` for a byte of a VisibleString outside the visible range (the
    // titles, the row labels, the identifier, the matrix name, the Seq-ids' strings, the
    // sequence letters) and posts an error for it as it reads it, so also for a PSSM that it then
    // fails to read (psiblast 2.17.0 `-in_pssm ps4_rnd_ctrl.asn`: `Error: (803.7) [psiblast]
    // PssmWithParameters.pssm.query.seq.descr.E.title` / `Bad char [0xC3] in string at line 12`
    // / the chunk, once per byte, after the formatting handler's warnings, before the PSSM
    // search's own messages, status 0; tblastn the same with `[tblastn]`).
    for bad in &bad_chars {
        let mut lines = Vec::new();
        bad.post(Some(reader.ncbi_name), &mut lines);
        diagnostics.write_all(&lines)?;
    }
    let read = match read {
        Ok(read) => read,
        // The reader's own limits (a query Seq-entry set, Seq-ids and Seq-data it does not type,
        // REAL numbers NCBI computes with undefined behaviour or in `long double`): LOSAT's
        // refusal, where NCBI reads the file.
        Err(error) if error.message.contains("is not supported by LOSAT") => {
            return Err(reader.refuse(format!(
                "a PSSM that LOSAT's reader does not read ({}, line {}: {})",
                path.display(),
                error.line,
                error.message
            )))
        }
        // NCBI's `CSerialException`.
        Err(_) => {
            return Err(crate::blastinput::app::options_error(&format!(
                "Unrecognized format for PSSM in {} (must be PssmWithParameters)",
                path.display()
            )))
        }
    };
    psi_search_and_object(read.pssm, read.skipped, reader)
}

/// `CPsiBlastValidate::Pssm`'s error for `pssm`, which NCBI raises (`CBlastException`, "BLAST
/// engine error") when the search is made, after the outfmt 0 prolog.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:569-625
/// ```c++
///     bool missing_scores(false);
///     if ( !pssm.GetPssm().CanGetFinalData() ||
///          !pssm.GetPssm().GetFinalData().CanGetScores() ||
///          pssm.GetPssm().GetFinalData().GetScores().empty() ) {
///         missing_scores = true;
///     }
///
///     bool missing_freq_ratios(false);
///     if ( !pssm.GetPssm().CanGetIntermediateData() ||
///          !pssm.GetPssm().GetIntermediateData().CanGetFreqRatios() ||
///          pssm.GetPssm().GetIntermediateData().GetFreqRatios().empty() ) {
///         missing_freq_ratios = true;
///     }
///
///     if (missing_freq_ratios && missing_scores) {
///         NCBI_THROW(CBlastException, eInvalidArgument,
///                    "PSSM data must contain either scores or frequency ratios");
///     }
///     ...
///     // Only unscaled PSSMs are supported
///     if (!missing_scores &&
///         pssm.GetPssm().GetFinalData().CanGetScalingFactor() &&
///         pssm.GetPssm().GetFinalData().GetScalingFactor() != 1) {
///         string msg("PSSM has a scaling factor of ");
///         msg += NStr::IntToString(pssm.GetPssm()
///                                  .GetFinalData()
///                                  .GetScalingFactor());
///         msg += ". PSI-BLAST does not accept scaled PSSMs";
///         NCBI_THROW(CBlastException, eInvalidArgument, msg);
///     }
///
///     if ( !pssm.HasQuery() ) {
///     ...
///     if ( !pssm.GetPssm().GetIsProtein() ) {
///         NCBI_THROW(CBlastException, eInvalidArgument,
///                    "PSSM does not represent protein scoring matrix");
///     }
/// ```
/// (The query checks between them are never reached from psiblast's application, which reads
/// the query of the PSSM before; LOSAT refuses those PSSMs, `psi_search_of`.)
fn validation_error(pssm: &PssmWithParameters) -> Option<String> {
    let scores = pssm.pssm.final_data.as_ref().map(|data| &data.scores);
    let missing_scores = scores.is_none_or(|scores| scores.is_empty());
    let missing_freq_ratios = pssm
        .pssm
        .intermediate_data
        .as_ref()
        .is_none_or(|data| data.freq_ratios.is_empty());
    if missing_freq_ratios && missing_scores {
        return Some("PSSM data must contain either scores or frequency ratios".to_string());
    }
    if !missing_scores {
        if let Some(factor) = pssm
            .pssm
            .final_data
            .as_ref()
            .and_then(|data| data.scaling_factor)
            .filter(|&factor| factor != 1)
        {
            return Some(format!(
                "PSSM has a scaling factor of {factor}. PSI-BLAST does not accept scaled PSSMs"
            ));
        }
    }
    if pssm.pssm.is_protein == Some(false) {
        return Some("PSSM does not represent protein scoring matrix".to_string());
    }
    None
}

/// The query's residues as NCBI's `GetSequenceProtein` reads them, and the 0-based places of
/// those that are no NCBIstdaa letter. NCBI's `CSeqVector` converts the Bioseq's coding to
/// NCBIstdaa letter by letter, with the tables of its own (a byte missing from a table becomes
/// 255), and `s_IsValidResidue` accepts a result below `BLASTAA_SIZE`. The returned letters have
/// one place per residue, an invalid one shown as `X` (NCBI leaves it out and stops, so nothing
/// reads it), the others the capital letters of the codes. The tables differ: `iupacaa` knows the 26 letters, either case, and neither `*`
/// nor `-`; `ncbieaa` knows the 26 capital letters, `*` (25) and `-` (0); `ncbistdaa` is its own
/// code (0 to 27).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:884-896
/// ```c++
///     for (i = 0; i < sv.size(); i++) {
///         // Change unsupported residues to X
///         if (sv[i] == AMINOACID_TO_NCBISTDAA[(int)'O']) {
///             replaced_residues.push_back(i);
///             *buf_var++ = AMINOACID_TO_NCBISTDAA[(int)'X'];
///         } else if (!s_IsValidResidue(sv[i])) {
///             invalid_residues.push_back(i);
///         } else {
///             *buf_var++ = sv[i];
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1783 (`CIupacaaToStdaa`,
/// `CEaaToStdaa`: 255 for every byte but the letters, `*` and `-` named above), and
/// sequtil_convert_imp.cpp:179-198 (the table each coding pair takes).
fn query_residues(seq_data: Option<&SeqData>) -> (Vec<u8>, Vec<usize>) {
    // The NCBIstdaa code of a capital letter, `*` or `-`: its place in the letter table.
    let code_of = |letter: u8| {
        NCBISTDAA_TO_AMINOACID
            .iter()
            .position(|&known| known == letter)
    };
    let codes: Vec<Option<u8>> = match seq_data {
        Some(SeqData::Ncbistdaa(codes)) => codes
            .iter()
            .map(|&code| (usize::from(code) < BLASTAA_SIZE).then_some(code))
            .collect(),
        Some(SeqData::Iupacaa(text)) => text
            .iter()
            .map(|&c| {
                c.is_ascii_alphabetic()
                    .then(|| code_of(c.to_ascii_uppercase()))
                    .flatten()
                    .map(|code| code as u8)
            })
            .collect(),
        Some(SeqData::Ncbieaa(text)) => text
            .iter()
            .map(|&c| {
                (c.is_ascii_uppercase() || c == b'*' || c == b'-')
                    .then(|| code_of(c))
                    .flatten()
                    .map(|code| code as u8)
            })
            .collect(),
        None => Vec::new(),
    };
    let letters = codes
        .iter()
        .map(|code| code.map_or(b'X', |code| NCBISTDAA_TO_AMINOACID[usize::from(code)]))
        .collect();
    let invalid = codes
        .iter()
        .enumerate()
        .filter(|(_, code)| code.is_none())
        .map(|(position, _)| position)
        .collect();
    (letters, invalid)
}

/// The PSSM search of `pssm`: its query record, rows and values, or LOSAT's refusal.
#[cfg(test)]
pub(crate) fn psi_search_of(pssm: PssmWithParameters) -> Result<PsiSearch> {
    psi_search_and_object(pssm, Vec::new(), PssmReader::PSIBLAST).map(|(search, _)| search)
}

/// The PSSM search of `pssm` and the PSSM as NCBI's `CPsiBlastImpl` leaves it before the first
/// round: a PSSM of frequency ratios only gets the scores and gapped values that the PSSM engine
/// computes for it (`with_computed_scores`). `skipped` names the members the reader
/// skipped (`blastinput::pssm_text::read_pssm_with_parameters_and_skips`).
pub(crate) fn psi_search_and_object(
    pssm: PssmWithParameters,
    skipped: Vec<&'static str>,
    reader: PssmReader<'_>,
) -> Result<(PsiSearch, InPssmFile)> {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:535-540
    // ```c++
    // 		_ASSERT(pssm->HasQuery());
    // 		_ASSERT(pssm->GetQuery().IsSeq());  // single query only!
    // 		scope->AddTopLevelSeqEntry(pssm->SetQuery());
    //
    // 		CRef<CSeq_loc> sl(new CSeq_loc());
    // 		sl->SetWhole().Assign(*pssm->GetQuery().GetSeq().GetFirstId());
    // ```
    // NCBI's application takes the query's Bioseq here (a PSSM without one ends with a C++
    // exception of the release build); the reader has refused a Seq-entry set.
    let Some(query) = pssm.pssm.query.as_ref() else {
        return Err(reader
            .refuse("a PSSM without a query sequence (NCBI BLAST+ stops with a C++ exception)"));
    };
    let validation_error = validation_error(&pssm);
    // A PSSM of frequency ratios only: NCBI computes its scores with the PSSM engine
    // (`x_CreatePssmScoresFromFrequencyRatios`, below).
    let scores: &[i32] = match pssm.pssm.final_data.as_ref() {
        Some(data) if !data.scores.is_empty() => data.scores.as_slice(),
        _ => &[],
    };
    let ratios_only = validation_error.is_none() && scores.is_empty();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:648-660
    // ```c++
    //     asn1_pssm.SetByRow(false);  // this is the default
    //     ...
    //     if (asn1_pssm.GetByRow() == false) {
    //         for (unsigned int i = 0; i < pssm->ncols; i++) {
    //             for (unsigned int j = 0; j < pssm->nrows; j++) {
    //                 asn1_pssm.SetFinalData().SetScores().
    //                     push_back(pssm->pssm[i][j]);
    // ```
    // The engine's scores are column by column, which NCBI's search then reads by the PSSM's
    // own `byRow` (`PsiBlastComputePssmScores` copies the scores, not `byRow`): a ratio-only
    // PSSM written row by row would be searched with its scores transposed.
    if ratios_only && pssm.pssm.by_row == Some(true) {
        return Err(reader.refuse(
            "a PSSM of frequency ratios only written by row (NCBI reads the scores it computes for it in the wrong order)",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537-540 (quoted above)
    // The query's Seq-ids: NCBI's scope holds the PSSM's Bioseq with all of them, and the query's
    // Seq-loc (so its Seq-aligns, the PSSMs of the later rounds and its warnings) takes the first.
    // One local id (`local str "Query_N"` in a psiblast PSSM) is shown as a FASTA query's local id
    // is: the reports show the title alone, and with `-parse_deflines` the id
    // (`PsiSearch::query_record`); the record's local id is the Seq-id's `GetSeqIdString()`, the
    // string or the integer's decimal digits (`CObject_id::AsString`). Any other id (`gi 4504347`
    // of a `-parse_deflines` PSSM, `swissprot`, `pdb`, ...), and several ids, are shown with or
    // without `-parse_deflines`, as `-parse_deflines` shows the ids of a FASTA query (NCBI BLAST+
    // 2.17.0: `Query= gi|4504347 <title>`, `qseqid` `gi|4504347`; `GetSeqIdString(cbs,
    // believe_local_id)` asks `believe_local_id` of a local id only, align_format_util.cpp).
    let ids = query.ids.as_slice();
    if ids.is_empty() {
        return Err(reader.refuse(
            "a PSSM query without a Seq-id (NCBI BLAST+ makes the query's Seq-loc of its first Seq-id)",
        ));
    }
    if let Some(refusal) = crate::blastinput::fasta_reader::believed::refused_id(
        ids,
        "query",
        "the PSSM",
        reader.program,
    ) {
        return Err(refusal);
    }
    // A local id among other ids: without `-parse_deflines` NCBI's reports show the title alone
    // and its tabular ids put the title's first word in the local id's place (`qseqid`
    // `gi|4504347|hemoglobin`, `s_ReplaceLocalId`), a mix that no FASTA query has (one without
    // `-parse_deflines` has its local id only); NCBI's psiblast writes the one Seq-id of the
    // query's Seq-aligns.
    if ids.len() > 1 && ids.iter().any(|id| matches!(id, ParsedSeqId::Local(_))) {
        return Err(reader.refuse(
            "a PSSM query with a local Seq-id among other Seq-ids (NCBI's reports show it as no FASTA query is shown)",
        ));
    }
    // A Textseq-id with neither an accession nor a name: NCBI's label of it is empty, and NCBI
    // BLAST+ 2.17.0 stops on it (`BLAST Database error: Accessions to Oids lookup error in tfp`,
    // exit 2, `-db`; line RG Q1b's review, S1).
    if ids.iter().any(|id| {
        id.textseq()
            .is_some_and(|text| text.accession.is_none() && text.name.is_none())
    }) {
        return Err(reader.refuse(
            "a PSSM query whose Seq-id has neither an accession nor a name (NCBI BLAST+ stops with a database error)",
        ));
    }
    // An empty local id: NCBI's reports show it as no id. (A byte outside the visible range is
    // the reader's, `read_in_pssm`.)
    if ids
        .iter()
        .any(|id| matches!(id, ParsedSeqId::Local(ObjectId::Str(name)) if name.is_empty()))
    {
        return Err(reader.refuse("a PSSM query whose Seq-id holds an empty local id"));
    }
    let local_id = match ids {
        [ParsedSeqId::Local(ObjectId::Str(id))] => Some(String::from_utf8_lossy(id).into_owned()),
        [ParsedSeqId::Local(ObjectId::Id(id))] => Some(id.to_string()),
        _ => None,
    };
    let titles: &[Vec<u8>] = query.titles.as_deref().unwrap_or(&[]);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:669-683
    // ```c++
    //         ITERATE(CBioseq::TDescr::Tdata, iter, data) {
    //             if((*iter)->IsTitle()) {
    //                 all_descr_str += (*iter)->GetTitle();
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
    // ```c++
    //     string all_id_str = GetSeqIdString(cbs, believe_query);
    //     all_id_str += " ";
    //     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
    // ```
    // NCBI's `Query=` and `# Query:` lines join every title of the Bioseq and trim the spaces at
    // the line's ends, `qseqid` reads the first title's first word (a title that starts with a
    // space has an empty first word, so `qseqid` is the local id), and the warnings read the first
    // title without its trailing `.` and spaces (`GetTitle`, below). A FASTA
    // defline gives one title without white space at its ends (psiblast's own PSSMs hold such a
    // title); the PSSM search reads that one only.
    if titles.len() > 1 {
        return Err(reader.refuse(format!(
            "a PSSM query with {} titles (NCBI's Query= line joins them, its tabular query id and warnings read the first)",
            titles.len()
        )));
    }
    let title: &[u8] = titles.first().map_or(&[], |title| title.as_slice());
    // `NStr::TruncateSpaces` (`eTrunc_Both`) of the `Query=` line would cut the title's ends.
    if crate::ncbi_str::truncate_spaces(title).len() != title.len() {
        return Err(reader.refuse(
            "a PSSM query whose title starts or ends with white space (NCBI's Query= line, tabular query id and warnings each read such a title their own way)",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/bioseq_extract_data_priv.cpp:345-372
    // ```c++
    // string
    // CBlastQuerySourceBioseqSet::GetTitle(int index) const
    // {
    //     ...
    //     ITERATE(CSeq_descr::Tdata, desc, descr) {
    //         if ((*desc)->Which() == CSeqdesc::e_Title && title == kEmptyStr) {
    //             title = (*desc)->GetTitle();
    //         }
    //         if ((*desc)->Which() == CSeqdesc::e_Molinfo) {
    //             has_molinfo = true;
    //         }
    //     }
    //     if (title != kEmptyStr && !has_molinfo) {
    //         while (NStr::EndsWith(title, ".") || NStr::EndsWith(title, " ")) {
    //             title.erase(title.end() - 1);
    //         }
    //         retval.assign(title);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:115-118 (quoted at
    // `CPsiBlastImpl::x_ExtractQueryFromPssm`: the PSSM's Bioseq goes through
    // `CObjMgrFree_QueryFactory`, whose query source is `CBlastQuerySourceBioseqSet`)
    // The warnings of a PSSM query (`blast_setup_cxx.cpp:533-543`, quoted at
    // `report::query_warnings::warning_query_id`) take this title, in every round (the query
    // factory is fixed for the run), while the reports take the Bioseq's title as it is.
    let warning_title: Vec<u8> = if skipped.contains(&"molinfo") {
        Vec::new()
    } else {
        let end = title
            .iter()
            .rposition(|&c| c != b'.' && c != b' ')
            .map_or(0, |last| last + 1);
        title[..end].to_vec()
    };
    if query.inst.mol.1 != 3 || query.inst.repr.1 != 2 {
        return Err(reader.refuse(format!(
            "a PSSM query whose Seq-inst is not a raw protein (repr {}, mol {})",
            query.inst.repr.0, query.inst.mol.0
        )));
    }
    let (letters, invalid_residues) = query_residues(query.inst.seq_data.as_ref());
    let length = query.inst.length.unwrap_or(-1);
    if letters.is_empty() || letters.len() as i64 != i64::from(length) {
        return Err(reader.refuse(
            "a PSSM query whose letters are missing or of another length than its Seq-inst (NCBI BLAST+ reads such a query past or short of its data)",
        ));
    }
    // A gap letter (`-` of ncbieaa, code 0 of ncbistdaa) is a valid residue to NCBI, which searches
    // with it as a letter of the query; LOSAT does not (the invalid residues, which NCBI reports
    // first, are the search's stop, `invalid_residues`).
    if invalid_residues.is_empty() && letters.contains(&b'-') {
        return Err(reader.refuse("a PSSM query with a gap letter"));
    }
    // A lower-case iupacaa letter converts as its capital, but NCBI's reports count the identities
    // of the Bioseq's own characters against the subject's (a lower-case query has none: probe
    // `lower_iupac`, 0 of 9 where LOSAT counts 8), which LOSAT's report does not follow.
    if invalid_residues.is_empty()
        && matches!(query.inst.seq_data.as_ref(),
            Some(SeqData::Iupacaa(text)) if text.iter().any(u8::is_ascii_lowercase))
    {
        return Err(reader.refuse(
            "a PSSM query with lower-case iupacaa letters (NCBI's identity counts compare the Bioseq's own characters)",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:240-261
    // ```c++
    //     const CPssm& pssm = pssm_asn.GetPssm();
    //     _ASSERT((size_t)pssm.GetFinalData().GetScores().size() ==
    //            (size_t)pssm.GetNumRows()*pssm_asn.GetPssm().GetNumColumns());
    //
    //     unique_ptr< CNcbiMatrix<int> > retval
    //         (new CNcbiMatrix<int>(BLASTAA_SIZE,
    //                               pssm.GetNumColumns(),
    //                               BLAST_SCORE_MIN));
    //
    //     Convert2Matrix(pssm.GetFinalData().GetScores(),
    //                    *retval, pssm.GetByRow(), pssm.GetNumRows(),
    //                    pssm.GetNumColumns());
    // ```
    // NCBI only asserts the sizes (its release build reads or writes past the matrix): LOSAT
    // refuses a PSSM whose sizes do not match its data, its query, or the 28 rows.
    let rows_n = pssm.pssm.num_rows;
    let columns_n = pssm.pssm.num_columns;
    let cells = i64::from(rows_n) * i64::from(columns_n);
    let ratios_len = pssm
        .pssm
        .intermediate_data
        .as_ref()
        .map_or(0, |data| data.freq_ratios.len());
    if validation_error.is_none()
        && (!(1..=BLASTAA_SIZE as i32).contains(&rows_n)
            || i64::from(columns_n) != letters.len() as i64
            || (!ratios_only && scores.len() as i64 != cells)
            || (ratios_len != 0 && ratios_len as i64 != cells))
    {
        return Err(reader.refuse(format!(
            "a PSSM whose dimensions ({rows_n} rows, {columns_n} columns) do not match its scores ({}), frequency ratios ({ratios_len}) or query ({} letters)",
            scores.len(),
            letters.len()
        )));
    }
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:121-122
    // ```c
    // #define BLAST_SCORE_MIN INT2_MIN   /**< minimum allowed score (for one letter comparison). */
    // #define BLAST_SCORE_MAX INT2_MAX   /**< maximum allowed score (for one letter comparison). */
    // ```
    // NCBI copies the scores into the matrix without a check and adds them in `Int4` (a score of
    // 40,000 runs; one near `INT4_MAX` would overflow, undefined behaviour). LOSAT keeps the scores
    // to one letter's range, so that the sums of a query that is not split (under 19,800 letters)
    // stay inside `i32`: a range rule of LOSAT's, not an NCBI error.
    if validation_error.is_none() {
        if let Some(score) = scores
            .iter()
            .find(|&&score| !(BLAST_SCORE_MIN..=BLAST_SCORE_MAX).contains(&score))
        {
            return Err(reader.refuse(format!(
                "a PSSM score of {score}, outside BLAST_SCORE_MIN..BLAST_SCORE_MAX ({BLAST_SCORE_MIN}..{BLAST_SCORE_MAX}, one letter's score range in NCBI BLAST+),"
            )));
        }
    }
    let by_row = pssm.pssm.by_row == Some(true);
    let mut rows = if validation_error.is_none() {
        convert_2_matrix(
            scores,
            rows_n as usize,
            columns_n as usize,
            by_row,
            BLAST_SCORE_MIN,
        )
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:264-284
    // ```c++
    //     if ( !pssm_asn.GetPssm().CanGetIntermediateData() ||
    //          !pssm_asn.GetPssm().GetIntermediateData().CanGetFreqRatios() ||
    //          pssm_asn.GetPssm().GetIntermediateData().GetFreqRatios().empty() ) {
    //         throw runtime_error("Cannot obtain frequency ratios from ASN.1 PSSM");
    //     }
    //     ...
    //     unique_ptr< CNcbiMatrix<double> > retval
    //         (new CNcbiMatrix<double>(BLASTAA_SIZE, pssm.GetNumColumns(), 0.0));
    //
    //     Convert2Matrix(pssm.GetIntermediateData().GetFreqRatios(),
    //                    *retval, pssm.GetByRow(), pssm.GetNumRows(),
    //                    pssm.GetNumColumns());
    // ```
    // The frequency ratios (`psi_matrix->freq_ratios`, 28 per position, 0.0 past the file's
    // rows), which the position-based composition-based statistics read; none when the PSSM has
    // none (NCBI's then stay zero).
    let freq_ratios = match pssm.pssm.intermediate_data.as_ref() {
        Some(data) if validation_error.is_none() && !data.freq_ratios.is_empty() => {
            convert_2_matrix(
                &data.freq_ratios,
                rows_n as usize,
                columns_n as usize,
                by_row,
                0.0,
            )
        }
        _ => Vec::new(),
    };
    let final_data = pssm.pssm.final_data.as_ref();
    let stat = |value: Option<f64>| value.filter(|&value| value != 0.0);
    let mut stats = PssmStats {
        lambda_ungapped: stat(final_data.and_then(|data| data.lambda_ungapped)),
        kappa_ungapped: stat(final_data.and_then(|data| data.kappa_ungapped)),
        h_ungapped: stat(final_data.and_then(|data| data.h_ungapped)),
        lambda: stat(final_data.map(|data| data.lambda)),
        kappa: stat(final_data.map(|data| data.kappa)),
        h: stat(final_data.map(|data| data.h)),
    };
    let mut computed = None;
    let mut validation_error = validation_error;
    if ratios_only {
        match scores_from_freq_ratios(&letters, &freq_ratios, reader)? {
            Ok((scores, gapped)) => {
                rows = scores;
                stats.lambda = stat(Some(gapped.lambda));
                stats.kappa = stat(Some(gapped.k));
                stats.h = stat(Some(gapped.h));
                computed = Some(gapped);
            }
            // The engine's error, raised where `CPsiBlastValidate::Pssm`'s are (after the
            // outfmt 0 prolog).
            Err(engine_error) => validation_error = Some(engine_error),
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-651
    // The query set-up (`SetupQueries_OMF`) comes after the checks and the scores above (the
    // search object's, `CPsiBlastImpl`), and stops a search with no valid query.
    // `-in_pssm` has the one query, so the search stops as the checks of the PSSM do, after
    // the outfmt 0 prolog (`report::query_warnings::invalid_residues_error`).
    if validation_error.is_none() {
        validation_error = invalid_residues_error(&invalid_residues);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:110-112
    // ```c++
    //     if (pssm->GetPssm().GetLambda() != CPssm::kInvalidStat) {
    //         score_blk->kbp_gap_psi[0]->Lambda = pssm->GetPssm().GetLambda();
    //     } else if (score_blk->kbp_gap_std[0]->Lambda > 0.0) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:43-50
    // ```c
    // static Boolean s_BlastKarlinBlkIsValid(const Blast_KarlinBlk* kbp)
    // {
    //     if ( !kbp ) {
    //         return FALSE;
    //     } else {
    //         return (kbp->Lambda > 0 && kbp->K > 0 && kbp->H > 0);
    //     }
    // }
    // ```
    // NCBI takes a negative value of the file into `kbp_psi` or `kbp_gap_psi` (0 is
    // `kInvalidStat`: the query's BLOSUM62 value). A `kbp_psi` that is not valid stops the search's
    // setup with NCBI's error of a query without Karlin-Altschul parameters
    // (`blastp::blast_engine::psiblast_stop_at_invalid_query`, line PS DECISIONS 30). A
    // `kbp_gap_psi` that is not valid does not stop NCBI: it runs with it (no gapped X-drop, a
    // cutoff of 1, E-values of -1), an output without meaning, refused as a PSSM whose scaling gives
    // no Karlin-Altschul parameters is (line PS DECISIONS 28), unless the setup stops first. NCBI's
    // check of the PSSM's file (`validation_error`) comes first.
    let negative = |values: [(&'static str, Option<f64>); 3]| {
        values
            .into_iter()
            .find(|(_, value)| value.is_some_and(|value| value < 0.0))
    };
    let ungapped_invalid = negative([
        ("lambdaUngapped", stats.lambda_ungapped),
        ("kappaUngapped", stats.kappa_ungapped),
        ("hUngapped", stats.h_ungapped),
    ])
    .is_some();
    let gapped_negative = negative([
        ("lambda", stats.lambda),
        ("kappa", stats.kappa),
        ("h", stats.h),
    ]);
    if let (Some((name, _)), None, false) =
        (gapped_negative, validation_error.as_ref(), ungapped_invalid)
    {
        return Err(reader.refuse(format!(
            "a PSSM whose gapped Karlin-Altschul parameter {name} is negative (NCBI BLAST+ 2.17.0 goes on with a Karlin-Altschul block that it holds invalid)"
        )));
    }
    let search = PsiSearch {
        pssm: QueryPssm { rows },
        freq_ratios,
        stats,
        query: match local_id {
            Some(local_id) => FastaRecord::new(local_id, title, &letters),
            None => FastaRecord::new(String::new(), title, &letters).believing(ids.to_vec()),
        },
        query_ids: ids.to_vec(),
        validation_error,
        setup_messages: Vec::new(),
        warning_title: Some(warning_title),
    };
    let mut pssm = pssm;
    if let Some(gapped) = computed {
        with_computed_scores(&mut pssm, &search.pssm.rows, gapped);
    }
    let file = InPssmFile {
        object: pssm,
        skipped,
        computed: computed.is_some(),
        valid: search.validation_error.is_none(),
    };
    Ok((search, file))
}

/// The scores of a PSSM of frequency ratios only, and the gapped values of its scaling:
/// `PsiBlastComputePssmScores` runs the PSSM engine on the ratios (padded to 28 rows) of the
/// query's letters as the search holds them (O as X), with the options' matrix (`reader`:
/// BLOSUM62 for psiblast, which runs no other matrix with `-in_pssm`) and the gap costs 11/1 (the
/// defaults of `CPsiBlastInputFreqRatios`), and keeps the scores and `lambda`, `kappa`, `h`; the
/// ungapped values stay absent (`kbp_psi` then takes `kbp_std`'s, `PsiSearch::ungapped_karlin`).
/// A matrix whose table has no 11/1 row stops the engine with NCBI's error (the inner `Err`,
/// `CBlastException` "BLAST engine error" when the search is made); a matrix NCBI does not know
/// is the options' error that follows the reading (`validate_protein_options`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:466-485
/// ```c++
/// void PsiBlastComputePssmScores(CRef<objects::CPssmWithParameters> pssm,
///                                const CBlastOptions& opts)
/// {
///     CConstRef<CBioseq> query(&pssm->GetQuery().GetSeq());
///     CRef<IQueryFactory> seq_fetcher(new CObjMgrFree_QueryFactory(query)); /* NCBI_FAKE_WARNING */
///
///     CRef<ILocalQueryData> query_data(seq_fetcher->MakeLocalQueryData(&opts));
///     BLAST_SequenceBlk* seqblk = query_data->GetSequenceBlk();
///     ...
///     unique_ptr< CNcbiMatrix<double> > freq_ratios
///         (CScorematPssmConverter::GetFreqRatios(*pssm));
///
///     CPsiBlastInputFreqRatios pssm_engine_input(seqblk->sequence,
///                                                seqblk->length,
///                                                *freq_ratios,
///                                                opts.GetMatrixName());
///     CPssmEngine pssm_engine(&pssm_engine_input);
///     CRef<CPssmWithParameters> pssm_with_scores(pssm_engine.Run());
/// ```
fn scores_from_freq_ratios(
    letters: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    reader: PssmReader<'_>,
) -> Result<std::result::Result<(Vec<[i32; BLASTAA_SIZE]>, crate::core::psi::KarlinBlk), String>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:979-985
    // ```c
    //         opts->matrix = strdup(matrix_name);
    //         /* Make it all upper case */
    //         for (i=0; i<strlen(opts->matrix); ++i)
    //             opts->matrix[i] = toupper((unsigned char) opts->matrix[i]);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3531-3561
    // ```c
    //    Int2 status =
    //       Blast_KarlinBlkGappedLoadFromTables(kbp, gap_open,
    //                                           gap_extend, matrix_name, FALSE);
    // ...
    //       else if (status == 2)
    //       {
    //          snprintf(buffer, sizeof(buffer), "Gap existence and extension values of %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, matrix_name);
    // ```
    // The engine's score block takes the matrix with the gap costs 11/1.
    let matrix_name = reader.matrix_name.to_ascii_uppercase();
    let Ok(matrix) = matrix_name.parse::<crate::config::ScoringMatrix>() else {
        return Ok(Err(format!("{matrix_name} is not a supported matrix")));
    };
    if crate::stats::protein_options::karlin_blk_gapped_load_from_tables(11, 1, &matrix_name, false)
        .is_err()
    {
        return Ok(Err(format!(
            "Gap existence and extension values of 11 and 1 not supported for {matrix_name}"
        )));
    }
    let query: Vec<u8> = letters
        .iter()
        .map(
            |&letter| match crate::utils::matrix::aa_char_to_ncbistdaa(letter) {
                crate::utils::matrix::ncbistdaa::O => crate::utils::matrix::ncbistdaa::X,
                code => code,
            },
        )
        .collect();
    let spec = crate::config::ProteinScoringSpec {
        matrix,
        gap_open: 11,
        gap_extend: 1,
    };
    let gapped = crate::stats::tables::lookup_protein_params(&spec);
    let sbp = score_block(
        spec.matrix,
        crate::core::psi::KarlinBlk {
            lambda: gapped.lambda,
            k: gapped.k,
            h: gapped.h,
        },
    )?;
    crate::core::psi::pssm_from_freq_ratios(&query, freq_ratios, &sbp)
        .map(Ok)
        .map_err(|_| {
            reader.refuse("a PSSM of frequency ratios from which NCBI's PSSM engine makes no scaled scores (NCBI BLAST+ 2.17.0 does not finish or stops with an error)")
        })
}

/// `GetScores` or `GetFreqRatios` + the copy of `PsiBlastSetupScoreBlock`:
/// `rows[position][residue]` of 28 cells, those past the file's rows `fill` (`BLAST_SCORE_MIN`
/// for the scores, 0.0 for the frequency ratios).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:215-237
/// ```c++
/// template <class T>
/// void Convert2Matrix(const list<T>& source, CNcbiMatrix<T>& dest,
///                     bool by_row, SIZE_TYPE num_rows, SIZE_TYPE num_columns)
/// {
///     typename list<T>::const_iterator itr = source.begin();
///     if (by_row == true) {
///         for (SIZE_TYPE r = 0; r < num_rows; r++) {
///             for (SIZE_TYPE c = 0; c < num_columns; c++) {
///                 dest(r, c) = *itr++;
///             }
///         }
///     } else {
///         for (SIZE_TYPE c = 0; c < num_columns; c++) {
///             for (SIZE_TYPE r = 0; r < num_rows; r++) {
///                 dest(r, c) = *itr++;
///             }
///         }
///     }
///     _ASSERT(itr == source.end());
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:140-144
/// ```c++
///         for (TSeqPos c = 0; c < scores->GetCols(); c++) {
///             for (TSeqPos r = 0; r < scores->GetRows(); r++) {
///                 score_blk->psi_matrix->pssm->data[c][r] = (*scores)(r, c);
///             }
///         }
/// ```
pub(crate) fn convert_2_matrix<T: Copy>(
    values: &[T],
    num_rows: usize,
    num_columns: usize,
    by_row: bool,
    fill: T,
) -> Vec<[T; BLASTAA_SIZE]> {
    let mut rows: Vec<[T; BLASTAA_SIZE]> = (0..num_columns).map(|_| [fill; BLASTAA_SIZE]).collect();
    let mut values = values.iter().copied();
    if by_row {
        for r in 0..num_rows {
            for row in rows.iter_mut() {
                row[r] = values.next().unwrap_or(fill);
            }
        }
    } else {
        for row in rows.iter_mut() {
            for cell in row.iter_mut().take(num_rows) {
                *cell = values.next().unwrap_or(fill);
            }
        }
    }
    rows
}

/// The PSSM file of `-in_pssm` as `read_in_pssm` reads it: the
/// object after `CPsiBlastImpl`'s constructor, the members the reader skipped, whether the scores
/// were computed from the frequency ratios, and whether NCBI's checks of the PSSM pass
/// (`CPsiBlastValidate::Pssm`; NCBI's constructor throws before any file is written otherwise).
pub(crate) struct InPssmFile {
    pub(crate) object: PssmWithParameters,
    pub(crate) skipped: Vec<&'static str>,
    pub(crate) computed: bool,
    pub(crate) valid: bool,
}

/// The scores and gapped values that `CPsiBlastImpl`'s constructor sets in a PSSM of frequency
/// ratios only (`rows[position][residue]` and `gapped`, the PSSM engine's), with the frequency
/// ratios padded to 28 rows when the file has fewer, and the gap costs of the search (BLOSUM62's
/// 11 and 1: LOSAT's psiblast runs no other, `scores_from_freq_ratios`) in its
/// parameters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:486-502
/// ```c++
///     if (pssm->GetPssm().GetNumRows() !=
///         pssm_with_scores->GetPssm().GetNumRows()) {
///         _ASSERT(pssm_with_scores->GetPssm().GetNumRows() == BLASTAA_SIZE);
///         s_AdjustFrequencyRatiosMatrixToMatchScoreMatrix(*pssm);
///     }
///     pssm->SetPssm().SetFinalData().SetScores() =
///         pssm_with_scores->GetPssm().GetFinalData().GetScores();
///     pssm->SetPssm().SetFinalData().SetLambda() =
///         pssm_with_scores->GetPssm().GetFinalData().GetLambda();
///     pssm->SetPssm().SetFinalData().SetKappa() =
///         pssm_with_scores->GetPssm().GetFinalData().GetKappa();
///     pssm->SetPssm().SetFinalData().SetH() =
///         pssm_with_scores->GetPssm().GetFinalData().GetH();
///
///     PsiBlastAddAncillaryPssmData(*pssm,
///                                   opts.GetGapOpeningCost(),
///                                   opts.GetGapExtensionCost());
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:440-462
/// ```c++
///     const size_t diff = (size_t)BLASTAA_SIZE - pssm.GetPssm().GetNumRows();
///     CPssmIntermediateData::TFreqRatios& freq_ratios =
///         pssm.SetPssm().SetIntermediateData().SetFreqRatios();
///
///     if (pssm.GetPssm().GetByRow() == true) {
///         freq_ratios.resize(pssm.GetPssm().GetNumColumns() * BLASTAA_SIZE, 0.0);
///     } else {
///         CPssmIntermediateData::TFreqRatios::iterator itr = freq_ratios.begin();
///         for (int c = 0; c < pssm.GetPssm().GetNumColumns(); c++) {
///             advance(itr, pssm.GetPssm().GetNumRows());
///             freq_ratios.insert(itr, diff, 0.0);
///         }
///     }
///
///     pssm.SetPssm().SetNumRows() = BLASTAA_SIZE;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:426-433
/// ```c++
/// PsiBlastAddAncillaryPssmData(objects::CPssmWithParameters& pssm,
///                          int gap_open,
///                          int gap_extend)
/// {
///     _ASSERT(pssm.GetParams().GetRpsdbparams().IsSetMatrixName());
///     pssm.SetParams().SetRpsdbparams().SetGapOpen(gap_open);
///     pssm.SetParams().SetRpsdbparams().SetGapExtend(gap_extend);
/// }
/// ```
/// The engine's scores are column by column (`x_PSIMatrix2Asn1`, quoted at
/// `psi_search_and_object`). A PSSM of frequency ratios written by row is LOSAT's refusal
/// (`psi_search_and_object`), so the ratios are padded column by column. A file without `rpsdbparams` gets
/// parameters without a matrix name, which NCBI cannot write: refused with `-out_pssm`
/// (`algorithm::psiblast::in_pssm_files::rounds_input`), and the object is left without them.
pub(crate) fn with_computed_scores(
    pssm: &mut PssmWithParameters,
    rows: &[[i32; BLASTAA_SIZE]],
    gapped: KarlinBlk,
) {
    let matrix = &mut pssm.pssm;
    let rows_n = usize::try_from(matrix.num_rows).unwrap_or(0);
    if rows_n < BLASTAA_SIZE {
        if let Some(data) = matrix.intermediate_data.as_mut() {
            let mut padded = Vec::with_capacity(rows.len() * BLASTAA_SIZE);
            for column in data.freq_ratios.chunks(rows_n.max(1)) {
                padded.extend_from_slice(column);
                padded.resize(padded.len() + BLASTAA_SIZE - rows_n, 0.0);
            }
            data.freq_ratios = padded;
        }
        matrix.num_rows = BLASTAA_SIZE as i32;
    }
    let scores: Vec<i32> = rows.iter().flatten().copied().collect();
    match matrix.final_data.as_mut() {
        Some(data) => {
            data.scores = scores;
            data.lambda = gapped.lambda;
            data.kappa = gapped.k;
            data.h = gapped.h;
        }
        None => {
            matrix.final_data = Some(PssmFinalData {
                scores,
                lambda: gapped.lambda,
                kappa: gapped.k,
                h: gapped.h,
                scaling_factor: None,
                lambda_ungapped: None,
                kappa_ungapped: None,
                h_ungapped: None,
                word_score_threshold: None,
            })
        }
    }
    if let Some(rps) = pssm
        .params
        .as_mut()
        .and_then(|params| params.rpsdbparams.as_mut())
    {
        rps.gap_open = Some(11);
        rps.gap_extend = Some(1);
    }
}

/// The score block of `CPssmEngine::x_InitializeScoreBlock` that the engine reads.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:525-577
// ```c++
//     const EBlastProgramType kProgramType = eBlastTypePsiBlast;
//     ...
//     BlastScoringOptionsSetMatrix(opts, matrix_name);
//     ...
//     status = BlastSetup_ScoreBlkInit(query_blk,
//                                      query_info,
//                                      opts,
//                                      kProgramType,
//                                      &retval,
//                                      kScaleFactor,
//                                      &errors,
//                                      &BlastFindMatrixPath);
//     ...
//     _ASSERT(retval->kbp_ideal);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2843-2848
// ```c
//    stdrfp = Blast_ResFreqNew(sbp);
//    Blast_ResFreqStdComp(sbp, stdrfp);
//    sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
//    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
//    sbp->kbp_ideal = Blast_KarlinBlkNew();
//    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1506-1528 (the matrix's
// `loscore` and `hiscore`: its smallest and largest scores between BLAST_SCORE_MIN and
// BLAST_SCORE_MAX)
// ```c
//          if (matrix[index1][index2] <= BLAST_SCORE_MIN ||
//              matrix[index1][index2] >= BLAST_SCORE_MAX)
//             continue;
//          if (sbp->loscore > matrix[index1][index2])
//             sbp->loscore = matrix[index1][index2];
//          if (sbp->hiscore < matrix[index1][index2])
//             sbp->hiscore = matrix[index1][index2];
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/matrix_freq_ratios.c:1661-1736
// ```c
//     if ( !strcasecmp(matrix_name, "BLOSUM62") ||
//          !strcasecmp(matrix_name, "BLOSUM62_20")) {
//         ...
//         retval->bit_scale_factor = 2;
//     ...
//     } else if ( !strcasecmp(matrix_name, "BLOSUM45") ) {
//         ...
//         retval->bit_scale_factor = 3;
//     } else if ( !strcasecmp(matrix_name, "BLOSUM80") ) {
//         ...
//         retval->bit_scale_factor = 2;
//     ... (BLOSUM50, BLOSUM90, PAM30, PAM70, PAM250: 2)
//     } else {
//         retval = _PSIMatrixFrequencyRatiosFree(retval);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:104-120
// ```c
//         } else {
//             retval =
//                 Blast_KarlinBlkGappedCalc(sbp->kbp_gap_std[index],
//                     scoring_options->gap_open, scoring_options->gap_extend,
//                     sbp->name, error_return);
//         }
//         ...
//         /* For right now, copy the contents from kbp_gap_std to
//          * kbp_gap_psi (as in old code - BLASTSetUpSearchInternalByLoc) */
//         if (program != eBlastTypeBlastn && program != eBlastTypeMapping) {
//             sbp->kbp_gap_psi[index] = Blast_KarlinBlkNew();
//             Blast_KarlinBlkCopy(sbp->kbp_gap_psi[index],
//                                 sbp->kbp_gap_std[index]);
//         }
// ```
// `gapped` is the gapped block of the matrix and gap costs (`Blast_KarlinBlkGappedCalc`), the
// search's: `kbp_gap_std[0]`, copied to `kbp_gap_psi[0]`.
pub(crate) fn score_block(matrix: ScoringMatrix, gapped: KarlinBlk) -> Result<PsiScoreBlk> {
    // psiblast's `Validate` refuses IDENTITY, the one matrix without frequency ratios
    // (`stats::protein_options`, the standard matrices only).
    let data = crate::core::composition_adjustment::adjust_scores::compo_matrix_data(matrix)
        .ok_or_else(|| anyhow::anyhow!("no frequency ratios for {}", matrix.ncbi_name()))?;
    let mut table = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];
    let (mut loscore, mut hiscore) = (i32::MAX, i32::MIN);
    for (i, row) in table.iter_mut().enumerate() {
        for (j, score) in row.iter_mut().enumerate() {
            *score = crate::utils::matrix::protein_score(matrix, i as u8, j as u8);
            if *score <= BLAST_SCORE_MIN || *score >= i16::MAX as i32 {
                continue;
            }
            loscore = loscore.min(*score);
            hiscore = hiscore.max(*score);
        }
    }
    let std_comp = crate::stats::karlin_calc::compute_std_aa_composition();
    let sfp = crate::stats::karlin_calc::compute_score_freq_profile_for_matrix(
        &std_comp, &std_comp, loscore, hiscore, matrix,
    );
    let ideal = crate::stats::karlin_calc::compute_karlin_params_ungapped(&sfp)
        .map_err(|error| crate::blastinput::app::engine_error(&error))?;
    Ok(PsiScoreBlk {
        matrix: table,
        ideal_lambda: ideal.lambda,
        ideal_k: ideal.k,
        gap_std_k: gapped.k,
        gap_psi: gapped,
        freq_ratios: data.freq_ratios,
        bit_scale_factor: if matrix == ScoringMatrix::Blosum45 {
            3
        } else {
            2
        },
        background: data.background,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_text_asn_test_is_ncbis() {
        assert!(text_asn(b"PssmWithParameters ::= {\n  pssm {\n"));
        assert!(text_asn(b"-- a comment\nSeqentry ::= set {\n"));
        // 27 of 34 bytes letters, digits, spaces or braces (the dash, the colons, the blank
        // line's end are not): under 80 %.
        assert!(!text_asn(b"-- a comment\n\nSeq-entry ::= set {\n"));
        assert!(!text_asn(b">tk\nMKTAYIAKQR\n"));
        assert!(!text_asn(b"PssmWithParameters {\n"));
        // Under 80 % letters, digits, spaces and braces.
        assert!(!text_asn(
            b"X ::= ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
        ));
        assert_eq!(guess_format(b""), PssmFileFormat::Unsupported);
        assert_eq!(
            guess_format(b"\x1f\x8b\x08\x00"),
            PssmFileFormat::Unsupported
        );
        assert_eq!(guess_format(b">q\nMKV\n"), PssmFileFormat::Unsupported);
        assert_eq!(guess_format(b"\x30\x80\xa0\x80"), PssmFileFormat::NotText);
        assert_eq!(
            guess_format(b"<?xml version=\"1.0\"?>\n"),
            PssmFileFormat::NotText
        );
    }

    fn fixture(name: &str) -> Vec<u8> {
        std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures")
                .join(name),
        )
        .unwrap()
    }

    fn gunzip(name: &str) -> Vec<u8> {
        use std::io::Read;
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(fixture(name).as_slice())
            .read_to_end(&mut out)
            .unwrap();
        out
    }

    // psiblast's own PSSM (the frozen one of `db.ps.asn.rnd.6`, 30 letters): the query record,
    // the 28 rows of each position and the six statistical values of its file.
    #[test]
    fn a_psiblast_pssm_becomes_a_search() {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
            "db_search/db.ps.asn.rnd.6.asn",
        ))
        .unwrap();
        let scores = pssm.pssm.final_data.as_ref().unwrap().scores.clone();
        let search = psi_search_of(pssm).unwrap();
        assert_eq!(search.query.local_id, "Query_1");
        assert_eq!(search.query.title(), b"rnd random");
        assert_eq!(search.query.seq().len(), 30);
        assert_eq!(search.pssm.rows.len(), 30);
        // byRow FALSE: position 1's 28 scores follow position 0's.
        assert_eq!(&search.pssm.rows[1][..], &scores[28..56]);
        assert!(search.validation_error.is_none());
        let stats = search.stats;
        assert!(stats.lambda_ungapped.is_some() && stats.kappa.is_some() && stats.h.is_some());
        // The gapped Lambda of psiblast's PSSMs is the table's 0.267 as NCBI reads `{ 267, 10, -3 }`.
        assert_eq!(stats.lambda, Some(267.0 * 1e-3));
    }

    // NCBI's unit-test PSSMs: one of frequency ratios only whose query has a gi and a swissprot
    // Seq-id, read (line RG Q1b; NCBI BLAST+ 2.17.0 psiblast: `Query= P01013 Gene X protein
    // (Ovalbumin-related)`, the warning label `129295`, the first id), and one without a query
    // (NCBI's application stops with a C++ exception), refused before the search.
    #[test]
    fn ncbis_unit_test_pssms_are_read_or_refused() {
        let read = |name: &str| {
            crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
        };
        let search = psi_search_of(read("pssm/pssm_freq_ratios.asn.gz")).unwrap();
        for parse_deflines in [false, true] {
            let query = search.query_record(parse_deflines);
            assert_eq!(
                &*query.query_line(),
                b"P01013 Gene X protein (Ovalbumin-related)"
            );
            assert_eq!(query.warning_id(), b"129295");
        }
        let error = psi_search_of(read("pssm/pssm_zero_fratios.asn.gz"))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("not supported by LOSAT's PSIBLAST"),
            "{error}"
        );
    }

    // A PSSM of frequency ratios only (the tk PSSM without its scores) gets NCBI's engine's
    // scores, which are the tk PSSM's own (NCBI 2.17.0 searches both to the same bytes, fixture
    // rows db.ps.inpssm.cbs.ratios.6 and db.ps.inpssm.cbs.tk.6), and the engine's gapped values;
    // its ungapped values stay absent. By row it is refused.
    #[test]
    fn a_pssm_of_frequency_ratios_only_gets_the_engines_scores() {
        let read = |name: &str| {
            crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
        };
        let tk = psi_search_of(read("pssm/ps4_tk_paja.asn.gz")).unwrap();
        let ratios = psi_search_of(read("pssm/ps4_tk_ratios.asn.gz")).unwrap();
        assert_eq!(ratios.pssm, tk.pssm);
        assert_eq!(ratios.freq_ratios, tk.freq_ratios);
        assert_eq!(ratios.stats.lambda, tk.stats.lambda);
        assert_eq!(ratios.stats.h, tk.stats.h);
        let (kappa, file_kappa) = (ratios.stats.kappa.unwrap(), tk.stats.kappa.unwrap());
        assert!((kappa - file_kappa).abs() < 1e-15, "{kappa} {file_kappa}");
        assert_eq!(ratios.stats.lambda_ungapped, None);
        assert_eq!(ratios.stats.kappa_ungapped, None);
        let error = psi_search_of(read("pssm/ps4_rnd_ratios_byrow.asn.gz"))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("a PSSM of frequency ratios only written by row"),
            "{error}"
        );
    }

    // `CPsiBlastValidate::Pssm`'s errors are kept for the search (after the prolog), in NCBI's
    // order; a PSSM whose sizes do not match is LOSAT's refusal.
    #[test]
    fn ncbis_checks_and_the_size_refusal() {
        let read = || {
            crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
                "db_search/db.ps.asn.rnd.6.asn",
            ))
            .unwrap()
        };
        let mut scaled = read();
        scaled.pssm.final_data.as_mut().unwrap().scaling_factor = Some(2);
        assert_eq!(
            psi_search_of(scaled).unwrap().validation_error.as_deref(),
            Some("PSSM has a scaling factor of 2. PSI-BLAST does not accept scaled PSSMs")
        );
        let mut nucleotide = read();
        nucleotide.pssm.is_protein = Some(false);
        assert_eq!(
            psi_search_of(nucleotide)
                .unwrap()
                .validation_error
                .as_deref(),
            Some("PSSM does not represent protein scoring matrix")
        );
        let mut empty = read();
        empty.pssm.final_data = None;
        empty.pssm.intermediate_data = None;
        assert_eq!(
            psi_search_of(empty).unwrap().validation_error.as_deref(),
            Some("PSSM data must contain either scores or frequency ratios")
        );
        let mut short = read();
        short.pssm.num_columns = 29;
        assert!(psi_search_of(short)
            .unwrap_err()
            .to_string()
            .contains("dimensions (28 rows, 29 columns)"));
    }

    fn rnd() -> PssmWithParameters {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
            "db_search/db.ps.asn.rnd.6.asn",
        ))
        .unwrap()
    }

    // The review of PS4a, S1: titles that NCBI's Query= line and its tabular query id read
    // differently (several titles; white space at a title's ends) are refused; one title is read.
    #[test]
    fn query_titles_that_ncbi_reads_two_ways_are_refused() {
        let with_titles = |titles: &[&[u8]]| {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().titles =
                Some(titles.iter().map(|title| title.to_vec()).collect());
            psi_search_of(pssm).map(|search| search.query.title)
        };
        assert_eq!(with_titles(&[b"rnd random"]).unwrap(), b"rnd random");
        assert_eq!(with_titles(&[b""]).unwrap(), b"");
        let error = with_titles(&[b"rnd", b"second"]).unwrap_err().to_string();
        assert!(error.contains("a PSSM query with 2 titles"), "{error}");
        for title in [&b" rnd"[..], b"rnd ", b"  sp1  two  words  "] {
            let error = with_titles(&[title]).unwrap_err().to_string();
            assert!(error.contains("starts or ends with white space"), "{error}");
        }
    }

    // The review of PS5b, B1: the warning label's title is `CBlastQuerySourceBioseqSet::GetTitle`'s
    // (bioseq_extract_data_priv.cpp:345-372, quoted at `psi_search_and_object`): the title without
    // its trailing `.` and spaces, none with a `molinfo` descriptor (NCBI 2.17.0 probes of
    // `PS/tmp/rev5b`, `d.*` and `m.molinfo*`); the record keeps the title for the reports.
    #[test]
    fn the_warning_title_is_the_query_sources() {
        let warning = |title: &[u8], skipped: Vec<&'static str>| {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().titles = Some(vec![title.to_vec()]);
            let (search, _) = psi_search_and_object(pssm, skipped, PssmReader::PSIBLAST).unwrap();
            assert_eq!(search.query.title, title);
            search.warning_title.unwrap()
        };
        assert_eq!(warning(b"rnd random", vec![]), b"rnd random");
        assert_eq!(warning(b"rnd random.", vec![]), b"rnd random");
        assert_eq!(warning(b"rnd. random.. . .", vec![]), b"rnd. random");
        assert_eq!(warning(b".", vec![]), b"");
        assert_eq!(warning(b"rnd random", vec!["Seqdesc"]), b"rnd random");
        assert_eq!(warning(b"rnd random", vec!["Seqdesc", "molinfo"]), b"");
    }

    // The review of PS4a, S2 and N3: the record's local id is the Seq-id's `GetSeqIdString()`
    // (NCBI 2.17.0: `local id 5` shows as `5`, `local id -3` as `-3`, in `qseqid` and the warnings).
    #[test]
    fn the_local_id_is_the_seq_id_string() {
        for (id, shown) in [
            (ObjectId::Str(b"myq".to_vec()), "myq"),
            (ObjectId::Id(5), "5"),
            (ObjectId::Id(-3), "-3"),
        ] {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().ids = vec![ParsedSeqId::Local(id)];
            let search = psi_search_of(pssm).unwrap();
            assert_eq!(search.query.local_id, shown);
            assert!(search.query.believed().is_none());
            assert!(search.query_record(true).believed().is_some());
        }
    }

    // Line RG Q1b (`Q/SURVEY_Q1b.md`, NCBI BLAST+ 2.17.0): a PSSM query whose Seq-id is not local
    // shows its ids with or without `-parse_deflines` (`Query= gi|4504347 ...`, `qseqid`
    // `gi|4504347`, the warnings' label `4504347`), and its Seq-aligns take its first id
    // (`m_gr`: of `gi 4504347, other NP_000549.1` the gi, while `Query=` shows `NP_000549.1`).
    #[test]
    fn a_query_with_other_seq_ids_shows_them() {
        use crate::seqid::model::TextseqId;
        let refseq = ParsedSeqId::Other(TextseqId {
            accession: Some(b"NP_000549".to_vec()),
            version: Some(1),
            ..TextseqId::default()
        });
        for (ids, line, shown, warning) in [
            (
                vec![ParsedSeqId::Gi(4504347)],
                "gi|4504347",
                "4504347",
                "4504347",
            ),
            (
                vec![ParsedSeqId::Gi(4504347), refseq.clone()],
                "NP_000549.1",
                "NP_000549.1",
                "4504347",
            ),
            (
                vec![refseq.clone(), ParsedSeqId::Gi(4504347)],
                "NP_000549.1",
                "NP_000549.1",
                "NP_000549",
            ),
        ] {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().ids = ids.clone();
            let search = psi_search_of(pssm).unwrap();
            for record in [search.query_record(false), search.query_record(true)] {
                let believed = record.believed().expect("the ids are shown");
                assert_eq!(record.local_id, line);
                assert_eq!(record.shown_id(), shown.as_bytes());
                assert_eq!(record.warning_id(), warning.as_bytes());
                assert_eq!(believed.best, ids[0]);
                assert_eq!(believed.ids, ids);
            }
            assert_eq!(search.query_seq_loc_id(), Some(&ids[0]));
        }
    }

    // `Q/SURVEY_Q1b.md` of line RG: a query without a Seq-id, a Geninfo import id (NCBI: `Failed to
    // resolve SeqId: gim|77`), a database ordinal, a local id among other ids (shown as no FASTA
    // query is), an empty local id and a Textseq-id with neither (a byte outside the visible range is
    // the reader's, `a_bad_byte_of_a_pssm_is_noted_as_ncbi_posts_it`)
    // an accession nor a name (NCBI: a database error) are refused.
    #[test]
    fn query_seq_ids_that_are_refused() {
        use crate::seqid::model::{Dbtag, GiimportId};
        for (ids, refusal) in [
            (vec![], "without a Seq-id"),
            (
                vec![ParsedSeqId::Giim(GiimportId {
                    id: 77,
                    db: None,
                    release: None,
                })],
                "is a Geninfo import id",
            ),
            (
                vec![ParsedSeqId::General(Dbtag {
                    db: b"BL_ORD_ID".to_vec(),
                    tag: ObjectId::Id(5),
                })],
                "names a database ordinal",
            ),
            (
                vec![
                    ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec())),
                    ParsedSeqId::Gi(4504347),
                ],
                "a local Seq-id among other Seq-ids",
            ),
            (
                vec![ParsedSeqId::Local(ObjectId::Str(Vec::new()))],
                "an empty local id",
            ),
            (
                vec![ParsedSeqId::Other(crate::seqid::model::TextseqId {
                    release: Some(b"x".to_vec()),
                    ..Default::default()
                })],
                "neither an accession nor a name",
            ),
        ] {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().ids = ids;
            let error = psi_search_of(pssm).unwrap_err().to_string();
            assert!(error.contains(refusal), "{error}");
        }
    }

    // The review of PS4a, N8: a score outside one letter pair's range (`Int2`) is refused;
    // `BLAST_SCORE_MIN` and `BLAST_SCORE_MAX` are scores.
    #[test]
    fn scores_outside_int2_are_refused() {
        let with_score = |score: i32| {
            let mut pssm = rnd();
            pssm.pssm.final_data.as_mut().unwrap().scores[5] = score;
            psi_search_of(pssm)
        };
        for score in [BLAST_SCORE_MIN, BLAST_SCORE_MAX] {
            assert_eq!(with_score(score).unwrap().pssm.rows[0][5], score);
        }
        for score in [BLAST_SCORE_MIN - 1, BLAST_SCORE_MAX + 1, i32::MAX] {
            let error = with_score(score).unwrap_err().to_string();
            assert!(
                error.contains(&format!("a PSSM score of {score}, outside")),
                "{error}"
            );
        }
        // A PSSM that NCBI's check refuses keeps NCBI's error (after the prolog).
        let mut scaled = rnd();
        let data = scaled.pssm.final_data.as_mut().unwrap();
        data.scaling_factor = Some(2);
        data.scores[5] = i32::MAX;
        assert!(psi_search_of(scaled).unwrap().validation_error.is_some());
    }

    // The invalid residues of the PSSM's query (review F2 of M4a-1b; NCBI 2.17.0 probes s19-s23 and
    // the `db.ps.inpssm.*` / `db.tbn.inpssm.*` rows): each coding converts with its own table, the
    // positions are 0-based, and the search stops with the engine error of the query set-up.
    #[test]
    fn residues_that_the_codings_do_not_have_stop_the_search() {
        let with_data = |data: SeqData| {
            let mut pssm = rnd();
            pssm.pssm.query.as_mut().unwrap().inst.seq_data = Some(data);
            psi_search_of(pssm).unwrap()
        };
        let text = |bytes: &[u8]| bytes.to_vec();
        let ok = b"ACDEFGHIKLMNPQRSTVWYACDEFGMAAA";
        let at26 = |byte: u8| {
            let mut letters = ok.to_vec();
            letters[26] = byte;
            letters
        };
        // `*` is no IUPACaa letter, but is NCBIeaa's (25) and NCBIstdaa's; `-` is NCBIeaa's gap
        // (a valid residue to NCBI, refused by LOSAT); `#` is none; lower case is IUPACaa's only.
        let error_of = |search: PsiSearch| search.validation_error.unwrap();
        let positions_26 = "Warning: Invalid residues found at positions 26 ";
        for bad in [b'*', b'#', b'\xE9', b'-'] {
            assert_eq!(
                error_of(with_data(SeqData::Iupacaa(at26(bad)))),
                positions_26,
                "iupacaa {bad:#x}"
            );
        }
        for bad in [b'#', b'\xE9', b'a'] {
            assert_eq!(
                error_of(with_data(SeqData::Ncbieaa(at26(bad)))),
                positions_26,
                "ncbieaa {bad:#x}"
            );
        }
        assert!(with_data(SeqData::Ncbieaa(at26(b'*')))
            .validation_error
            .is_none());
        // Lower case converts as capitals (the search would run) but is refused, below.
        let mut lower = rnd();
        lower.pssm.query.as_mut().unwrap().inst.seq_data =
            Some(SeqData::Iupacaa(text(&ok.to_ascii_lowercase())));
        let error = psi_search_of(lower).unwrap_err().to_string();
        assert!(error.contains("lower-case iupacaa"), "{error}");
        // NCBIstdaa: codes below 28 are letters, anything else is invalid.
        let mut codes: Vec<u8> = ok
            .iter()
            .map(|&c| NCBISTDAA_TO_AMINOACID.iter().position(|&l| l == c).unwrap() as u8)
            .collect();
        assert!(with_data(SeqData::Ncbistdaa(codes.clone()))
            .validation_error
            .is_none());
        for code in [28, 30, 255] {
            codes[26] = code;
            assert_eq!(
                error_of(with_data(SeqData::Ncbistdaa(codes.clone()))),
                positions_26
            );
        }
        // NCBI lists the first 20 and says so.
        let many = with_data(SeqData::Ncbieaa(vec![b'#'; 30]));
        assert_eq!(
            error_of(many),
            "Warning: Invalid residues found at positions 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,... (only first 20 shown) "
        );
    }

    #[test]
    fn scores_fill_28_rows_by_column_or_by_row() {
        // Two rows (A, B) of three positions, column by column.
        let rows = convert_2_matrix(&[1, 2, 3, 4, 5, 6], 2, 3, false, BLAST_SCORE_MIN);
        assert_eq!(rows.len(), 3);
        assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
        assert_eq!(rows[1][2], BLAST_SCORE_MIN);
        // The same values row by row.
        let rows = convert_2_matrix(&[1, 3, 5, 2, 4, 6], 2, 3, true, BLAST_SCORE_MIN);
        assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
        assert_eq!(rows[2][27], BLAST_SCORE_MIN);
    }
}
