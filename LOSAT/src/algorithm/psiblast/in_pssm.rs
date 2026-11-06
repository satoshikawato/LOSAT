//! psiblast `-in_pssm F` (series PS, item PS4a): the PSSM file read as NCBI's `CPsiBlastArgs`
//! reads it, NCBI's errors for files it cannot read, NCBI's checks of the PSSM
//! (`CPsiBlastValidate::Pssm`), the query record made from the PSSM's Bioseq, the score rows
//! (`CScorematPssmConverter::GetScores`) and the statistical values of `PsiBlastSetupScoreBlock`,
//! and LOSAT's refusals of what the PSSM search does not run yet (DECISIONS 11 of line PS).
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
use std::sync::Arc;

use anyhow::{bail, Result};

use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::pssm_search::{PsiSearch, PssmStats};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::blast_util::{BLASTAA_SIZE, BLAST_SCORE_MIN};
use crate::report::serial::pssm::{ObjectId, PssmSeqId, PssmWithParameters, SeqData};

/// NCBISTDAA's letters (`NCBISTDAA_TO_AMINOACID`, blast_encoding.c:115-118).
const NCBISTDAA_LETTERS: [u8; BLASTAA_SIZE] = *b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

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

/// The refusal of the PSSM search's limits ("... is not supported by LOSAT's PSIBLAST").
fn refuse(what: impl std::fmt::Display) -> anyhow::Error {
    anyhow::anyhow!("{what} is not supported by LOSAT's PSIBLAST")
}

/// Reads `-in_pssm path`: NCBI's errors for a file it cannot read or guess, LOSAT's refusals of
/// what it does not read, then the PSSM search's query, rows and values.
pub(crate) fn read_in_pssm(path: &Path) -> Result<PsiSearch> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
    // ```c++
    // void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
    // {
    //     if ( !m_Ios ) {
    //         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
    //             "File is not accessible",AsString()));
    // ```
    let bytes = std::fs::read(path).map_err(|_| crate::cli::inaccessible("in_pssm", path))?;
    match guess_format(&bytes) {
        PssmFileFormat::TextAsn => {}
        PssmFileFormat::Unsupported => {
            return Err(crate::blastinput::app::options_error(
                "Unsupported format for PSSM",
            ))
        }
        PssmFileFormat::NotText => {
            return Err(refuse(format!(
                "a PSSM file that is not ASN.1 text (binary ASN.1, XML or another format; {})",
                path.display()
            )))
        }
    }
    let pssm = match crate::blastinput::pssm_text::read_pssm_with_parameters(&bytes) {
        Ok(pssm) => pssm,
        // The reader's own limits (a query Seq-entry set, Seq-ids and Seq-data it does not type,
        // REAL numbers NCBI computes with undefined behaviour or in `long double`): LOSAT's
        // refusal, where NCBI reads the file.
        Err(error) if error.message.contains("is not supported by LOSAT") => {
            return Err(refuse(format!(
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
    psi_search_of(pssm)
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

/// A string that NCBI's text reader would change: a byte outside the visible range
/// (0x20..0x7E), which NCBI replaces by `#` with an error message (`FixVisibleChars`).
fn visible(bytes: &[u8]) -> bool {
    bytes.iter().all(|&c| (0x20..=0x7e).contains(&c))
}

/// The PSSM search of `pssm`: its query record, rows and values, or LOSAT's refusal.
pub(crate) fn psi_search_of(pssm: PssmWithParameters) -> Result<PsiSearch> {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:534-539
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
        return Err(refuse(
            "a PSSM without a query sequence (NCBI BLAST+'s psiblast stops with a C++ exception)",
        ));
    };
    let validation_error = validation_error(&pssm);
    // A PSSM of frequency ratios only: NCBI computes its scores with the PSSM engine
    // (`x_CreatePssmScoresFromFrequencyRatios`), which item PS4b ports.
    let scores = match pssm.pssm.final_data.as_ref() {
        Some(data) if !data.scores.is_empty() => data.scores.as_slice(),
        _ if validation_error.is_some() => &[],
        _ => {
            return Err(refuse(
                "a PSSM without scores (frequency ratios only, from which NCBI computes the scores) yet",
            ))
        }
    };
    // The query: one local string Seq-id (`Query_N` of a psiblast PSSM), its first title, its
    // letters.
    let local_id = match query.ids.as_slice() {
        [PssmSeqId::Local(ObjectId::Str(id))] if visible(id) && !id.is_empty() => {
            String::from_utf8_lossy(id).into_owned()
        }
        _ => {
            return Err(refuse(
                "a PSSM query whose Seq-id is not one local string id (an accession, gi or other id is shown through the parsed Seq-ids of the query)",
            ))
        }
    };
    let title: &[u8] = query
        .titles
        .as_ref()
        .and_then(|titles| titles.first())
        .map_or(&[], |title| title.as_slice());
    if !visible(title)
        || query.titles.iter().flatten().any(|title| !visible(title))
        || pssm
            .pssm
            .row_labels
            .iter()
            .flatten()
            .any(|label| !visible(label))
    {
        return Err(refuse(
            "a PSSM whose strings hold a byte outside the visible ASCII range (NCBI reads it as '#' with an error message)",
        ));
    }
    if query.inst.mol.1 != 3 || query.inst.repr.1 != 2 {
        return Err(refuse(format!(
            "a PSSM query whose Seq-inst is not a raw protein (repr {}, mol {})",
            query.inst.repr.0, query.inst.mol.0
        )));
    }
    let letters: Vec<u8> = match query.inst.seq_data.as_ref() {
        Some(SeqData::Ncbistdaa(codes)) => codes
            .iter()
            .map(|&code| NCBISTDAA_LETTERS.get(code as usize).copied())
            .collect::<Option<Vec<u8>>>()
            .unwrap_or_default(),
        Some(SeqData::Iupacaa(text)) | Some(SeqData::Ncbieaa(text)) => text
            .iter()
            .map(|&c| NCBISTDAA_LETTERS.contains(&c).then_some(c))
            .collect::<Option<Vec<u8>>>()
            .unwrap_or_default(),
        None => Vec::new(),
    };
    let length = query.inst.length.unwrap_or(-1);
    if letters.is_empty() || letters.len() as i64 != i64::from(length) || letters.contains(&b'-') {
        return Err(refuse(
            "a PSSM query whose letters are missing, of another length than its Seq-inst, or with a gap or a byte outside NCBIstdaa",
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
            || scores.len() as i64 != cells
            || (ratios_len != 0 && ratios_len as i64 != cells))
    {
        return Err(refuse(format!(
            "a PSSM whose dimensions ({rows_n} rows, {columns_n} columns) do not match its scores ({}), frequency ratios ({ratios_len}) or query ({} letters)",
            scores.len(),
            letters.len()
        )));
    }
    let rows = if validation_error.is_none() {
        convert_scores(
            scores,
            rows_n as usize,
            columns_n as usize,
            pssm.pssm.by_row == Some(true),
        )
    } else {
        Vec::new()
    };
    let final_data = pssm.pssm.final_data.as_ref();
    let stat = |value: Option<f64>| value.filter(|&value| value != 0.0);
    let stats = PssmStats {
        lambda_ungapped: stat(final_data.and_then(|data| data.lambda_ungapped)),
        kappa_ungapped: stat(final_data.and_then(|data| data.kappa_ungapped)),
        h_ungapped: stat(final_data.and_then(|data| data.h_ungapped)),
        lambda: stat(final_data.map(|data| data.lambda)),
        kappa: stat(final_data.map(|data| data.kappa)),
        h: stat(final_data.map(|data| data.h)),
    };
    Ok(PsiSearch {
        pssm: QueryPssm { rows },
        stats,
        query: FastaRecord::new(local_id, title, &letters),
        validation_error,
    })
}

/// `GetScores` + the copy of `PsiBlastSetupScoreBlock`: `rows[position][residue]` of 28 cells,
/// those past the file's rows `BLAST_SCORE_MIN`.
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
fn convert_scores(
    scores: &[i32],
    num_rows: usize,
    num_columns: usize,
    by_row: bool,
) -> Vec<[i32; BLASTAA_SIZE]> {
    let mut rows: Vec<[i32; BLASTAA_SIZE]> = (0..num_columns)
        .map(|_| [BLAST_SCORE_MIN; BLASTAA_SIZE])
        .collect();
    let mut values = scores.iter().copied();
    if by_row {
        for r in 0..num_rows {
            for row in rows.iter_mut() {
                row[r] = values.next().unwrap_or(BLAST_SCORE_MIN);
            }
        }
    } else {
        for row in rows.iter_mut() {
            for cell in row.iter_mut().take(num_rows) {
                *cell = values.next().unwrap_or(BLAST_SCORE_MIN);
            }
        }
    }
    rows
}

/// The refusal of psiblast's PSSM file options with `-in_pssm` (item PS5), made when the PSSM has
/// been read.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248
/// ```c++
///             while (itr) {
///
///                 SavePssmToFile(pssm, &itr);
/// ```
/// With `-in_pssm` NCBI writes the input PSSM before the first round (`-out_pssm` alone too,
/// `-save_each_pssm` as `F.1`), and the rounds' PSSMs after it.
pub(crate) fn check_pssm_file_options(
    options: &crate::algorithm::psiblast::pssm::PsiblastPssmOptions,
) -> Result<()> {
    if options.out_pssm.is_some()
        || options.out_ascii_pssm.is_some()
        || options.save_each_pssm
        || options.save_pssm_after_last_round
    {
        bail!(
            "the PSSM files (-out_pssm, -out_ascii_pssm, -save_each_pssm, -save_pssm_after_last_round) with -in_pssm are not supported by LOSAT's PSIBLAST yet"
        );
    }
    Ok(())
}

/// LOSAT's refusals of a psiblast `-in_pssm` search that the PSSM search does not run yet
/// (DECISIONS 11 of line PS), checked when the options have been processed, before the search.
pub(crate) fn check_in_pssm_options(args: &ResolvedBlastpArgs, psi: &Arc<PsiSearch>) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1497
    // ```c
    //     if (ext_params->options->compositionBasedStats > 0 ||
    //         ext_params->options->eTbackExt == eSmithWatermanTbck) {
    //         status =
    //             Blast_RedoAlignmentCore_MT(program_number, ...);
    // ```
    // The position-based redo (modes 1, 2 and 3, which NCBI runs as 1) is item PS4b's.
    if args.comp_based_stats.mode != BlastpCompositionMode::NoCompositionBasedStats {
        bail!(
            "-comp_based_stats {} with -in_pssm (the position-based composition-based statistics) is not supported by LOSAT's PSIBLAST yet",
            args.comp_based_stats.to_ncbi_cli_string()
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:191-198 (the SEG mask
    // of the query, then `PsiBlastSetupScoreBlock`)
    // The query's SEG mask changes the PSSM search's hits (DESIGN_ps4 §1.6).
    if args.seg.params().is_some() || args.soft_masking {
        bail!("-seg and -soft_masking true with -in_pssm are not supported by LOSAT's PSIBLAST");
    }
    // The outfmt 0, 6 and 7 reports (the others are item PS15b's).
    let choice = crate::blastinput::app::parse_formatting_string(&args.outfmt)?;
    if !matches!(choice.number, 0 | 6 | 7) {
        bail!(
            "output format {} with -in_pssm is not supported by LOSAT's PSIBLAST yet",
            choice.number
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:91-94
    // ```c++
    //     else if (Blast_ProgramIsPsiBlast(program))
    //         retval = 10000;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196 (the
    // chunks of a split query are set up without the PSSM)
    // A PSSM query that NCBI splits is searched with BLOSUM62 in its preliminary stage.
    let sizes = crate::blastinput::app::protein_query_split_sizes(args.app.name(), 10000)?;
    if crate::algorithm::common::protein_query_split::split_protein_batch(
        &[psi.query.seq().len()],
        sizes,
    )
    .is_some()
    {
        bail!(
            "a PSSM query of {} letters, which NCBI BLAST+ splits into chunks searched without the PSSM, is not supported by LOSAT's PSIBLAST",
            psi.query.seq().len()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_text_asn_test_is_ncbis() {
        assert!(text_asn(b"PssmWithParameters ::= {\n  pssm {\n"));
        assert!(text_asn(b"-- a comment\n\nSeq-entry ::= set {\n"));
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

    // NCBI's unit-test PSSMs: one of frequency ratios only (PS4b) and one without a query (NCBI's
    // application stops with a C++ exception); both refused before the search.
    #[test]
    fn ncbis_unit_test_pssms_are_refused() {
        let read = |name: &str| {
            crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
        };
        let error = psi_search_of(read("pssm/pssm_freq_ratios.asn.gz"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("a PSSM without scores"), "{error}");
        let error = psi_search_of(read("pssm/pssm_zero_fratios.asn.gz"))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("not supported by LOSAT's PSIBLAST"),
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

    #[test]
    fn scores_fill_28_rows_by_column_or_by_row() {
        // Two rows (A, B) of three positions, column by column.
        let rows = convert_scores(&[1, 2, 3, 4, 5, 6], 2, 3, false);
        assert_eq!(rows.len(), 3);
        assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
        assert_eq!(rows[1][2], BLAST_SCORE_MIN);
        // The same values row by row.
        let rows = convert_scores(&[1, 3, 5, 2, 4, 6], 2, 3, true);
        assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
        assert_eq!(rows[2][27], BLAST_SCORE_MIN);
    }
}
