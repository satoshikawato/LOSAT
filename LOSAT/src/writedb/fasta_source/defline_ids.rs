//! The first parse of the ids of a definition line with `-parse_seqids`: the id branch of
//! `CFastaDeflineReader::ParseDefline` with makeblastdb's flags (`fRequireID |
//! fParseRawID | fQuickIDCheck`, `fAssumeProt` or `fAssumeNuc`): the range suffix
//! (`ParseRange`), the ids (`x_ProcessIDs`, with `CSeq_id::ParseIDs`) and their checks
//! (`CSeqIdCheck`, `CFastaIdValidate`). Errors end the reading (`CFastaBioseqSource` gives
//! the reader no listener, so `s_PostError` throws); warnings are posted as `LOG_POST`
//! lines (`c++/src/objtools/readers/fasta_reader_utils.cpp`, `seqid_validate.cpp`).

use super::reader::ParseError;
use super::Warn;
use crate::seqid::accguide::{parse_flags, SeqIdWarning};
use crate::seqid::label::content;
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::seqid::parse::{
    has_invalid_local_char, is_blank, parse_ids, MAX_ACCESSION_LENGTH, MAX_GENERAL_TAG_LENGTH,
    MAX_LOCAL_ID_LENGTH,
};

/// A diagnostic of NCBI's Seq-id code as makeblastdb prints it (`ERR_POST_X` of
/// Seq_id.cpp: `NCBI_ERRCODE_X` 1306 and the subcode; the application's name after the
/// severity).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2676-2679
/// ```c++
///                 if ((flags & fParse_PartialOK) != 0) {
///                     ERR_POST_X(7, Warning << e.what());
/// ```
/// makeblastdb 2.17.0 printed `Warning: (1306.7) [makeblastdb] NCBI C++ Exception:` and the
/// report for `>xyz|ABC|DEF` (survey B fixture p16_unsupported_type).
pub(crate) fn post_seq_id_warning(warning: SeqIdWarning, warn: Warn<'_>) {
    warn(
        format!(
            "Warning: (1306.{}) [makeblastdb] {}",
            warning.subcode, warning.message
        )
        .as_bytes(),
    );
}

/// `CFastaDeflineReader::ParseRange(s, start, end)`: the length of the range suffix of
/// the id token (`:<from>-<to>` or `:c<from>-<to>`), 0 for none, with its ends.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:229-273
/// ```c++
/// TSeqPos CFastaDeflineReader::ParseRange(
///     const CTempString& s,
///     TSeqPos& start,
///     TSeqPos& end,
///     ILineErrorListener * pMessageListener)
/// {
///
///     if (s.empty()) {
///         return 0;
///     }
///
///     bool    on_start = false;
///     bool    negative = false;
///     TSeqPos mult = 1;
///     size_t  pos;
///     start = end = 0;
///     for (pos = s.length() - 1;  pos > 0;  --pos) {
///         unsigned char c = s[pos];
///         if (c >= '0'  &&  c <= '9') {
///             if (on_start) {
///                 start += (c - '0') * mult;
///             } else {
///                 end += (c - '0') * mult;
///             }
///             mult *= 10;
///         } else if (c == '-'  &&  !on_start  &&  mult > 1) {
///             on_start = true;
///             mult = 1;
///         } else if (c == ':'  &&  on_start  &&  mult > 1) {
///             break;
///         } else if (c == 'c'  &&  pos > 0  &&  s[--pos] == ':'
///                    &&  on_start  &&  mult > 1) {
///             negative = true;
///             break;
///         } else {
///             return 0; // syntax error
///         }
///     }
///     if ((negative ? (end > start) : (start > end))  ||  s[pos] != ':') {
///         return 0;
///     }
///     --start;
///     --end;
///     return TSeqPos(s.length() - pos);
/// }
/// ```
/// `TSeqPos` is a 32-bit unsigned integer: the sums and `mult` wrap as NCBI's do.
pub(super) fn parse_range(s: &[u8]) -> (usize, u32, u32) {
    if s.is_empty() {
        return (0, 0, 0);
    }
    let mut on_start = false;
    let mut negative = false;
    let mut mult: u32 = 1;
    let (mut start, mut end) = (0u32, 0u32);
    let mut pos = s.len() - 1;
    while pos > 0 {
        let c = s[pos];
        if c.is_ascii_digit() {
            let value = u32::from(c - b'0').wrapping_mul(mult);
            if on_start {
                start = start.wrapping_add(value);
            } else {
                end = end.wrapping_add(value);
            }
            mult = mult.wrapping_mul(10);
        } else if c == b'-' && !on_start && mult > 1 {
            on_start = true;
            mult = 1;
        } else if c == b':' && on_start && mult > 1 {
            break;
        } else if c == b'c'
            && pos > 0
            && {
                pos -= 1;
                s[pos] == b':'
            }
            && on_start
            && mult > 1
        {
            negative = true;
            break;
        } else {
            return (0, 0, 0);
        }
        pos -= 1;
    }
    if (if negative { end > start } else { start > end }) || s[pos] != b':' {
        return (0, 0, 0);
    }
    (s.len() - pos, start.wrapping_sub(1), end.wrapping_sub(1))
}

/// `GetSeqIdString()` (without the version) of an id; empty for an id left unset.
fn seq_id_string(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    content(id, false, &mut out);
    out
}

/// The message of a `CObjReaderParseException` that `s_PostError` throws.
fn error(message: String) -> ParseError {
    ParseError { message }
}

/// `CFastaDeflineReader::x_ProcessIDs(id_string, info, ids, listener, CSeqIdCheck)` with
/// makeblastdb's flags: the ids of the token, checked. `protein` is `fAssumeProt` (else
/// `fAssumeNuc`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:333-404
/// ```c++
/// void CFastaDeflineReader::x_ProcessIDs(
///     const CTempString& id_string,
///     const SDeflineParseInfo& info,
///     TIds& ids,
///     ILineErrorListener* pMessageListener,
///     FIdCheck f_id_check
///     )
/// {
///     if (info.fBaseFlags & CReaderBase::fAllIdsAsLocal)
///     {
///         ...
///     }
///
///     CSeq_id::TParseFlags flags =
///         CSeq_id::fParse_PartialOK |
///         CSeq_id::fParse_AnyLocal;
///
///     if (info.fFastaFlags & CFastaReader::fParseRawID) {
///         flags |= CSeq_id::fParse_RawText;
///     }
///
///     string local_copy;
///     auto to_parse = id_string;
///     if (id_string.find(',') != NPOS &&
///         id_string.find('|') == NPOS) {
///         const string err_message =
///             "Near line " + NStr::NumericToString(info.lineNumber)
///             + ", the sequence id string contains 'comma' symbol, which has been replaced with 'underscore' "
///             + "symbol. Please correct the sequence id string.";
///
///         s_PostWarning(pMessageListener,
///             info.lineNumber,
///             id_string,
///             err_message,
///             ILineError::eProblem_GeneralParsingError,
///             CObjReaderParseException::eFormat);
///
///         local_copy = id_string;
///         for (auto& rit : local_copy)
///             if (rit == ',')
///                 rit = '_';
///
///         to_parse = local_copy;
///     }
///
///     try {
///         CSeq_id::ParseIDs(ids, to_parse, flags);
///         ids.remove_if([](CRef<CSeq_id> id_ref)
///                 { return NStr::IsBlank(id_ref->GetSeqIdString()); });
///     }
///     catch(...) {
///         ids.clear();
///     }
///
///     if (ids.empty()) {
///         s_PostError(pMessageListener,
///                 info.lineNumber,
///                 id_string,
///                 "Could not construct seq-id from '" + id_string + "'",
///                 ILineError::eProblem_GeneralParsingError,
///                 CObjReaderParseException::eNoIDs);
///
///         ids.push_back(Ref(new CSeq_id(CSeq_id::e_Local, id_string)));
///         return;
///     }
///     // Convert anything that looks like a GI to a local id
///     if ( info.fBaseFlags & CReaderBase::fNumericIdsAsLocal ) {
///         x_ConvertNumericToLocal(ids);
///     }
///
///     f_id_check(ids, info, pMessageListener);
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:85,101-104
/// ```c++
/// static void s_PostWarning(ILineErrorListener* pMessageListener,
///     ...
///     if (!pMessageListener) {
///         LOG_POST_X(1, Warning << pLineExpt->Message());
///         return;
///     }
/// ```
/// Neither `fAllIdsAsLocal` nor `fNumericIdsAsLocal` is set (base flags 0), and the
/// registry's `MAX_SEQID_LENGTH` is not read (makedb-v5 DECISIONS 12), so the default
/// limits apply.
pub(super) fn process_ids(
    id_string: &[u8],
    protein: bool,
    line_number: u64,
    warn: Warn<'_>,
) -> Result<Vec<ParsedSeqId>, ParseError> {
    let flags = parse_flags::PARTIAL_OK | parse_flags::ANY_LOCAL | parse_flags::RAW_TEXT;
    let replaced;
    let mut to_parse = id_string;
    if id_string.contains(&b',') && !id_string.contains(&b'|') {
        warn(
            format!(
                "Near line {line_number}, the sequence id string contains 'comma' symbol, which \
                 has been replaced with 'underscore' symbol. Please correct the sequence id \
                 string."
            )
            .as_bytes(),
        );
        replaced = id_string
            .iter()
            .map(|&c| if c == b',' { b'_' } else { c })
            .collect::<Vec<u8>>();
        to_parse = &replaced;
    }
    let ids: Vec<ParsedSeqId> =
        match parse_ids(to_parse, flags, &mut |w| post_seq_id_warning(w, &mut *warn)) {
            Ok(ids) => ids
                .into_iter()
                .flatten()
                .filter(|id| !is_blank(&seq_id_string(id)))
                .collect(),
            Err(_) => Vec::new(),
        };
    if ids.is_empty() {
        return Err(error(format!(
            "Could not construct seq-id from '{}'",
            String::from_utf8_lossy(id_string)
        )));
    }
    validate_ids(&ids, protein, line_number, warn)?;
    Ok(ids)
}

/// `CFastaIdValidate(flags)(ids, lineNum, CIdErrorReporter(nullptr))` with `fQuickIDCheck`
/// and `fAssumeProt` or `fAssumeNuc`: an error throws, a warning is a `LOG_POST` line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:107-133
/// ```c++
/// void CFastaIdValidate::operator()(const TIds& ids,
///         int lineNum,
///         FReportError fReportError)
/// {
///     if (ids.empty()) {
///         return;
///     }
///
///     if (!(m_Flags&CFastaReader::fAssumeProt)) {
///         CheckForExcessiveNucData(*ids.back(), lineNum, fReportError);
///     }
///
///     if (!(m_Flags&CFastaReader::fAssumeNuc)) {
///         CheckForExcessiveProtData(*ids.back(), lineNum, fReportError);
///     }
///
///
///     for (const auto& pId : ids) {
///         if (pId->IsLocal() &&
///             !IsValidLocalID(*pId)) {
///             const auto& idString = pId->GetSeqIdString();
///             string msg = "'" + idString + "' is not a valid local ID";
///             fReportError(eDiag_Error, lineNum, idString, CFastaIdValidate::eBadLocalID, msg);
///         }
///         CheckIDLength(*pId, lineNum, fReportError);
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:215-237
/// ```c++
/// bool CFastaIdValidate::IsValidLocalID(const CSeq_id& id) const
/// {
///     if (id.IsLocal()) {
///         if (id.GetLocal().IsId()) {
///             return true;
///         }
///         if (id.GetLocal().IsStr()) {
///             return IsValidLocalString(id.GetLocal().GetStr());
///         }
///     }
///     return false;
/// }
///
///
/// bool CFastaIdValidate::IsValidLocalString(const CTempString& idString) const
/// {
///     const CTempString& checkString =
///         (m_Flags & CFastaReader::fQuickIDCheck) ?
///         idString.substr(0,1) :
///         idString;
///
///     return !(CSeq_id::CheckLocalID(checkString)&CSeq_id::fInvalidChar);
/// }
/// ```
/// The ids come from the first-parse token, which ends at the first byte of 0x20 or below
/// or of 0x80 or more (decision 28), so the first byte is never one that
/// `CSeq_id_find_pred` rejects for being a control or non-ASCII byte; only ` >[]|"` can
/// decide here.
fn validate_ids(
    ids: &[ParsedSeqId],
    protein: bool,
    line_number: u64,
    warn: Warn<'_>,
) -> Result<(), ParseError> {
    let Some(last) = ids.last() else {
        return Ok(());
    };
    if protein {
        check_for_excessive_prot_data(last, warn);
    } else {
        check_for_excessive_nuc_data(last, warn)?;
    }
    for id in ids {
        if let ParsedSeqId::Local(ObjectId::Str(text)) = id {
            if has_invalid_local_char(&text[..text.len().min(1)]) {
                return Err(error(format!(
                    "'{}' is not a valid local ID",
                    String::from_utf8_lossy(text)
                )));
            }
        }
        check_id_length(id, line_number)?;
    }
    Ok(())
}

/// `s_GetIDLengthErrorString(length, idType, maxAllowedLength, lineNum)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:154-166
/// ```c++
/// static string s_GetIDLengthErrorString(size_t length,
///         const string& idType,
///         size_t maxAllowedLength,
///         int lineNum)
/// {
///     string err_message =
///         "Near line " + NStr::NumericToString(lineNum) +
///         + ", the " + idType + " is too long.  Its length is " + NStr::NumericToString(length)
///         + " but the maximum allowed " + idType + " length is "+  NStr::NumericToString(maxAllowedLength)
///         + ".  Please find and correct all " + idType + "s that are too long.";
///
///     return err_message;
/// }
/// ```
fn id_length_error(length: usize, id_type: &str, max: usize, line_number: u64) -> ParseError {
    error(format!(
        "Near line {line_number}, the {id_type} is too long.  Its length is {length} but the \
         maximum allowed {id_type} length is {max}.  Please find and correct all {id_type}s that \
         are too long."
    ))
}

/// `CFastaIdValidate::CheckIDLength(id, lineNum, fReportError)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:169-213
/// ```c++
/// void CFastaIdValidate::CheckIDLength(const CSeq_id& id, int lineNum, FReportError fReportError) const
/// {
///     if (id.IsLocal()) {
///         if (id.GetLocal().IsStr() &&
///             id.GetLocal().GetStr().length() > kMaxLocalIDLength) {
///             const auto& msg =
///                 s_GetIDLengthErrorString(id.GetLocal().GetStr().length(),
///                         "local id",
///                         kMaxLocalIDLength,
///                         lineNum);
///             fReportError(eDiag_Error, lineNum, id.GetSeqIdString(), CFastaIdValidate::eIDTooLong, msg);
///         }
///         return;
///     }
///
///     if (id.IsGeneral()) {
///         if (id.GetGeneral().IsSetTag() &&
///             id.GetGeneral().GetTag().IsStr()) {
///             const auto length = id.GetGeneral().GetTag().GetStr().length();
///             if (length > kMaxGeneralTagLength) {
///                 const auto& msg =
///                     s_GetIDLengthErrorString(id.GetGeneral().GetTag().GetStr().length(),
///                             "general id string",
///                             kMaxGeneralTagLength,
///                             lineNum);
///                 fReportError(eDiag_Error, lineNum, id.GetSeqIdString(), CFastaIdValidate::eIDTooLong, msg);
///             }
///         }
///         return;
///     }
///
///     auto pTextId = id.GetTextseq_Id();
///     if (pTextId && pTextId->IsSetAccession()) {
///         const auto length = pTextId->GetAccession().length();
///         if (length > kMaxAccessionLength) {
///             const auto& msg =
///                 s_GetIDLengthErrorString(length,
///                                 "accession",
///                                 kMaxAccessionLength,
///                                 lineNum);
///             fReportError(eDiag_Error, lineNum, id.GetSeqIdString(), CFastaIdValidate::eIDTooLong, msg);
///         }
///     }
/// }
/// ```
fn check_id_length(id: &ParsedSeqId, line_number: u64) -> Result<(), ParseError> {
    match id {
        ParsedSeqId::Local(tag) => {
            if let ObjectId::Str(text) = tag {
                if text.len() > MAX_LOCAL_ID_LENGTH {
                    return Err(id_length_error(
                        text.len(),
                        "local id",
                        MAX_LOCAL_ID_LENGTH,
                        line_number,
                    ));
                }
            }
        }
        ParsedSeqId::General(dbtag) => {
            if let ObjectId::Str(text) = &dbtag.tag {
                if text.len() > MAX_GENERAL_TAG_LENGTH {
                    return Err(id_length_error(
                        text.len(),
                        "general id string",
                        MAX_GENERAL_TAG_LENGTH,
                        line_number,
                    ));
                }
            }
        }
        _ => {
            if let Some(accession) = id.textseq().and_then(|text| text.accession.as_ref()) {
                if accession.len() > MAX_ACCESSION_LENGTH {
                    return Err(id_length_error(
                        accession.len(),
                        "accession",
                        MAX_ACCESSION_LENGTH,
                        line_number,
                    ));
                }
            }
        }
    }
    Ok(())
}

/// `kWarnNumNucCharsAtEnd`, `kErrNumNucCharsAtEnd`, `kWarnNumAminoAcidCharsAtEnd`.
///
/// NCBI reference (598d8ae6): c++/include/objtools/readers/seqid_validate.hpp:122-124
/// ```c++
///     size_t  kWarnNumNucCharsAtEnd = 20;
///     size_t  kErrNumNucCharsAtEnd = 25;
///     size_t  kWarnNumAminoAcidCharsAtEnd = 50;
/// ```
const WARN_NUM_NUC_CHARS_AT_END: usize = 20;
const ERR_NUM_NUC_CHARS_AT_END: usize = 25;
const WARN_NUM_AMINO_ACID_CHARS_AT_END: usize = 50;

/// `CFastaIdValidate::CheckForExcessiveNucData(id, lineNum, fReportError)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:240-292
/// ```c++
/// void CFastaIdValidate::CheckForExcessiveNucData(
///     const CSeq_id& id,
///     int lineNum,
///     FReportError fReportError
/// ) const
/// {
///     const auto& idString = id.GetSeqIdString();
///     if (idString.length() > kWarnNumNucCharsAtEnd) {
///         auto numNucChars = CountPossibleNucResidues(idString);
///         if (numNucChars > kWarnNumNucCharsAtEnd) {
///             const string err_message =
///             "Fasta Reader: sequence id ends with " +
///             NStr::NumericToString(numNucChars) +
///             " valid nucleotide characters. " +
///             " Was the sequence accidentally placed in the definition line?";
///
///             auto severity = (numNucChars > kErrNumNucCharsAtEnd) ?
///                             eDiag_Error :
///                             eDiag_Warning;
///
///             fReportError(severity, lineNum, idString, eUnexpectedNucResidues, err_message);
///             return;
///         }
///     }
/// }
///
///
/// static bool s_IsPossibleNuc(unsigned char c)
/// {
///     switch( c ) {
///     case 'N':
///     case 'A':
///     case 'C':
///     case 'G':
///     case 'T':
///     case 'a':
///     case 'c':
///     case 'g':
///     case 't':
///         return true;
///     default:
///         return false;
///     }
/// }
///
///
/// size_t CFastaIdValidate::CountPossibleNucResidues(
///         const string& idString)
/// {
///     const auto first_it = rbegin(idString);
///     const auto it = find_if_not(first_it, rend(idString), s_IsPossibleNuc);
///
///     return static_cast<size_t>(distance(first_it, it));
/// }
/// ```
fn check_for_excessive_nuc_data(id: &ParsedSeqId, warn: Warn<'_>) -> Result<(), ParseError> {
    let id_string = seq_id_string(id);
    if id_string.len() > WARN_NUM_NUC_CHARS_AT_END {
        let count = id_string
            .iter()
            .rev()
            .take_while(|c| b"NACGTacgt".contains(c))
            .count();
        if count > WARN_NUM_NUC_CHARS_AT_END {
            let message = format!(
                "Fasta Reader: sequence id ends with {count} valid nucleotide characters.  Was \
                 the sequence accidentally placed in the definition line?"
            );
            if count > ERR_NUM_NUC_CHARS_AT_END {
                return Err(error(message));
            }
            warn(message.as_bytes());
        }
    }
    Ok(())
}

/// `CFastaIdValidate::CheckForExcessiveProtData(id, lineNum, fReportError)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:295-326
/// ```c++
/// void CFastaIdValidate::CheckForExcessiveProtData(
///         const CSeq_id& id,
///         int lineNum,
///         FReportError fReportError) const
/// {
///     const auto& idString = id.GetSeqIdString();
///
///     // Check for Aa sequence
///     if (idString.length() > kWarnNumAminoAcidCharsAtEnd) {
///         const auto numAaChars = CountPossibleAminoAcids(idString);
///         if (numAaChars > kWarnNumAminoAcidCharsAtEnd) {
///             const string err_message =
///             "Fasta Reader: sequence id ends with " +
///             NStr::NumericToString(numAaChars) +
///             " valid amino-acid characters. " +
///             " Was the sequence accidentally placed in the definition line?";
///             fReportError(eDiag_Warning, lineNum, idString, eUnexpectedAminoAcids, err_message);
///         }
///     }
/// }
///
///
/// size_t CFastaIdValidate::CountPossibleAminoAcids(
///         const string& idString)
/// {
///     const auto first_it = rbegin(idString);
///     const auto it = find_if_not(first_it, rend(idString),
///             [](char c) { return (c >= 'A' && c <= 'Z') ||
///                                 (c >= 'a' && c <= 'z'); });
///
///     return static_cast<size_t>(distance(first_it, it));
/// }
/// ```
fn check_for_excessive_prot_data(id: &ParsedSeqId, warn: Warn<'_>) {
    let id_string = seq_id_string(id);
    if id_string.len() > WARN_NUM_AMINO_ACID_CHARS_AT_END {
        let count = id_string
            .iter()
            .rev()
            .take_while(|c| c.is_ascii_alphabetic())
            .count();
        if count > WARN_NUM_AMINO_ACID_CHARS_AT_END {
            warn(
                format!(
                    "Fasta Reader: sequence id ends with {count} valid amino-acid characters.  \
                     Was the sequence accidentally placed in the definition line?"
                )
                .as_bytes(),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_of_the_id_token() {
        assert_eq!(parse_range(b"seqA:1-10"), (5, 0, 9));
        assert_eq!(parse_range(b"seqB:c10-1"), (6, 9, 0));
        assert_eq!(parse_range(b":1-10"), (5, 0, 9));
        assert_eq!(parse_range(b"x:10-1"), (0, 0, 0));
        assert_eq!(parse_range(b"x:c1-10"), (0, 0, 0));
        assert_eq!(parse_range(b"x:1-"), (0, 0, 0));
        assert_eq!(parse_range(b"x:-1"), (0, 0, 0));
        assert_eq!(parse_range(b"1-10"), (0, 0, 0));
        assert_eq!(parse_range(b"NC_000914.1:1-50"), (5, 0, 49));
        assert_eq!(parse_range(b"a:5-5"), (4, 4, 4));
        assert_eq!(parse_range(b""), (0, 0, 0));
    }

    fn check(id: &str, protein: bool) -> (Result<Vec<ParsedSeqId>, String>, Vec<String>) {
        let mut warnings = Vec::new();
        let result = process_ids(id.as_bytes(), protein, 7, &mut |line: &[u8]| {
            warnings.push(String::from_utf8_lossy(line).into_owned())
        });
        (result.map_err(|e| e.message), warnings)
    }

    // Messages of makeblastdb 2.17.0 (survey B fixtures p05_raw, p06*, p07, p08, n03, n04).
    #[test]
    fn checks_and_messages_of_the_first_parse() {
        let (ids, warnings) = check("a,b", true);
        assert_eq!(
            ids.unwrap(),
            [ParsedSeqId::Local(ObjectId::Str(b"a_b".to_vec()))]
        );
        assert_eq!(
            warnings,
            ["Near line 7, the sequence id string contains 'comma' symbol, which has been replaced with 'underscore' symbol. Please correct the sequence id string."]
        );
        let acc31 = format!("ref|NP_{}|", "1".repeat(28));
        assert_eq!(
            check(&acc31, true).0.unwrap_err(),
            "Near line 7, the accession is too long.  Its length is 31 but the maximum allowed accession length is 30.  Please find and correct all accessions that are too long."
        );
        assert!(check(&format!("ref|NP_{}|", "1".repeat(27)), true)
            .0
            .is_ok());
        let local51 = format!("lcl|{}", "a".repeat(51));
        let (result, warnings) = check(&local51, true);
        assert_eq!(
            result.unwrap_err(),
            "Near line 7, the local id is too long.  Its length is 51 but the maximum allowed local id length is 50.  Please find and correct all local ids that are too long."
        );
        assert_eq!(
            warnings,
            ["Fasta Reader: sequence id ends with 51 valid amino-acid characters.  Was the sequence accidentally placed in the definition line?"]
        );
        assert_eq!(
            check("lcl|[abc]", true).0.unwrap_err(),
            "'[abc]' is not a valid local ID"
        );
        assert!(check("lcl|a[bc]", true).0.is_ok());
        let nuc21 = format!("lcl|x{}", "ACGT".repeat(5) + "A");
        let (result, warnings) = check(&nuc21, false);
        assert!(result.is_ok());
        assert_eq!(
            warnings,
            ["Fasta Reader: sequence id ends with 21 valid nucleotide characters.  Was the sequence accidentally placed in the definition line?"]
        );
        let nuc26 = format!("lcl|x{}", "acgt".repeat(6) + "AC");
        assert_eq!(
            check(&nuc26, false).0.unwrap_err(),
            "Fasta Reader: sequence id ends with 26 valid nucleotide characters.  Was the sequence accidentally placed in the definition line?"
        );
        assert_eq!(
            check("gi|abc|", true).0.unwrap_err(),
            "Could not construct seq-id from 'gi|abc|'"
        );
    }
}
