//! The first parse of the ids of a definition line with `-parse_seqids`: the id branch of
//! `CFastaDeflineReader::ParseDefline` with the flags of makeblastdb (`fRequireID`,
//! `fParseRawID` and `fQuickIDCheck`) or of the mask reader (`fParseRawID`), `fAssumeProt` or
//! `fAssumeNuc`: the range suffix (`ParseRange`), the ids (`x_ProcessIDs`, with
//! `CSeq_id::ParseIDs`) and their checks (`CSeqIdCheck`, `CFastaIdValidate`). Errors end the
//! reading (neither `CFastaBioseqSource` nor `CMaskFastaReader` gives the reader a listener,
//! so `s_PostError` throws); warnings are posted as `LOG_POST` lines
//! (`c++/src/objtools/readers/fasta_reader_utils.cpp`, `seqid_validate.cpp`).

use super::reader::{ErrCode, ParseError, ReaderWarn, ReaderWarning, Thrower};
use super::Warn;
use crate::ncbi_str::is_blank;
use crate::seqid::accguide::{parse_flags, SeqIdWarning};
use crate::seqid::label::content;
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::seqid::parse::{
    has_invalid_local_char, parse_ids, MAX_ACCESSION_LENGTH, MAX_GENERAL_TAG_LENGTH,
    MAX_LOCAL_ID_LENGTH,
};

/// A diagnostic of NCBI's Seq-id code as makeblastdb prints it (`ERR_POST_X` of
/// Seq_id.cpp: `NCBI_ERRCODE_X` 1306 and the subcode; the application's name after the
/// severity).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2676-2679
/// makeblastdb 2.17.0 printed `Warning: (1306.7) [makeblastdb] NCBI C++ Exception:` and the
/// report for `>xyz|ABC|DEF` (survey B fixture p16_unsupported_type).
pub(crate) fn post_seq_id_warning(warning: SeqIdWarning, warn: Warn<'_>) {
    post_seq_id_warning_of("makeblastdb", warning, warn);
}

/// The diagnostic of `post_seq_id_warning` as an application of the prefix `name` writes
/// it (`SetDiagPostPrefix`; the BLAST programs' `-parse_deflines` name theirs).
pub(crate) fn post_seq_id_warning_of(name: &str, warning: SeqIdWarning, warn: Warn<'_>) {
    warn(&seq_id_warning_line(Some(name), warning));
}

/// The line of `post_seq_id_warning_of` without its line feed, after the prefix `[name] ` or none
/// where the application sets none (`common::database::list_warning_line`).
pub(crate) fn seq_id_warning_line(name: Option<&str>, warning: SeqIdWarning) -> Vec<u8> {
    let tag = name.map(|name| format!("[{name}] ")).unwrap_or_default();
    [
        format!("Warning: (1306.{}) {tag}", warning.subcode).as_bytes(),
        &warning.message,
    ]
    .concat()
}

/// `CFastaDeflineReader::ParseRange` of `s` with its start and end: the length of the range suffix
/// of the id token (`:<from>-<to>` or `:c<from>-<to>`), 0 for none, with its ends.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:229-273
/// `TSeqPos` is a 32-bit unsigned integer: the sums and `mult` wrap as NCBI's do.
pub(crate) fn parse_range(s: &[u8]) -> (usize, u32, u32) {
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

/// The `CObjReaderParseException` that `s_PostError` throws, with no listener, the line number,
/// the error message, the problem and the error code (cited at `reader::Thrower::PostError`).
fn error(code: ErrCode, message: String, line_number: u64) -> ParseError {
    ParseError::new(Thrower::PostError, code, message, line_number)
}

/// `CFastaDeflineReader::x_ProcessIDs` of `id_string`, with `info`, `ids`, a listener and
/// `CSeqIdCheck`: the ids of the token, checked. `protein` is `fAssumeProt` (the excessive data
/// check of the nucleotides is made without it, also for convert2blastmask's reader, which sets
/// neither; NCBI's reader without `fAssumeNuc` also makes the check of the proteins, which only
/// posts a warning, and convert2blastmask writes no warning, so LOSAT does not make it there),
/// `quick_id_check` `fQuickIDCheck` (makeblastdb, the BLAST programs' `-parse_deflines`; not the
/// mask readers), `parse_raw_id` `fParseRawID` (set by every caller but convert2blastmask's).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:333-404
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:85,101-104
/// Neither `fAllIdsAsLocal` nor `fNumericIdsAsLocal` is set (base flags 0), and the
/// registry's `MAX_SEQID_LENGTH` is not read (makedb-v5 DECISIONS 12), so the default
/// limits apply.
pub(crate) fn process_ids(
    id_string: &[u8],
    protein: bool,
    quick_id_check: bool,
    parse_raw_id: bool,
    line_number: u64,
    warn: ReaderWarn<'_>,
) -> Result<Vec<ParsedSeqId>, ParseError> {
    let mut flags = parse_flags::PARTIAL_OK | parse_flags::ANY_LOCAL;
    if parse_raw_id {
        flags |= parse_flags::RAW_TEXT;
    }
    let replaced;
    let mut to_parse = id_string;
    if id_string.contains(&b',') && !id_string.contains(&b'|') {
        warn(ReaderWarning::Line(
            format!(
                "Near line {line_number}, the sequence id string contains 'comma' symbol, which \
                 has been replaced with 'underscore' symbol. Please correct the sequence id \
                 string."
            )
            .as_bytes(),
        ));
        replaced = id_string
            .iter()
            .map(|&c| if c == b',' { b'_' } else { c })
            .collect::<Vec<u8>>();
        to_parse = &replaced;
    }
    let ids: Vec<ParsedSeqId> =
        match parse_ids(to_parse, flags, &mut |w| warn(ReaderWarning::SeqId(w))) {
            Ok(ids) => ids
                .into_iter()
                .flatten()
                .filter(|id| !is_blank(&seq_id_string(id)))
                .collect(),
            Err(_) => Vec::new(),
        };
    if ids.is_empty() {
        return Err(error(
            ErrCode::NoIds,
            format!(
                "Could not construct seq-id from '{}'",
                String::from_utf8_lossy(id_string)
            ),
            line_number,
        ));
    }
    validate_ids(&ids, protein, quick_id_check, line_number, warn)?;
    Ok(ids)
}

/// The call of a `CFastaIdValidate` of `flags` on `ids`, the line number and a
/// `CIdErrorReporter` without a listener, with `fAssumeProt` or
/// `fAssumeNuc`, with or without `fQuickIDCheck`: an error throws, a warning is a
/// `LOG_POST` line. The exception's code is the one of `CIdErrorReporter`'s map:
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:305-311
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:107-133
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:215-237
/// The ids come from the first-parse token, which ends at the first byte of 0x20 or below
/// or of 0x80 or more (decision 28), so no byte of a local id is one that
/// `CSeq_id_find_pred` rejects for being a control or non-ASCII byte; only ` >[]|"` can
/// decide here, in the first byte (`fQuickIDCheck`) or in any.
fn validate_ids(
    ids: &[ParsedSeqId],
    protein: bool,
    quick_id_check: bool,
    line_number: u64,
    warn: ReaderWarn<'_>,
) -> Result<(), ParseError> {
    let Some(last) = ids.last() else {
        return Ok(());
    };
    if protein {
        check_for_excessive_prot_data(last, warn);
    } else {
        check_for_excessive_nuc_data(last, line_number, warn)?;
    }
    for id in ids {
        if let ParsedSeqId::Local(ObjectId::Str(text)) = id {
            let check = if quick_id_check {
                &text[..text.len().min(1)]
            } else {
                &text[..]
            };
            if has_invalid_local_char(check) {
                return Err(error(
                    ErrCode::InvalidId,
                    format!(
                        "'{}' is not a valid local ID",
                        String::from_utf8_lossy(text)
                    ),
                    line_number,
                ));
            }
        }
        check_id_length(id, line_number)?;
    }
    Ok(())
}

/// `s_GetIDLengthErrorString` of the length, the id type, the maximum allowed length and the
/// line number.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:154-166
fn id_length_error(length: usize, id_type: &str, max: usize, line_number: u64) -> ParseError {
    error(
        ErrCode::IdTooLong,
        format!(
            "Near line {line_number}, the {id_type} is too long.  Its length is {length} but the \
             maximum allowed {id_type} length is {max}.  Please find and correct all {id_type}s \
             that are too long."
        ),
        line_number,
    )
}

/// `CFastaIdValidate::CheckIDLength` of `id` and the line number, with `fReportError`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:169-213
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
const WARN_NUM_NUC_CHARS_AT_END: usize = 20;
const ERR_NUM_NUC_CHARS_AT_END: usize = 25;
const WARN_NUM_AMINO_ACID_CHARS_AT_END: usize = 50;

/// `CFastaIdValidate::CheckForExcessiveNucData` of `id` and the line number, with
/// `fReportError`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:240-292
fn check_for_excessive_nuc_data(
    id: &ParsedSeqId,
    line_number: u64,
    warn: ReaderWarn<'_>,
) -> Result<(), ParseError> {
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
                return Err(error(ErrCode::Format, message, line_number));
            }
            warn(ReaderWarning::Line(message.as_bytes()));
        }
    }
    Ok(())
}

/// `CFastaIdValidate::CheckForExcessiveProtData` of `id` and the line number, with
/// `fReportError`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/seqid_validate.cpp:295-326
fn check_for_excessive_prot_data(id: &ParsedSeqId, warn: ReaderWarn<'_>) {
    let id_string = seq_id_string(id);
    if id_string.len() > WARN_NUM_AMINO_ACID_CHARS_AT_END {
        let count = id_string
            .iter()
            .rev()
            .take_while(|c| c.is_ascii_alphabetic())
            .count();
        if count > WARN_NUM_AMINO_ACID_CHARS_AT_END {
            warn(ReaderWarning::Line(
                format!(
                    "Fasta Reader: sequence id ends with {count} valid amino-acid characters.  \
                     Was the sequence accidentally placed in the definition line?"
                )
                .as_bytes(),
            ));
        }
    }
}

/// Whether a `CFastaIdValidate` with flags 0 (the alignment reader's `CDefaultIdValidate`:
/// neither `fAssumeProt` nor `fAssumeNuc`, no `fQuickIDCheck`) reports anything for `ids`: an id
/// that ends in more than 20 nucleotide letters or more than 50 letters (warnings, or an error from
/// 25 nucleotide letters), a local id with a byte that `CheckLocalID` rejects, an id over the
/// length limits. The alignment reader has no message listener, so every report, a warning too,
/// throws (`CAlnErrorReporter::Report`, cited at `blastinput::aln_reader`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:162-168
pub(crate) fn alignment_ids_are_reported(ids: &[ParsedSeqId]) -> bool {
    let mut reported = false;
    let mut warn = |_: ReaderWarning<'_>| reported = true;
    // The nucleotide check and the protein check of `CFastaIdValidate::operator()` are the two
    // `validate_ids` calls (the checks of the local ids and of the lengths are the same in both).
    let nuc = validate_ids(ids, false, false, 0, &mut warn);
    let prot = validate_ids(ids, true, false, 0, &mut warn);
    nuc.is_err() || prot.is_err() || reported
}

#[cfg(test)]
#[path = "defline_ids_tests.rs"]
mod tests;
