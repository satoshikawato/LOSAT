//! What every scanner of NCBI's alignment reader shares (`CAlnScanner`, `AlnUtil` and the
//! `CSequenceInfo` the scanner adjusts): the scan result, the failure kinds, the sequence info,
//! the base verification of the rows and the id and line helpers. Each format's scanner is one
//! file next to this one.

use crate::ncbi_str::truncate_spaces;
use crate::utils::format_guess::split_in_two;

/// The alphabet of `CAlnReader::GetAlphabetLetters(eAlpha_Protein)`, which `psiblast -in_msa`
/// sets (`SetClustal(eAlpha_Protein)`) and no scanner changes.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:217-219
pub(crate) const PROTEIN_ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*";

/// Why a scanner gives no alignment.
#[derive(Debug)]
pub(crate) enum ScanFailure {
    /// The scanner or a report threw `SShowStopper`/`CObjReaderParseException`: `x_ReadAsciiMsa`
    /// swallows it ([`super::AlnError::NotAvailable`]).
    Failed,
    /// NCBI's scanner does something undefined here (it crashes, or reads or writes memory that
    /// the input does not give): LOSAT refuses, with this cause.
    Undefined(&'static str),
    /// An exception of NCBI's toolkit that no `catch` of `CAlnReader` or `x_ReadAsciiMsa` takes
    /// (`CStringException` of `NStr::StringToInt`): it ends the run as any uncaught exception.
    Exception(crate::ncbi_exception::NcbiException),
    /// A `std::exception` that no `catch` of the reader takes: psiblast's `CATCH_ALL` writes
    /// `Error: ` and its text, exit status 255.
    StdException(String),
}

pub(crate) type Scan<T> = Result<T, ScanFailure>;

/// The sequence info of the read (`CSequenceInfo`), as far as the scanners of `-in_msa` change
/// it. The alphabet is always [`PROTEIN_ALPHABET`]. The three gap sets of NCBI's info (beginning,
/// middle, end) are set together everywhere, so there is one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SeqInfo {
    pub(crate) gap: Vec<u8>,
    pub(crate) match_chars: Vec<u8>,
    pub(crate) missing: Vec<u8>,
}

impl SeqInfo {
    /// What `CAlnReader::SetClustal(eAlpha_Protein)` sets: the gap `-`, no match or missing
    /// character.
    pub(crate) fn protein() -> Self {
        SeqInfo {
            gap: b"-".to_vec(),
            match_chars: Vec::new(),
            missing: Vec::new(),
        }
    }
}

/// `SAlignmentFile` as the scanners fill it.
pub(crate) struct Scanned {
    pub(crate) ids: Vec<Vec<u8>>,
    pub(crate) rows: Vec<Vec<u8>>,
    /// `mDeflines` (the text after the id of each defline, as the scanner keeps it).
    pub(crate) deflines: Vec<Vec<u8>>,
    pub(crate) info: SeqInfo,
}

/// `xGetExistingSeqIdInfo`: how `seq_id` compares with the ids so far.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner.cpp:205-223
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IdMatch {
    Identical,
    DiffersByCase,
    Different,
}

pub(crate) fn existing_id(ids: &[Vec<u8>], seq_id: &[u8]) -> IdMatch {
    if ids.iter().any(|id| id.as_slice() == seq_id) {
        IdMatch::Identical
    } else if ids.iter().any(|id| id.eq_ignore_ascii_case(seq_id)) {
        IdMatch::DiffersByCase
    } else {
        IdMatch::Different
    }
}

/// `NStr::TruncateSpaces` then `NStr::Split(.., " \t", fSplit_MergeDelimiters)` and `Join(.., "")`
/// of `AlnUtil::StripBlanks`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_util.cpp:222-231
pub(crate) fn strip_blanks(line: &[u8]) -> Vec<u8> {
    truncate_spaces(line)
        .iter()
        .copied()
        .filter(|c| !matches!(c, b' ' | b'\t'))
        .collect()
}

/// `NStr::Split(str, delims, fSplit_MergeDelimiters)`: runs of delimiters make one boundary, a
/// leading delimiter leaves an empty first token and a trailing one an empty last token; an empty
/// string has no token.
pub(crate) fn split_merge<'a>(str: &'a [u8], delims: &[u8]) -> Vec<&'a [u8]> {
    if str.is_empty() {
        return Vec::new();
    }
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i < str.len() {
        if delims.contains(&str[i]) {
            tokens.push(&str[start..i]);
            while i < str.len() && delims.contains(&str[i]) {
                i += 1;
            }
            start = i;
        } else {
            i += 1;
        }
    }
    tokens.push(&str[start..]);
    tokens
}

/// `AlnUtil::ProcessDataLine` (no offset): the first token and the rest joined; fewer than two
/// tokens is a failure.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_util.cpp:195-211
pub(crate) fn process_data_line(line: &[u8]) -> Scan<(Vec<u8>, Vec<u8>)> {
    let tokens = split_merge(line, b" \t");
    if tokens.len() < 2 {
        return Err(ScanFailure::Failed);
    }
    Ok((tokens[0].to_vec(), tokens[1..].concat()))
}

/// What `AlnUtil::ProcessDefline` hands back: the id part (empty when the defline starts with
/// `[`) and the rest.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_util.cpp:140-166
pub(crate) fn process_defline(line: &[u8]) -> Scan<(Vec<u8>, Vec<u8>)> {
    if !line.starts_with(b">") {
        return Err(ScanFailure::Failed);
    }
    let Some(data_start) = line[1..]
        .iter()
        .position(|c| !matches!(c, b' ' | b'\t'))
        .map(|at| at + 1)
    else {
        return Err(ScanFailure::Failed);
    };
    let def_line = &line[data_start..];
    if def_line.starts_with(b"[") {
        return Ok((Vec::new(), def_line.to_vec()));
    }
    // NCBI takes `defLine.substr(dataStart)` of the string it already cut at `dataStart`: the
    // start is past the end (`std::out_of_range`) when the id part is shorter than that.
    if data_start > def_line.len() {
        return Err(ScanFailure::StdException(format!(
            "basic_string::substr: __pos (which is {data_start}) > this->size() (which is {})",
            def_line.len()
        )));
    }
    let (id, info, _) = split_in_two(&def_line[data_start..], b" \t", true);
    Ok((id.to_vec(), info.to_vec()))
}

/// `AlnUtil::CheckId`: the id of a data line against the ids in order (`id_count` is the row
/// the line belongs to), in the first block or after it. Every branch that does not return throws.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_util.cpp:44-135
pub(crate) fn check_id(
    seq_id: &[u8],
    ordered_ids: &[Vec<u8>],
    id_count: usize,
    first_block: bool,
) -> Scan<()> {
    if ordered_ids.len() > id_count && seq_id == ordered_ids[id_count].as_slice() {
        return Ok(());
    }
    let mut found = None;
    for (at, id) in ordered_ids.iter().enumerate() {
        if id.as_slice() == seq_id {
            found = Some((at, true));
            break;
        }
        if id.eq_ignore_ascii_case(seq_id) {
            found = Some((at, false));
            break;
        }
    }
    if first_block {
        return match found {
            Some(_) => Err(ScanFailure::Failed),
            None => Ok(()),
        };
    }
    // After the first block every other outcome is one of the throws ([`ScanFailure::Failed`]).
    Err(ScanFailure::Failed)
}

/// `CAlnScanner::xVerifySingleSequenceData` of every row (the base verification): the row may
/// start with gap and missing characters, its first other character must be a letter of the
/// alphabet, and every character after that one of the alphabet, the gap, missing or match
/// characters (the tail state could only hold the end gap and the missing characters, which the
/// body allows already).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner.cpp:112-169
pub(crate) fn verify_rows(info: &SeqInfo, rows: &[Vec<u8>]) -> Scan<()> {
    for row in rows {
        let mut head = true;
        for &c in row {
            if head {
                if info.gap.contains(&c) || info.missing.contains(&c) {
                    continue;
                }
                head = false;
                if !PROTEIN_ALPHABET.contains(&c) {
                    return Err(ScanFailure::Failed);
                }
                continue;
            }
            let legal = PROTEIN_ALPHABET.contains(&c)
                || info.gap.contains(&c)
                || info.missing.contains(&c)
                || info.match_chars.contains(&c);
            if !legal {
                return Err(ScanFailure::Failed);
            }
        }
    }
    Ok(())
}
