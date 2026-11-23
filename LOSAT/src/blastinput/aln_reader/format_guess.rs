//! NCBI's guess of an alignment file's format (`CAlnFormatGuesser`), ported whole: LOSAT reads the
//! Clustal files that NCBI reads and names the other formats that NCBI would read as such.

use super::base::split_merge;
use super::lines::PeekAheadStream;
use crate::ncbi_str::{string_to_int, string_to_uint, truncate_spaces};

/// `EAlignFormat` of NCBI's guess.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlignFormat {
    Nexus,
    Clustal,
    FastaGap,
    Phylip,
    Sequin,
    MultAlign,
    Unknown,
}

/// `CAlnFormatGuesser::SAMPLE_SIZE`.
const SAMPLE_SIZE: usize = 10;

/// `CAlnFormatGuesser::GetFormat` of the stream `iStr`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:50-74
/// `xInitSample`: the first `SAMPLE_SIZE` lines, each with the white space cut from both ends
/// (`NStr::TruncateSpacesInPlace`), as many as the stream has.
pub(crate) fn guess_format(input: &mut PeekAheadStream<'_>) -> AlignFormat {
    let mut sample: Vec<Vec<u8>> = Vec::new();
    for _ in 0..SAMPLE_SIZE {
        let Some(line) = input.peek_line() else {
            break;
        };
        sample.push(truncate_spaces(line).to_vec());
    }
    if sample.is_empty() {
        return AlignFormat::Unknown;
    }
    if sample_is_nexus(&sample) {
        return AlignFormat::Nexus;
    }
    if sample_is_clustal(&mut sample, input) {
        return AlignFormat::Clustal;
    }
    if sample_is_fasta_gap(&sample) {
        return AlignFormat::FastaGap;
    }
    if sample_is_phylip(&sample) {
        return AlignFormat::Phylip;
    }
    if sample_is_sequin(&sample) {
        return AlignFormat::Sequin;
    }
    if sample_is_mult_align(&sample) {
        return AlignFormat::MultAlign;
    }
    AlignFormat::Unknown
}

/// `NStr::StartsWith` of `str` with `prefix`.
fn starts_with(str: &[u8], prefix: &[u8]) -> bool {
    str.starts_with(prefix)
}

/// `xSampleIsNexus`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:106-108
fn sample_is_nexus(sample: &[Vec<u8>]) -> bool {
    starts_with(&sample[0].to_ascii_lowercase(), b"#nexus")
}

/// `xSampleIsClustal`: a first line `clustalw` or `clustal w` (a `CLUSTAL O` header of Clustal
/// Omega is neither), or three blocks of lines each ended by an empty line that follows a
/// conservation line. The lines past the sample come from the stream as it is, not cut, and the
/// sample grows with them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:123-161
fn sample_is_clustal(sample: &mut Vec<Vec<u8>>, input: &mut PeekAheadStream<'_>) -> bool {
    let first = sample[0].to_ascii_lowercase();
    if starts_with(&first, b"clustalw") || starts_with(&first, b"clustal w") {
        return true;
    }
    let mut line_count = 0;
    let mut block_count = 0;
    while block_count < 3 {
        let line: Vec<u8> = if line_count < sample.len() {
            sample[line_count].clone()
        } else {
            let line = input.peek_line().map(<[u8]>::to_vec).unwrap_or_default();
            sample.push(line.clone());
            line
        };
        if line_count > 0 && line.is_empty() {
            let maybe_consensus = &sample[line_count - 1];
            if !maybe_consensus.iter().any(|c| b".:*".contains(c)) {
                return false;
            }
            if maybe_consensus.iter().any(|c| !b" .:*".contains(c)) {
                return false;
            }
            block_count += 1;
        }
        line_count += 1;
    }
    true
}

/// `xSampleIsFastaGap`: the first line that does not start with `;` starts with `>`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:175-185
fn sample_is_fasta_gap(sample: &[Vec<u8>]) -> bool {
    match sample.iter().find(|line| !starts_with(line, b";")) {
        Some(line) => starts_with(line, b">"),
        None => false,
    }
}

/// `xSampleIsPhylip`: a first line of two numbers.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:199-211
fn sample_is_phylip(sample: &[Vec<u8>]) -> bool {
    let tokens = split_merge(&sample[0], b" \t");
    if tokens.len() != 2 {
        return false;
    }
    let digits = |token: &[u8]| token.iter().all(u8::is_ascii_digit);
    digits(tokens[0]) && digits(tokens[1])
}

/// `xSampleIsSequin`: an empty first line and a second line of the offsets 10, 20, 30, ...
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:225-240
/// (The sample has a second line here: `sample_is_clustal` has grown it to two lines at least
/// when the first line is empty.)
fn sample_is_sequin(sample: &[Vec<u8>]) -> bool {
    if !sample[0].is_empty() {
        return false;
    }
    let Some(second) = sample.get(1) else {
        return false;
    };
    let tokens = split_merge(second, b" \t");
    if tokens.is_empty() {
        return false;
    }
    for (index, token) in tokens.iter().enumerate() {
        // `fConvErr_NoThrow`: 0 for a token that is no number.
        let offset = string_to_uint(token).unwrap_or(0);
        if u64::from(offset) != 10 + 10 * index as u64 {
            return false;
        }
    }
    true
}

/// `xSampleIsMultAlign`: an optional `//` line, an empty line, a line `1 N` (N up to 50), and a
/// line with the name and the N letters of the first block.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_formatguess.cpp:254-295
fn sample_is_mult_align(sample: &[Vec<u8>]) -> bool {
    let mut line_index = 0;
    if starts_with(&sample[0], b"//") {
        line_index += 1;
    }
    if sample.len() < line_index + 4 {
        return false;
    }
    if !sample[line_index].is_empty() {
        return false;
    }
    line_index += 1;
    let tokens = split_merge(&sample[line_index], b" \t");
    if tokens.len() != 2 {
        return false;
    }
    let (Ok(start_offset), Ok(end_offset)) =
        (string_to_int(tokens[0], 0), string_to_int(tokens[1], 0))
    else {
        return false;
    };
    if start_offset != 1 {
        return false;
    }
    if end_offset > 50 {
        return false;
    }
    line_index += 1;
    let data_columns = split_merge(&sample[line_index], b" \t");
    if data_columns.len() < 2 {
        return false;
    }
    let data_size: usize = data_columns[1..].iter().map(|column| column.len()).sum();
    i64::try_from(data_size)
        .is_ok_and(|size| size == i64::from(end_offset) - i64::from(start_offset) + 1)
}
