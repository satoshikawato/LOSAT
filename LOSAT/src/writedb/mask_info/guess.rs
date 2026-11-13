//! The format guess of `makeblastdb -mask_data`: `CFormatGuess::GuessFormat` with the preferred
//! formats `eBinaryASN` and `eTextASN` and every other format disabled.
//!
//! `GuessFormat` reads a test buffer (shared with the FASTA guess of `format_guess.rs`), then tries
//! the formats in the order of `sm_CheckOrder`; of the two that are enabled `eTextASN` comes
//! first, `eBinaryASN` last. A file that is neither (or that has no test buffer: empty, or
//! made of comments that end before the buffer is full) is `eUnknown`. A file with no test buffer
//! gets one more test, `TestFormatNewick`, whose answer `eNewick` is just as unusable for
//! `s_ReadObject` as `eUnknown`: "Unknown encoding".

use super::MaskFileFormat;
use crate::writedb::format_guess::{c_isspace, ensure_test_buffer, split_tokenize};

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:453-509
// ```c++
// CFormatGuess::EFormat
// CFormatGuess::GuessFormat( EOnError onerror )
// {
//     ...
//     if (!EnsureTestBuffer()) {
//         ...
//         if(TestFormatNewick(eQuick)) {
//             return CFormatGuess::eNewick;
//         }
//         return CFormatGuess::eUnknown;
//     }
//
//     EMode mode = eQuick;
//
//     // First, try to use hints
//     if ( !m_Hints.IsEmpty() ) {
//         for (size_t f = 0; f < sm_CheckOrder_Size; ++f) {
//             EFormat fmt = EFormat( sm_CheckOrder[f] );
//             if (m_Hints.IsPreferred(fmt)  &&  x_TestFormat(fmt, mode)) {
//                 return fmt;
//             }
//         }
//     }
//
//     // Check other formats, skip the ones that are disabled through hints
//     for (size_t f = 0; f < sm_CheckOrder_Size; ++f) {
//         EFormat fmt = EFormat( sm_CheckOrder[f] );
//         if ( ! m_Hints.IsDisabled(fmt)  &&  x_TestFormat(fmt, mode) ) {
//             return fmt;
//         }
//     }
//     return eUnknown;
// }
// ```
// NCBI reference (598d8ae6): c++/include/util/format_guess.hpp:420-424
// ```c++
// CFormatGuess::CFormatHints::DisableAllNonpreferred(void)
// {
//     m_Disabled = ~m_Preferred;
//     return *this;
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:46-88
// ```c++
// static const CFormatGuess::EFormat sm_CheckOrder[] =
// {
//     ...
//     CFormatGuess::eFasta,
//     CFormatGuess::eTextASN,
//     CFormatGuess::eAlignment,
//     CFormatGuess::eTaxplot,
//     CFormatGuess::eTable,
//     CFormatGuess::eBinaryASN,
//     ...
// };
// ```
/// The format of a mask file, `None` for `eUnknown` (and `eNewick`).
pub(super) fn guess_format(data: &[u8]) -> Option<MaskFileFormat> {
    let mut input = data;
    let sample = ensure_test_buffer(&mut input)?;
    if test_format_text_asn(&sample) {
        return Some(MaskFileFormat::TextAsn);
    }
    if test_format_binary_asn(&sample) {
        return Some(MaskFileFormat::BinaryAsn);
    }
    None
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1161-1185
// ```c++
// bool
// CFormatGuess::TestFormatBinaryAsn(
//     EMode /* not used */ )
// {
//     if ( ! EnsureTestBuffer() ) {
//         return false;
//     }
//
//     //
//     //  Criterion: Presence of any non-printing characters
//     //
//     EConfidence conf = eNo;
//     for (int i = 0;  i < m_iTestDataSize;  ++i) {
//         if ( !isgraph((unsigned char) m_pTestBuffer[i])  &&
//              !isspace((unsigned char) m_pTestBuffer[i]) )
//         {
//             if (m_pTestBuffer[i] == '\1') {
//                 conf = eMaybe;
//             } else {
//                 return true;
//             }
//         }
//     }
//     return (conf == eYes);
// }
// ```
/// `isgraph` and `isspace` of the "C" locale: a byte outside `0x21..=0x7e` and the six white
/// space bytes (NUL, other control bytes, DEL and every byte from 0x80 on) is a "non-printing
/// character"; only `0x01` does not count. `conf` never becomes `eYes`.
fn test_format_binary_asn(sample: &[u8]) -> bool {
    sample
        .iter()
        .any(|&c| !(0x21..=0x7e).contains(&c) && !c_isspace(c) && c != 1)
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1678-1717
// ```c++
// bool
// CFormatGuess::TestFormatTextAsn(
//     EMode /* not used */ )
// {
//     if ( ! EnsureStats() ) {
//         return false;
//     }
//
//     // reject obvious misfits:
//     if ( m_iTestDataSize == 0 || m_pTestBuffer[0] == '>' ) {
//         return false;
//     }
//
//     // criteria:
//     // at least 80% text-ish,
//     // 1st field of the first non-blank not comment line must start with letter.
//     // "::=" as the 2nd field of the first non-blank non comment line.
//     //
//     double dAlNumFraction = (double)(m_iStatsCountAlNumChars+m_iStatsCountBraces) /
//                             (double)m_iTestDataSize;
//     if ( dAlNumFraction < 0.80 ) {
//         return false;
//     }
//
//     string strBuffer(m_pTestBuffer, m_iTestDataSize);
//     CNcbiIstrstream TestBuffer(strBuffer);
//     string strLine;
//
//     while ( ! TestBuffer.fail() ) {
//         vector<string> Fields;
//         NcbiGetline(TestBuffer, strLine, "\n\r");
//         NStr::Split(strLine, " \t", Fields, NStr::fSplit_Tokenize);
//         if ( IsAsnComment( Fields  ) ) {
//             continue;
//         }
//         return ( Fields.size() >= 2 && Fields[1] == "::=" && isalpha(Fields[0][0]));
//     }
//     return false;
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3680-3689
// ```c++
// bool
// CFormatGuess::IsAsnComment(
//     const vector<string>& Fields )
// {
//     if ( Fields.size() == 0 ) {
//         return true;
//     }
//     return ( NStr::StartsWith( Fields[0], "--" ) );
// }
// ```
fn test_format_text_asn(sample: &[u8]) -> bool {
    if sample.is_empty() || sample[0] == b'>' {
        return false;
    }
    let (al_num, braces) = ensure_stats(sample);
    if ((al_num + braces) as f64) / (sample.len() as f64) < 0.80 {
        return false;
    }
    // NcbiGetline with "\n\r": the lines, an empty line gives no fields and counts as a comment.
    for line in split_tokenize(sample, |c| c == b'\r' || c == b'\n') {
        let fields = split_tokenize(line, |c| c == b' ' || c == b'\t');
        if fields.is_empty() || fields[0].starts_with(b"--") {
            continue;
        }
        return fields.len() >= 2 && fields[1] == b"::=" && fields[0][0].is_ascii_alphabetic();
    }
    false
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:687-749
// ```c++
// bool
// CFormatGuess::EnsureStats()
// //  ----------------------------------------------------------------------------
// {
//     ...
//     while ( ! TestBuffer.fail() ) {
//         NcbiGetline( TestBuffer, strLine, "\r\n" );
// // code in CFormatGuess::Format counts line ends
// // so, we will count them here as well
//         if (!strLine.empty()) {
//             strLine += '\n';
//         }
//         size_t size = strLine.size();
//         bool is_header = size > 0 && strLine[0] == '>';
//         for ( size_t i=0; i < size; ++i ) {
//             unsigned char c = strLine[i];
//             unsigned char type = symbol_type_table[c];
//
//             if ( type & (fAlpha | fDigit | fSpace) ) {
//                 ++m_iStatsCountAlNumChars;
//             }
//             else if (c == '{'  ||  c == '}') {
//                 ++m_iStatsCountBraces;
//             }
//             ...
// ```
/// `m_iStatsCountAlNumChars` and `m_iStatsCountBraces` of the sample (`format_guess.rs` counts
/// the first but not the second, which only the ASN test reads).
fn ensure_stats(sample: &[u8]) -> (usize, usize) {
    let (mut al_num, mut braces) = (0usize, 0usize);
    for line in split_tokenize(sample, |c| c == b'\r' || c == b'\n') {
        for &c in line.iter().chain(std::iter::once(&b'\n')) {
            if c != 0 && (c.is_ascii_alphanumeric() || c_isspace(c)) {
                al_num += 1;
            } else if c == b'{' || c == b'}' {
                braces += 1;
            }
        }
    }
    (al_num, braces)
}
