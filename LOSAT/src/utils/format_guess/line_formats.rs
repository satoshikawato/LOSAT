//! The tests of `CFormatGuess` that read the split lines (`m_TestLines`), in the order of `sm_CheckOrder`, with the
//! line tests they call.

use super::{
    c_isprint, c_string, contains, find, split_in_two, tokens, truncate_spaces, FormatGuess,
    Source, Stop, Tokens, F_DNA_MAIN_ALPHABET, SYMBOL_TYPE_TABLE,
};
use crate::ncbi_str::{is_space, string_to_int, string_to_non_negative_int, string_to_uint};

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:196-245
// ```c++
// static bool s_IsTokenPosInt(
//     const string& strToken )
//     //  ----------------------------------------------------------------------------
// {
//     size_t tokenSize = strToken.size();
//     if (tokenSize == 0) {
//         return false;
//     }
//     if (tokenSize == 1  &&  strToken[0] == '0') {
//         return true;
//     }
//     if (strToken[0] < '1'  ||  '9' < strToken[0]) {
//         return false;
//     }
//     for (size_t i=1; i<tokenSize; ++i) {
//         if (strToken[i] < '0'  ||  '9' < strToken[i]) {
//             return false;
//         }
//     }
//     return true;
// }
//
// //  ----------------------------------------------------------------------------
// static bool s_IsTokenInteger(
//     const string& strToken )
// //  ----------------------------------------------------------------------------
// {
//     if ( ! strToken.empty() && (strToken[0] == '-'  ||  strToken[0] == '+')) {
//         return s_IsTokenPosInt( strToken.substr( 1 ) );
//     }
//     return s_IsTokenPosInt( strToken );
// }
//
// //  ----------------------------------------------------------------------------
// static bool s_IsTokenDouble(
//     const string& strToken )
// {
//     string token( strToken );
//     NStr::ReplaceInPlace( token, ".", "1", 0, 1 );
//     if ( token.size() > 1 && token[0] == '-' ) {
//         token[0] = '1';
//     }
//     if (token.size() > 1 && token[0] == '0') {
//         token[0] = '1';
//     }
//     return s_IsTokenPosInt(token);
// }
// ```
/// `s_IsTokenPosInt`.
fn is_token_pos_int(token: &[u8]) -> bool {
    match token {
        [] => false,
        [b'0'] => true,
        [first, rest @ ..] => (b'1'..=b'9').contains(first) && rest.iter().all(u8::is_ascii_digit),
    }
}

/// `s_IsTokenInteger`.
fn is_token_integer(token: &[u8]) -> bool {
    match token {
        [b'-' | b'+', rest @ ..] => is_token_pos_int(rest),
        _ => is_token_pos_int(token),
    }
}

/// `s_IsTokenDouble`: the first `.` (`ReplaceInPlace(token, ".", "1", 0, 1)`), then a first `-` and a first `0`
/// of a token of two bytes or more, count as `1`.
fn is_token_double(token: &[u8]) -> bool {
    let dot = token.iter().position(|&c| c == b'.');
    // The byte at `i` after the three replacements.
    let at = |i: usize| {
        let mut c = if dot == Some(i) { b'1' } else { token[i] };
        if i == 0 && token.len() > 1 && (c == b'-' || c == b'0') {
            c = b'1';
        }
        c
    };
    match token.len() {
        0 => false,
        1 if at(0) == b'0' => true,
        len => (b'1'..=b'9').contains(&at(0)) && (1..len).all(|i| at(i).is_ascii_digit()),
    }
}

/// `token.size() == 1 && NPOS != token.find_first_of(set)`.
fn one_of(token: &[u8], set: &[u8]) -> bool {
    token.len() == 1 && set.contains(&token[0])
}

impl<S: Source> FormatGuess<S> {
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2085-2122
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     bool ignoreFirstColumn = false;
    //     unsigned int uPslLineCount = 0;
    //     list<string>::iterator it = m_TestLines.begin();
    //     while (it != m_TestLines.end()  &&  NStr::StartsWith(*it, "#")) {
    //         it++;
    //     }
    //     if (it == m_TestLines.end()) {
    //         return false;
    //     }
    //     if (!IsLinePsl(*it, ignoreFirstColumn)) {
    //         ignoreFirstColumn = true;
    //         if (!IsLinePsl(*it, ignoreFirstColumn)) {
    //             return false;
    //         }
    //     }
    //     uPslLineCount++;
    //     it++;
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         if ( ! IsLinePsl(*it, ignoreFirstColumn) ) {
    //             return false;
    //         }
    //         uPslLineCount++;
    //     }
    //     return (uPslLineCount != 0);
    // ```
    /// `TestFormatPsl`.
    pub(super) fn test_format_psl(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut lines = self.lines().skip_while(|line| line.starts_with(b"#"));
        let Some(first) = lines.next() else {
            return false;
        };
        let ignore_first_column = if is_line_psl(first, false) {
            false
        } else if is_line_psl(first, true) {
            true
        } else {
            return false;
        };
        lines.all(|line| is_line_psl(line, ignore_first_column))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:764-773
    // ```c++
    // bool CFormatGuess::TestFormatRepeatMasker(
    //     EMode /* not used */ )
    // {
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     return IsInputRepeatMaskerWithHeader() ||
    //         IsInputRepeatMaskerWithoutHeader();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2851-2943
    // ```c++
    // bool CFormatGuess::IsInputRepeatMaskerWithHeader()
    // {
    //     ...
    //     string labels_1st_line[] = { "SW", "perc", "query", "position", "matching", "" };
    //     string labels_2nd_line[] = { "score", "div.", "del.", "ins.", "sequence", "" };
    //
    //     //
    //     //  Purge junk lines:
    //     //
    //     list<string>::iterator it = m_TestLines.begin();
    //     for  ( ; it != m_TestLines.end(); ++it ) {
    //         NStr::TruncateSpacesInPlace( *it );
    //         if ( *it != "" ) {
    //             break;
    //         }
    //     }
    //
    //     if ( it == m_TestLines.end() ) {
    //         return false;
    //     }
    //
    //     //
    //     //  Verify first line of labels:
    //     //
    //     size_t current_offset = 0;
    //     for ( size_t i=0; labels_1st_line[i] != ""; ++i ) {
    //         current_offset = NStr::FindCase( *it, labels_1st_line[i], current_offset );
    //         if ( current_offset == NPOS ) {
    //             return false;
    //         }
    //     }
    //
    //     //
    //     //  Verify second line of labels:
    //     //
    //     ++it;
    //     if ( it == m_TestLines.end() ) {
    //         return false;
    //     }
    //     current_offset = 0;
    //     for ( size_t j=0; labels_2nd_line[j] != ""; ++j ) {
    //         current_offset = NStr::FindCase( *it, labels_2nd_line[j], current_offset );
    //         if ( current_offset == NPOS ) {
    //             return false;
    //         }
    //     }
    //
    //     //
    //     //  Should have at least one extra line:
    //     //
    //     ++it;
    //     if ( it == m_TestLines.end() ) {
    //         return false;
    //     }
    //
    //     return true;
    // }
    //
    // bool CFormatGuess::IsInputRepeatMaskerWithoutHeader()
    // {
    //     ...
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         string str = NStr::TruncateSpaces( *it );
    //         if ( str == "" ) {
    //             continue;
    //         }
    //         if ( ! IsLineRmo( str ) ) {
    //             return false;
    //         }
    //     }
    //
    //     return true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5503-5510
    // ```c++
    // SIZE_TYPE NStr::FindCase(const CTempString str, const CTempString pattern, SIZE_TYPE start)
    // {
    //     SIZE_TYPE pos = Find(CTempString(str, start), pattern, eCase);
    //     if (pos == NPOS) {
    //         return NPOS;
    //     }
    //     return pos + start;
    // }
    // ```
    /// `TestFormatRepeatMasker`. The header test trims the lines up to the first that is not blank in place: the
    /// later tests read them so.
    pub(super) fn test_format_repeat_masker(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        self.repeat_masker_with_header() || self.repeat_masker_without_header()
    }

    fn repeat_masker_with_header(&mut self) -> bool {
        const LABELS_1ST_LINE: [&[u8]; 5] = [b"SW", b"perc", b"query", b"position", b"matching"];
        const LABELS_2ND_LINE: [&[u8]; 5] = [b"score", b"div.", b"del.", b"ins.", b"sequence"];
        // The lines up to the first that is not blank are trimmed (`extra`, added by a later test, is never here).
        debug_assert!(self.lines.extra.is_none());
        let count = self.lines.count;
        let at = self
            .lines()
            .position(|line| !truncate_spaces(line).is_empty())
            .unwrap_or(count);
        self.lines.trimmed = (at + 1).min(count);
        let labelled = |line: &[u8], labels: &[&[u8]]| {
            let mut current_offset = 0;
            labels.iter().all(|label| {
                match find(line.get(current_offset..).unwrap_or(&[]), label) {
                    Some(pos) => {
                        current_offset += pos;
                        true
                    }
                    None => false,
                }
            })
        };
        let mut lines = self.lines().skip(at);
        lines
            .next()
            .is_some_and(|line| labelled(line, &LABELS_1ST_LINE))
            && lines
                .next()
                .is_some_and(|line| labelled(line, &LABELS_2ND_LINE))
            && lines.next().is_some()
    }

    fn repeat_masker_without_header(&self) -> bool {
        self.lines().all(|line| {
            let str = truncate_spaces(line);
            str.is_empty() || is_line_rmo(str)
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2500-2518
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         if (NStr::StartsWith(*it, "##fileformat=VCFv")) {
    //             return true;
    //         }
    //     }
    //     return false;
    // ```
    /// `TestFormatVcf`.
    pub(super) fn test_format_vcf(&mut self) -> bool {
        self.ensure_stats_and_lines()
            && self.lines_may_hold(b"##fileformat=VCFv")
            && self
                .lines()
                .any(|line| line.starts_with(b"##fileformat=VCFv"))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:863-898
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     unsigned int uGvfLineCount = 0;
    //     list<string>::iterator it = m_TestLines.begin();
    //
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         //
    //         //  Make sure to ignore any UCSC track and browser lines prior to the
    //         //  start of data
    //         //
    //         if ( it->empty() || (*it)[0] == '#' ) {
    //             if (NStr::StartsWith(*it, "##gvf-version")) {
    //                 return true;
    //             }
    //             continue;
    //         }
    //         if ( !uGvfLineCount && NStr::StartsWith( *it, "browser " ) ) {
    //             continue;
    //         }
    //         if ( !uGvfLineCount && NStr::StartsWith( *it, "track " ) ) {
    //             continue;
    //         }
    //         if ( ! IsLineGvf( *it ) ) {
    //             return false;
    //         }
    //         ++uGvfLineCount;
    //     }
    //     return (uGvfLineCount != 0);
    // ```
    /// `TestFormatGvf`.
    pub(super) fn test_format_gvf(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() {
            return Ok(false);
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                if line.starts_with(b"##gvf-version") {
                    return Ok(true);
                }
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line_gvf(line)? {
                return Ok(false);
            }
            count += 1;
        }
        Ok(count != 0)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:901-936
    // ```c++
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         //
    //         //  Make sure to ignore any UCSC track and browser lines prior to the
    //         //  start of data
    //         //
    //         if (!uGffLineCount && NStr::StartsWith(*it, "##gff-version")) {
    //             return NStr::StartsWith(*it, "##gff-version 3");
    //         }
    //         if ( it->empty() || (*it)[0] == '#' ) {
    //             continue;
    //         }
    //         if ( !uGffLineCount && NStr::StartsWith( *it, "browser " ) ) {
    //             continue;
    //         }
    //         if ( !uGffLineCount && NStr::StartsWith( *it, "track " ) ) {
    //             continue;
    //         }
    //         if ( ! IsLineGff3( *it ) ) {
    //             return false;
    //         }
    //         ++uGffLineCount;
    //     }
    //     return (uGffLineCount != 0);
    // ```
    /// `TestFormatGff3` (after `EnsureTestBuffer() && EnsureSplitLines()`).
    pub(super) fn test_format_gff3(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if count == 0 && line.starts_with(b"##gff-version") {
                return line.starts_with(b"##gff-version 3");
            }
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line_gff3(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:829-860 (`TestFormatGtf`), :977-1009
    // (`TestFormatGff2`, the same loop with `IsLineGff2`)
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     unsigned int uGtfLineCount = 0;
    //     list<string>::iterator it = m_TestLines.begin();
    //
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         //
    //         //  Make sure to ignore any UCSC track and browser lines prior to the
    //         //  start of data
    //         //
    //         if ( it->empty() || (*it)[0] == '#' ) {
    //             continue;
    //         }
    //         if ( !uGtfLineCount && NStr::StartsWith( *it, "browser " ) ) {
    //             continue;
    //         }
    //         if ( !uGtfLineCount && NStr::StartsWith( *it, "track " ) ) {
    //             continue;
    //         }
    //         if ( ! IsLineGtf( *it ) ) {
    //             return false;
    //         }
    //         ++uGtfLineCount;
    //     }
    //     return (uGtfLineCount != 0);
    // ```
    /// The loop of `TestFormatGtf` and `TestFormatGff2` with its line test.
    fn gff_lines(&mut self, is_line: fn(&[u8]) -> bool) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    /// `TestFormatGtf`.
    pub(super) fn test_format_gtf(&mut self) -> bool {
        self.gff_lines(is_line_gtf)
    }

    /// `TestFormatGff2`.
    pub(super) fn test_format_gff2(&mut self) -> bool {
        self.gff_lines(is_line_gff2)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:939-974
    // ```c++
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         //
    //         //  Make sure to ignore any UCSC track and browser lines prior to the
    //         //  start of data
    //         //
    //         if (!uGffLineCount && NStr::StartsWith(*it, "##gff-version 3")) {
    //             return false;
    //         }
    //         if ( it->empty() || (*it)[0] == '#' ) {
    //             continue;
    //         }
    //         if ( !uGffLineCount && NStr::StartsWith( *it, "browser " ) ) {
    //             return false;
    //         }
    //         if ( !uGffLineCount && NStr::StartsWith( *it, "track " ) ) {
    //             return false;
    //         }
    //         if ( !IsLineAugustus( *it ) ) {
    //             return false;
    //         }
    //         ++uGffLineCount;
    //     }
    //     return (uGffLineCount != 0);
    // ```
    /// `TestFormatAugustus` (after `EnsureTestBuffer() && EnsureSplitLines()`).
    pub(super) fn test_format_augustus(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if count == 0 && line.starts_with(b"##gff-version 3") {
                return false;
            }
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                return false;
            }
            if !is_line_augustus(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1012-1038
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     /// first line should be a FASTA defline
    //     list<string>::iterator it = m_TestLines.begin();
    //     if (it->empty()  ||  (*it)[0] != '>') {
    //         return false;
    //     }
    //
    //     /// there should be additional data lines, and they should be easily parseable,
    //     ///  with five columns
    //     ++it;
    //     if (it == m_TestLines.end()) {
    //         return false;
    //     }
    //     for ( /**/;  it != m_TestLines.end();  ++it) {
    //         if ( !IsLineGlimmer3( *it ) ) {
    //             return false;
    //         }
    //     }
    //     return true;
    // ```
    /// `TestFormatGlimmer3`: `Err` when a frame is beyond `int` (`NStr::StringToInt` throws).
    pub(super) fn test_format_glimmer3(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() {
            return Ok(false);
        }
        if self.line_at(0).and_then(<[u8]>::first) != Some(&b'>') || self.lines.len() < 2 {
            return Ok(false);
        }
        for line in self.lines().skip(1) {
            if !is_line_glimmer3(line)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1041-1058
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         try {
    //             if ( !IsLineAgp( *it ) ) {
    //                 return false;
    //             }
    //         } catch(...) {
    //             return false;
    //         }
    //     }
    //     return true;
    // ```
    /// `TestFormatAgp` (`IsLineAgp` throws nothing).
    pub(super) fn test_format_agp(&mut self) -> bool {
        self.ensure_lines() && self.lines().all(is_line_agp)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1861-1890
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         if ( NStr::StartsWith( *it, "track" ) ) {
    //             if ( NStr::Find( *it, "type=wiggle_0" ) != NPOS ) {
    //                 return true;
    //             }
    //             if ( NStr::Find( *it, "type=bedGraph" ) != NPOS ) {
    //                 return true;
    //             }
    //         }
    //         if ( NStr::StartsWith(*it, "fixedStep") ) { /* MSS-140 */
    //             if ( NStr::Find(*it, "chrom=")  &&  NStr::Find(*it, "start=") ) {
    //                 return true;
    //             }
    //         }
    //         if ( NStr::StartsWith(*it, "variableStep") ) { /* MSS-140 */
    //             if ( NStr::Find(*it, "chrom=") ) {
    //                 return true;
    //             }
    //             return true;
    //         }
    //     }
    //     return false;
    // ```
    /// `TestFormatWiggle`. `NStr::Find(...)` as a condition is true unless it is 0 (`NPOS` is true too).
    pub(super) fn test_format_wiggle(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        if ![&b"track"[..], b"fixedStep", b"variableStep"]
            .iter()
            .any(|start| self.lines_may_hold(start))
        {
            return false;
        }
        let found = |line: &[u8], pattern: &[u8]| find(line, pattern) != Some(0);
        self.lines().any(|line| {
            (line.starts_with(b"track")
                && (contains(line, b"type=wiggle_0") || contains(line, b"type=bedGraph")))
                || (line.starts_with(b"fixedStep")
                    && found(line, b"chrom=")
                    && found(line, b"start="))
                || line.starts_with(b"variableStep")
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1747-1807
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     bool bTrackLineFound( false );
    //     bool bHasStartAndStop ( false );
    //     size_t columncount = 0;
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         string str = NStr::TruncateSpaces( *it );
    //         if ( str.empty() ) {
    //             continue;
    //         }
    //
    //         // 'chr 8' fixup, the bedreader does this too
    //         if (str.find("chr ") == 0 ||
    //             str.find("Chr ") == 0 ||
    //             str.find("CHR ") == 0)
    //             str.erase(3, 1);
    //
    //         //
    //         //  while occurrence of the following decorations _is_ a good sign, they could
    //         //  also be indicator for a variety of other UCSC data formats
    //         //
    //         if ( NStr::StartsWith( str, "track" ) ) {
    //             bTrackLineFound = true;
    //             continue;
    //         }
    //         if ( NStr::StartsWith( str, "browser" ) ) {
    //             continue;
    //         }
    //         if ( NStr::StartsWith( str, "#" ) ) {
    //             continue;
    //         }
    //
    //         vector<string> columns;
    //         NStr::Split(str, " \t", columns, NStr::fSplit_Tokenize);
    //         if (columns.size() < 3 || columns.size() > 12) {
    //             return false;
    //         }
    //         if ( columns.size() != columncount ) {
    //             if ( columncount == 0 ) {
    //                 columncount = columns.size();
    //             }
    //             else {
    //                 return false;
    //             }
    //         }
    //         if(columns.size() >= 3) {
    //             if (s_IsTokenPosInt(columns[1]) &&
    //                 s_IsTokenPosInt(columns[2])) {
    //                 bHasStartAndStop = true;
    //             }
    //         }
    //     }
    //
    //     return (bHasStartAndStop || bTrackLineFound);
    // ```
    /// `TestFormatBed`.
    pub(super) fn test_format_bed(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        let mut track_line_found = false;
        let mut has_start_and_stop = false;
        let mut columncount = 0;
        for line in self.lines() {
            let mut str = std::borrow::Cow::Borrowed(truncate_spaces(line));
            if str.is_empty() {
                continue;
            }
            if str.starts_with(b"chr ") || str.starts_with(b"Chr ") || str.starts_with(b"CHR ") {
                str.to_mut().remove(3);
            }
            if str.starts_with(b"track") {
                track_line_found = true;
                continue;
            }
            if str.starts_with(b"browser") || str.starts_with(b"#") {
                continue;
            }
            let columns = tokens(&str, b" \t");
            if columns.len() < 3 || columns.len() > 12 {
                return false;
            }
            if columns.len() != columncount {
                if columncount == 0 {
                    columncount = columns.len();
                } else {
                    return false;
                }
            }
            if is_token_pos_int(columns.at(1)) && is_token_pos_int(columns.at(2)) {
                has_start_and_stop = true;
            }
        }
        has_start_and_stop || track_line_found
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1810-1858
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     bool LineFound = false;
    //     size_t columncount = 15;
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         if ( NStr::TruncateSpaces( *it ).empty() ) {
    //             continue;
    //         }
    //         //
    //         //  while occurrence of the following decorations _is_ a good sign, they could
    //         //  also be indicator for a variety of other UCSC data formats
    //         //
    //         if ( NStr::StartsWith( *it, "track" ) ) {
    //             continue;
    //         }
    //         if ( NStr::StartsWith( *it, "browser" ) ) {
    //             continue;
    //         }
    //         if ( NStr::StartsWith( *it, "#" ) ) {
    //             continue;
    //         }
    //
    //         vector<string> columns;
    //         NStr::Split(*it, " \t", columns, NStr::fSplit_Tokenize);
    //         if ( columns.size() != columncount ) {
    //             return false;
    //         } else {
    //             if (!s_IsTokenPosInt(columns[1]) ||   //chr start
    //                 !s_IsTokenPosInt(columns[2]) ||   //chr end
    //                 !s_IsTokenPosInt(columns[4]) ||   //score
    //                 !s_IsTokenPosInt(columns[6]) ||   //thick draw start
    //                 !s_IsTokenPosInt(columns[7]))     //thick draw end
    //                     return false;
    //             string strand = NStr::TruncateSpaces(columns[5]);
    //
    //             if (strand != "+" && strand != "-")
    //                 return false;
    //
    //             LineFound = true;
    //         }
    //     }
    //     return LineFound;
    // ```
    /// `TestFormatBed15`.
    pub(super) fn test_format_bed15(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        let mut line_found = false;
        for line in self.lines() {
            if truncate_spaces(line).is_empty()
                || line.starts_with(b"track")
                || line.starts_with(b"browser")
                || line.starts_with(b"#")
            {
                continue;
            }
            let columns = tokens(line, b" \t");
            if columns.len() != 15
                || ![1, 2, 4, 6, 7]
                    .iter()
                    .all(|&i| is_token_pos_int(columns.at(i)))
                || !matches!(truncate_spaces(columns.at(5)), b"+" | b"-")
            {
                return false;
            }
            line_found = true;
        }
        line_found
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1893-1925
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         const int BUFFSIZE = 1024;
    //         if (m_pTestBuffer) {
    //             delete [] m_pTestBuffer;
    //         }
    //         m_pTestBuffer = new char[BUFFSIZE+1];
    //         m_Stream.read( m_pTestBuffer, BUFFSIZE );
    //         m_iTestDataSize = m_Stream.gcount();
    //         m_pTestBuffer[m_iTestDataSize] = 0;
    //         m_Stream.clear();  // in case we reached eof
    //         CStreamUtils::Stepback( m_Stream, m_pTestBuffer, m_iTestDataSize );
    //         m_TestLines.push_back(m_pTestBuffer);
    //     }
    //
    //     unsigned int uHgvsLineCount = 0;
    //     list<string>::iterator it = m_TestLines.begin();
    //
    //     for ( ;  it != m_TestLines.end();  ++it) {
    //         if ( it->empty() || (*it)[0] == '#' ) {
    //             continue;
    //         }
    //         if ( ! IsLineHgvs( *it ) ) {
    //             return false;
    //         }
    //         ++uHgvsLineCount;
    //     }
    //     return (uHgvsLineCount != 0);
    // ```
    /// `TestFormatHgvs` (its own buffer, as `TestFormatNewick`'s, is never read in `GuessFormat`, whose
    /// `TestFormatNewick` before it leaves lines).
    pub(super) fn test_format_hgvs(&mut self) -> bool {
        if !self.ensure_stats() || !self.ensure_split_lines() {
            const BUFFSIZE: usize = 1024;
            let read = self.stream.read(BUFFSIZE);
            self.stream.clear(); // in case we reached eof
            self.stream.stepback(&read);
            self.buffer = Some(read);
            self.terminated = true;
            let line = self.c_string_line();
            self.push_line(line);
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if !is_line_hgvs(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1189-1234
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     //
    //     // criteria are odd:
    //     //
    //     list<string>::const_iterator iter = m_TestLines.begin();
    //     list<string> toks;
    //
    //     /// first line: one token, one number
    //     NStr::Split(*iter++, "\t ", toks, NStr::fSplit_Tokenize);
    //     if (toks.size() != 1  ||
    //         toks.front().find_first_not_of("0123456789") != string::npos) {
    //         return false;
    //     }
    //
    //     // now, for remaining ones, we expect an alphanumeric item first,
    //     // followed by a set of floating-point values.  Unless we are at the last
    //     // line, the number of values should increase monotonically
    //     for (size_t i = 1;  iter != m_TestLines.end();  ++i, ++iter) {
    //         toks.clear();
    //         NStr::Split(*iter, "\t ", toks, NStr::fSplit_Tokenize);
    //         if (toks.size() != i) {
    //             /// we can ignore the last line ; it may be truncated
    //             list<string>::const_iterator it = iter;
    //             ++it;
    //             if (it != m_TestLines.end()) {
    //                 return false;
    //             }
    //         }
    //
    //         list<string>::const_iterator it = toks.begin();
    //         for (++it;  it != toks.end();  ++it) {
    //             if ( ! s_IsTokenDouble( *it ) ) {
    //                 return false;
    //             }
    //         }
    //     }
    //
    //     return true;
    // ```
    /// `TestFormatDistanceMatrix`. A last line without a token: libstdc++'s `++` of the end of an empty
    /// `std::list` is the end again (the list is a ring through its sentinel), so its value loop does not run.
    pub(super) fn test_format_distance_matrix(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let first = tokens(self.line_at(0).unwrap_or_default(), b"\t ");
        if first.len() != 1 || !first.at(0).iter().all(u8::is_ascii_digit) {
            return false;
        }
        let rest = self.lines.len() - 1;
        for (i, line) in self.lines().skip(1).enumerate() {
            let mut toks = line
                .split(|&c| c == b'\t' || c == b' ')
                .filter(|token| !token.is_empty());
            let count = toks.clone().count();
            if count != i + 1 && i + 1 != rest {
                return false;
            }
            if !toks.by_ref().skip(1).all(is_token_double) {
                return false;
            }
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1237-1251
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     ITERATE (list<string>, it, m_TestLines) {
    //         if ( !IsLineFlatFileSequence( *it ) ) {
    //             return false;
    //         }
    //     }
    //     return true;
    // ```
    /// `TestFormatFlatFileSequence`.
    pub(super) fn test_format_flat_file_sequence(&mut self) -> bool {
        self.ensure_lines() && self.lines().all(is_line_flat_file_sequence)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2400-2497
    // ```c++
    // bool
    // UniProtGetLineData(
    //     list<string>::iterator& lineIt,
    //     list<string>::iterator endIt,
    //     string& lineCode,
    //     string& lineData)
    // //  ----------------------------------------------------------------------------
    // {
    //     if (lineIt == endIt) {
    //         return false;
    //     }
    //     try {
    //         NStr::SplitInTwo(
    //             *lineIt, " ", lineCode, lineData, NStr::fSplit_MergeDelimiters);
    //     }
    //     catch(CException&) {
    //         lineCode = *lineIt;
    //         lineData = "";
    //     }
    //     lineIt++;
    //     return true;
    // }
    //
    // //  ----------------------------------------------------------------------------
    // bool CFormatGuess::TestFormatFlatFileUniProt(
    //     EMode /*unused*/)
    // {
    //     // see: https://web.expasy.org/docs/userman.html#genstruc
    //
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     // smell test:
    //     // note: sample size at least 8000 characters, line length soft limited to
    //     //  75 characters
    //     if (m_TestLines.size() < 15) { // number of required records
    //         return false;
    //     }
    //
    //     ...
    //     string lineCode, lineData, lookingFor;
    //     auto recordIt = m_TestLines.begin();
    //     auto endIt = m_TestLines.end();
    //     NStr::SplitInTwo(
    //         *recordIt, " ", lineCode, lineData, NStr::fSplit_MergeDelimiters);
    //
    //     lookingFor = "ID"; // excatly one
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     recordIt++;
    //
    //     lookingFor = "AC"; // one or more
    //     if (!UniProtGetLineData(recordIt, endIt, lineCode, lineData)) {
    //         return false;
    //     }
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     while (lineCode == lookingFor) {
    //         if (!UniProtGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "DT"; // three (first hard difference from UniProt)
    //     for (int i = 0; i < 3; ++i) {
    //         if (lineCode != lookingFor) {
    //             return false;
    //         }
    //         if (!UniProtGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return false;
    //         }
    //     }
    //
    //
    //     lookingFor = "DE"; // one or more
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     while (lineCode == lookingFor) {
    //         if (!UniProtGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return true;
    //         }
    //     }
    //
    //     // optional "GN" line or first "OS" line
    //     if (lineCode != "GN"  &&  lineCode != "OS") {
    //         return false;
    //     }
    //
    //     //  once here it's UniProt or someone is messing with me
    //     return true;
    // }
    // ```
    /// `TestFormatFlatFileUniProt` (`SplitInTwo` throws nothing).
    pub(super) fn test_format_flat_file_uniprot(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 15 {
            return false;
        }
        let mut records = self.lines().map(|line| split_in_two(line, b" ", true).0);
        if records.next() != Some(b"ID".as_slice()) {
            return false;
        }
        let Some(mut line_code) = records.next() else {
            return false;
        };
        if line_code != b"AC" {
            return false;
        }
        while line_code == b"AC" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for _ in 0..3 {
            if line_code != b"DT" {
                return false;
            }
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        if line_code != b"DE" {
            return false;
        }
        while line_code == b"DE" {
            let Some(next) = records.next() else {
                return true;
            };
            line_code = next;
        }
        line_code == b"GN" || line_code == b"OS"
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2270-2397
    // ```c++
    // bool
    // EnaGetLineData(
    //     list<string>::iterator& lineIt,
    //     list<string>::iterator endIt,
    //     string& lineCode,
    //     string& lineData)
    // //  ----------------------------------------------------------------------------
    // {
    //     while (lineIt != endIt  &&  NStr::StartsWith(*lineIt, "XX")) {
    //         lineIt++;
    //     }
    //     if (lineIt == endIt) {
    //         return false;
    //     }
    //     try {
    //         NStr::SplitInTwo(
    //             *lineIt, " ", lineCode, lineData, NStr::fSplit_MergeDelimiters);
    //     }
    //     catch(CException&) {
    //         lineCode = *lineIt;
    //         lineData = "";
    //     }
    //     lineIt++;
    //     return true;
    // }
    //
    // //  ----------------------------------------------------------------------------
    // bool CFormatGuess::TestFormatFlatFileEna(
    //     EMode /*unused*/)
    // {
    //     ...
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     ...
    //     if (m_TestLines.size() < 19) { // number of required records
    //         return false;
    //     }
    //
    //     string lineCode, lineData, lookingFor;
    //     auto recordIt = m_TestLines.begin();
    //     auto endIt = m_TestLines.end();
    //     NStr::SplitInTwo(
    //         *recordIt, " ", lineCode, lineData, NStr::fSplit_MergeDelimiters);
    //
    //     lookingFor = "ID"; // excactly one
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     recordIt++;
    //
    //     lookingFor = "AC"; // one or more
    //     if (!EnaGetLineData(recordIt, endIt, lineCode, lineData)) {
    //         return false;
    //     }
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     while (lineCode == lookingFor) {
    //         if (!EnaGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "PR"; // zero or more
    //     while (lineCode == lookingFor) {
    //         if (!EnaGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "DT"; // two (first hard difference from UniProt)
    //     for (int i = 0; i < 2; ++i) {
    //         if (lineCode != lookingFor) {
    //             return false;
    //         }
    //         if (!EnaGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "DE"; // one or more
    //     if (lineCode != lookingFor) {
    //         return false;
    //     }
    //     while (lineCode == lookingFor) {
    //         if (!EnaGetLineData(recordIt, endIt, lineCode, lineData)) {
    //             return true;
    //         }
    //     }
    //
    //     lookingFor = "KW"; // one or more
    //     ... (the same as "DE", then "OS" and "OC")
    //
    //     //  once here it's Ena or someone is messing with me
    //     return true;
    // }
    // ```
    /// `TestFormatFlatFileEna`.
    pub(super) fn test_format_flat_file_ena(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 19 {
            return false;
        }
        if split_in_two(self.line_at(0).unwrap_or_default(), b" ", true).0 != b"ID" {
            return false;
        }
        let mut records = self
            .lines()
            .skip(1)
            .filter(|line| !line.starts_with(b"XX"))
            .map(|line| split_in_two(line, b" ", true).0);
        let Some(mut line_code) = records.next() else {
            return false;
        };
        if line_code != b"AC" {
            return false;
        }
        while line_code == b"AC" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        while line_code == b"PR" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for _ in 0..2 {
            if line_code != b"DT" {
                return false;
            }
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for looking_for in [&b"DE"[..], b"KW", b"OS", b"OC"] {
            if line_code != looking_for {
                return false;
            }
            while line_code == looking_for {
                let Some(next) = records.next() else {
                    return true;
                };
                line_code = next;
            }
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2125-2165
    // ```c++
    // bool
    // GenbankGetKeywordLine(
    //     list<string>::iterator& lineIt,
    //     list<string>::iterator endIt,
    //     string& keyword,
    //     string& data)
    // //  ----------------------------------------------------------------------------
    // {
    //     if (lineIt == endIt) {
    //         return false;
    //     }
    //     if (lineIt->size() > 79) {
    //         return false;
    //     }
    //
    //     vector<int> validIndents = {0, 2, 3, 5, 12, 21};
    //     auto firstNotBlank = lineIt->find_first_not_of(" ");
    //     while (firstNotBlank != 0) {
    //         if (std::find(validIndents.begin(), validIndents.end(), firstNotBlank) ==
    //                 validIndents.end()) {
    //             auto firstNotBlankOrDigit = lineIt->find_first_not_of(" 1234567890");
    //             if (firstNotBlankOrDigit != 10) {
    //                 return false;
    //             }
    //         }
    //         lineIt++;
    //         if (lineIt == endIt) {
    //             return false;
    //         }
    //         firstNotBlank = lineIt->find_first_not_of(" ");
    //     }
    //     try {
    //         NStr::SplitInTwo(
    //             *lineIt, " ", keyword, data, NStr::fSplit_MergeDelimiters);
    //     }
    //     catch (CException&) {
    //         return false;
    //     }
    //     lineIt++;
    //     return true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2168-2267
    // ```c++
    //     if ( ! EnsureStats() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     ...
    //     if (m_TestLines.size() < 9) { // number of required records
    //         return false;
    //     }
    //
    //     string keyword, data, lookingFor;
    //     auto recordIt = m_TestLines.begin();
    //     auto endIt = m_TestLines.end();
    //     NStr::SplitInTwo(
    //         *recordIt, " ", keyword, data, NStr::fSplit_MergeDelimiters);
    //
    //     lookingFor = "LOCUS"; // excactly one
    //     if (keyword != lookingFor) {
    //         return false;
    //     }
    //     recordIt++;
    //     if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //         return false;
    //     }
    //
    //     lookingFor = "DEFINITION"; // one or more
    //     if (keyword != lookingFor) {
    //         return false;
    //     }
    //     while (keyword == lookingFor) {
    //         if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "ACCESSION"; // one or more
    //     ... (the same as "DEFINITION")
    //
    //     bool nidSeen = false;
    //     lookingFor = "NID"; // zero or one, can come before or after VERSION
    //     if (keyword == lookingFor) {
    //         nidSeen = true;
    //         if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "VERSION"; // exactly one
    //     if (keyword != lookingFor) {
    //         return false;
    //     }
    //     if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //         return false;
    //     }
    //
    //     if (!nidSeen) {
    //         lookingFor = "NID"; // zero or one
    //         if (keyword == lookingFor) {
    //             if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //                 return false;
    //             }
    //         }
    //     }
    //
    //     lookingFor = "PROJECT"; // zero or more
    //     while (keyword == lookingFor) {
    //         if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "DBLINK"; // zero or more
    //     while (keyword == lookingFor) {
    //         if (!GenbankGetKeywordLine(recordIt, endIt, keyword, data)) {
    //             return false;
    //         }
    //     }
    //
    //     lookingFor = "KEYWORDS"; // one or more
    //     if (keyword != lookingFor) {
    //         return false;
    //     }
    //
    //     // I am convinced now. There may be flaws farther down but this input
    //     //  definitely wants to be a Genbank flat file.
    //     return true;
    // ```
    /// `TestFormatFlatFileGenbank`.
    pub(super) fn test_format_flat_file_genbank(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 9 {
            return false;
        }
        if split_in_two(self.line_at(0).unwrap_or_default(), b" ", true).0 != b"LOCUS" {
            return false;
        }
        let mut records = self.lines().skip(1);
        // `GenbankGetKeywordLine`: the keyword of the next line (after the continuation lines), the lines past it.
        let mut next_keyword = || -> Option<&[u8]> {
            let mut current = records.next()?;
            if current.len() > 79 {
                return None;
            }
            const VALID_INDENTS: [usize; 6] = [0, 2, 3, 5, 12, 21];
            let first_not = |line: &[u8], set: &[u8]| line.iter().position(|c| !set.contains(c));
            let mut first_not_blank = first_not(current, b" ");
            while first_not_blank != Some(0) {
                if !first_not_blank.is_some_and(|indent| VALID_INDENTS.contains(&indent))
                    && first_not(current, b" 1234567890") != Some(10)
                {
                    return None;
                }
                current = records.next()?;
                first_not_blank = first_not(current, b" ");
            }
            Some(split_in_two(current, b" ", true).0)
        };
        let Some(mut keyword) = next_keyword() else {
            return false;
        };
        for looking_for in [&b"DEFINITION"[..], b"ACCESSION"] {
            if keyword != looking_for {
                return false;
            }
            while keyword == looking_for {
                let Some(next) = next_keyword() else {
                    return false;
                };
                keyword = next;
            }
        }
        let nid_seen = keyword == b"NID";
        if nid_seen {
            let Some(next) = next_keyword() else {
                return false;
            };
            keyword = next;
        }
        if keyword != b"VERSION" {
            return false;
        }
        let Some(next) = next_keyword() else {
            return false;
        };
        keyword = next;
        if !nid_seen && keyword == b"NID" {
            let Some(next) = next_keyword() else {
                return false;
            };
            keyword = next;
        }
        for looking_for in [&b"PROJECT"[..], b"DBLINK"] {
            while keyword == looking_for {
                let Some(next) = next_keyword() else {
                    return false;
                };
                keyword = next;
            }
        }
        keyword == b"KEYWORDS"
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1254-1274
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         if (it->empty()) {
    //             continue;
    //         }
    //
    //         if (it->find(">Feature ") != 0 && it->find(">Features ") != 0) {
    //             return false;
    //         }
    //         break;
    //     }
    //
    //     return true;
    // ```
    /// `TestFormatFiveColFeatureTable`: the first line that is not empty starts with `>Feature ` or `>Features `
    /// (true when there is none: the empty line `TestFormatNewick` adds for a buffer that starts with NUL).
    pub(super) fn test_format_five_col_feature_table(&mut self) -> bool {
        self.ensure_lines()
            && self
                .lines()
                .find(|line| !line.is_empty())
                .is_none_or(|line| {
                    line.starts_with(b">Feature ") || line.starts_with(b">Features ")
                })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1727-1744
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         string str = *it;
    //         int rsid, chr, pos, numMatched;
    //         numMatched = sscanf( it->c_str(), "rs%d\t%d\t%d", &rsid, &chr, &pos);
    //         if ( numMatched == 3) {
    //             return true;
    //         }
    //     }
    //     return false;
    // ```
    /// `TestFormatSnpMarkers`.
    pub(super) fn test_format_snp_markers(&mut self) -> bool {
        self.ensure_lines()
            && self.lines_may_hold(b"rs")
            && self
                .lines()
                .any(|line| line.starts_with(b"rs") && scans_snp_marker(c_string(line)))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1316-1345
    // ```c++
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //
    //     if (TestFormatCLUSTAL()) {
    //         return true;
    //     }
    //     ...
    //     ITERATE( list<string>, it, m_TestLines ) {
    //         if ( NPOS != it->find( "#NEXUS" ) ) {
    //             return true;
    //         }
    //     }
    //     return false;
    // ```
    /// `TestFormatAlignment`.
    pub(super) fn test_format_alignment(&mut self) -> bool {
        self.ensure_lines()
            && (self.test_format_clustal()
                || (self.lines_may_hold(b"#NEXUS")
                    && self.lines().any(|line| contains(line, b"#NEXUS"))))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1434-1576
    // ```c++
    // bool
    // CFormatGuess::TestFormatCLUSTAL()
    // {
    //
    //     if (!EnsureTestBuffer()) {
    //         return false;
    //     }
    //
    //     string strBuffer(m_pTestBuffer, m_iTestDataSize);
    //     CNcbiIstrstream TestBuffer(strBuffer);
    //     string strLine;
    //
    //     SClustalBlockInfo block_info;
    //
    //     bool has_valid_block = false;
    //     size_t seg_length = 0;
    //     size_t seg_length_prev = 0;
    //
    //
    //     const bool buffer_full = m_iTestDataSize == m_iTestBufferSize;
    //
    //     while ( !TestBuffer.eof() ) {
    //         NcbiGetline(TestBuffer, strLine, "\r\n");
    //
    //         if (buffer_full &&
    //             TestBuffer.eof()) { // Skip last line if buffer is full
    //             break;              // to avoid misidentification due to line truncation
    //         }
    //
    //         if (TestBuffer.fail()) {
    //             break;
    //         }
    //
    //         if (NStr::StartsWith(strLine, "CLUSTAL")) {
    //             continue;
    //         }
    //
    //         if (NStr::IsBlank(strLine)) {
    //             if (block_info.m_InBlock) {
    //                 if (block_info.m_Size < 2) {
    //                     return false;
    //                 }
    //                 block_info.Reset();
    //             }
    //             continue;
    //         }
    //
    //         if (x_LooksLikeCLUSTALConservedInfo(strLine)) {
    //             if (! block_info.m_InBlock || block_info.m_Size<2) {
    //                 return false;
    //             }
    //             block_info.Reset();
    //             continue;
    //         }
    //
    //         string seq_id;
    //         if (!x_TryProcessCLUSTALSeqData(strLine, seq_id, seg_length)) {
    //             return false;
    //         }
    //
    //         if (seg_length > 60) {
    //             return false;
    //         }
    //         if (block_info.m_InBlock) {
    //             if(seg_length != seg_length_prev) {
    //                 return false;
    //             }
    //             has_valid_block = true;
    //         }
    //
    //         if (block_info.m_Ids.find(seq_id) != block_info.m_Ids.end()) {
    //             return false;
    //         }
    //         block_info.m_Ids.insert(seq_id);
    //
    //         seg_length_prev = seg_length;
    //         block_info.m_InBlock = true;
    //         ++(block_info.m_Size);
    //     }
    //
    //     return has_valid_block;
    // }
    // ```
    /// `TestFormatCLUSTAL`. `NcbiGetline(..., "\r\n")` (quoted at `stats_of`) takes `\r\n` as one line end; a last
    /// line without an end sets `eofbit` (with `failbit` too when it is empty), and the read after a last line end
    /// finds nothing (`eofbit | failbit`).
    fn test_format_clustal(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        let buffer_full = data.len() == self.buffer_size;
        let mut pos = 0;
        let mut eof = false;
        let mut in_block = false;
        let mut block_size = 0usize;
        let mut ids: Vec<&[u8]> = Vec::new();
        let mut has_valid_block = false;
        let mut seg_length_prev = 0usize;
        while !eof {
            let start = pos;
            while pos < data.len() && data[pos] != b'\r' && data[pos] != b'\n' {
                pos += 1;
            }
            let line = &data[start..pos];
            let fail;
            if pos == data.len() {
                eof = true;
                fail = line.is_empty();
            } else {
                fail = false;
                if data[pos] == b'\r' && data.get(pos + 1) == Some(&b'\n') {
                    pos += 1;
                }
                pos += 1;
            }
            if (buffer_full && eof) || fail {
                break;
            }
            if line.starts_with(b"CLUSTAL") {
                continue;
            }
            if line.iter().all(|&c| is_space(c)) {
                if in_block {
                    if block_size < 2 {
                        return false;
                    }
                    in_block = false;
                    block_size = 0;
                    ids.clear();
                }
                continue;
            }
            if looks_like_clustal_conserved_info(line) {
                if !in_block || block_size < 2 {
                    return false;
                }
                in_block = false;
                block_size = 0;
                ids.clear();
                continue;
            }
            let Some((seq_id, seg_length)) = try_process_clustal_seq_data(line) else {
                return false;
            };
            if seg_length > 60 {
                return false;
            }
            if in_block {
                if seg_length != seg_length_prev {
                    return false;
                }
                has_valid_block = true;
            }
            if ids.contains(&seq_id) {
                return false;
            }
            ids.push(seq_id);
            seg_length_prev = seg_length;
            in_block = true;
            block_size += 1;
        }
        has_valid_block
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1579-1621
    // ```c++
    // bool
    // CFormatGuess::TestFormatTable(
    //     EMode /* not used */ )
    // {
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //     if ( ! IsAsciiText()) {//gp-13007:  "table" means "ascii table"
    //         return false;
    //     }
    //     ...
    //     //' ' ' \t' '\t' ',' '|'
    //     if (x_TestTableDelimiter(" "))
    //         return true;
    //     else if (x_TestTableDelimiter(" \t"))
    //         return true;
    //     else if (x_TestTableDelimiter("\t"))
    //         return true;
    //     else if (x_TestTableDelimiter(","))
    //         return true;
    //     else if (x_TestTableDelimiter("|"))
    //         return true;
    //
    //     return false;
    // }
    // ```
    /// `TestFormatTable`.
    pub(super) fn test_format_table(&mut self) -> bool {
        if !self.ensure_lines() || !self.is_ascii_text() {
            return false;
        }
        [&b" "[..], b" \t", b"\t", b",", b"|"]
            .iter()
            .any(|delims| self.test_table_delimiter(delims))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1521-1576 (`x_TestTableDelimiter`)
    // ```c++
    //     list<string>::const_iterator iter = m_TestLines.begin();
    //     list<string> toks;
    //
    //     // Skip initial lines since not all headers start with comments like # or ;:
    //     // Don't skip though if file is very short - add up to 3, 1 for each line
    //     // over 5:
    //     for (size_t i=5; i<7; ++i)
    //         if (m_TestLines.size() > i)
    //             ++iter;
    //
    //     /// determine the number of observed columns
    //     size_t ncols = 0;
    //     for ( ;  iter != m_TestLines.end();  ++iter) {
    //         if (iter->empty()  ||  (*iter)[0] == '#'  ||  (*iter)[0] == ';') {
    //             continue;
    //         }
    //
    //         toks.clear();
    //         NStr::Split(*iter, delims, toks, NStr::fSplit_Tokenize);
    //         ncols = toks.size();
    //         break;
    //     }
    //     if ( ncols < 2 ) {
    //         return false;
    //     }
    //
    //     size_t nlines = 1;
    //     // verify that columns all have the same size
    //     // we can add an exception for the last line
    //     for ( ;  iter != m_TestLines.end();  ++iter) {
    //         if (iter->empty()  ||  (*iter)[0] == '#'  ||  (*iter)[0] == ';') {
    //             continue;
    //         }
    //
    //         toks.clear();
    //         NStr::Split(*iter, delims, toks, NStr::fSplit_Tokenize);
    //         if (toks.size() != ncols) {
    //             list<string>::const_iterator it = iter;
    //             ++it;
    //             if (it != m_TestLines.end() || (m_iTestDataSize < m_iTestBufferSize) ) {
    //                 return false;
    //             }
    //         } else {
    //             ++nlines;
    //         }
    //         // Tokens should only contain printable characters
    //         for (const auto& token : toks) {
    //             auto it = find_if(token.begin(), token.end(),
    //                     [](unsigned char c){ return !isprint(c); });
    //             if (it != token.end()) {
    //                 return false;
    //             }
    //         }
    //     }
    //     return ( nlines >= 3 );
    // ```
    /// `x_TestTableDelimiter(delims)`: the count starts again at the line that gave the columns, so that line counts
    /// twice.
    fn test_table_delimiter(&self, delims: &[u8]) -> bool {
        let skip = (5..7).filter(|&i| self.lines.len() > i).count();
        let commented = |line: &[u8]| line.is_empty() || line[0] == b'#' || line[0] == b';';
        let Some((start, first)) = self
            .lines()
            .enumerate()
            .skip(skip)
            .find(|(_, line)| !commented(line))
        else {
            return false;
        };
        let ncols = tokens(first, delims).len();
        if ncols < 2 {
            return false;
        }
        let data_size = self.data().len();
        let mut nlines = 1;
        for (i, line) in self.lines().enumerate().skip(start) {
            if commented(line) {
                continue;
            }
            if tokens(line, delims).len() != ncols {
                if i + 1 != self.lines.len() || data_size < self.buffer_size {
                    return false;
                }
            } else {
                nlines += 1;
            }
            // Every token printable: every byte that is not a delimiter.
            if line.iter().any(|c| !delims.contains(c) && !c_isprint(*c)) {
                return false;
            }
        }
        nlines >= 3
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:801-826
    // ```c++
    // bool
    // CFormatGuess::TestFormatPhrapAce(
    //     EMode /* not used */ )
    // {
    //     if ( ! EnsureTestBuffer() || ! EnsureSplitLines() ) {
    //         return false;
    //     }
    //
    //     if (memchr(m_pTestBuffer, 0, m_iTestDataSize)) { // Cannot contain NuLL bytes
    //         return false;                                // RW-1102
    //     }
    //
    //     bool foundId = false;
    //     for (const auto& line : m_TestLines) {
    //         if (foundId) {
    //             if (s_LooksLikeNucSeqData(line)) {
    //                 return true;
    //             }
    //         }
    //         else if (IsLinePhrapId(line)) {
    //             foundId = true;
    //         }
    //     }
    //     return false;
    // }
    // ```
    /// `TestFormatPhrapAce`.
    pub(super) fn test_format_phrap_ace(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() || self.data().contains(&0) {
            return Ok(false);
        }
        let mut found_id = false;
        for line in self.lines() {
            if found_id {
                if looks_like_nuc_seq_data(line) {
                    return Ok(true);
                }
            } else if is_line_phrap_id(line)? {
                found_id = true;
            }
        }
        Ok(false)
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3615-3677
// ```c++
// bool
// CFormatGuess::IsLinePsl(
//     const string& line,
//     bool ignoreFirstLine)
// //  ----------------------------------------------------------------------------
// {
//     vector<string> tokens;
//     int firstColumn = (ignoreFirstLine ? 1 : 0);
//     NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize);
//     if (tokens.size() - firstColumn != 21) {
//         return false;
//     }
//     // first 8 columns are positive integers:
//     for (auto column = firstColumn; column < firstColumn + 8; ++column) {
//         if (!s_IsTokenPosInt(tokens[column]) ) {
//             return false;
//         }
//     }
//     // next is one or two "+" or "-":
//     const string& token = tokens[firstColumn + 8];
//     if (token.empty()  ||  token.size() > 2) {
//         return false;
//     }
//     if (token.find_first_not_of("-+") != string::npos) {
//         return false;
//     }
//     // columns 11 - 13 are positive integers:
//     for (auto column = firstColumn + 10; column < firstColumn + 13; ++column) {
//         if (!s_IsTokenPosInt(tokens[column]) ) {
//             return false;
//         }
//     }
//     // columns 15 - 18 are positive integers:
//     for (auto column = firstColumn + 14; column < firstColumn + 18; ++column) {
//         if (!s_IsTokenPosInt(tokens[column]) ) {
//             return false;
//         }
//     }
//     ...
//     return true;
// }
// ```
/// `IsLinePsl(line, ignoreFirstLine)`.
fn is_line_psl(line: &[u8], ignore_first_column: bool) -> bool {
    let tokens = tokens(line, b" \t");
    let first = usize::from(ignore_first_column);
    if tokens.len().wrapping_sub(first) != 21 {
        return false;
    }
    let token = tokens.at(first + 8);
    (first..first + 8).all(|column| is_token_pos_int(tokens.at(column)))
        && !token.is_empty()
        && token.len() <= 2
        && token.iter().all(|&c| c == b'-' || c == b'+')
        && (first + 10..first + 13).all(|column| is_token_pos_int(tokens.at(column)))
        && (first + 14..first + 18).all(|column| is_token_pos_int(tokens.at(column)))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3542-3612
// ```c++
//     const size_t MIN_VALUES_PER_RECORD = 14;
//
//     //
//     //  Make sure there is enough stuff on that line:
//     //
//     list<string> values;
//     if ( NStr::Split(line, " \t", values, NStr::fSplit_Tokenize).size() < MIN_VALUES_PER_RECORD ) {
//         return false;
//     }
//
//     //
//     //  Look at specific values and make sure they are of the correct type:
//     //
//
//     //  1: positive integer:
//     list<string>::iterator it = values.begin();
//     if ( ! s_IsTokenPosInt( *it ) ) {
//         return false;
//     }
//
//     //  2: float:
//     ++it;
//     if ( ! s_IsTokenDouble( *it ) ) {
//         return false;
//     }
//
//     //  3: float:
//     ++it;
//     if ( ! s_IsTokenDouble( *it ) ) {
//         return false;
//     }
//
//     //  4: float:
//     ++it;
//     if ( ! s_IsTokenDouble( *it ) ) {
//         return false;
//     }
//
//     //  5: string, not checked
//     ++it;
//
//     //  6: positive integer:
//     ++it;
//     if ( ! s_IsTokenPosInt( *it ) ) {
//         return false;
//     }
//
//     //  7: positive integer:
//     ++it;
//     if ( ! s_IsTokenPosInt( *it ) ) {
//         return false;
//     }
//
//     //  8: positive integer, likely in paretheses, not checked:
//     ++it;
//
//     //  9: '+' or 'C':
//     ++it;
//     if ( *it != "+" && *it != "C" ) {
//         return false;
//     }
//     ...
//     return true;
// ```
/// `IsLineRmo(line)`.
fn is_line_rmo(line: &[u8]) -> bool {
    let values = tokens(line, b" \t");
    values.len() >= 14
        && is_token_pos_int(values.at(0))
        && (1..4).all(|i| is_token_double(values.at(i)))
        && is_token_pos_int(values.at(5))
        && is_token_pos_int(values.at(6))
        && matches!(values.at(8), b"+" | b"C")
}

/// The columns 3 to 7 of the GFF family (`tokens[3]` to `tokens[7]` of `IsLineGtf`, `IsLineGff3`, `IsLineGff2`):
/// two positions, a score, a strand of `strands`, a phase.
fn gff_columns(tokens: &Tokens<'_>, strands: &[u8]) -> bool {
    is_token_pos_int(tokens.at(3))
        && is_token_pos_int(tokens.at(4))
        && is_token_double(tokens.at(5))
        && one_of(tokens.at(6), strands)
        && one_of(tokens.at(7), b".0123")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3242-3270
// ```c++
//     vector<string> tokens;
//     if ( NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize).size() < 8 ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[3] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[4] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenDouble( tokens[5] ) ) {
//         return false;
//     }
//     if ( tokens[6].size() != 1 || NPOS == tokens[6].find_first_of( ".+-" ) ) {
//         return false;
//     }
//     if ( tokens[7].size() != 1 || NPOS == tokens[7].find_first_of( ".0123" ) ) {
//         return false;
//     }
//     if ( tokens.size() < 9 ||
//          (NPOS == tokens[8].find( "gene_id" ) && NPOS == tokens[8].find( "transcript_id" ) ) ) {
//         return false;
//     }
//     return true;
// ```
/// `IsLineGtf(line)`.
fn is_line_gtf(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    tokens.len() >= 9
        && gff_columns(&tokens, b".+-")
        && (contains(tokens.at(8), b"gene_id") || contains(tokens.at(8), b"transcript_id"))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3273-3346
// ```c++
//     vector<string> tokens;
//     if ( NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize).size() < 8 ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[3] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[4] ) ) {
//         return false;
//     }
//
//     //make sure that "type" is a GVF admissible value:
//     {{
//         bool typeOk = false;
//         list<string> terms;
//         terms.push_back("snv");
//         ... ("cnv" through "sequence_alteration", 25 terms)
//         ITERATE(list<string>, termiter, terms) {
//             if(NStr::EqualNocase(*termiter, tokens[2])) {
//                 typeOk = true;
//                 break;
//             }
//         }
//         if (!typeOk) {
//             return false;
//         }
//     }}
//     if ( ! s_IsTokenDouble( tokens[5] ) ) {
//         return false;
//     }
//     if ( tokens[6].size() != 1 || NPOS == tokens[6].find_first_of( ".+-" ) ) {
//         return false;
//     }
//     if ( tokens[7].size() != 1 || NPOS == tokens[7].find_first_of( ".0123" ) ) {
//         return false;
//     }
//
//     //make sure all the mandatory attributes are present:
//     string attrs = tokens[8];
//     if (string::npos == attrs.find("ID="))
//         return false;
//     if (string::npos == attrs.find("Variant_seq=")) {
//         return false;
//     }
//     return true;
// ```
/// `IsLineGvf(line)`: a line of eight columns that passes the column tests reads `tokens[8]` past the end of the
/// vector (`Stop::Undecided`).
fn is_line_gvf(line: &[u8]) -> Result<bool, Stop> {
    const TERMS: [&[u8]; 25] = [
        b"snv",
        b"cnv",
        b"copy_number_variation",
        b"gain",
        b"copy_number_gain",
        b"loss",
        b"copy_number_loss",
        b"loss_of_heterozygosity",
        b"complex",
        b"complex_substitution",
        b"complex_sequence_alteration",
        b"indel",
        b"insertion",
        b"inversion",
        b"substitution",
        b"deletion",
        b"duplication",
        b"translocation",
        b"upd",
        b"uniparental_disomy",
        b"maternal_uniparental_disomy",
        b"paternal_uniparental_disomy",
        b"tandom_duplication",
        b"structural_variation",
        b"sequence_alteration",
    ];
    let tokens = tokens(line, b" \t");
    if tokens.len() < 8
        || !is_token_pos_int(tokens.at(3))
        || !is_token_pos_int(tokens.at(4))
        || !TERMS
            .iter()
            .any(|term| term.eq_ignore_ascii_case(tokens.at(2)))
        || !is_token_double(tokens.at(5))
        || !one_of(tokens.at(6), b".+-")
        || !one_of(tokens.at(7), b".0123")
    {
        return Ok(false);
    }
    let Some(attrs) = tokens.get(8) else {
        return Err(Stop::Undecided(
            "TestFormatGvf reads the ninth column of a GVF line of eight",
        ));
    };
    Ok(contains(attrs, b"ID=") && contains(attrs, b"Variant_seq="))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3349-3390
// ```c++
//     vector<string> tokens;
//     if ( NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize).size() < 8 ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[3] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[4] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenDouble( tokens[5] ) ) {
//         return false;
//     }
//     if ( tokens[6].size() != 1 || NPOS == tokens[6].find_first_of( ".+-?" ) ) {
//         return false;
//     }
//     if ( tokens[7].size() != 1 || NPOS == tokens[7].find_first_of( ".0123" ) ) {
//         return false;
//     }
//     if ( tokens.size() < 9 || tokens[8].empty()) {
//         return false;
//     }
//     if ( tokens.size() >= 9 && tokens[8].size() > 1) {
//         const string& col9 = tokens[8];
//         if ( NPOS == NStr::Find(col9, "ID") &&
//              NPOS == NStr::Find(col9, "Parent") &&
//              NPOS == NStr::Find(col9, "Target") &&
//              NPOS == NStr::Find(col9, "Name") &&
//              NPOS == NStr::Find(col9, "Alias") &&
//              NPOS == NStr::Find(col9, "Note") &&
//              NPOS == NStr::Find(col9, "Dbxref") &&
//              NPOS == NStr::Find(col9, "Xref") ) {
//             return false;
//         }
//     }
//
//     return true;
// ```
/// `IsLineGff3(line)`.
fn is_line_gff3(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    if tokens.len() < 9 || !gff_columns(&tokens, b".+-?") {
        return false;
    }
    let col9 = tokens.at(8);
    col9.len() <= 1
        || [
            &b"ID"[..],
            b"Parent",
            b"Target",
            b"Name",
            b"Alias",
            b"Note",
            b"Dbxref",
            b"Xref",
        ]
        .iter()
        .any(|key| contains(col9, key))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3486-3511
// ```c++
//     vector<string> tokens;
//     const size_t num_cols = NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize).size();
//     if ( num_cols < 8 ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[3] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenPosInt( tokens[4] ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenDouble( tokens[5] ) ) {
//         return false;
//     }
//     if ( tokens[6].size() != 1 || NPOS == tokens[6].find_first_of( ".+-" ) ) {
//         return false;
//     }
//     if ( tokens[7].size() != 1 || NPOS == tokens[7].find_first_of( ".0123" ) ) {
//         return false;
//     }
//     return true;
// ```
/// `IsLineGff2(line)`.
fn is_line_gff2(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    tokens.len() >= 8 && gff_columns(&tokens, b".+-")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3393-3483
// ```c++
//     vector<string> tokens;
//     string remaining(line), head, tail;
//
//     //column 0: ID, string
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 1: method, most likely "AUGUSTUS" but don't want to rely on this
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 2: feature type, controlled vocabulary
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)) {
//         return false;
//     }
//     remaining = tail;
//     string featureType = head;
//
//     //column 3: start, integer
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)  ||  !s_IsTokenPosInt(head)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 4: stop, integer
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)  ||  !s_IsTokenPosInt(head)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 5: score, double
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)  ||  !s_IsTokenDouble(head)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 6: strand, one in "+-.?"
//     const string legalStrands{"+-.?"};
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)  ||  head.size() != 1  ||
//             string::npos == legalStrands.find(head)) {
//         return false;
//     }
//     remaining = tail;
//
//     //column 7: phase, one in ".0123"
//     const string legalPhases{".0123"};
//     if (!NStr::SplitInTwo(remaining, " \t", head, tail)  ||  head.size() != 1  ||
//             string::npos == legalPhases.find(head)) {
//         return false;
//     }
//     remaining = tail;
//
//     //everything else: attributes, format depends on featureType
//     if (remaining.empty()) {
//         return false;
//     }
//
//     if (featureType == "gene") {
//         if (NPOS != NStr::Find(remaining, ";")) {
//             return false;
//         }
//         if (NPOS != NStr::Find(remaining, " ")) {
//             return false;
//         }
//         return true;
//     }
//     if (featureType == "transcript") {
//         ... (the same as "gene")
//     }
//     if (NPOS == NStr::Find(remaining, "transcript_id")) {
//         return false;
//     }
//     if (NPOS == NStr::Find(remaining, "gene_id")) {
//         return false;
//     }
//     return true;
// ```
/// `IsLineAugustus(line)`: `SplitInTwo` without flags, so a run of delimiters gives empty columns.
fn is_line_augustus(line: &[u8]) -> bool {
    // The columns 0 to 7: an id, a method, the feature type, start, stop, score, strand and phase.
    let tests: [fn(&[u8]) -> bool; 8] = [
        |_| true,
        |_| true,
        |_| true,
        is_token_pos_int,
        is_token_pos_int,
        is_token_double,
        |head| one_of(head, b"+-.?"),
        |head| one_of(head, b".0123"),
    ];
    let mut remaining = line;
    let mut feature_type: &[u8] = &[];
    for (column, test) in tests.iter().enumerate() {
        let (head, tail, found) = split_in_two(remaining, b" \t", false);
        if !found || !test(head) {
            return false;
        }
        if column == 2 {
            feature_type = head;
        }
        remaining = tail;
    }
    if remaining.is_empty() {
        return false;
    }
    if feature_type == b"gene" || feature_type == b"transcript" {
        return !contains(remaining, b";") && !contains(remaining, b" ");
    }
    contains(remaining, b"transcript_id") && contains(remaining, b"gene_id")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3201-3239
// ```c++
//     list<string> toks;
//     NStr::Split(line, "\t ", toks, NStr::fSplit_Tokenize);
//     if (toks.size() != 5) {
//         return false;
//     }
//
//     list<string>::iterator i = toks.begin();
//
//     /// first column: skip (ascii identifier)
//     ++i;
//
//     /// second, third columns: both ints
//     if ( ! s_IsTokenInteger( *i++ ) ) {
//         return false;
//     }
//     if ( ! s_IsTokenInteger( *i++ ) ) {
//         return false;
//     }
//
//     /// fourth column: int in the range of -3...3
//     if ( ! s_IsTokenInteger( *i ) ) {
//         return false;
//     }
//     int frame = NStr::StringToInt( *i++ );
//     if (frame < -3  ||  frame > 3) {
//         return false;
//     }
//
//     /// fifth column: score; double
//     if ( ! s_IsTokenDouble( *i ) ) {
//         return false;
//     }
//
//     return true;
// ```
/// `IsLineGlimmer3(line)`: `NStr::StringToInt` with its default flags throws for a frame beyond `int`
/// (`CStringException`, NCBI BLAST+ 2.17.0 psiblast `-in_pssm`: "Cannot convert string '99999999999' to int,
/// overflow", exit status 255).
fn is_line_glimmer3(line: &[u8]) -> Result<bool, Stop> {
    let toks = tokens(line, b"\t ");
    if toks.len() != 5
        || !is_token_integer(toks.at(1))
        || !is_token_integer(toks.at(2))
        || !is_token_integer(toks.at(3))
    {
        return Ok(false);
    }
    let frame = string_to_int(toks.at(3), 0).map_err(Stop::Exception)?;
    Ok((-3..=3).contains(&frame) && is_token_double(toks.at(4)))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3129-3198
// ```c++
//     string line( strLine );
//     size_t uCommentStart = NStr::Find( line, "#" );
//
//     if ( NPOS != uCommentStart ) {
//         line = line.substr( 0, uCommentStart );
//     }
//     NStr::TruncateSpacesInPlace( line );
//     if ( line.empty() ) {
//         return true;
//     }
//
//     vector<string> tokens;
//     if ( NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize).size() < 8 ) {
//         return false;
//     }
//
//     if ( tokens[1].size() > 1 && tokens[1][0] == '-' ) {
//         tokens[1][0] = '1';
//     }
//     if ( -1 == NStr::StringToNonNegativeInt( tokens[1] ) ) {
//         return false;
//     }
//
//     if ( tokens[2].size() > 1 && tokens[2][0] == '-' ) {
//         tokens[2][0] = '1';
//     }
//     if ( -1 == NStr::StringToNonNegativeInt( tokens[2] ) ) {
//         return false;
//     }
//
//     if ( tokens[3].size() > 1 && tokens[3][0] == '-' ) {
//         tokens[3][0] = '1';
//     }
//     if ( -1 == NStr::StringToNonNegativeInt( tokens[3] ) ) {
//         return false;
//     }
//
//     if ( tokens[4].size() != 1 || NPOS == tokens[4].find_first_of( "ADFGPNOW" ) ) {
//         return false;
//     }
//     if ( tokens[4] == "N" ) {
//         if ( -1 == NStr::StringToNonNegativeInt( tokens[5] ) ) {
//             return false;
//         }
//     }
//     else {
//         if ( -1 == NStr::StringToNonNegativeInt( tokens[6] ) ) {
//             return false;
//         }
//         if ( -1 == NStr::StringToNonNegativeInt( tokens[7] ) ) {
//             return false;
//         }
//         if ( tokens.size() != 9 ) {
//             return false;
//         }
//         if ( tokens[8].size() != 1 || NPOS == tokens[8].find_first_of( "+-" ) ) {
//             return false;
//         }
//     }
//
//     return true;
// ```
/// `IsLineAgp(strLine)`.
fn is_line_agp(line: &[u8]) -> bool {
    let line = match find(line, b"#") {
        Some(comment) => &line[..comment],
        None => line,
    };
    let line = truncate_spaces(line);
    if line.is_empty() {
        return true;
    }
    let tokens = tokens(line, b" \t");
    if tokens.len() < 8 {
        return false;
    }
    // A first `-` counts as `1`: a number of 11 bytes or fewer is converted from a copy (`1` and ten more digits
    // are 10^10 or more, past `INT_MAX`, as any longer number).
    let non_negative_after_dash = |token: &[u8]| {
        if token.len() > 1 && token[0] == b'-' {
            if token.len() > 11 {
                return false;
            }
            let mut copy = [0u8; 11];
            copy[..token.len()].copy_from_slice(token);
            copy[0] = b'1';
            return string_to_non_negative_int(&copy[..token.len()]) != -1;
        }
        string_to_non_negative_int(token) != -1
    };
    if !(1..4).all(|i| non_negative_after_dash(tokens.at(i))) || !one_of(tokens.at(4), b"ADFGPNOW")
    {
        return false;
    }
    if tokens.at(4) == b"N" {
        string_to_non_negative_int(tokens.at(5)) != -1
    } else {
        string_to_non_negative_int(tokens.at(6)) != -1
            && string_to_non_negative_int(tokens.at(7)) != -1
            && tokens.len() == 9
            && one_of(tokens.at(8), b"+-")
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3795-3843
// ```c++
// bool CFormatGuess::IsLineHgvs(
//     const string& line )
// {
//     // This simple check can mistake Newwick, so Newwick is checked first
//     //  /[:alnum:]+:(g|c|r|p|m|mt|n)\.[:alnum:]+/  as in NC_000001.9:g.1234567C>T
//     int State = 0;
//     ITERATE(string, Iter, line) {
//         char Char = *Iter;
//         char Next = '\0';
//         string::const_iterator NextI = Iter;
//         ++NextI;
//         if(NextI != line.end())
//             Next = *NextI;
//
//         if(State == 0) {
//             if( isalnum(Char) )
//                 State++;
//         } else if(State == 1) {
//             if(Char == ':')
//                 State++;
//         } else if(State == 2) {
//             if (Char == 'g' ||
//                 Char == 'c' ||
//                 Char == 'r' ||
//                 Char == 'p' ||
//                 Char == 'n' ||
//                 Char == 'm' ) {
//                 State++;
//                 if (Char=='m' && Next == 't') {
//                     ++Iter;
//                 }
//             } else {
//                 return false;
//             }
//         } else if(State == 3) {
//             if(Char == '.')
//                 State++;
//             else
//                 return false;
//         } else if(State == 4) {
//             if( isalnum(Char) )
//                 State++;
//         }
//     }
//
//     return (State == 5);
// }
// ```
/// `IsLineHgvs(line)`.
fn is_line_hgvs(line: &[u8]) -> bool {
    let mut state = 0;
    let mut i = 0;
    while i < line.len() {
        let c = line[i];
        let next = line.get(i + 1).copied().unwrap_or(0);
        match state {
            0 if c.is_ascii_alphanumeric() => state = 1,
            1 if c == b':' => state = 2,
            2 => {
                if !matches!(c, b'g' | b'c' | b'r' | b'p' | b'n' | b'm') {
                    return false;
                }
                state = 3;
                if c == b'm' && next == b't' {
                    i += 1;
                }
            }
            3 => {
                if c != b'.' {
                    return false;
                }
                state = 4;
            }
            4 if c.is_ascii_alphanumeric() => state = 5,
            _ => {}
        }
        i += 1;
    }
    state == 5
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3074-3100
// ```c++
//     // blocks of ten residues (or permitted punctuation characters)
//     // with a count at the start or end; require at least four
//     // (normally six)
//     SIZE_TYPE pos = line.find_first_not_of("0123456789 \t");
//     if (pos == NPOS  ||  pos + 45 >= line.size()) {
//         return false;
//     }
//
//     for (SIZE_TYPE i = 0;  i < 45;  ++i) {
//         char c = line[pos + i];
//         if (i % 11 == 10) {
//             if ( !isspace(c) ) {
//                 return false;
//             }
//         } else {
//             if ( !isalpha(c)  &&  c != '-'  &&  c != '*') {
//                 return false;
//             }
//         }
//     }
//
//     return true;
// ```
/// `IsLineFlatFileSequence(line)`.
fn is_line_flat_file_sequence(line: &[u8]) -> bool {
    let Some(pos) = line.iter().position(|c| !b"0123456789 \t".contains(c)) else {
        return false;
    };
    if pos + 45 >= line.len() {
        return false;
    }
    line[pos..pos + 45].iter().enumerate().all(|(i, &c)| {
        if i % 11 == 10 {
            is_space(c)
        } else {
            c.is_ascii_alphabetic() || c == b'-' || c == b'*'
        }
    })
}

/// `sscanf(line, "rs%d\t%d\t%d", ...) == 3` of glibc on a C string: `rs`, then three conversions of `%d` (white
/// space, a sign, one digit or more; the value does not matter), each `\t` of the format any white space, none
/// included.
fn scans_snp_marker(line: &[u8]) -> bool {
    let Some(mut rest) = line.strip_prefix(b"rs") else {
        return false;
    };
    for i in 0..3 {
        if i > 0 {
            rest = &rest[rest.iter().take_while(|&&c| is_space(c)).count()..];
        }
        rest = &rest[rest.iter().take_while(|&&c| is_space(c)).count()..];
        if let [b'+' | b'-', after @ ..] = rest {
            rest = after;
        }
        let digits = rest.iter().take_while(|c| c.is_ascii_digit()).count();
        if digits == 0 {
            return false;
        }
        rest = &rest[digits..];
    }
    true
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1348-1362
// ```c++
// bool CFormatGuess::x_LooksLikeCLUSTALConservedInfo(const string& line) const
// {
//     for (auto c : line) {
//         if ( isspace(c)) {
//             continue;
//         }
//
//         if (c != ':' &&
//             c != '*' &&
//             c != '.') {
//             return false;
//         }
//     }
//     return true;
// }
// ```
/// `x_LooksLikeCLUSTALConservedInfo(line)`.
fn looks_like_clustal_conserved_info(line: &[u8]) -> bool {
    line.iter()
        .all(|&c| is_space(c) || matches!(c, b':' | b'*' | b'.'))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1366-1406
// ```c++
//     vector<string> toks;
//     NStr::Split(line, " \t", toks, NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
//     const size_t num_toks = toks.size();
//
//     if (num_toks != 2 &&
//         num_toks != 3) {
//         return false;
//     }
//
//     const string& seqdata = toks[1];
//
//
//     unsigned int cumulated_res = 0;
//     if (num_toks == 3) {
//         cumulated_res = NStr::StringToUInt(toks[2], NStr::fConvErr_NoThrow);
//         if (cumulated_res == 0) {
//             return false;
//         }
//     }
//
//     // Check sequence data
//     ESequenceType seqtype =
//         SequenceType(seqdata.c_str(), static_cast<unsigned int>(seqdata.size()), eST_Strict);
//
//     if (seqtype == eUndefined) {
//         return false;
//     }
//
//     if (num_toks == 3) {
//         size_t num_gaps = count(seqdata.begin(), seqdata.end(), '-');
//         if (((seqdata.size() - num_gaps) > cumulated_res)) {
//             return false;
//         }
//     }
//
//
//     id = toks[0];
//     seg_length = seqdata.size();
//
//     return true;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:569-590 (`S2N_CONVERT_ERROR`: with
// `fConvErr_NoThrow` the conversion returns 0)
/// `x_TryProcessCLUSTALSeqData(line, id, seg_length)`: the id and the length.
fn try_process_clustal_seq_data(line: &[u8]) -> Option<(&[u8], usize)> {
    let toks = tokens(line, b" \t");
    if toks.len() != 2 && toks.len() != 3 {
        return None;
    }
    let seqdata = toks.at(1);
    let mut cumulated_res = 0;
    if toks.len() == 3 {
        cumulated_res = string_to_uint(toks.at(2)).unwrap_or(0);
        if cumulated_res == 0 {
            return None;
        }
    }
    if !strict_sequence_type(seqdata) {
        return None;
    }
    if toks.len() == 3 {
        let num_gaps = seqdata.iter().filter(|&&c| c == b'-').count();
        if (seqdata.len() - num_gaps) as u64 > u64::from(cumulated_res) {
            return None;
        }
    }
    Some((toks.at(0), seqdata.len()))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:307-374 (`CFormatGuess::SequenceType`, from `case
// eST_Strict`)
// ```c++
//     for (unsigned i = 0; i < length; ++i) {
//         unsigned char c = str[i];
//         unsigned char type = symbol_type_table[c];
//         if ( type & fDNA_Main_Alphabet ) {
//             ++main_nuc_content;
//         } else if ( type & fDNA_Ambig_Alphabet ) {
//             ++ambig_content;
//         } else if ( !(type & (fSpace | fDigit)) ) {
//             ++bad_nuc_content;
//         }
//
//         if ( type & fProtein_Alphabet ) {
//             ++amino_acid_content;
//         } else if ( type & fAlpha ) {
//             ++exotic_aa_content;
//         } else if ( !(type & (fSpace | fDigit)) ) {
//             ++bad_aa_content;
//         }
//     }
//
//     switch (strictness) {
//     ...
//     case eST_Strict: // Must be 100% encodable
//         if (bad_nuc_content == 0  &&  ambig_content <= main_nuc_content / 3) {
//             return eNucleotide;
//         } else if (bad_aa_content == 0
//                    &&  exotic_aa_content <= amino_acid_content / 9) {
//             return eProtein;
//         }
//     }
//
//     return eUndefined;
// ```
/// `SequenceType(str, size, eST_Strict) != eUndefined`.
fn strict_sequence_type(str: &[u8]) -> bool {
    use super::{F_ALPHA, F_DIGIT, F_DNA_AMBIG_ALPHABET, F_PROTEIN_ALPHABET, F_SPACE};
    let (mut main_nuc, mut ambig, mut bad_nuc, mut amino_acid, mut exotic_aa, mut bad_aa) =
        (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    for &c in str {
        let type_ = SYMBOL_TYPE_TABLE[c as usize];
        if type_ & F_DNA_MAIN_ALPHABET != 0 {
            main_nuc += 1;
        } else if type_ & F_DNA_AMBIG_ALPHABET != 0 {
            ambig += 1;
        } else if type_ & (F_SPACE | F_DIGIT) == 0 {
            bad_nuc += 1;
        }
        if type_ & F_PROTEIN_ALPHABET != 0 {
            amino_acid += 1;
        } else if type_ & F_ALPHA != 0 {
            exotic_aa += 1;
        } else if type_ & (F_SPACE | F_DIGIT) == 0 {
            bad_aa += 1;
        }
    }
    (bad_nuc == 0 && ambig <= main_nuc / 3) || (bad_aa == 0 && exotic_aa <= amino_acid / 9)
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3514-3539
// ```c++
//     vector<string> values;
//     if ( NStr::Split(line, " \t", values, NStr::fSplit_Tokenize).empty() ) {
//         return false;
//     }
//
//     //
//     //  Old style: "^DNA \\w+ "
//     //
//     if ( values[0] == "DNA" ) {
//         return true;
//     }
//
//     //
//     //  New style: "^AS [0-9]+ [0-9]+"
//     //
//     if ( values[0] == "AS" ) {
//         return ( 0 <= NStr::StringToNonNegativeInt( values[1] ) &&
//           0 <= NStr::StringToNonNegativeInt( values[2] ) );
//     }
//
//     return false;
// ```
/// `IsLinePhrapId(line)`: an `AS` line of fewer than three tokens reads past the end of the vector
/// (`Stop::Undecided`) where `&&` does not stop first.
fn is_line_phrap_id(line: &[u8]) -> Result<bool, Stop> {
    let values = tokens(line, b" \t");
    match values.get(0) {
        Some(b"DNA") => Ok(true),
        Some(b"AS") => {
            const PAST: &str = "TestFormatPhrapAce reads past the tokens of an AS line";
            let Some(second) = values.get(1) else {
                return Err(Stop::Undecided(PAST));
            };
            if string_to_non_negative_int(second) < 0 {
                return Ok(false);
            }
            let Some(third) = values.get(2) else {
                return Err(Stop::Undecided(PAST));
            };
            Ok(string_to_non_negative_int(third) >= 0)
        }
        _ => Ok(false),
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:778-798
// ```c++
// static bool s_LooksLikeNucSeqData(const string& line, size_t minLength=10) {
//     if (line.size()<minLength) {
//         return false;
//     }
//
//     int nucCount=0;
//     for (auto c : line) {
//         if (isalpha(c)) {
//             auto index = static_cast<int>(c);
//             if (symbol_type_table[index] & fDNA_Main_Alphabet) {
//                 ++nucCount;
//             }
//             continue;
//         }
//
//         if (!isspace(c)) {
//             return false;
//         }
//     }
//
//     return (nucCount/line.size() > 0.9);
// }
// ```
/// `s_LooksLikeNucSeqData(line)`: `nucCount/line.size()` divides integers, so every byte is a letter of ACGNTU.
fn looks_like_nuc_seq_data(line: &[u8]) -> bool {
    if line.len() < 10 {
        return false;
    }
    let mut nuc_count = 0usize;
    for &c in line {
        if c.is_ascii_alphabetic() {
            if SYMBOL_TYPE_TABLE[c as usize] & F_DNA_MAIN_ALPHABET != 0 {
                nuc_count += 1;
            }
            continue;
        }
        if !is_space(c) {
            return false;
        }
    }
    (nuc_count / line.len()) as f64 > 0.9
}
