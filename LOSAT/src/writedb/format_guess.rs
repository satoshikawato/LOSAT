//! Port of the decision `makeblastdb` takes in `CMakeBlastDBApp::x_VerifyInputFilesType` (NCBI BLAST+ 2.17.0) when it prints
//! "Input file X does NOT appear to be FASTA (processing anyway)": whether `CFormatGuess::GuessFormat` returns `eFasta` for
//! the input stream. Only the decision is ported, not the warning.
//!
//! The rule, as read from the source:
//!
//! * `makeblastdb` prefers `eBinaryASN`, `eTextASN` and `eFasta` and disables every other format. `GuessFormat` then walks
//!   `sm_CheckOrder`, where `eFasta` comes before `eTextASN` and `eBinaryASN`, and returns the first preferred format whose
//!   test is true. `eFasta` is the first of the three, so the result is `eFasta` exactly when `TestFormatFasta` is true
//!   (a later format can only claim an input that `TestFormatFasta` rejected). The disabled formats are never tried.
//! * `GuessFormat` needs a test buffer first (`EnsureTestBuffer`). Without one it answers `eNewick` or `eUnknown`,
//!   never `eFasta`.
//! * `TestFormatFasta` works on that sample only: the first 8096 bytes of the stream, or more when the sample is all comment
//!   (see `ensure_test_buffer`).
//!
//! C `isalpha`, `isdigit`, `isspace` and `isprint` are the ones of the "C" locale (ASCII only). The C++ works on bytes, so
//! this port works on bytes.
//!
//! Not ported because they do not influence the decision: the `eNewick` test `GuessFormat` makes when there is no test
//! buffer (its answer is never `eFasta`), the brace count `EnsureStats` keeps for `TestFormatTextAsn`, the AGP/other
//! format tests, `CStreamUtils::Stepback` (it only puts the sample back into the stream; the sample is always the first
//! `Multiplier * 8096` bytes of the stream) and the `m_Stream.good()` / `x_TestInput` checks (the stream of a file that
//! exists is good).

use std::io::Read;

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:174-183
// ```c++
// enum ESymbolType {
//     fDNA_Main_Alphabet  = 1<<0, ///< Just ACGTUN-.
//     fDNA_Ambig_Alphabet = 1<<1, ///< Anything else representable in ncbi4na.
//     fProtein_Alphabet   = 1<<2, ///< Allows BZX*-, but not JOU.
//     fLineEnd            = 1<<3,
//     fAlpha              = 1<<4,
//     fDigit              = 1<<5,
//     fSpace              = 1<<6,
//     fInvalid            = 1<<7
// };
// ```
const F_DNA_MAIN_ALPHABET: u8 = 1 << 0;
const F_DNA_AMBIG_ALPHABET: u8 = 1 << 1;
const F_PROTEIN_ALPHABET: u8 = 1 << 2;
const F_LINE_END: u8 = 1 << 3;
const F_ALPHA: u8 = 1 << 4;
const F_DIGIT: u8 = 1 << 5;
const F_SPACE: u8 = 1 << 6;
const F_INVALID: u8 = 1 << 7;

/// C `isspace` in the "C" locale: space, `\t`, `\n`, `\v`, `\f`, `\r`.
fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// C `isprint` in the "C" locale.
fn c_isprint(c: u8) -> bool {
    (0x20..=0x7e).contains(&c)
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:249-286
// ```c++
// static void init_symbol_type_table(void)
// {
//     if ( symbol_type_table[0] == 0 ) {
//         for ( const char* s = "ACGNTU"; *s; ++s ) {
//             int c = *s;
//             symbol_type_table[c] |= fDNA_Main_Alphabet;
//             c = tolower(c);
//             symbol_type_table[c] |= fDNA_Main_Alphabet;
//         }
//         for ( const char* s = "BDHKMRSVWY"; *s; ++s ) {
//             int c = *s;
//             symbol_type_table[c] |= fDNA_Ambig_Alphabet;
//             c = tolower(c);
//             symbol_type_table[c] |= fDNA_Ambig_Alphabet;
//         }
//         for ( const char* s = "ACDEFGHIKLMNPQRSTVWYBZX"; *s; ++s ) {
//             int c = *s;
//             symbol_type_table[c] |= fProtein_Alphabet;
//             c = tolower(c);
//             symbol_type_table[c] |= fProtein_Alphabet;
//         }
//         symbol_type_table[(int)'-'] |= fDNA_Main_Alphabet | fProtein_Alphabet;
//         symbol_type_table[(int)'*'] |= fProtein_Alphabet;
//         for ( const char* s = "\r\n"; *s; ++s ) {
//             int c = *s;
//             symbol_type_table[c] |= fLineEnd;
//         }
//         for ( int c = 1; c < 256; ++c ) {
//             if ( isalpha((unsigned char)c) )
//                 symbol_type_table[c] |= fAlpha;
//             if ( isdigit((unsigned char)c) )
//                 symbol_type_table[c] |= fDigit;
//             if ( isspace((unsigned char)c) )
//                 symbol_type_table[c] |= fSpace;
//         }
//         symbol_type_table[0] |= fInvalid;
//     }
// }
// ```
fn init_symbol_type_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    for &c in b"ACGNTU" {
        table[c as usize] |= F_DNA_MAIN_ALPHABET;
        table[c.to_ascii_lowercase() as usize] |= F_DNA_MAIN_ALPHABET;
    }
    for &c in b"BDHKMRSVWY" {
        table[c as usize] |= F_DNA_AMBIG_ALPHABET;
        table[c.to_ascii_lowercase() as usize] |= F_DNA_AMBIG_ALPHABET;
    }
    for &c in b"ACDEFGHIKLMNPQRSTVWYBZX" {
        table[c as usize] |= F_PROTEIN_ALPHABET;
        table[c.to_ascii_lowercase() as usize] |= F_PROTEIN_ALPHABET;
    }
    table[b'-' as usize] |= F_DNA_MAIN_ALPHABET | F_PROTEIN_ALPHABET;
    table[b'*' as usize] |= F_PROTEIN_ALPHABET;
    for &c in b"\r\n" {
        table[c as usize] |= F_LINE_END;
    }
    for (c, entry) in table.iter_mut().enumerate().skip(1) {
        let b = c as u8;
        if b.is_ascii_alphabetic() {
            *entry |= F_ALPHA;
        }
        if b.is_ascii_digit() {
            *entry |= F_DIGIT;
        }
        if c_isspace(b) {
            *entry |= F_SPACE;
        }
    }
    table[0] |= F_INVALID;
    table
}

/// Reads from `input` until `data` holds `want` bytes, the input ends, or a read fails. A read error counts as end of data.
/// `ended` remembers that the input is at its end (NCBI's stream stays at EOF as well).
fn fill_to<R: Read>(input: &mut R, data: &mut Vec<u8>, want: usize, ended: &mut bool) {
    if *ended || data.len() >= want {
        return;
    }
    let before = data.len();
    let missing = (want - before) as u64;
    // `read_to_end` retries interrupted reads and keeps the bytes read before an error.
    let _ = input.by_ref().take(missing).read_to_end(data);
    if data.len() < want {
        *ended = true;
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3747-3764
// ```c++
// bool
// CFormatGuess::IsAsciiText()
// {
//     const double REQUIRED_ASCII_RATIO = 0.9;
//
//     // first stab - are we text?  comments are only valid if we are text
//     size_t count = 0;
//     size_t count_print = 0;
//     for (int i = 0;  i < m_iTestDataSize;  ++i, ++count) {
//         if (isprint((unsigned char) m_pTestBuffer[i])) {
//             ++count_print;
//         }
//     }
//     if (count_print < (double)count * REQUIRED_ASCII_RATIO) {
//         return false;
//     }
//     return true;
// }
// ```
fn is_ascii_text(data: &[u8]) -> bool {
    const REQUIRED_ASCII_RATIO: f64 = 0.9;
    let count = data.len();
    let count_print = data.iter().filter(|&&c| c_isprint(c)).count();
    if (count_print as f64) < (count as f64) * REQUIRED_ASCII_RATIO {
        return false;
    }
    true
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:2509-2510
// ```c++
//         /// All delimiters are merged and trimmed, to get non-empty tokens only
//         fSplit_Tokenize        = fSplit_MergeDelimiters | fSplit_Truncate
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7439-7545
// ```c++
// bool CStrTokenizeBase::Advance(CTempStringList* part_collector, SIZE_TYPE* ptr_part_start, SIZE_TYPE* ptr_delim_pos)
// {
//     SIZE_TYPE pos, part_start, delim_pos = 0, quote_pos = 0;
//     bool      found_text = false, done;
//     char      active_quote = '\0';
//
//     // Skip leading delimiters.
//     // NOTE: We cannot process
//     if (!m_Pos  &&  (m_Flags & NStr::fSplit_Truncate_Begin) != 0) {
//         SkipDelims();
//     }
//     pos = part_start = m_Pos;
//     done = (pos == NPOS);
//     // save part start position
//     if (ptr_part_start) {
//         *ptr_part_start = part_start;
//     }
//
//     // Checks
//     if (pos >= m_Str.size()) {
//         pos = NPOS;
//         done = true;
//     }
//     if (ptr_delim_pos) {
//         *ptr_delim_pos = NPOS;
//     }
//
//     // Each chunk covers the half-open interval [part_start, delim_pos).
//
//     while ( !done  &&
//             ((delim_pos = m_Str.find_first_of(m_InternalDelim, pos)) != NPOS)) {
//
//         SIZE_TYPE next_start = pos = delim_pos + 1;
//         bool      handled    = false;
//         char      c          = m_Str[delim_pos];
//
//         if ((m_Flags & NStr::fSplit_CanEscape) != 0  &&  c == '\\') {
//             // treat the following character literally
//             if (++pos > m_Str.size()) {
//                 NCBI_THROW2(CStringException, eFormat, "Unescaped trailing \\", delim_pos);
//             }
//             handled = true;
//
//         } else if ((m_Flags & NStr::fSplit_CanQuote) != 0) {
//             if (active_quote != '\0') {
//                 if (c == active_quote) {
//                     if (pos < m_Str.size()  &&  m_Str[pos] == active_quote) {
//                         // count a doubled quote as one literal occurrence
//                         ++pos;
//                     } else {
//                         active_quote = '\0';
//                     }
//                 } else {
//                     continue; // not actually a boundary
//                 }
//                 handled = true;
//             } else if (((m_Flags & NStr::fSplit_CanSingleQuote) != 0  &&  c == '\'') ||
//                        ((m_Flags & NStr::fSplit_CanDoubleQuote) != 0  &&  c == '"')) {
//                 active_quote = c;
//                 quote_pos    = delim_pos;
//                 handled = true;
//             }
//         }
//
//         if ( !handled ) {
//             if ((m_Flags & NStr::fSplit_ByPattern) != 0) {
//                 if (delim_pos + m_Delim.size() <= m_Str.size()
//                     &&  (memcmp(m_Delim.data() + 1, m_Str.data() + pos,
//                                 m_Delim.size() - 1) == 0)) {
//                     done = true;
//                     next_start = pos = delim_pos + m_Delim.size();
//                 } else {
//                     continue;
//                 }
//             } else {
//                 done = true;
//             }
//             // save delimiter position
//             if (ptr_delim_pos) {
//                 *ptr_delim_pos = delim_pos;
//             }
//         }
//
//         if (delim_pos > part_start) {
//             found_text = true;
//             if (part_collector != NULL) {
//                 part_collector->Add(m_Str.substr(part_start, delim_pos - part_start));
//             }
//         }
//         part_start = next_start;
//     }
//
//     if (active_quote != '\0') {
//         NCBI_THROW2(CStringException, eFormat, string("Unbalanced ") + active_quote, quote_pos);
//     }
//
//     if (delim_pos == NPOS) {
//         found_text = true;
//         if (part_collector != NULL) {
//             part_collector->Add(m_Str.substr(part_start));
//         }
//         m_Pos = NPOS;
//     } else {
//         m_Pos = pos;
//         MergeDelims();
//     }
//     return found_text;
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-321
// ```c++
//     void Do(TContainer&     target,
//             TPosContainer&  token_pos,
//             const TString&  empty_str = TString())
//     {
//         auto target_initial_size = target.size();
//
//         // Special cases
//         if (m_Str.empty()) {
//             return;
//         } else if (m_Delim.empty()) {
//             target.push_back(m_Str);
//             token_pos.push_back(0);
//             return;
//         }
//         // Do target space reservation (if applicable)
//         TReserveTrait::Reserve(*this, target, token_pos);
//
//         // Tokenization
//
//         CTempStringList part_collector(m_Storage);
//         SIZE_TYPE prev_pos;
//         SIZE_TYPE delim_pos = NPOS;
//         m_Pos = 0;
//         do {
//             Advance(&part_collector, &prev_pos, &delim_pos);
//             target.push_back(empty_str); // reserve space for value added next line
//             part_collector.Join(&target.back());
//             part_collector.Clear();
//             token_pos.push_back(prev_pos);
//         } while ( !AtEnd() );
//
//         if ( (m_Flags & NStr::fSplit_Truncate_End) == 0 ) {
//             // account training delimiter
//             if (delim_pos != NPOS) {
//                 // add empty token after last delimiter
//                 target.push_back(empty_str);
//                 token_pos.push_back(delim_pos + 1);
//             }
//         }
//         else {
//             // truncate trailing delimiters
//
//             // number of newly added items
//             SIZE_TYPE num_new = target.size() - target_initial_size;
//             // number of empty trailing items to delete
//             SIZE_TYPE num_del = 0;
//             for (auto i = target.rbegin(); i != target.rend() && num_new--; ++i) {
//                 if (!i->empty()) {
//                     break;
//                 }
//                 num_del++;
//             }
//             if (num_del) {
//                 target.resize(target.size()-num_del);
//                 token_pos.resize(token_pos.size()-num_del);
//             }
//         }
//     }
// ```
/// `NStr::Split(data, delims, lines, NStr::fSplit_Tokenize)` with `delims` a set of single characters (`fSplit_ByPattern` is
/// not set, so every character of the delimiter string is a delimiter): the non-empty runs of bytes between delimiters.
/// Leading, trailing and repeated delimiters give no empty token.
fn split_tokenize<F: Fn(u8) -> bool>(data: &[u8], is_delim: F) -> Vec<&[u8]> {
    data.split(|&c| is_delim(c))
        .filter(|token| !token.is_empty())
        .collect()
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3692-3744
// ```c++
// bool
// CFormatGuess::EnsureSplitLines()
// //  ----------------------------------------------------------------------------
// {
//     if ( m_bSplitDone ) {
//         return !m_TestLines.empty();
//     }
//     m_bSplitDone = true;
//
//     //
//     //  Make sure the given data is ASCII before checking potential line breaks:
//     //
//     const size_t MIN_HIGH_RATIO = 20;
//     size_t high_count = 0;
//     for ( streamsize i=0; i < m_iTestDataSize; ++i ) {
//         if ( 0x80 & m_pTestBuffer[i] ) {
//             ++high_count;
//         }
//     }
//     if ( 0 < high_count && m_iTestDataSize / high_count < MIN_HIGH_RATIO ) {
//         return false;
//     }
//
//     //
//     //  Let's expect at least one line break in the given data:
//     //
//     string data( m_pTestBuffer, (size_t)m_iTestDataSize );
//     m_TestLines.clear();
//
//     if ( string::npos != data.find("\r\n") ) {
//         NStr::Split(data, "\r\n", m_TestLines, NStr::fSplit_Tokenize);
//     }
//     else if ( string::npos != data.find("\n") ) {
//         NStr::Split(data, "\n", m_TestLines, NStr::fSplit_Tokenize);
//     }
//     else if ( string::npos != data.find("\r") ) {
//         NStr::Split(data, "\r", m_TestLines, NStr::fSplit_Tokenize);
//     }
//     else if ( m_iTestDataSize == m_iTestBufferSize) {
//         //most likely single truncated line
//         return false;
//     }
//     else {
//         //test buffer contains the entire file
//         m_TestLines.push_back(data);
//     }
//
//     if ( m_iTestDataSize == m_iTestBufferSize   &&  m_TestLines.size() > 1 ) {
//         //multiple lines, last likely truncated
//         m_TestLines.pop_back();
//     }
//     return !m_TestLines.empty();
// }
// ```
/// `CFormatGuess::EnsureSplitLines()` called right after `m_TestLines.clear()` and `m_bSplitDone = false` (as
/// `IsAllComment` does). Returns `m_TestLines`; the C++ return value (`!m_TestLines.empty()`) is not used by the only
/// caller, and an early `return false` leaves `m_TestLines` empty, which is an empty vector here.
fn ensure_split_lines(data: &[u8], buffer_size: usize) -> Vec<&[u8]> {
    const MIN_HIGH_RATIO: usize = 20;
    let high_count = data.iter().filter(|&&c| 0x80 & c != 0).count();
    if 0 < high_count && data.len() / high_count < MIN_HIGH_RATIO {
        return Vec::new();
    }
    let mut lines: Vec<&[u8]>;
    if data.windows(2).any(|w| w == b"\r\n") {
        lines = split_tokenize(data, |c| c == b'\r' || c == b'\n');
    } else if data.contains(&b'\n') {
        lines = split_tokenize(data, |c| c == b'\n');
    } else if data.contains(&b'\r') {
        lines = split_tokenize(data, |c| c == b'\r');
    } else if data.len() == buffer_size {
        // most likely single truncated line
        return Vec::new();
    } else {
        // test buffer contains the entire file
        lines = vec![data];
    }
    if data.len() == buffer_size && lines.len() > 1 {
        // multiple lines, last likely truncated
        lines.pop();
    }
    lines
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3767-3792
// ```c++
// bool
// CFormatGuess::IsAllComment()
// {
//     if (!IsAsciiText()) {
//         return false;
//     }
//
//     m_bSplitDone = false;
//     m_TestLines.clear();
//     EnsureSplitLines();
//
//     ITERATE(list<string>, it, m_TestLines) {
//         if(it->empty()) {
//             continue;
//         }
//         if (NStr::StartsWith(*it, "#")) {
//             continue;
//         }
//         if(NStr::StartsWith(*it, "--")) {
//             continue;
//         }
//         return false;
//     }
//
//     return true;
// }
// ```
fn is_all_comment(data: &[u8], buffer_size: usize) -> bool {
    if !is_ascii_text(data) {
        return false;
    }
    // Lines are not trimmed: " # x" is not a comment. When EnsureSplitLines gives no line (it returned false early),
    // the loop does not run and the sample counts as all comment.
    for line in ensure_split_lines(data, buffer_size) {
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b"#") {
            continue;
        }
        if line.starts_with(b"--") {
            continue;
        }
        return false;
    }
    true
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:625-684
// ```c++
// bool
// CFormatGuess::EnsureTestBuffer()
// //  ----------------------------------------------------------------------------
// {
//     if ( m_pTestBuffer ) {
//         return true;
//     }
//     if ( ! m_Stream.good() ) {
//         return false;
//     }
//
//     // Fix to the all-comment problem.
//     // Read a test buffer,
//     // Test it for being all comment
//     // If its all comment, read a twice as long buffer
//     // Stop when its no longer all comment, end of the stream,
//     //   or Multiplier hits 1024
//
//     const streamsize k_TestBufferGranularity = 8096;
//
//     int Multiplier = 1;
//
//     while(true) {
//         m_iTestBufferSize = Multiplier * k_TestBufferGranularity;
//         m_pTestBuffer = new char[ m_iTestBufferSize ];
//         m_Stream.read( m_pTestBuffer, m_iTestBufferSize );
//         m_iTestDataSize = m_Stream.gcount();
//         if (m_iTestDataSize == 0) {
//             delete[] m_pTestBuffer;
//             m_pTestBuffer = 0;
//             m_iTestBufferSize = 0;
//             return false; //empty file
//         }
//         m_Stream.clear();  // in case we reached eof
//         CStreamUtils::Stepback( m_Stream, m_pTestBuffer, m_iTestDataSize );
//
//         if (IsAllComment()) {
//             if (Multiplier >= 1024)  {
//                 // this is how far we will go and no further.
//                 // if it's indeed all comments then none of the format specific
//                 //  tests will assert.
//                 // if something was misidentified as a comment then the relevant
//                 //  format specific test may still have a good sample to work with.
//                 // so it does not hurt to at least try.
//                 return true;
//             }
//             Multiplier *= 2;
//             delete [] m_pTestBuffer;
//             m_pTestBuffer = NULL;
//             if (m_iTestDataSize < m_iTestBufferSize)  {
//                 return false;
//             }
//             continue;
//         } else {
//             break;
//         }
//     }
//
//     return true;
// }
// ```
/// The sample (`m_pTestBuffer`, `m_iTestDataSize` bytes), or `None` when `EnsureTestBuffer` returns false (empty stream, or
/// an all-comment sample that ended before the buffer was full). Each pass reads the first `Multiplier * 8096` bytes of the stream again (NCBI pushes the bytes back with
/// `CStreamUtils::Stepback`), so the sample is a prefix of the input and bytes are read once here.
fn ensure_test_buffer<R: Read>(input: &mut R) -> Option<Vec<u8>> {
    const TEST_BUFFER_GRANULARITY: usize = 8096;
    let mut data: Vec<u8> = Vec::new();
    let mut ended = false;
    let mut multiplier: usize = 1;
    loop {
        let buffer_size = multiplier * TEST_BUFFER_GRANULARITY;
        fill_to(input, &mut data, buffer_size, &mut ended);
        let data_size = data.len();
        if data_size == 0 {
            return None; // empty file
        }
        if is_all_comment(&data, buffer_size) {
            if multiplier >= 1024 {
                return Some(data);
            }
            multiplier *= 2;
            if data_size < buffer_size {
                return None;
            }
            continue;
        }
        return Some(data);
    }
}

/// What `CFormatGuess::EnsureStats` counts and `TestFormatFasta` reads.
struct Stats {
    /// `m_iStatsCountAlNumChars`
    al_num_chars: usize,
    /// `m_iStatsCountData`
    data: usize,
    /// `m_iStatsCountDnaChars`
    dna_chars: usize,
    /// `m_iStatsCountAaChars`
    aa_chars: usize,
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:687-749
// ```c++
// bool
// CFormatGuess::EnsureStats()
// //  ----------------------------------------------------------------------------
// {
//     if ( m_bStatsAreValid ) {
//         return true;
//     }
//     if ( ! EnsureTestBuffer() ) {
//         return false;
//     }
//
//     string strBuffer(m_pTestBuffer, m_iTestDataSize);
//     CNcbiIstrstream TestBuffer(strBuffer);
//     string strLine;
//
//     init_symbol_type_table();
//     // Things we keep track of:
//     //   m_iStatsCountAlNumChars: number of characters that are letters or
//     //     digits
//     //   m_iStatsCountData: number of characters not part of a line starting
//     //     with '>', ignoring whitespace
//     //   m_iStatsCountDnaChars: number of characters counted in m_iStatsCountData
//     //     from the DNA alphabet
//     //   m_iStatsCountAaChars: number of characters counted in m_iStatsCountData
//     //     from the AA alphabet
//     //  m_iStatsCountBraces: Opening { and closing } braces
//     //
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
//             if ( !is_header ) {
//                 if ( !(type & fSpace) ) {
//                     ++m_iStatsCountData;
//                 }
//
//                 if ( type & fDNA_Main_Alphabet ) {
//                     ++m_iStatsCountDnaChars;
//                 }
//                 if ( type & fProtein_Alphabet ) {
//                     ++m_iStatsCountAaChars;
//                 }
//             }
//         }
//     }
//     m_bStatsAreValid = true;
//     return true;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:55-134
// ```c++
// CNcbiIstream& NcbiGetline(CNcbiIstream& is, string& str, const string& delims,
//                           SIZE_TYPE* count)
// {
//     str.erase();
//
//     IOS_BASE::fmtflags f = is.flags();
//     is.unsetf(IOS_BASE::skipws);
// #ifdef NO_PUBSYNC
//     if ( !is.ipfx(1) ) {
//         is.flags(f);
//         is.setstate(NcbiFailbit);
//         return is;
//     }
// #else
//     CNcbiIstream::sentry s(is);
//     if ( !s ) {
//         is.flags(f);
//         is.setstate(NcbiFailbit);
//         return is;
//     }
// #endif //NO_PUBSYNC
//     _ASSERT( is.good() );
//
//     char buf[1024];
//     SIZE_TYPE pos = 0;
//     SIZE_TYPE size = 0;
//     SIZE_TYPE max_size = str.max_size();
//     SIZE_TYPE delim_count = 0;
//     IOS_BASE::iostate iostate = NcbiGoodbit/*0*/;
//     for (;;) {
//         CT_INT_TYPE ch = is.rdbuf()->sbumpc();
//         if ( CT_EQ_INT_TYPE(ch, CT_EOF) ) {
//             iostate = NcbiEofbit;
//             break;
//         }
//         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
//         if (delim_pos != NPOS) {
//             // Special case -- if two different delimiters are back to
//             // back and in the same order as in delims, treat them as
//             // a single delimiter (necessary for correct handling of
//             // DOS/MAC-style CR/LF endings).
//             ch = is.rdbuf()->sgetc();
//             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
//                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
//                 is.rdbuf()->sbumpc();
//                 delim_count = 2;
//             } else {
//                 delim_count = 1;
//             }
//             break;
//         }
//         if (size == max_size) {
//             CT_INT_TYPE bk = is.rdbuf()->sungetc();
//             iostate = CT_EQ_INT_TYPE(bk, ch) ? NcbiFailbit : NcbiBadbit;
//             break;
//         }
//
//         buf[pos++] = CT_TO_CHAR_TYPE(ch);
//         if (pos == sizeof(buf)) {
//             str.append(buf, pos);
//             pos  = 0;
//         }
//         size++;
//     }
//     if (pos > 0)
//         str.append(buf, pos);
//     if (count != NULL)
//         *count = size + delim_count;
//
// #ifdef NO_PUBSYNC
//     is.isfx();
// #endif //NO_PUBSYNC
//     is.flags(f);
//     if (iostate) {
//         if (iostate == NcbiEofbit  &&  str.empty())
//             iostate |= NcbiFailbit;
//         is.clear(iostate);
//     }
//     return is;
// }
// ```
/// `EnsureStats` on the sample. `NcbiGetline(..., "\r\n")` ends a line at `\r` or `\n`; a `\r` followed by `\n` is one
/// delimiter and every other pair (`\n\r`, `\r\r`, `\n\n`) is two, so an empty line only adds an empty `strLine`, which
/// `EnsureStats` does not count. The lines are therefore the non-empty runs between `\r`/`\n` bytes, each with one `\n`
/// appended; the loop ends at the first read that finds nothing left (`fail()`), after the last line.
fn ensure_stats(data: &[u8]) -> Stats {
    let symbol_type_table = init_symbol_type_table();
    let mut stats = Stats {
        al_num_chars: 0,
        data: 0,
        dna_chars: 0,
        aa_chars: 0,
    };
    for line in split_tokenize(data, |c| c == b'\r' || c == b'\n') {
        let is_header = line[0] == b'>';
        // `strLine += '\n'`
        for &c in line.iter().chain(std::iter::once(&b'\n')) {
            let type_ = symbol_type_table[c as usize];
            if type_ & (F_ALPHA | F_DIGIT | F_SPACE) != 0 {
                stats.al_num_chars += 1;
            }
            if !is_header {
                if type_ & F_SPACE == 0 {
                    stats.data += 1;
                }
                if type_ & F_DNA_MAIN_ALPHABET != 0 {
                    stats.dna_chars += 1;
                }
                if type_ & F_PROTEIN_ALPHABET != 0 {
                    stats.aa_chars += 1;
                }
            }
        }
    }
    stats
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3145-3180
// ```c++
// template <class TStr>
// TStr s_TruncateSpaces(const TStr& str, NStr::ETrunc where,
//                       const TStr& empty_str)
// {
//     SIZE_TYPE length = str.length();
//     if (length == 0) {
//         return empty_str;
//     }
//     SIZE_TYPE beg = 0;
//     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
//         _ASSERT(beg < length);
//         while ( isspace((unsigned char) str[beg]) ) {
//             if (++beg == length) {
//                 return empty_str;
//             }
//         }
//     }
//     SIZE_TYPE end = length;
//     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
//         _ASSERT(beg < end);
//         while (isspace((unsigned char) str[--end])) {
//             if (beg == end) {
//                 return empty_str;
//             }
//         }
//         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
//         ++end;
//     }
//     _ASSERT(beg < end  &&  end <= length);
//     if ( beg | (end - length) ) { // if either beg != 0 or end != length
//         return str.substr(beg, end - beg);
//     }
//     else {
//         return str;
//     }
// }
// ```
/// `NStr::TruncateSpaces_Unsafe(text, NStr::eTrunc_Begin)`: drops leading `isspace` bytes (all of them gives the empty string).
fn truncate_spaces_begin(text: &[u8]) -> &[u8] {
    let skip = text.iter().take_while(|&&c| c_isspace(c)).count();
    &text[skip..]
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1624-1637
// ```c++
// void SkipCommentAndBlank(CTempString &text)
// {
//     const CTempString COMMENT_SYMBOLS(";#!");
//     const CTempString NEW_LINE_SYMBOLS("\r\n");
//     while (true) {
//         text = NStr::TruncateSpaces_Unsafe(text, NStr::eTrunc_Begin);
//         if ( COMMENT_SYMBOLS.find(text[0]) != CTempString::npos ) {
//             CTempString::size_type pos = text.find_first_of(NEW_LINE_SYMBOLS, 1);
//             text = text.substr(pos);
//         } else {
//             break;
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:340-347
// ```c++
// inline
// char CTempString::operator[](size_type pos) const
// {
//     if ( pos < m_Length ) {
//         return m_String[pos];
//     }
//     return '\0';
// }
// ```
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:537-549
// ```c++
// inline
// CTempString::size_type CTempString::find_first_of(const CTempString match,
//                                                   size_type pos) const
// {
//     if (match.length()  &&  pos < length()) {
//         const_iterator it = std::find_first_of(begin() + pos, end(),
//                                                match.begin(), match.end());
//         if (it != end()) {
//             return it - begin();
//         }
//     }
//     return npos;
// }
// ```
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:362-373
// ```c++
// inline
// void CTempString::x_Init(const char* str, size_type str_len,
//                          size_type pos)
// {
//     if ( pos >= str_len ) {
//         clear();
//     }
//     else {
//         m_String = str + pos;
//         m_Length = str_len - pos;
//     }
// }
// ```
/// `SkipCommentAndBlank`. `text[0]` of an empty text is NUL, which is not in `;#!`. A comment line without a following
/// line end (`find_first_of` gives `npos`) makes `substr(npos)` the empty string.
fn skip_comment_and_blank(mut text: &[u8]) -> &[u8] {
    loop {
        text = truncate_spaces_begin(text);
        let first = text.first().copied().unwrap_or(0);
        if matches!(first, b';' | b'#' | b'!') {
            // `find_first_of(NEW_LINE_SYMBOLS, 1)`: search from index 1
            let pos = text
                .iter()
                .skip(1)
                .position(|&c| c == b'\r' || c == b'\n')
                .map(|p| p + 1);
            text = match pos {
                Some(p) => &text[p..],
                None => &[],
            };
        } else {
            break;
        }
    }
    text
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1640-1675
// ```c++
// CFormatGuess::TestFormatFasta( EMode /* not used */ )
// {
//     if ( ! EnsureStats() ) {
//         return false;
//     }
//
//     // reject obvious misfits:
//     CTempString header(m_pTestBuffer, m_iTestDataSize);
//     SkipCommentAndBlank(header);
//     if ( m_iTestDataSize == 0 || header.length() == 0 || header[0] != '>' ) {
//         return false;
//     }
//     if ( m_iStatsCountData == 0 ) {
//         if (0.75 > double(m_iStatsCountAlNumChars)/double(m_iTestDataSize) ) {
//             return false;
//         }
//         return ( NStr::Find( m_pTestBuffer, "|" ) <= 10 );
//     }
//
//     // remaining decision based on text stats:
//     double dAlNumFraction = (double)m_iStatsCountAlNumChars / (double)m_iTestDataSize;
//     double dDnaFraction   = (double)m_iStatsCountDnaChars / (double)m_iStatsCountData;
//     double dAaFraction    = (double)m_iStatsCountAaChars / (double)m_iStatsCountData;
//
//     // want at least 80% text-ish overall:
//     if ( dAlNumFraction < 0.8 ) {
//         return false;
//     }
//
//     // want more than 91 percent of either DNA content or AA content in what we
//     // presume is data:
//     if ( dDnaFraction > 0.91 || dAaFraction > 0.91 ) {
//         return true;
//     }
//     return false;
// }
// ```
fn test_format_fasta(sample: &[u8], stats: &Stats) -> bool {
    // reject obvious misfits:
    let header = skip_comment_and_blank(sample);
    if sample.is_empty() || header.is_empty() || header[0] != b'>' {
        return false;
    }
    if stats.data == 0 {
        if 0.75 > (stats.al_num_chars as f64) / (sample.len() as f64) {
            return false;
        }
        // NCBI calls `NStr::Find(m_pTestBuffer, "|")` with a `char*` that is not NUL-terminated at the end of the sample,
        // so the C string runs past the sample into memory that was never written (undefined behaviour). This port does not
        // reproduce that: it searches the sample's own bytes, stopping at the first NUL byte or at the end of the sample.
        let searched = match sample.iter().position(|&c| c == 0) {
            Some(nul) => &sample[..nul],
            None => sample,
        };
        return searched
            .iter()
            .position(|&c| c == b'|')
            .is_some_and(|i| i <= 10);
    }

    // remaining decision based on text stats:
    let al_num_fraction = stats.al_num_chars as f64 / sample.len() as f64;
    let dna_fraction = stats.dna_chars as f64 / stats.data as f64;
    let aa_fraction = stats.aa_chars as f64 / stats.data as f64;

    // want at least 80% text-ish overall:
    if al_num_fraction < 0.8 {
        return false;
    }

    // want more than 91 percent of either DNA content or AA content in what we
    // presume is data:
    if dna_fraction > 0.91 || aa_fraction > 0.91 {
        return true;
    }
    false
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:841-850
// ```c++
// CMakeBlastDBApp::TFormat
// CMakeBlastDBApp::x_GuessFileType(CNcbiIstream & input)
// {
//     CFormatGuess fg(input);
//     fg.GetFormatHints().AddPreferredFormat(CFormatGuess::eBinaryASN);
//     fg.GetFormatHints().AddPreferredFormat(CFormatGuess::eTextASN);
//     fg.GetFormatHints().AddPreferredFormat(CFormatGuess::eFasta);
//     fg.GetFormatHints().DisableAllNonpreferred();
//     return fg.GuessFormat();
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:452-491
// ```c++
// CFormatGuess::EFormat
// CFormatGuess::GuessFormat( EOnError onerror )
// {
//     //sqd-4036:
//     // make sure we got something to work with
//     //
//     if (!x_TestInput(m_Stream, onerror)) {
//         return eUnknown;
//     }
//     if (!EnsureTestBuffer()) {
//         //one condition that won't allow us to get a good test buffer is an ascii
//         // file without any line breaks. so before giving up, let's specifically
//         // try any formats that would allow for that:
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
/// `true` when makeblastdb takes the input for FASTA (no warning): `x_GuessFileType` returns `CFormatGuess::eFasta`.
///
/// The preferred formats are tried in the order of `sm_CheckOrder`, where `eFasta` precedes `eTextASN` and `eBinaryASN`,
/// and all other formats are disabled, so the answer is `TestFormatFasta` on the test buffer (false when there is none).
/// A read error counts as end of data.
pub fn appears_fasta<R: Read>(mut input: R) -> bool {
    let Some(buffer) = ensure_test_buffer(&mut input) else {
        return false;
    };
    let stats = ensure_stats(&buffer);
    test_format_fasta(&buffer, &stats)
}

#[cfg(test)]
mod tests {
    use super::appears_fasta;
    use std::io::{self, Read};

    // Every expected value below is the answer of makeblastdb 2.17.0 (`-dbtype prot` and `nucl` agree): it printed
    // "does NOT appear to be FASTA" for `false` and not for `true`. `name` is the case of
    // db2/fg/cases/<name>.fa that holds the same bytes (made by scripts/db2_fg_gen.py).

    fn check(name: &str, input: &[u8], want: bool) {
        assert_eq!(
            appears_fasta(input),
            want,
            "case {name} ({} bytes)",
            input.len()
        );
    }

    #[test]
    fn small_inputs() {
        let cases: &[(&str, &[u8], bool)] = &[
            ("two_records", b">a\nACGT\n>b\nTTTT\n", true),
            ("protein", b">p\nMKVLAAGIVGLLLAQPAMA\nMKV\n", true),
            ("no_final_newline", b">s\nACGTACGTACGT", true),
            ("crlf", b">s\r\nACGTACGT\r\nACGTACGT\r\n", true),
            ("cr_only", b">s\rACGTACGT\rACGTACGT\r", true),
            ("lfcr", b">s\n\rACGTACGT\n\rACGTACGT\n\r", true),
            ("mixed_ends", b">s\rACGTACGT\r\nACGTACGT\nACGTACGT\n", true),
            ("seq_without_header", b"ACGTACGTACGT\nACGTACGT\n", false),
            ("plain_text", b"hello world, this is not a sequence file\nsecond line\n", false),
            ("fastq", b"@r1\nACGT\n+\nIIII\n", false),
            ("gff", b"##gff-version 3\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n", false),
            ("text_asn", b"Seq-entry ::= seq {\n  id { local str \"a\" } ,\n  inst { repr raw, mol na, length 4 } }\n", false),
            ("leading_space_header_inline", b" >s\nACGTACGT\n", false),
            ("leading_space_header", b"  \n\t>s\nACGTACGT\n", false),
            ("hash_comment_then_fasta", b"#ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
            ("semicolon_comment_then_fasta", b";ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
            ("bang_comment_then_fasta", b"!ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
            ("comment_two_then_fasta", b"#ACGTACGTACGT\n;ACGTACGTACGT\n   >s\nACGT\n", false),
            ("dash_comment_then_fasta", b"--ACGTACGTACGTACGTACGT\n>s\nACGTACGT\n", false),
            ("comment_nonl_only", b"#ACGTACGT", false),
            ("nul_in_seq", b">s\nACGT\x00ACGT\nACGTACGT\n", true),
            ("nul_first", b"\x00>s\nACGT\nACGT\n", false),
            ("bom_fasta", b"\xef\xbb\xbf>s\nACGT\nACGT\n", false),
            ("latin1_header_fasta", b">caf\xe9 s\nACGTACGTACGT\nACGTACGT\n", true),
            ("fasta_header_gt_in_seq", b">s\nAC>GT\nACGT\n", false),
            ("only_gt", b">", false),
            ("gt_bar", b">|", false),
            ("tiny_abc", b">abc", false),
            ("lonely_nl", b"\n", false),
            ("nul_only_file", &[0u8; 100], false),
            // defline-only files: the '|' must be at index <= 10 of the sample
            ("defline_bar_at_9", b">aaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", true),
            ("defline_bar_at_10", b">aaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", true),
            ("defline_bar_at_11", b">aaaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", false),
            ("defline_bar_at_10_crlf", b">aaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\r\n", true),
            ("defline_bar_at_11_nonl", b">aaaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", false),
            ("defline_no_bar", b">aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n", false),
            ("defline_ten_chars_bar", b">sp|P12345|FOO_HUMAN something\n", true),
            ("defline_bar_second_line_early", b">aaaaaaaaaaaaaaa\n>b|bb\n", false),
            ("nul_defline_only_before_bar", b">ab\x00cd|ef\n", false),
            ("empty_seq_spaces", b">s\n \n", false),
        ];
        assert!(cases.len() >= 30);
        for &(name, input, want) in cases {
            check(name, input, want);
        }
    }

    /// There is no test buffer for an empty stream, so `GuessFormat` cannot answer `eFasta`. (`makeblastdb` stops with
    /// "is empty" before it asks, so there is no oracle answer for this.)
    #[test]
    fn empty_input_is_not_fasta() {
        assert!(!appears_fasta(&b""[..]));
    }

    // The generators below build the same bytes as scripts/db2_fg_gen.py.

    fn repeat_cut(unit: &[u8], n: usize) -> Vec<u8> {
        let mut v = unit.repeat(n / unit.len() + 1);
        v.truncate(n);
        v
    }

    /// `n` bytes of `#AAAA...` comment lines (DNA-looking, so a later FASTA part still passes the 0.91 test).
    fn comment_lines(n: usize, sep: &[u8]) -> Vec<u8> {
        let mut line = b"#".to_vec();
        line.extend_from_slice(&[b'A'; 79]);
        line.extend_from_slice(sep);
        repeat_cut(&line, n)
    }

    fn dna_lines(n: usize, sep: &[u8]) -> Vec<u8> {
        let mut line = b"ACGT".repeat(15);
        line.extend_from_slice(sep);
        repeat_cut(&line, n)
    }

    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    fn fasta_of_size(n: usize) -> Vec<u8> {
        cat(&[b">s\n", &dna_lines(n - 3, b"\n")])
    }

    /// `comment_lines(n)` with `k` evenly spread body bytes (not the '#', not the line end) replaced by `byte`.
    fn comment_marked(n: usize, k: usize, byte: u8) -> Vec<u8> {
        let mut b = comment_lines(n, b"\n");
        let pos: Vec<usize> = (0..b.len())
            .filter(|&i| b[i] != b'\n' && b[i] != b'#')
            .collect();
        for j in 0..k {
            b[pos[j * pos.len() / k]] = byte;
        }
        b
    }

    /// FASTA of `n_total` bytes where `count` evenly spread sequence bytes are replaced by `junk`.
    fn fasta_with(n_total: usize, junk: u8, count: usize) -> Vec<u8> {
        let mut body = cat(&[b">s\n", &dna_lines(n_total - 3, b"\n")]);
        let pos: Vec<usize> = (3..body.len()).filter(|&i| body[i] != b'\n').collect();
        for j in 0..count {
            body[pos[j * pos.len() / count]] = junk;
        }
        body
    }

    /// A header and 1000 data characters of which `good` are `letter` and the others the digit `1`, 50 to a line.
    fn frac_case(good: usize, letter: u8) -> Vec<u8> {
        let mut chars = vec![letter; good];
        chars.extend_from_slice(&vec![b'1'; 1000 - good]);
        let lines: Vec<&[u8]> = chars.chunks(50).collect();
        cat(&[b">s\n", &lines.join(&b'\n'), b"\n"])
    }

    /// A header of `%` and `>` and `n - header` bytes of DNA: `alnum` alphanumeric-or-space characters in `n` bytes.
    fn alnum_case(alnum: usize) -> Vec<u8> {
        let n = 1000;
        let mut header = b">".to_vec();
        header.extend_from_slice(&vec![b'%'; n - alnum - 1]);
        header.push(b'\n');
        cat(&[&header, &dna_lines(n - header.len(), b"\n")])
    }

    /// A defline-only file of 400 bytes with `alnum` alphanumeric-or-space characters (the '|' is at index 1).
    fn defline_alnum_case(alnum: usize) -> Vec<u8> {
        let mut v = b">|".to_vec();
        v.extend_from_slice(&vec![b'%'; 400 - alnum - 2]);
        v.extend_from_slice(&vec![b'x'; alnum - 1]);
        v.push(b'\n');
        v
    }

    #[test]
    fn sample_size_rule() {
        // the sample is the first 8096 bytes (more only for all-comment samples)
        for n in [8095, 8096, 8097, 8098, 16191, 16192, 16193] {
            check(&format!("size_{n}"), &fasta_of_size(n), true);
        }
        let junk = b"%%%%%%%%%%%%\n".repeat(5000);
        check(
            "head_fasta_tail_junk",
            &cat(&[&fasta_of_size(8096), &junk]),
            true,
        );
        check(
            "head_fasta_tail_junk_8000",
            &cat(&[&fasta_of_size(8000), &junk]),
            true,
        );
        let junk_head = b"%%%%%%%%%%%%\n".repeat(700);
        check(
            "head_junk_tail_fasta",
            &cat(&[&junk_head, &fasta_of_size(100000)]),
            false,
        );
        let binary: Vec<u8> = (0..=255u8).collect::<Vec<_>>().repeat(200);
        check(
            "head_fasta_tail_binary",
            &cat(&[&fasta_of_size(8096), &binary]),
            true,
        );
        let binary_short: Vec<u8> = (1..=255u8).collect::<Vec<_>>().repeat(3);
        check(
            "head_fasta_tail_binary_short",
            &cat(&[&fasta_of_size(4000), &binary_short]),
            false,
        );
    }

    #[test]
    fn all_comment_sample_doubles() {
        // only comments: false at every size (the buffer doubles until the data ends before the buffer does)
        for n in [
            8095, 8096, 8097, 9000, 16191, 16192, 16193, 17000, 33000, 70000,
        ] {
            check(&format!("allcomment_{n}"), &comment_lines(n, b"\n"), false);
        }
        check(
            "allcomment_crlf_20000",
            &comment_lines(20000, b"\r\n"),
            false,
        );
        check("allcomment_cr_20000", &comment_lines(20000, b"\r"), false);
        let dashes = (b"--".iter().chain(&[b'A'; 77]).chain(b"\n"))
            .copied()
            .collect::<Vec<u8>>();
        check("allcomment_dash_20000", &dashes.repeat(250), false);
        // comments, then FASTA: found as soon as the doubled buffer reaches it
        let tail = cat(&[b"\n>s\n", &dna_lines(30000, b"\n")]);
        for (c, want) in [
            (8000, true),
            (8100, true),
            (9000, true),
            (16300, true),
            (17000, true),
            (40000, true),
            (4145152, true),
            (8000000, true),
            (8290300, true),
            (8290400, false),
            (8400000, false),
        ] {
            check(
                &format!("comment_{c}_then_fasta"),
                &cat(&[&comment_lines(c, b"\n"), &tail]),
                want,
            );
        }
        // the sample ends inside the comments: the cap of 1024 * 8096 bytes
        check("comment_4145152", &comment_lines(4145152, b"\n"), false);
        check("comment_8300000", &comment_lines(8300000, b"\n"), false);
        // CR and CR LF comment lines
        let cr_tail = cat(&[b">s\r", &dna_lines(30000, b"\r")]);
        check(
            "comment_9000_then_fasta_cr",
            &cat(&[&comment_lines(9000, b"\r"), &cr_tail]),
            false,
        );
        let crlf_tail = cat(&[b">s\r\n", &dna_lines(30000, b"\r\n")]);
        check(
            "comment_9000_then_fasta_crlf",
            &cat(&[&comment_lines(9000, b"\r\n"), &crlf_tail]),
            false,
        );
        check(
            "comment_cr_then_fasta",
            b"#ACGTACGTACGT\r>s\rACGT\rACGT\r",
            true,
        );
    }

    #[test]
    fn printable_and_high_byte_limits() {
        let tail = cat(&[b"\n>s\n", &dna_lines(100000, b"\n")]);
        // at most 10% of the sample may be non-printable (the line ends count) for the comment test to apply
        for (c, want) in [
            (600, true),
            (709, true),
            (710, true),
            (711, false),
            (800, false),
        ] {
            check(
                &format!("ctrl_comment8096_{c}_then_fasta"),
                &cat(&[&comment_marked(8096, c, 1), &tail]),
                want,
            );
        }
        for (c, want) in [(1000, true), (1421, true), (1422, false), (1700, false)] {
            let second = comment_marked(8096, c, 1);
            check(
                &format!("ctrl_comment_2nd_{c}_then_fasta"),
                &cat(&[&comment_lines(8096, b"\n"), &second, &tail]),
                want,
            );
        }
        // more than 1 byte in 20 >= 0x80 makes EnsureSplitLines give up, and the sample counts as all comment
        for (h, want) in [
            (0, true),
            (200, true),
            (404, true),
            (405, true),
            (450, true),
            (451, false),
            (500, false),
        ] {
            check(
                &format!("hi_fasta_9000_{h}"),
                &fasta_with(9000, 0xE9, h),
                want,
            );
        }
        for (h, want) in [(350, true), (404, true), (405, false), (810, false)] {
            check(
                &format!("hi_fasta_8096_{h}"),
                &fasta_with(8096, 0xE9, h),
                want,
            );
        }
        for (h, want) in [(100, true), (404, false), (500, false)] {
            check(
                &format!("hi_fasta_5000_{h}"),
                &fasta_with(5000, 0xE9, h),
                want,
            );
        }
        check(
            "high_in_seq_few",
            &cat(&[
                b">s\n",
                &b"ACGT".repeat(100),
                b"\xff",
                &b"ACGT".repeat(100),
                b"\n",
            ]),
            true,
        );
    }

    #[test]
    fn fraction_limits() {
        // alphanumeric-or-space characters / sample size must be >= 0.8
        for (a, want) in [
            (780, false),
            (798, false),
            (799, true),
            (800, true),
            (810, true),
        ] {
            check(&format!("alnum_{a}_of_1000"), &alnum_case(a), want);
        }
        // defline-only files: the same ratio must be >= 0.75
        for (a, want) in [(290, false), (299, false), (300, true), (301, true)] {
            check(
                &format!("defline_alnum_{a}_of_400"),
                &defline_alnum_case(a),
                want,
            );
        }
        // more than 91% of the data characters of DNA (ACGNTU-) or of protein letters
        for (g, want) in [
            (880, false),
            (909, false),
            (910, false),
            (911, true),
            (920, true),
        ] {
            check(&format!("dna_{g}_of_1000"), &frac_case(g, b'G'), want);
            check(&format!("aa_{g}_of_1000"), &frac_case(g, b'L'), want);
        }
        // J and O are not protein letters; '*' and '-' are not alphanumeric, so the 0.8 test fails first
        for (letter, name) in [(b'J', "junkletter"), (b'O', "o"), (b'*', "star")] {
            check(
                &format!("{name}_911_of_1000"),
                &frac_case(911, letter),
                false,
            );
        }
        check("n_911_of_1000", &frac_case(911, b'N'), true);
        check("x_911_of_1000", &frac_case(911, b'X'), true);
        check("dash_911_of_1000", &frac_case(911, b'-'), false);
    }

    /// A reader that gives its data in pieces of `chunk` bytes, then fails.
    struct Chunked<'a> {
        data: &'a [u8],
        chunk: usize,
        fail_after_data: bool,
        interrupt_next: bool,
    }

    impl Read for Chunked<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.interrupt_next {
                self.interrupt_next = false;
                return Err(io::ErrorKind::Interrupted.into());
            }
            self.interrupt_next = true;
            if self.data.is_empty() {
                return if self.fail_after_data {
                    Err(io::Error::other("broken"))
                } else {
                    Ok(0)
                };
            }
            let n = self.data.len().min(self.chunk).min(buf.len());
            buf[..n].copy_from_slice(&self.data[..n]);
            self.data = &self.data[n..];
            Ok(n)
        }
    }

    #[test]
    fn reader_behaviour_does_not_matter() {
        let data = cat(&[
            &comment_lines(9000, b"\n"),
            b"\n>s\n",
            &dna_lines(3000, b"\n"),
        ]);
        for chunk in [1, 7, 4096, 100000] {
            for fail in [false, true] {
                let reader = Chunked {
                    data: &data,
                    chunk,
                    fail_after_data: fail,
                    interrupt_next: true,
                };
                assert!(appears_fasta(reader), "chunk {chunk} fail {fail}");
            }
        }
        // a read error is the end of the data: the same answer as for the shorter input
        let short = &data[..9004];
        let reader = Chunked {
            data: short,
            chunk: 100,
            fail_after_data: true,
            interrupt_next: false,
        };
        assert_eq!(appears_fasta(reader), appears_fasta(short));
    }

    /// Reads every file of the TSV named by `LOSAT_FORMAT_GUESS_CHECK` (`path<TAB>fasta|other`, `other` = makeblastdb
    /// warned) and compares. The TSV is built from makeblastdb 2.17.0's answers by `scripts/db2_fg_check.py`.
    #[test]
    #[ignore]
    fn agrees_with_makeblastdb_on_listed_files() {
        let list = std::env::var("LOSAT_FORMAT_GUESS_CHECK")
            .expect("LOSAT_FORMAT_GUESS_CHECK must name a TSV file (path<TAB>fasta|other)");
        let text = std::fs::read_to_string(&list).expect("read the TSV");
        let mut disagree: Vec<String> = Vec::new();
        let mut total = 0usize;
        for line in text.lines().filter(|l| !l.is_empty()) {
            let (path, answer) = line.split_once('\t').expect("path<TAB>fasta|other");
            let want = match answer {
                "fasta" => true,
                "other" => false,
                other => panic!("bad answer {other:?} in {line:?}"),
            };
            let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
            let got = appears_fasta(io::BufReader::new(file));
            total += 1;
            if got != want {
                disagree.push(format!(
                    "{path}: makeblastdb {answer}, port {}",
                    if got { "fasta" } else { "other" }
                ));
            }
        }
        println!("{total} files, {} disagreements", disagree.len());
        for d in &disagree {
            println!("DISAGREE {d}");
        }
        assert!(disagree.is_empty());
    }
}
