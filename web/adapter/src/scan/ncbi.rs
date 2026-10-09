//! Parser kinds 1 and 2 of the index scan (docs/web/abi_v2.md §9): the records that NCBI
//! BLAST+'s FASTA input source reads, as the engine reads them
//! (`LOSAT::blastinput::fasta_reader`, the port of `CBlastFastaInputSource`,
//! `CFastaReader` and `CStreamLineReader`), from the bytes of a file
//! (`FastaInputSource::from_bytes`), with the data loaders on (the CLI's default).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-345
//! ```c++
//!     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//!                                     CFastaReader::fParseRawID:
//!                                     (CFastaReader::fNoParseID |
//!                                      CFastaReader::fDLOptional);
//!     flags += (m_ReadProteins
//!               ? CFastaReader::fAssumeProt
//!               : CFastaReader::fAssumeNuc);
//!     const char* env_var = getenv("BLASTINPUT_GEN_DELTA_SEQ");
//!     if (env_var == NULL || (env_var && string(env_var) == kEmptyStr)) {
//!         flags += CFastaReader::fNoSplit;
//!     }
//!     flags+= CFastaReader::fHyphensIgnoreAndWarn;
//!     flags+= CFastaReader::fDisableNoResidues;
//!     flags+= CFastaReader::fQuickIDCheck;
//!     if (m_Config.GetDataLoaderConfig().UseDataLoaders()) {
//! ```
//! Kind 1 has the flags of nucleotide input (`fAssumeNuc`: BLASTN, TBLASTX, the TBLASTN
//! subject), kind 2 those of protein input (`fAssumeProt`: BLASTP, the TBLASTN query).
//!
//! The scan takes the input in chunks of any size and holds only a header's first word,
//! the first line of an input whose first line may be a Seq-id, at most 70 bytes of a
//! record's first data line and at most 257 bytes whose loss is not yet decided. It writes
//! no messages. Each record is reported with every record that the reader returns
//! (records without residues too): the ID is the title up to its first space (empty without
//! a title), the length is the number of stored residues, and `residue_counts` counts the
//! stored residues upper-cased (`U` is `T` for kind 1, as the reader stores it).
//!
//! The scan fails with the reader's error (`CheckDataLine`) and rejects explicitly, with
//! `not supported by LOSAT Web`:
//! - a first line that NCBI's data loaders would read as a Seq-id (the engine's test);
//! - a `>?` gap line (maintainer decision 3: gap residues have no bytes to index);
//! - a record longer than 2147483647 letters (the engine's limit);
//! - a record whose residues neither follow the uniform layout nor the read-forward rules
//!   (`Forward`), which happens only where NCBI's line reader joins two lines of the file
//!   (a CR line end inside a file of LF line ends, or the reverse).

use super::{LineLayout, ScanRecord, CHECKPOINT_EVERY};

/// Bytes of a regular file that an `ifstream` buffers at once (`FastaStream::new` in
/// `LOSAT/src/blastinput/fasta_reader/stream.rs`, checked against an independent C++ oracle).
const FILE_WINDOW: u64 = 8191;

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:117
// ```c++
// const size_t CPushback_Streambuf::kMinBufSize = 4096;
// ```
const MIN_BUF_SIZE: u64 = 4096;

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-417
// ```c++
//         if (how == ePushback_Stepback
//             ||  (how == ePushback_Copy
//                  &&  buf_size <= (del_ptr
//                                   ? CPushback_Streambuf::kMinBufSize
//                                   : CPushback_Streambuf::kMinBufSize >> 4))) {
// ```
const STEP_BACK_MAX: usize = (MIN_BUF_SIZE >> 4) as usize;

/// The longest record (`check_length` of the engine's reader).
const MAX_LETTERS: u64 = i32::MAX as u64;

/// `CheckDataLine`'s window.
const CHECK_WINDOW: usize = 70;

/// What a byte of a data line is (`ParseDataLine`'s character types for the flags).
#[derive(Clone, Copy, PartialEq, Eq)]
enum CharType {
    Residue,
    Comment,
    Other,
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:854-935
// ```c++
//         case 'A': case 'B': case 'C': case 'D':
//         case 'G': case 'H':
//         case 'K':
//         case 'M':
//         case 'R': case 'S': case 'T': case 'U': case 'V': case 'W':
//         case 'Y':
//         case 'a': case 'b': case 'c': case 'd':
//         case 'g': case 'h':
//         case 'k':
//         case 'm':
//         case 'r': case 's': case 't': case 'u': case 'v': case 'w':
//         case 'y':
//             char_type = eCharType_MaskedNonGap;
//         case 'E': case 'F':
//         case 'I': case 'J':
//         case 'L':
//         case 'O': case 'P': case 'Q':
//         case 'Z':
//         case '*':
//             if( bIsNuc ) {
//                 char_type = eCharType_Bad;
//         case 'e': case 'f':
//         case 'i': case 'j':
//         case 'l':
//         case 'o': case 'p': case 'q':
//         case 'z':
//             char_type = (bIsNuc ? eCharType_Bad : eCharType_MaskedNonGap );
//         case 'N':
//             char_type = ( bIsNuc && bAllowLetterGaps ?
//                      eCharType_Gap : eCharType_NormalNonGap );
//         case 'X':
//             char_type = ( bIsNuc ? eCharType_Bad :
//                      bAllowLetterGaps ? eCharType_Gap :
//                      eCharType_NormalNonGap);
//         case ';':
//             char_type = eCharType_Comment;
//         default:
//             char_type = eCharType_Bad;
// ```
// `bIsNuc` is `fAssumeNuc`; `bAllowLetterGaps` is false (no `fParseGaps`); `-`, white space
// and bad bytes are not stored (`fHyphensIgnoreAndWarn`), so they are `Other` here.
const fn char_types(protein: bool) -> [CharType; 256] {
    let mut table = [CharType::Other; 256];
    let nucleotide = b"ABCDGHKMRSTUVWYN";
    let mut index = 0;
    while index < nucleotide.len() {
        table[nucleotide[index] as usize] = CharType::Residue;
        table[nucleotide[index].to_ascii_lowercase() as usize] = CharType::Residue;
        index += 1;
    }
    if protein {
        let protein_only = b"EFIJLOPQZX";
        let mut index = 0;
        while index < protein_only.len() {
            table[protein_only[index] as usize] = CharType::Residue;
            table[protein_only[index].to_ascii_lowercase() as usize] = CharType::Residue;
            index += 1;
        }
        table[b'*' as usize] = CharType::Residue;
    }
    table[b';' as usize] = CharType::Comment;
    table
}

const NUCLEOTIDE_TYPES: [CharType; 256] = char_types(false);
const PROTEIN_TYPES: [CharType; 256] = char_types(true);

/// `isspace` of the C locale (`TruncateSpaces_Unsafe`, `ParseDefline`); CR and LF never
/// reach a line's content.
const fn space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// What a byte of the input is to NCBI's line reader.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Role {
    /// A byte of a line.
    Content,
    /// A byte of a line's end of line.
    Delimiter,
    /// An end-of-line byte that the line reader consumes without ending a line (the LF
    /// after a CR that it pushed back, the CR after an LF that it pushed back).
    Dropped,
    /// A byte that the line reader never returns (the tail lost to the pushback's
    /// end-of-file state).
    Lost,
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:167-173
// ```c++
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
// ```
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Eol {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}

/// The pushback of a line's tail after a CR in the CRLF style (the first pushback of a
/// file that can be followed by a second one).
struct FirstPush {
    /// Offset of the first pushed-back byte.
    start: u64,
    /// Offset of the LF that ended the physical line (consumed, not pushed back).
    lf: Option<u64>,
}

/// NCBI's `CStreamLineReader` over the bytes of a file, pushed in chunks: each byte goes to
/// the record reader with its role, in the order of the input. The line reader's states
/// are those of `FastaStream` and `LineReader` in the engine; a pushed-back tail is read
/// again in the input's order, so the pushback is the end-of-line byte that it skips
/// (`Role::Dropped`) and the style that reads the tail.
struct Lines {
    eol: Eol,
    /// Offset of the next byte.
    offset: u64,
    /// A CR ended the line (`Unknown`, `Cr`, `Mixed`): an LF next is part of its end of line.
    cr_end: bool,
    /// The line has a CR (`Lf`, `CrLf`): an LF next makes it a CRLF.
    lf_cr: bool,
    drop_lf: bool,
    drop_cr: bool,
    first_push: Option<FirstPush>,
    /// The bytes after an LF of a CR-style line read through the pushback buffer (the
    /// second pushback, whose bytes may be lost) and the offset of the first.
    held: Option<(u64, Vec<u8>)>,
}

impl Lines {
    fn new() -> Self {
        Self {
            eol: Eol::Unknown,
            offset: 0,
            cr_end: false,
            lf_cr: false,
            drop_lf: false,
            drop_cr: false,
            first_push: None,
            held: None,
        }
    }

    fn byte(&mut self, byte: u8, records: &mut Records) {
        let offset = self.offset;
        self.offset += 1;
        if let Some((_, held)) = self.held.as_mut() {
            if byte != b'\r' {
                held.push(byte);
                if held.len() > STEP_BACK_MAX {
                    // More than 256 bytes are pushed back into a new buffer, which clears
                    // the end-of-file state: nothing is lost.
                    self.release(records);
                }
                return;
            }
            // The CR that ends the line: the second pushback is not at the end of the input.
            self.release(records);
        }
        self.visible(byte, offset, records);
    }

    fn release(&mut self, records: &mut Records) {
        if let Some((start, held)) = self.held.take() {
            for (index, &byte) in held.iter().enumerate() {
                self.visible(byte, start + index as u64, records);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-243
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLUnknown(void)
    // {
    //     _ASSERT(m_AutoEOL);
    //     NcbiGetline(*m_Stream, m_Line, "\r\n", &m_LastReadSize);
    //     m_Stream->unget();
    //     CT_INT_TYPE eol = m_Stream->get();
    //     if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\r'))) {
    //         m_EOLStyle = eEOL_cr;
    //     } else if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\n'))) {
    //         m_EOLStyle = eEOL_crlf;
    //     }
    //     return m_EOLStyle;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-268
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLSimple(char eol,
    //                                                                    char alt_eol)
    // {
    //     SIZE_TYPE pos;
    //     NcbiGetline(*m_Stream, m_Line, eol, &m_LastReadSize);
    //     if (m_AutoEOL  &&  (pos = m_Line.find(alt_eol)) != NPOS) {
    //         ++pos;
    //         if (eol != '\n'  ||  pos != m_Line.size()) {
    //             // an *immediately* preceding CR is quite all right
    //             CStreamUtils::Pushback(*m_Stream, m_Line.data() + pos,
    //                                    m_Line.size() - pos);
    //             m_EOLStyle = eEOL_mixed;
    //         }
    //         m_Line.resize(pos - 1);
    //         m_LastReadSize = pos;
    //         return (m_EOLStyle == eEOL_mixed) ? m_EOLStyle : eEOL_crlf;
    //     } else if (m_AutoEOL  &&  eol == '\r'  &&
    //                CT_EQ_INT_TYPE(m_Stream->peek(), CT_TO_INT_TYPE(alt_eol))) {
    //         m_Stream->get();
    //         ++m_LastReadSize;
    //         return eEOL_crlf;
    //     }
    //     return (eol == '\r') ? eEOL_cr : eEOL_lf;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLCRLF(void)
    // {
    //     if (m_AutoEOL) {
    //         EEOLStyle style = x_AdvanceEOLSimple('\n', '\r');
    //         if (style == eEOL_mixed) {
    //             // found an embedded CR
    //             m_EOLStyle = eEOL_cr;
    //         } else if (style != eEOL_crlf) {
    //             m_EOLStyle = eEOL_lf;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:90-102
    // ```c++
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
    // ```
    /// One byte that the line reader reads (not held). A CR in the LF or CRLF style that
    /// is not followed by its LF ends the line and pushes the rest of the physical line
    /// back: the LF that ended the physical line was consumed (`drop_lf`), and the rest is
    /// read in the mixed style (from LF) or in the CR style (from CRLF). An LF inside a
    /// CR-style line ends the line and pushes the rest back: the CR that ended the read was
    /// consumed (`drop_cr`), and the rest is read in the mixed style.
    fn visible(&mut self, byte: u8, offset: u64, records: &mut Records) {
        if byte == b'\n' && self.drop_lf {
            self.drop_lf = false;
            if let Some(push) = self.first_push.as_mut() {
                push.lf = Some(offset);
            }
            records.byte(byte, offset, Role::Dropped);
            return;
        }
        if byte == b'\r' && self.drop_cr {
            self.drop_cr = false;
            records.byte(byte, offset, Role::Dropped);
            return;
        }
        if std::mem::take(&mut self.cr_end) {
            if self.eol == Eol::Unknown {
                self.eol = if byte == b'\n' { Eol::CrLf } else { Eol::Cr };
            }
            if byte == b'\n' {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                return;
            }
            records.line_end();
        }
        if std::mem::take(&mut self.lf_cr) {
            if byte == b'\n' {
                // CR LF: the style does not change.
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                return;
            }
            // A CR inside the line: it ends the line; this byte starts the pushed-back rest.
            records.line_end();
            self.drop_lf = true;
            if self.eol == Eol::CrLf {
                self.eol = Eol::Cr;
                self.first_push = Some(FirstPush {
                    start: offset,
                    lf: None,
                });
            } else {
                self.eol = Eol::Mixed;
            }
        }
        match (self.eol, byte) {
            (Eol::Unknown | Eol::Mixed | Eol::Cr, b'\r') => {
                records.byte(byte, offset, Role::Delimiter);
                self.cr_end = true;
            }
            (Eol::Unknown | Eol::Mixed, b'\n') => {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                if self.eol == Eol::Unknown {
                    self.eol = Eol::CrLf;
                }
            }
            (Eol::Cr, b'\n') => {
                // The first LF of a CR-style line.
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                self.eol = Eol::Mixed;
                self.drop_cr = true;
                if self.first_push.is_some() {
                    // A second pushback: its bytes go back into the first one's buffer when
                    // they fit, and the end-of-file state stays if the read ended there.
                    self.held = Some((offset + 1, Vec::new()));
                }
            }
            (Eol::Lf | Eol::CrLf, b'\r') => {
                records.byte(byte, offset, Role::Delimiter);
                self.lf_cr = true;
            }
            (Eol::Lf | Eol::CrLf, b'\n') => {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                if self.eol == Eol::CrLf {
                    self.eol = Eol::Lf;
                }
            }
            _ => records.byte(byte, offset, Role::Content),
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:418-434
    // ```c++
    //             CT_CHAR_TYPE* bp = sb->gptr();
    //             size_t avail = bp - sb->m_Buf;
    //             size_t take  = avail < buf_size ? avail : buf_size;
    //             if (take) {
    //                 bp -= take;
    //                 buf_size -= take;
    //                 if (how != ePushback_Stepback  &&  bp != buf + buf_size) {
    //                     memmove(bp, buf + buf_size, take);
    //                 }
    //                 sb->setg(bp, bp, sb->egptr());
    //             }
    //         }
    //     }
    //
    //     if ( !buf_size ) {
    //         delete[] (CT_CHAR_TYPE*) del_ptr;
    //         return;
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    /// The end of the input. A second pushback at the end of the input (the CR-style line
    /// ended there) whose bytes all step back into the current pushback buffer keeps the
    /// stream's end-of-file state: `AtEOF` is true and those bytes are lost. A new buffer
    /// (more than 256 bytes, or more than the buffer has read) clears the state.
    fn finish(&mut self, records: &mut Records) {
        if let Some((start, held)) = self.held.take() {
            let lost = match &self.first_push {
                Some(FirstPush {
                    start: first,
                    lf: Some(lf),
                }) => {
                    !held.is_empty()
                        && last_refill(lf + 1, lf - first, self.offset) >= held.len() as u64
                }
                _ => false,
            };
            if lost {
                for (index, &byte) in held.iter().enumerate() {
                    records.byte(byte, start + index as u64, Role::Lost);
                }
            } else {
                self.held = Some((start, held));
                self.release(records);
            }
        }
        // A CR at the end ends the line (in the LF and CRLF styles too: `pos ==
        // m_Line.size()`, nothing is pushed back).
        self.cr_end = false;
        self.lf_cr = false;
        records.line_end();
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
// ```c++
// CT_INT_TYPE CPushback_Streambuf::underflow(void)
// {
//     // we are here because there is no more data in the pushback buffer
//     _ASSERT(gptr()  &&  gptr() >= egptr());
//
// #ifdef NCBI_COMPILER_MIPSPRO
//     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
//         return CT_EOF;
//     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
// #endif //NCBI_COMPILER_MIPSPRO
//
//     x_FillBuffer((size_t) m_Sb->in_avail());
//     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-327
// ```c++
// void CPushback_Streambuf::x_FillBuffer(size_t max_size)
// {
//     _ASSERT(m_Sb);
//     if ( !max_size ) {
//         ++max_size;
//     }
//
//     CPushback_Streambuf* sb = dynamic_cast<CPushback_Streambuf*> (m_Sb);
//     if ( !sb ) {
//         CT_CHAR_TYPE* bp = 0;
//         size_t buf_size = m_DelPtr
//             ? (size_t)(m_Buf - (CT_CHAR_TYPE*) m_DelPtr) + m_BufSize : 0;
//         if (buf_size < kMinBufSize) {
//             buf_size = kMinBufSize;
//             bp = new CT_CHAR_TYPE[buf_size];
//         }
//         streamsize r = (streamsize)(buf_size < max_size ? buf_size : max_size);
//         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
//         if (n <= 0) {
//             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
//             delete[] bp;
//             return;
// ```
/// The length of the last refill of the first pushback buffer before the end of a file of
/// `total` bytes: the buffer holds the `pushed` bytes of the first pushback, which ended
/// with the LF before `resume`, and then refills with `in_avail()` bytes, at most its size
/// (at least 4096). `in_avail()` of the file's stream is the rest of its 8191-byte window,
/// or the rest of the file when the window is used up; a read of a whole window or more
/// with the window used up bypasses it (`FastaStream::raw_peek` in the engine, whose
/// windows start at multiples of 8191 until the first pushback).
fn last_refill(resume: u64, pushed: u64, total: u64) -> u64 {
    let capacity = pushed.max(MIN_BUF_SIZE);
    let mut position = resume;
    let mut window = if resume.is_multiple_of(FILE_WINDOW) {
        0
    } else {
        ((resume / FILE_WINDOW + 1) * FILE_WINDOW).min(total) - resume
    };
    let mut last = pushed;
    loop {
        let read = if window > 0 {
            let read = capacity.min(window);
            window -= read;
            read
        } else {
            let rest = total - position;
            if rest == 0 {
                return last;
            }
            let read = capacity.min(rest);
            if read < FILE_WINDOW {
                window = FILE_WINDOW.min(rest) - read;
            }
            read
        };
        position += read;
        last = read;
    }
}

/// How a line that starts with `>` begins (at most 4 bytes, with their offsets).
#[derive(Clone, Copy)]
struct Prefix {
    bytes: [(u8, u64); 4],
    len: usize,
}

impl Prefix {
    fn new(byte: u8, offset: u64) -> Self {
        Self {
            bytes: [(byte, offset); 4],
            len: 1,
        }
    }

    fn push(&mut self, byte: u8, offset: u64) {
        self.bytes[self.len] = (byte, offset);
        self.len += 1;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-374
    // ```c++
    //         if (c == '>' ) {
    //             CTempString next_line = *++GetLineReader();
    //             string strmodified;
    //             if( NStr::StartsWith(next_line, ">?_") ) {
    //                 CTempString tmp = next_line.substr(3);
    //                 strmodified = ">";
    //                 strmodified.append(tmp.data(), tmp.length());
    //                 next_line = strmodified;
    //             }
    //             if( NStr::StartsWith(next_line, ">?") ) {
    //                 // This is actually a data line. an assembly gap, in particular, which
    //                 // we handle farther below
    //                 GetLineReader().UngetLine();
    //             } else {
    //                 if(need_defline) {
    //                     ParseDefLine(next_line, pMessageListener);
    //                     need_defline = false;
    //                     continue;
    //                 } else {
    //                     GetLineReader().UngetLine();
    //                     // start of the next sequence
    //                     break;
    //                 }
    //             }
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:402-409
    // ```c++
    //                 string strmodified;
    //                 if( NStr::StartsWith(line, ">?_") ) {
    //                     CTempString tmp = line.substr(3);
    //                     strmodified = ">";
    //                     strmodified.append(tmp.data(), tmp.length());
    //                     line = strmodified;
    //                 }
    //                 ParseDataLine(line, pMessageListener);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-779
    // ```c++
    // void CFastaReader::ParseDataLine(
    //     const TStr& s, ILineErrorListener * pMessageListener)
    // {
    //     if( NStr::StartsWith(s, ">?") ) {
    //         ParseGapLine(s, pMessageListener);
    //         return;
    // ```
    /// What the line is once enough of it is known (`at_end`: the line has no more bytes):
    /// `Some(None)` a gap line, `Some(Some(skip))` not a gap line whose `>?_` (when
    /// `skip` is 3) becomes `>`, `None` not known yet.
    fn kind(&self, at_end: bool) -> Option<Option<usize>> {
        let byte = |index: usize| (index < self.len).then(|| self.bytes[index].0);
        match (byte(1), byte(2), byte(3)) {
            (None, ..) => at_end.then_some(Some(1)),
            (Some(b'?'), None, _) => at_end.then_some(None),
            (Some(b'?'), Some(b'_'), None) => at_end.then_some(Some(3)),
            (Some(b'?'), Some(b'_'), Some(b'?')) => Some(None),
            (Some(b'?'), Some(b'_'), Some(_)) => Some(Some(3)),
            (Some(b'?'), Some(_), _) => Some(None),
            (Some(_), ..) => Some(Some(1)),
        }
    }
}

/// The ID of a defline being read: the title's first word.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-165
/// ```c++
///     const size_t len = defline.length();
///     if (len <= 1 ||
///         NStr::IsBlank(defline.substr(1))) {
///         return;
///     }
///
///     if (defline[0] != '>') {
///         NCBI_THROW2(CObjReaderParseException, eFormat,
///             "Invalid defline. First character is not '>'", 0);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:168-180
/// ```c++
///     // ignore spaces between '>' and the sequence ID
///     size_t start;
///     for(start = 1 ; start < len; ++start ) {
///         if( ! isspace(defline[start]) ) {
///             break;
///         }
///     }
///
///     size_t pos;
///     size_t title_start = NPOS;
///     if ((fFastaFlags & CFastaReader::fNoParseID)) {
///         title_start = start;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-225
/// ```c++
///     // trim leading whitespace from title (is this appropriate?)
///     while (title_start < len
///         &&  isspace((unsigned char)defline[title_start])) {
///         ++title_start;
///     }
///
///     if (title_start < len) {
///         for (pos = title_start + 1;  pos < len;  ++pos) {
///             if ((unsigned char)defline[pos] < ' ') {
///             break;
///             }
///         }
///         // Parse the title elsewhere - after the molecule has been deduced
///         data.titles.push_back(
///             SLineTextAndLoc(
///                 defline.substr(title_start, pos - title_start), lineNumber));
///     }
/// ```
/// `x_ApplyMods` then drops the title's white space at both ends (its first byte is not
/// white space, and a space ends the first word anyway). The ID is the title up to its
/// first space (`FastaRecord::shown_id`, plan §1).
#[derive(Default)]
struct Title {
    started: bool,
    word_done: bool,
    title_done: bool,
}

impl Title {
    fn byte(&mut self, byte: u8, id: &mut Vec<u8>) {
        if self.title_done {
            return;
        }
        if !self.started {
            if !space(byte) {
                self.started = true;
                id.push(byte);
            }
            return;
        }
        if byte < b' ' {
            self.title_done = true;
        } else if !self.word_done {
            if byte == b' ' {
                self.word_done = true;
            } else {
                id.push(byte);
            }
        }
    }
}

/// `CheckDataLine` of a record's first data lines (no residue stored yet): the first 70
/// bytes of the line (after `TruncateSpaces_Unsafe` and the `>?_` change) and how long
/// the line is up to its last byte other than white space.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:710-760
/// ```c++
/// void CFastaReader::CheckDataLine(
///     const TStr& s, ILineErrorListener * pMessageListener)
/// {
///     // make sure the first data line has at least SOME resemblance to
///     // actual sequence data.
///     if (TestFlag(fSkipCheck)  ||  ! m_SeqData.empty() ) {
///         return;
///     }
///     const bool bIgnoreHyphens = TestFlag(fHyphensIgnoreAndWarn);
///     size_t good = 0, bad = 0;
///     // in case the data has huge sequences all on the first line we do need
///     // a cutoff and "70" seems reasonable since it's the default width of
///     // CFastaOstream (as of 2017-03-09)
///     size_t len_to_check = min(s.length(),
///                               static_cast<size_t>(70));
///     const bool bIsNuc = (
///         ( TestFlag(fAssumeNuc) && TestFlag(fForceType) ) ||
///         ( m_CurrentSeq && m_CurrentSeq->IsSetInst() &&
///           m_CurrentSeq->GetInst().IsSetMol() &&  m_CurrentSeq->IsNa() ) );
///     size_t ambig_nuc = 0;
///     for (size_t pos = 0;  pos < len_to_check;  ++pos) {
///         unsigned char c = s[pos];
///         if (s_ASCII_IsAlpha(c) ||  c == '*') {
///             ++good;
///             if( bIsNuc && s_ASCII_IsAmbigNuc(c) ) {
///                 ++ambig_nuc;
///             }
///         } else if( c == '-' ) {
///             if( ! bIgnoreHyphens ) {
///                 ++good;
///             }
///             // if bIgnoreHyphens == true, the "hyphens are ignored" warning
///             // will be triggered elsewhere
///         } else if (isspace(c)  ||  (c >= '0' && c <= '9')) {
///             // treat whitespace and digits as neutral
///         } else if (c == ';') {
///             break; // comment -- ignore rest of line
///         } else {
///             ++bad;
///         }
///     }
///     if (bad >= good / 3  &&
///         (len_to_check > 3  ||  good == 0  ||  bad > good))
///     {
///         FASTA_ERROR( LineNumber(),
///             "CFastaReader: Near line " << LineNumber()
///             << ", there's a line that doesn't look like plausible data, "
///             "but it's not marked as defline or comment.",
///             CObjReaderParseException::eFormat);
/// ```
/// The warning about ambiguous nucleotides is never written (`bIsNuc` is false before
/// `AssembleSeq`).
struct Check {
    bytes: [u8; CHECK_WINDOW],
    len: usize,
    /// One past the last byte other than white space, at most 71.
    end: usize,
}

impl Check {
    fn new() -> Self {
        Self {
            bytes: [0; CHECK_WINDOW],
            len: 0,
            end: 0,
        }
    }

    fn byte(&mut self, byte: u8) {
        let index = self.len;
        if index < CHECK_WINDOW {
            self.bytes[index] = byte;
            self.len += 1;
            if !space(byte) {
                self.end = index + 1;
            }
        } else if !space(byte) {
            self.end = CHECK_WINDOW + 1;
        }
    }

    fn plausible(&self) -> bool {
        let len_to_check = self.end.min(CHECK_WINDOW);
        let (mut good, mut bad) = (0usize, 0usize);
        for &c in &self.bytes[..len_to_check] {
            if c.is_ascii_alphabetic() || c == b'*' {
                good += 1;
            } else if c == b'-' || space(c) || c.is_ascii_digit() {
            } else if c == b';' {
                break;
            } else {
                bad += 1;
            }
        }
        !(bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good))
    }
}

/// The state of the line being read.
enum LineState {
    /// No byte of the line yet.
    Fresh,
    /// The line starts with `>`: a defline or a gap line.
    Angle(Prefix),
    /// The defline of the current record.
    Title(Title),
    /// White space at the start of a line that does not start with `>`.
    Leading,
    /// A data line whose first byte other than white space is `>` (a gap line or not),
    /// with its `CheckDataLine`.
    DataAngle(Prefix, Option<Check>),
    /// A data line: whether a `;` ended its data, and its `CheckDataLine`.
    Data {
        semicolon: bool,
        check: Option<Check>,
    },
    /// A comment line.
    Skip,
}

/// The rules by which an application reads the residues of a record forward from a
/// checkpoint (docs/web/abi_v2.md §9, kinds 1 and 2), applied to every byte of the input
/// from the record's first residue (checkpoint 0) on. The byte at a checkpoint is a
/// residue. Lines end at every CR and every LF. A line whose first byte other than space,
/// tab, VT and FF is `!`, `#` or `;` is skipped; in other lines, a `;` skips the rest of
/// the line, and every byte that the kind stores is a residue (kind 1: `ABCDGHKMNRSTUVWY`
/// in either case; kind 2: every ASCII letter and `*`); every other byte (white space,
/// `-`, `>`, digits, other bytes) is skipped.
///
/// The scan checks that these rules find the reader's residues: from the first residue,
/// the residues they find must begin with the reader's residues, at the same offsets
/// (then reading from any checkpoint finds them too: the rules are in the same state at
/// every residue that they find). They differ only where the line reader joins two lines
/// (`Role::Dropped`): the joined line's comment or `;` covers the second line, or the
/// second line's comment mark is read as data.
#[derive(Default)]
struct Forward {
    active: bool,
    state: ForwardState,
    /// The rules found a residue that the reader did not store.
    ahead: bool,
    /// The rules have found every residue that the reader stored, at its offset.
    readable: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum ForwardState {
    #[default]
    LineStart,
    Data,
    Skip,
}

impl Forward {
    /// The rules from the record's first residue (a checkpoint: in a data line).
    fn started() -> Self {
        Self {
            active: true,
            state: ForwardState::Data,
            ahead: false,
            readable: true,
        }
    }

    /// Whether `byte` is a residue under the rules.
    fn step(&mut self, byte: u8, types: &[CharType; 256]) -> bool {
        if byte == b'\r' || byte == b'\n' {
            self.state = ForwardState::LineStart;
            return false;
        }
        if self.state == ForwardState::LineStart {
            if matches!(byte, b' ' | b'\t' | 0x0b | 0x0c) {
                return false;
            }
            if matches!(byte, b'!' | b'#' | b';') {
                self.state = ForwardState::Skip;
                return false;
            }
            self.state = ForwardState::Data;
        }
        if self.state == ForwardState::Skip {
            return false;
        }
        match types[byte as usize] {
            CharType::Residue => true,
            CharType::Comment => {
                self.state = ForwardState::Skip;
                false
            }
            CharType::Other => false,
        }
    }

    /// Compares the reader (`stored`) and the rules (`found`) at one byte.
    fn compare(&mut self, stored: bool, found: bool) {
        if stored {
            if self.ahead || !found {
                self.readable = false;
            }
        } else if found {
            self.ahead = true;
        }
    }
}

/// A record being read.
struct Builder {
    id: Vec<u8>,
    header_offset: u64,
    sequence_offset: u64,
    length: u64,
    residue_counts: Box<[u64; 256]>,
    checkpoints: Vec<u64>,
    uniform: bool,
    width: u64,
    eol: u64,
    previous: u64,
}

impl Builder {
    fn new(header_offset: u64, sequence_offset: u64) -> Self {
        Self {
            id: Vec::new(),
            header_offset,
            sequence_offset,
            length: 0,
            residue_counts: Box::new([0; 256]),
            checkpoints: Vec::new(),
            uniform: true,
            width: 0,
            eol: 1,
            previous: 0,
        }
    }

    /// Stores a residue at `offset` (upper-cased; a nucleotide `U` as `T`). The layout is
    /// uniform while every residue is where the formula of the first line's width and the
    /// gap after it puts it.
    fn store(&mut self, residue: u8, offset: u64) -> Result<(), String> {
        let index = self.length;
        if index.is_multiple_of(CHECKPOINT_EVERY) {
            self.checkpoints.push(offset);
        }
        if index == 0 {
            self.uniform = offset == self.sequence_offset;
        } else if self.uniform {
            if self.width == 0 {
                if offset != self.previous + 1 {
                    self.width = index;
                    self.eol = offset - self.previous - 1;
                }
            } else if offset
                != self.sequence_offset
                    + index / self.width * (self.width + self.eol)
                    + index % self.width
            {
                self.uniform = false;
            }
        }
        self.previous = offset;
        self.residue_counts[residue as usize] += 1;
        self.length += 1;
        // NCBI reference (598d8ae6): c++/include/corelib/ncbimisc.hpp:879
        // ```c++
        // typedef unsigned int TSeqPos;
        // ```
        if self.length > MAX_LETTERS {
            return Err(format!(
                "a record longer than {MAX_LETTERS} letters is not supported by LOSAT Web"
            ));
        }
        Ok(())
    }

    fn finish(
        self,
        end_offset: u64,
        forward: &Forward,
        index: usize,
    ) -> Result<ScanRecord, String> {
        let id = String::from_utf8_lossy(&self.id).into_owned();
        let layout = if self.uniform {
            let width = if self.width == 0 {
                self.length
            } else {
                self.width
            };
            LineLayout::Uniform {
                width,
                eol: self.eol,
            }
        } else if forward.readable {
            LineLayout::Checkpoints {
                offsets: self.checkpoints,
            }
        } else {
            return Err(format!(
                "record {} ({id:?}): NCBI BLAST+ reads two of its lines as one (a CR line end in a file of LF line ends, or an LF in a file of CR line ends), so LOSAT Web's index cannot locate its residues; such input is not supported by LOSAT Web (use one kind of line end)",
                index + 1
            ));
        };
        Ok(ScanRecord {
            id,
            header_offset: self.header_offset,
            sequence_offset: self.sequence_offset,
            end_offset,
            length: self.length,
            layout,
            residue_counts: self.residue_counts,
        })
    }
}

/// The record reader (`CFastaReader::ReadOneSeq` as `CBlastInputReader` calls it), fed
/// with the line reader's bytes.
struct Records {
    protein: bool,
    types: &'static [CharType; 256],
    line_number: u64,
    in_line: bool,
    line_start: u64,
    state: LineState,
    /// The first line from its first byte other than white space, held while it may be a
    /// Seq-id.
    first_line: Option<Vec<u8>>,
    /// No record has started (`need_defline` of the first `ReadOneSeq`).
    need_defline: bool,
    current: Option<Builder>,
    /// The current record's `sequence_offset` is the next byte that the line reader
    /// reads after its defline.
    awaiting_sequence: bool,
    forward: Forward,
    records: Vec<ScanRecord>,
    error: Option<String>,
}

impl Records {
    fn new(protein: bool) -> Self {
        Self {
            protein,
            types: if protein {
                &PROTEIN_TYPES
            } else {
                &NUCLEOTIDE_TYPES
            },
            line_number: 0,
            in_line: false,
            line_start: 0,
            state: LineState::Fresh,
            first_line: None,
            need_defline: true,
            current: None,
            awaiting_sequence: false,
            forward: Forward::default(),
            records: Vec::new(),
            error: None,
        }
    }

    fn fail(&mut self, message: String) {
        if self.error.is_none() {
            self.error = Some(message);
        }
    }

    /// One byte of the input, in the input's order, with its role.
    fn byte(&mut self, byte: u8, offset: u64, role: Role) {
        if self.error.is_some() {
            return;
        }
        if matches!(role, Role::Content | Role::Delimiter) {
            if std::mem::take(&mut self.awaiting_sequence) {
                if let Some(current) = self.current.as_mut() {
                    current.sequence_offset = offset;
                }
            }
            if !self.in_line {
                // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-161
                // ```c++
                // CStreamLineReader& CStreamLineReader::operator++(void)
                // {
                //     /* If at EOF - noop */
                //     if (AtEOF()) {
                //         m_Line = string();
                //         return *this;
                //     }
                //     ++m_LineNumber;
                // ```
                self.in_line = true;
                self.line_number += 1;
                self.line_start = offset;
                self.state = LineState::Fresh;
            }
        }
        let stored = role == Role::Content && self.content(byte, offset);
        if self.forward.active {
            let found = self.forward.step(byte, self.types);
            self.forward.compare(stored, found);
        } else if stored {
            // The record's first residue.
            self.forward = Forward::started();
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:339-340
    // ```c++
    //     while ( !GetLineReader().AtEOF() ) {
    //         char c = GetLineReader().PeekChar();
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:376-389
    // ```c++
    //         CTempString line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    //
    //         if (line.empty()) {
    //             continue; // ignore lines containing only whitespace
    //         }
    //         c = line[0];
    //
    //         if (c == '!'  ||  c == '#' || c == ';') {
    //             // no content, just a comment or blank line
    //             continue;
    //         } else if (need_defline) {
    //             if (TestFlag(fDLOptional)) {
    //                 ParseDefLine(">", pMessageListener);
    //                 need_defline = false;
    // ```
    /// A byte of the current line; whether it is a stored residue.
    fn content(&mut self, byte: u8, offset: u64) -> bool {
        if let Some(line) = self.first_line.as_mut() {
            line.push(byte);
        }
        match std::mem::replace(&mut self.state, LineState::Skip) {
            LineState::Fresh if byte == b'>' => {
                self.state = LineState::Angle(Prefix::new(byte, offset));
                false
            }
            LineState::Fresh | LineState::Leading => self.leading(byte, offset),
            LineState::Angle(mut prefix) => {
                prefix.push(byte, offset);
                self.angle(prefix, false);
                false
            }
            LineState::Title(mut title) => {
                if let Some(current) = self.current.as_mut() {
                    title.byte(byte, &mut current.id);
                }
                self.state = LineState::Title(title);
                false
            }
            LineState::DataAngle(mut prefix, check) => {
                prefix.push(byte, offset);
                self.data_angle(prefix, check, false)
            }
            LineState::Data { semicolon, check } => {
                self.state = LineState::Data { semicolon, check };
                if let LineState::Data {
                    check: Some(check), ..
                } = &mut self.state
                {
                    check.byte(byte);
                }
                self.residue(byte, offset)
            }
            LineState::Skip => false,
        }
    }

    /// The first byte other than white space of a line that does not start with `>`.
    fn leading(&mut self, byte: u8, offset: u64) -> bool {
        if space(byte) {
            self.state = LineState::Leading;
            return false;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:132-133
        // ```c++
        //         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
        //         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
        // ```
        if self.line_number == 1 && byte.is_ascii_alphanumeric() {
            self.first_line = Some(vec![byte]);
        }
        if matches!(byte, b'!' | b'#' | b';') {
            return false;
        }
        if self.need_defline {
            // `ParseDefLine(">")`: the input's first record, without a title.
            self.need_defline = false;
            self.current = Some(Builder::new(0, 0));
        }
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:715-717
        // ```c++
        //     if (TestFlag(fSkipCheck)  ||  ! m_SeqData.empty() ) {
        //         return;
        //     }
        // ```
        let check = self
            .current
            .as_ref()
            .is_some_and(|current| current.length == 0)
            .then(Check::new);
        if byte == b'>' {
            self.state = LineState::DataAngle(Prefix::new(byte, offset), check);
            return false;
        }
        let mut check = check;
        if let Some(check) = check.as_mut() {
            check.byte(byte);
        }
        self.state = LineState::Data {
            semicolon: false,
            check,
        };
        self.residue(byte, offset)
    }

    /// A line that starts with `>`, once its kind is known: a defline starts a record
    /// (and ends the current one), a gap line is rejected.
    fn angle(&mut self, prefix: Prefix, at_end: bool) {
        match prefix.kind(at_end) {
            None => self.state = LineState::Angle(prefix),
            Some(None) => self.gap(),
            Some(Some(skip)) => {
                if let Some(current) = self.current.take() {
                    self.close(current, self.line_start);
                }
                self.need_defline = false;
                let mut current = Builder::new(self.line_start, self.line_start);
                let mut title = Title::default();
                for &(byte, _) in &prefix.bytes[skip..prefix.len] {
                    title.byte(byte, &mut current.id);
                }
                self.current = Some(current);
                self.forward = Forward::default();
                self.state = LineState::Title(title);
            }
        }
    }

    /// A data line whose first byte other than white space is `>`, once its kind is known.
    fn data_angle(&mut self, prefix: Prefix, check: Option<Check>, at_end: bool) -> bool {
        match prefix.kind(at_end) {
            None => {
                self.state = LineState::DataAngle(prefix, check);
                false
            }
            Some(None) => {
                self.gap();
                false
            }
            Some(Some(skip)) => {
                let mut check = check;
                if let Some(check) = check.as_mut() {
                    check.byte(prefix.bytes[0].0);
                    for &(byte, _) in &prefix.bytes[skip..prefix.len] {
                        check.byte(byte);
                    }
                }
                self.state = LineState::Data {
                    semicolon: false,
                    check,
                };
                // `>`, `?` and `_` are not residues: only the line's last byte, when this
                // byte decided the kind, can be one.
                match (at_end, prefix.len) {
                    (false, len) if len > 1 => {
                        let (byte, offset) = prefix.bytes[len - 1];
                        self.residue(byte, offset)
                    }
                    _ => false,
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:938-950
    // ```c++
    //         switch(char_type) {
    //         case eCharType_NormalNonGap:
    //             CloseGap(pos == 0);
    //             m_SeqData[m_CurrentPos] = c;
    //             CloseMask();
    //             ++m_CurrentPos;
    //             break;
    //         case eCharType_MaskedNonGap:
    //             CloseGap(pos == 0);
    //             m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
    //             OpenMask();
    //             ++m_CurrentPos;
    //             break;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:970-973
    // ```c++
    //         case eCharType_Comment:
    //             // artificially advance pos to the end to break the pos loop
    //             pos = s_len;
    //             break;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447
    // ```c++
    //         CRef<CSeq_data> data(new CSeq_data(m_SeqData, format));
    //         if ( !TestFlag(fLeaveAsText) ) {
    //             CSeqportUtil::Pack(data, inst.GetLength());
    //         }
    // ```
    /// A byte of a data line's data: a residue is stored upper-cased (the engine stores a
    /// nucleotide `U` as `T`; the lowercase mask only changes case); `;` ends the data.
    fn residue(&mut self, byte: u8, offset: u64) -> bool {
        let LineState::Data { semicolon, .. } = &mut self.state else {
            return false;
        };
        if *semicolon {
            return false;
        }
        match self.types[byte as usize] {
            CharType::Residue => {}
            CharType::Comment => {
                *semicolon = true;
                return false;
            }
            CharType::Other => return false,
        }
        let mut residue = byte.to_ascii_uppercase();
        if !self.protein && residue == b'U' {
            residue = b'T';
        }
        if let Some(Err(error)) = self
            .current
            .as_mut()
            .map(|current| current.store(residue, offset))
        {
            self.fail(error);
        }
        true
    }

    /// A `>?` gap line (maintainer decision 3: `register` and `scan` reject them).
    fn gap(&mut self) {
        let line = self.line_number;
        self.fail(format!(
            "line {line} is a gap line ('>?'), whose residues have no bytes in the input for LOSAT Web's index to locate; gap lines are not supported by LOSAT Web"
        ));
    }

    /// The end of the current line, if one is open.
    fn line_end(&mut self) {
        if !std::mem::take(&mut self.in_line) || self.error.is_some() {
            return;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-162
        // ```c++
        //     virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener * pMessageListener) {
        //
        //         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
        //         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
        //             try {
        //                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
        // 							CSeq_id::fParse_ValidLocal)));
        // 		if (id->IsLocal()  &&  !NStr::StartsWith(line, "lcl|") ) {
        //                     // Expected to throw an exception.
        //                     id.Reset(new CSeq_id(line));
        // 		}
        //                 CRef<CBioseq> bioseq(x_CreateBioseq(id));
        //                 CRef<CSeq_entry> retval(new CSeq_entry());
        //                 retval->SetSeq(*bioseq);
        //                 return retval;
        //             } catch (const CSeqIdException& e) {
        //                 if (NStr::Find(e.GetMsg(), "Malformatted ID") != NPOS) {
        //                     // This is probably just plain fasta, so just
        //                     // defer to CFastaReader
        //                 } else {
        //                     throw;
        //                 }
        //             } catch (const exception&) {
        //                 throw;
        //             } catch (...) {
        //                 // in case of other exceptions, just defer to CFastaReader
        //             }
        //         } // end if ( !line.empty() ...
        //
        //         // If all fails, fall back to parent's implementation
        //         GetLineReader().UngetLine();
        //         return CFastaReader::ReadOneSeq(pMessageListener);
        //     }
        // ```
        // Every later record starts at the `>` line that ended the one before, so only the
        // input's first line is tried.
        if let Some(line) = self.first_line.take() {
            if seq_id_line(&line, self.protein) {
                let line = String::from_utf8_lossy(&line);
                self.fail(format!(
                    "the first line ({:?}) is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT Web (start the input with a '>' defline)",
                    line.trim_end()
                ));
                return;
            }
        }
        match std::mem::replace(&mut self.state, LineState::Fresh) {
            LineState::Angle(prefix) => self.angle(prefix, true),
            LineState::DataAngle(prefix, check) => {
                self.data_angle(prefix, check, true);
            }
            state => self.state = state,
        }
        match std::mem::replace(&mut self.state, LineState::Fresh) {
            LineState::Title(_) => self.awaiting_sequence = true,
            LineState::Data {
                check: Some(check), ..
            } if !check.plausible() => {
                let line = self.line_number;
                self.fail(format!(
                    "CFastaReader: Near line {line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
                ));
            }
            _ => {}
        }
    }

    fn close(&mut self, current: Builder, end_offset: u64) {
        match current.finish(end_offset, &self.forward, self.records.len()) {
            Ok(record) => self.records.push(record),
            Err(error) => self.fail(error),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:428-432
    // ```c++
    //     if (need_defline  &&  GetLineReader().AtEOF()) {
    //         FASTA_ERROR(LineNumber(),
    //             "CFastaReader: Expected defline around line " << LineNumber(),
    //             CObjReaderParseException::eEOF);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:203-208
    // ```c++
    //     while (!End()) {
    //         try { retval->AddQuery(m_Source->GetNextSequence(scope)); }
    //         catch (const CObjReaderParseException& e) {
    //             auto err = e.GetErrCode();
    //             if (err == CObjReaderParseException::eEOF) {
    //                 break;
    // ```
    /// The end of the input (`total` bytes): the current record ends there; an input with
    /// only white space and comment lines has no record.
    fn finish(mut self, total: u64) -> Result<Vec<ScanRecord>, String> {
        if self.error.is_none() {
            if let Some(mut current) = self.current.take() {
                if self.awaiting_sequence {
                    current.sequence_offset = total;
                }
                self.close(current, total);
            }
        }
        match self.error {
            Some(error) => Err(error),
            None => Ok(self.records),
        }
    }
}

/// Whether NCBI BLAST+, with its data loaders on, reads `line` (the first line of an input,
/// without its white space at the start) as a Seq-id: the engine's reader decides
/// (`ReadError::Unsupported` for such a line, `seq_id.rs`).
fn seq_id_line(line: &[u8], protein: bool) -> bool {
    use LOSAT::blastinput::fasta_reader::{FastaInputSource, ReadError, ReaderConfig};
    let mut bytes = Vec::with_capacity(line.len() + 1);
    bytes.extend_from_slice(line);
    bytes.push(b'\n');
    let mut source =
        FastaInputSource::from_bytes(&bytes, ReaderConfig::query("LOSAT Web", protein, true));
    matches!(
        source.next_sequence(&mut |_| Ok(())),
        Err(ReadError::Unsupported(_))
    )
}

/// A streaming scan with parser kind 1 (`protein` false) or 2 (`protein` true).
pub struct NcbiScanner {
    lines: Lines,
    records: Records,
}

impl NcbiScanner {
    pub fn new(protein: bool) -> Self {
        Self {
            lines: Lines::new(),
            records: Records::new(protein),
        }
    }

    /// Scans the next chunk of the input.
    pub fn feed(&mut self, chunk: &[u8]) {
        for &byte in chunk {
            if self.records.error.is_some() {
                return;
            }
            self.lines.byte(byte, &mut self.records);
        }
    }

    /// Ends the input and returns the records, or the scan's error.
    pub fn finish(mut self) -> Result<Vec<ScanRecord>, String> {
        if self.records.error.is_none() {
            self.lines.finish(&mut self.records);
        }
        let total = self.lines.offset;
        self.records.finish(total)
    }
}
