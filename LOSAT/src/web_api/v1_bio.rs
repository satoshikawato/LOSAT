//! ABI v1's `bio` FASTA records (plan TD-1, which freezes its accepted inputs and
//! messages): the checks that ABI v1 makes of them, and the records of NCBI's reader made
//! from those that pass (`from_bio`), which its searches take. LOSAT's programs and ABI v2
//! read their inputs with NCBI's reader (`fasta_reader`) and make none of these checks
//! (session SFc, S10: moved from `blastinput/bio_checks.rs` and `FastaRecord::from_bio`).
//!
//! ABI v1 parses FASTA with `bio`, while NCBI BLAST+ reads it with `CFastaReader`. For the
//! inputs accepted here the two agree, except that `U` is read as `T`, which `from_bio`
//! applies. The others are rejected, because the reader of NCBI would change them:
//!
//! - an empty defline, or one that starts with white space or has a control character or
//!   a non-ASCII byte (NCBI names an empty one `Query_1`, skips the white space, ends the
//!   title at the control character, and
//!   splits the ID and wraps the title by bytes, where `bio` splits by Unicode white space
//!   and the report wraps by characters): `check_deflines`, on the bytes of the file,
//!   because `bio` drops the character that ends the ID;
//! - a record without residues, or a residue that is not an IUPAC nucleotide letter (NCBI
//!   warns that the sequence contains no data, ignores white space and hyphens, ends the
//!   line at `;`, and removes other characters with a warning): `check_residues_of`.
//!
//! BLASTP and TBLASTX accept the deflines that ABI v1 accepted (plan TD-1), and reject only
//! those whose bytes the shared report now makes neither ABI v1's nor NCBI's
//! (`check_bio_deflines_of`, session SFd): a defline to which `bio` gives an empty ID when
//! the record that `from_bio` makes of it is not NCBI's reader's (a non-ASCII white space
//! character that `bio` drops and NCBI keeps, a carriage return that ends NCBI's line, a
//! control character that ends NCBI's title); a record without a title after a `>?` line,
//! where its local ID shows; and in BLASTP's outfmt 0, a title that `bio` reads otherwise
//! and that ends with a non-ASCII character.

use anyhow::{bail, Result};
use bio::io::fasta;

use crate::blastinput::fasta_reader::FastaRecord;

/// The IUPAC nucleotide letters, both cases.
const IUPAC_NUCLEOTIDE: [bool; 256] = {
    let letters = b"ACGTUMRWSYKVHDBNacgtumrwsykvhdbn";
    let mut table = [false; 256];
    let mut index = 0;
    while index < letters.len() {
        table[letters[index] as usize] = true;
        index += 1;
    }
    table
};

/// Rejects the deflines that NCBI reads differently (see the module), in the bytes of a
/// FASTA file split into lines as `bio` splits them. `bio` trims the white space at the
/// end of a line, and so does NCBI.
///
/// NCBI reference: c++/src/objtools/readers/fasta_reader_utils.cpp:168-225
/// ```c
///     // ignore spaces between '>' and the sequence ID
///     size_t start;
///     for(start = 1 ; start < len; ++start ) {
///         if( ! isspace(defline[start]) ) {
///             break;
///         }
///     }
///     ...
///         while (pos < len && defline[pos] > ' ') {
///             pos++;
///         }
///     ...
///     if (title_start < len) {
///         for (pos = title_start + 1;  pos < len;  ++pos) {
///             if ((unsigned char)defline[pos] < ' ') {
///             break;
///             }
///         }
/// ```
/// `role` is `query` or `subject`.
pub fn check_deflines(bytes: &[u8], role: &str) -> Result<()> {
    check_deflines_of(bytes, role, "BLASTN")
}

/// `check_deflines` for the FASTA input of `program` (as named in the message).
pub fn check_deflines_of(bytes: &[u8], role: &str, program: &str) -> Result<()> {
    check_deflines_with(bytes, role, program)
}

fn check_deflines_with(bytes: &[u8], role: &str, program: &str) -> Result<()> {
    let mut record = 0;
    for line in bytes.split(|&byte| byte == b'\n') {
        let Some(defline) = line.strip_prefix(b">") else {
            continue;
        };
        record += 1;
        // NCBI keeps the white space at the end of a line (it strips the carriage return
        // of CRLF), which `bio` drops; it matters only to the title warning.
        let raw = defline.strip_suffix(b"\r").unwrap_or(defline);
        let defline = defline.trim_ascii_end();
        let problem = if defline.first() == Some(&b'?') {
            "starts with '?' (NCBI BLAST+ reads '>?' as a gap in the sequence, and '>?_' as a defline without the prefix)".to_string()
        } else if raw.len() != defline.len() && ends_with_nucleotides(defline) {
            "ends with white space after 20 nucleotide letters (NCBI BLAST+'s warning about the letters depends on that white space, which LOSAT's reader drops)".to_string()
        } else if defline.is_empty() {
            "is empty".to_string()
        } else if defline.first().is_some_and(u8::is_ascii_whitespace) {
            "starts with white space".to_string()
        } else if let Some(&byte) = defline.iter().find(|&&byte| byte < b' ') {
            format!("has the control character 0x{byte:02x}")
        } else if !defline.is_ascii() {
            "has a non-ASCII byte".to_string()
        } else {
            continue;
        };
        bail!(
            "{role} record {record} has a defline that {problem}; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's {program} (use ASCII deflines without control characters)"
        );
    }
    Ok(())
}

/// `isspace` of the C locale, which NCBI's defline parser and `NStr` use.
const fn c_isspace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// What an ABI v1 search shows of its records (`check_bio_deflines_of`).
#[derive(Clone, Copy, Debug, Default)]
pub struct Shown {
    /// The names of records without a title (`Query_N`, `Subject_N`, `unnamed`): the TBLASTX
    /// report, BLASTP's subjects, and the query IDs of BLASTP's tabular formats. (BLASTP's
    /// outfmt 0 and 7 query lines show the title only, which such a record lacks.)
    pub names: bool,
    /// The local IDs (`Query_N`, `Subject_N`) of records without a title: the TBLASTX
    /// report, and the query IDs of BLASTP's tabular formats (`qseqid`, `qacc`,
    /// `qaccver`). (A BLASTP subject without a title is `unnamed`; the outfmt 0 and 7
    /// query lines show the title only.)
    pub local_ids: bool,
    /// The titles of BLASTP's outfmt 0 subjects (`CDeflineGenerator::GenerateDefline`).
    pub outfmt0_titles: bool,
}

/// ABI v1's BLASTP and TBLASTX: rejects the deflines for which the bytes of the search are
/// now neither ABI v1's (S11) nor NCBI's (see the module), in the bytes of a FASTA file split
/// into lines as `bio` splits them. `bio` 1.6 (`src/io/fasta.rs:331-333`) trims the Unicode
/// white space at the end of the line and ends the ID at the first `char::is_whitespace`.
///
/// - A defline that starts with such a character, or has nothing else, has an empty ID, and
///   `from_bio`'s title is the rest of the line without that character and the C white space
///   after it. ABI v1 printed `unknown` for such a record until it searched the records of
///   NCBI's reader (session SFc, S7 and S8). It is rejected when `from_bio`'s record is not
///   NCBI's (its title, or the lines that a carriage return starts), except a record left
///   without a title where its name does not show (`shown.names`; ABI v1's bytes, as
///   before). The others are NCBI's: an empty defline, one of C white space (with carriage
///   returns), and a title after C white space.
/// - Where the search shows the local IDs (`shown.local_ids`), a record without a title after
///   a `>?` line is rejected: NCBI reads that line as a gap in the record before it, not as a
///   record, so it numbers the record lower than `bio` (a `>?` first line opens a record
///   without a title for NCBI, so it shifts nothing).
/// - In BLASTP's outfmt 0 (`shown.outfmt0_titles`), a defline with an ID whose `from_bio`
///   record is not NCBI's (a `>?` line is a gap for NCBI, and NCBI drops a `?_` prefix) is
///   rejected when its title ends with a non-ASCII character: NCBI's
///   `x_CleanAndCompress` drops such a last byte (a signed `char`), which ABI v1's report did
///   not, so the shared report's bytes for `bio`'s title are neither. (Other titles that `bio`
///   reads otherwise keep ABI v1's bytes.)
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-225
/// ```c++
///     const size_t len = defline.length();
///     if (len <= 1 ||
///         NStr::IsBlank(defline.substr(1))) {
///         return;
///     }
///     ...
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
///     ...
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-322
/// ```c++
///     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
///                                     CFastaReader::fParseRawID:
///                                     (CFastaReader::fNoParseID |
///                                      CFastaReader::fDLOptional);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038
/// ```c++
///     NStr::TruncateSpacesInPlace(processed_title);
/// ```
/// NCBI's line reader ends a line at a carriage return as at a line feed (a CRLF is one
/// line end), and the reader skips a line of white space.
///
/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-222
/// ```c++
/// CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLUnknown(void)
/// {
///     _ASSERT(m_AutoEOL);
///     NcbiGetline(*m_Stream, m_Line, "\r\n", &m_LastReadSize);
/// ```
/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:249-258
/// ```c++
///     NcbiGetline(*m_Stream, m_Line, eol, &m_LastReadSize);
///     if (m_AutoEOL  &&  (pos = m_Line.find(alt_eol)) != NPOS) {
///         ++pos;
///         if (eol != '\n'  ||  pos != m_Line.size()) {
///             // an *immediately* preceding CR is quite all right
///             CStreamUtils::Pushback(*m_Stream, m_Line.data() + pos,
///                                    m_Line.size() - pos);
///             m_EOLStyle = eEOL_mixed;
///         }
///         m_Line.resize(pos - 1);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:376-380
/// ```c++
///         CTempString line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
///
///         if (line.empty()) {
///             continue; // ignore lines containing only whitespace
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-362
/// ```c++
///         if (c == '>' ) {
///             CTempString next_line = *++GetLineReader();
///             string strmodified;
///             if( NStr::StartsWith(next_line, ">?_") ) {
///                 CTempString tmp = next_line.substr(3);
///                 strmodified = ">";
///                 strmodified.append(tmp.data(), tmp.length());
///                 next_line = strmodified;
///             }
///             if( NStr::StartsWith(next_line, ">?") ) {
///                 // This is actually a data line. an assembly gap, in particular, which
///                 // we handle farther below
///                 GetLineReader().UngetLine();
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:386-389
/// ```c++
///         } else if (need_defline) {
///             if (TestFlag(fDLOptional)) {
///                 ParseDefLine(">", pMessageListener);
///                 need_defline = false;
/// ```
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3955-3958
/// ```c++
///         size_t pos = m_MainTitle.find_last_not_of (".,;~ ");
///         if (pos != NPOS) {
///             m_MainTitle.erase (pos + 1);
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4070-4073
/// ```c++
///     size_t pos = decoded.find_last_not_of (",;~ ");
///     if (pos != NPOS) {
///         decoded.erase (pos + 1);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:304-306
/// ```c++
///     if (curr > 0 && curr != ' ') {
///         *out++ = curr;
///     }
/// ```
/// `role` is `query` or `subject`; `program` is named in the message.
pub fn check_bio_deflines_of(bytes: &[u8], role: &str, program: &str, shown: Shown) -> Result<()> {
    let mut record = 0;
    let mut gaps = 0;
    for line in bytes.split(|&byte| byte == b'\n') {
        let Some(defline) = line.strip_prefix(b">") else {
            continue;
        };
        record += 1;
        // `bio` reads UTF-8 only, so ABI v1 searches no other bytes.
        let Ok(text) = std::str::from_utf8(defline) else {
            continue;
        };
        let reading = BioReading::of(text);
        // A record without a title changes only where its name shows.
        let shows = shown.names || !reading.title.is_empty();
        if let Some(problem) = reading.empty_id_problem().filter(|_| shows) {
            bail!(
                "{role} record {record} has a defline that starts with white space and {problem}; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's {program} (begin the defline with a character that is not white space)"
            );
        }
        if shown.local_ids && reading.id_empty && reading.title.is_empty() && gaps > 0 {
            bail!(
                "{role} record {record} has a defline without a title after a line that starts with '>?', which NCBI BLAST+ reads as a gap, not as a record (NCBI BLAST+ numbers this record {}); this is not supported by LOSAT's {program} (give the record a title)",
                record - gaps
            );
        }
        if shown.outfmt0_titles {
            if let Some(problem) = reading.outfmt0_title_problem() {
                bail!(
                    "{role} record {record} has a defline that NCBI BLAST+ reads with another title (it {problem}) and that ends with a non-ASCII character, whose last byte NCBI BLAST+'s outfmt 0 drops; this is not supported by LOSAT's {program} outfmt 0 (use ASCII deflines without control characters)"
                );
            }
        }
        // NCBI's gap lines after the first line (`>?`, and `>?_?` read as `>?`).
        let gap = match defline.strip_prefix(b"?_") {
            Some(rest) => rest.starts_with(b"?"),
            None => defline.starts_with(b"?"),
        };
        if gap && record > 1 {
            gaps += 1;
        }
    }
    Ok(())
}

/// A defline (the line after `>` as `bio` splits the lines) as `bio`, `from_bio` and NCBI's
/// reader read it (`check_bio_deflines_of`).
struct BioReading<'a> {
    /// The line without the Unicode white space at its end (`bio`'s header).
    header: &'a str,
    /// Whether `bio`'s ID is empty.
    id_empty: bool,
    /// `from_bio`'s title.
    title: Vec<u8>,
    /// NCBI's line: `bio`'s line up to its first carriage return.
    ncbi_line: &'a str,
    /// NCBI's defline after `>`: the line without a `?_` prefix, `None` for a gap line
    /// (`>?`, also after a `?_` prefix).
    ncbi_defline: Option<&'a str>,
    /// Whether the rest of `bio`'s line, which NCBI reads as further lines, is C white
    /// space (lines that add nothing to the record).
    rest_is_blank: bool,
}

impl<'a> BioReading<'a> {
    fn of(defline: &'a str) -> Self {
        // `bio`'s ID and description (bio-1.6.0/src/io/fasta.rs:331-333).
        let header = defline.trim_end();
        let mut fields = header.splitn(2, char::is_whitespace);
        let id = fields.next().unwrap_or("");
        let title = bio_title(id, fields.next());
        let (ncbi_line, rest) = match defline.find('\r') {
            Some(at) => (&defline[..at], &defline[at + 1..]),
            None => (defline, ""),
        };
        let ncbi_defline = match ncbi_line.strip_prefix("?_") {
            Some(text) => (!text.starts_with('?')).then_some(text),
            None => (!ncbi_line.starts_with('?')).then_some(ncbi_line),
        };
        BioReading {
            header,
            id_empty: id.is_empty(),
            title,
            ncbi_line,
            ncbi_defline,
            rest_is_blank: rest.bytes().all(c_isspace),
        }
    }

    /// Whether `from_bio`'s record of the defline is NCBI's reader's record.
    fn alike(&self) -> bool {
        self.rest_is_blank
            && self
                .ncbi_defline
                .is_some_and(|defline| self.title == ncbi_title(defline.as_bytes()))
    }

    /// Why NCBI's reader makes another record of a defline to which `bio` gives an empty ID;
    /// `None` when the ID is not empty or the records are alike.
    fn empty_id_problem(&self) -> Option<String> {
        if !self.id_empty || self.alike() {
            return None;
        }
        // A carriage return inside `bio`'s header leaves the header's last character, which
        // is not white space, in the rest.
        if !self.rest_is_blank {
            return Some("has a carriage return before the end of its line".to_string());
        }
        // `bio` drops the character that ends the empty ID, and the whole of a line of white
        // space; NCBI keeps a non-ASCII one (its bytes are not C white space).
        let start = match self.header.chars().next() {
            Some(first) if first.is_ascii() => None,
            Some(first) => non_ascii_space(first.encode_utf8(&mut [0; 4])),
            None => non_ascii_space(self.ncbi_line),
        };
        if start.is_some() {
            return start;
        }
        // NCBI's title ends at a control character after its first byte.
        if let Some(&byte) = self.title.iter().skip(1).find(|&&byte| byte < b' ') {
            return Some(format!("has the control character 0x{byte:02x}"));
        }
        // `bio` drops the Unicode white space at the end of the line; NCBI keeps a non-ASCII
        // one before a control character. (The carriage return, if any, is after the header.)
        Some(
            non_ascii_space(&self.ncbi_line[self.header.len().min(self.ncbi_line.len())..])
                .unwrap_or_else(|| "is read with another title".to_string()),
        )
    }

    /// Why BLASTP's outfmt 0 title of a defline with an ID is neither ABI v1's nor NCBI's:
    /// `from_bio`'s record is not NCBI's, and its title (without the periods, commas,
    /// semicolons, tildes and spaces at its end that `GenerateDefline` strips) ends with a
    /// non-ASCII byte, which `x_CleanAndCompress` drops; `None` otherwise.
    fn outfmt0_title_problem(&self) -> Option<String> {
        if self.id_empty || self.alike() {
            return None;
        }
        let last = self
            .title
            .iter()
            .rposition(|byte| !b".,;~ ".contains(byte))
            .map(|at| self.title[at]);
        if !last.is_some_and(|byte| byte >= 0x80) {
            return None;
        }
        match self.ncbi_defline {
            None => return Some("starts with '?', a gap in the sequence".to_string()),
            Some(defline) if defline.len() != self.ncbi_line.len() => {
                return Some("starts with '?_', which NCBI BLAST+ drops".to_string())
            }
            Some(_) => {}
        }
        if !self.rest_is_blank {
            return Some("has a carriage return before the end of its line".to_string());
        }
        // NCBI's title ends at the first byte below a space (its first byte is the ID's);
        // where that is inside `bio`'s header, a control character ends it, and otherwise
        // `bio` split the ID or trimmed the end at a non-ASCII white space character.
        let line = self.ncbi_line.as_bytes();
        match line.iter().skip(1).position(|&byte| byte < b' ') {
            Some(at) if at + 1 < self.header.len() => {
                Some(format!("has the control character 0x{:02x}", line[at + 1]))
            }
            _ => Some(
                non_ascii_space(self.ncbi_line)
                    .unwrap_or_else(|| "is read with another title".to_string()),
            ),
        }
    }
}

/// The first non-ASCII white space character of `text`, as a problem of a defline.
fn non_ascii_space(text: &str) -> Option<String> {
    text.chars()
        .find(|&c| !c.is_ascii() && c.is_whitespace())
        .map(|c| {
            format!(
                "has the non-ASCII white space character U+{:04X}",
                u32::from(c)
            )
        })
}

/// The title of NCBI's reader for a defline's line (after `>`, up to its line end): none
/// for a line of C white space; otherwise from the first byte that is not C white space up
/// to the first byte below a space after it, without the white space at its end
/// (`check_bio_deflines_of`).
fn ncbi_title(line: &[u8]) -> Vec<u8> {
    let Some(start) = line.iter().position(|&byte| !c_isspace(byte)) else {
        return Vec::new();
    };
    let end = line[start + 1..]
        .iter()
        .position(|&byte| byte < b' ')
        .map_or(line.len(), |at| start + 1 + at);
    let title = &line[start..end];
    let kept = title
        .iter()
        .rposition(|&byte| !c_isspace(byte))
        .map_or(0, |at| at + 1);
    title[..kept].to_vec()
}

/// `from_bio`'s title for `bio`'s ID and description: the ID, a space and the
/// description, without the C white space at its start, which NCBI's defline parser skips
/// (fasta_reader_utils.cpp:209-213, quoted at `from_bio`).
fn bio_title(id: &str, desc: Option<&str>) -> Vec<u8> {
    let mut title = Vec::with_capacity(id.len() + 1);
    title.extend_from_slice(id.as_bytes());
    if let Some(desc) = desc {
        title.push(b' ');
        title.extend_from_slice(desc.as_bytes());
    }
    let start = title
        .iter()
        .position(|&byte| !c_isspace(byte))
        .unwrap_or(title.len());
    title.drain(..start);
    title
}

/// Rejects a sequence line with a non-ASCII byte. `bio` drops Unicode white space at the
/// end of a line (U+00A0 and others) and keeps other bytes; NCBI trims only ASCII white
/// space, so it reads such bytes as invalid residues (a warning, or an error on the first
/// data line of a record).
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:1004-1006
/// ```c
///             stringstream warn_strm;
///             warn_strm << "FASTA-Reader: Ignoring invalid " << x_NucOrProt()
///                 << "residues at position(s): ";
/// ```
/// Text before the first defline (a byte order mark, for example) is not a sequence line;
/// `bio` cannot read it, and the reader names it. `role` is `query` or `subject`.
pub fn check_sequence_lines(bytes: &[u8], role: &str) -> Result<()> {
    check_sequence_lines_of(bytes, role, "BLASTN")
}

/// `check_sequence_lines` for the FASTA input of `program`.
pub fn check_sequence_lines_of(bytes: &[u8], role: &str, program: &str) -> Result<()> {
    let mut record = 0;
    for line in bytes.split(|&byte| byte == b'\n') {
        if line.starts_with(b">") {
            record += 1;
        } else if record > 0 && !line.is_ascii() {
            bail!(
                "{role} record {record} has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's {program} (use IUPAC nucleotide letters)"
            );
        }
    }
    Ok(())
}

/// Whether the text of a defline (after `>`) makes NCBI's reader warn that the title ends
/// with nucleotides: it is longer than 20 bytes and its last 20 are unambiguous letters.
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:1624-1643
/// ```c
///     const static size_t kWarnNumNucCharsAtEnd = 20;
///     const static size_t kWarnAminoAcidCharsAtEnd = 50;
///
///     const size_t length = sLineText.length();
///     SIZE_TYPE pos_to_check = length-1;
///
///     if((length > kWarnNumNucCharsAtEnd) && !TestFlag(fAssumeProt)) {
///         // find last non-nuc character, within the last kWarnNumNucCharsAtEnd characters
///         const SIZE_TYPE last_pos_to_check_for_nuc = (sLineText.length() - kWarnNumNucCharsAtEnd);
///         for( ; pos_to_check >= last_pos_to_check_for_nuc; --pos_to_check ) {
///             if( ! s_ASCII_IsUnAmbigNuc(sLineText[pos_to_check]) ) {
///                 // found a character which is not an unambiguous nucleotide
///                 break;
///             }
///         }
///         if( pos_to_check < last_pos_to_check_for_nuc ) {
///             FASTA_WARNING(iLineNum,
///                 "FASTA-Reader: Title ends with at least " << kWarnNumNucCharsAtEnd
///                 << " valid nucleotide characters.  Was the sequence "
///                 << "accidentally put in the title line?",
/// ```
fn ends_with_nucleotides(text: &[u8]) -> bool {
    text.len() > 20
        && text[text.len() - 20..]
            .iter()
            .all(|byte| b"ACGTacgt".contains(byte))
}

/// Rejects the records with a residue that NCBI reads differently (see the module).
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935
/// ```c
///         case '-':
///             char_type = (
///                 bHyphensAreGaps ? eCharType_Gap :
///                 bHyphensIgnoreAndWarn ? eCharType_HyphenToIgnoreAndWarn :
///                 eCharType_NormalNonGap );
///             break;
///         case ';':
///             char_type = eCharType_Comment;
///             break;
///
///         case '\t': case '\n': case '\v': case '\f': case '\r': case ' ':
///             continue;
///
///         default:
///             char_type = eCharType_Bad;
///             break;
/// ```
/// `role` is `query` or `subject`; `program` is named in the message.
pub fn check_residues_of(records: &[fasta::Record], role: &str, program: &str) -> Result<()> {
    first_problem(records, |index, record| {
        invalid_residue(index, record, role, program)
    })
}

/// The residue check of `check_protein_input_of` for the `bio` records already read.
pub fn check_protein_residues_of(
    records: &[fasta::Record],
    role: &str,
    program: &str,
) -> Result<()> {
    first_problem(records, |index, record| {
        let position = record
            .seq()
            .iter()
            .position(|&byte| !(byte.is_ascii_alphabetic() || byte == b'*'))?;
        let byte = record.seq()[position];
        let shown = if byte.is_ascii_graphic() {
            format!("'{}'", byte as char)
        } else {
            format!("byte 0x{byte:02x}")
        };
        Some(format!(
            "{role} record {} ({}) has the residue {shown} at position {}; NCBI BLAST+ removes such a character from a protein sequence, which is not supported by LOSAT's {program} (use letters and '*')",
            index + 1,
            record_id(record),
            position + 1
        ))
    })
}

/// Rejects a record without residues: NCBI reads it without a message and reports it when
/// it sets up the search ("Sequence contains no data"), which LOSAT does not reproduce.
/// `role` is `query` or `subject`; `program` is named in the message.
pub fn check_records_have_residues_of(
    records: &[fasta::Record],
    role: &str,
    program: &str,
) -> Result<()> {
    first_problem(records, |index, record| {
        no_residues(index, record, role, program)
    })
}

/// Both checks, record by record (the order of ABI v1, which plan TD-1 freezes).
pub fn check_records(records: &[fasta::Record], role: &str) -> Result<()> {
    first_problem(records, |index, record| {
        no_residues(index, record, role, "BLASTN")
            .or_else(|| invalid_residue(index, record, role, "BLASTN"))
    })
}

fn first_problem<R>(records: &[R], problem: impl Fn(usize, &R) -> Option<String>) -> Result<()> {
    match records
        .iter()
        .enumerate()
        .find_map(|(index, record)| problem(index, record))
    {
        Some(message) => bail!("{message}"),
        None => Ok(()),
    }
}

fn no_residues(index: usize, record: &fasta::Record, role: &str, program: &str) -> Option<String> {
    record.seq().is_empty().then(|| {
        format!(
            "{role} record {} ({}) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's {program}",
            index + 1,
            record_id(record)
        )
    })
}

/// The ID that LOSAT's messages name a record by: its `bio` ID.
fn record_id(record: &fasta::Record) -> String {
    record.id().to_string()
}

fn invalid_residue(
    index: usize,
    record: &fasta::Record,
    role: &str,
    program: &str,
) -> Option<String> {
    let position = record
        .seq()
        .iter()
        .position(|&byte| !IUPAC_NUCLEOTIDE[byte as usize])?;
    let byte = record.seq()[position];
    let shown = if byte.is_ascii_graphic() {
        format!("'{}'", byte as char)
    } else {
        format!("0x{byte:02x}")
    };
    Some(format!(
        "{role} record {} ({}) has {shown} at residue {}, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not supported by LOSAT's {program} (use IUPAC nucleotide letters)",
        index + 1,
        record.id(),
        position + 1
    ))
}

/// The record of a `bio` record (the bridge of the port plan, §1): only where LOSAT's
/// checks guarantee that `bio` reads the input as NCBI's reader does (ABI v1, plan TD-1;
/// ABI v2 read its inputs so until session SFc, S10). `n` is the `N` of the local ID (`Query_N`,
/// `Subject_N`, 1-based) and `prefix` its prefix; `protein` is the molecule of the input
/// (`fAssumeProt`, else `fAssumeNuc`).
///
/// The title is the defline as `bio` splits it (the ID, a space and the description),
/// without the white space at its start, which NCBI's defline parser skips (`bio` gives
/// such a defline an empty ID); that is NCBI's title for the deflines those checks
/// accept. The messages are the title warning that NCBI writes at the end of such a
/// record.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-213
/// ```c++
///     // trim leading whitespace from title (is this appropriate?)
///     while (title_start < len
///         &&  isspace((unsigned char)defline[title_start])) {
///         ++title_start;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038-2043
/// ```c++
///     NStr::TruncateSpacesInPlace(processed_title);
///     if (!processed_title.empty()) {
///         auto pDesc = Ref(new CSeqdesc());
///         pDesc->SetTitle() = processed_title;
///         bioseq.SetDescr().Set().push_back(std::move(pDesc));
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:469-482
/// ```c++
///     auto n = m_Counter.load();
///     if (advance)
///         m_Counter++;
///
///     if (m_Prefix.empty()  &&  m_Suffix.empty()) {
///         seq_id->SetLocal().SetId(n);
///     } else {
///         string& id = seq_id->SetLocal().SetStr();
///         id.reserve(128);
///         id += m_Prefix;
///         id += NStr::IntToString(n);
///         id += m_Suffix;
///     }
///     return seq_id;
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447
/// ```c++
///         CRef<CSeq_data> data(new CSeq_data(m_SeqData, format));
///         if ( !TestFlag(fLeaveAsText) ) {
///             CSeqportUtil::Pack(data, inst.GetLength());
///         }
/// ```
/// A nucleotide `U` is stored as `T` (`u` as `t` inside the lowercase mask), as the
/// reader stores it (`reader.rs`, `assemble_seq`).
pub(crate) fn from_bio(
    record: &fasta::Record,
    n: usize,
    prefix: &str,
    protein: bool,
) -> FastaRecord {
    let title = bio_title(record.id(), record.desc());
    let mut sequence = record.seq().to_vec();
    if !protein {
        for residue in sequence.iter_mut() {
            *residue = match *residue {
                b'U' => b'T',
                b'u' => b't',
                other => other,
            };
        }
    }
    let mut warnings = Vec::new();
    if let Some(warning) =
        crate::blastinput::fasta_reader::seq_data_in_title_warning(&title, protein)
    {
        warnings.extend_from_slice(warning);
        warnings.push(b'\n');
    }
    FastaRecord {
        local_id: format!("{prefix}{n}"),
        title,
        sequence,
        warnings,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};

    /// A search that shows the names of records without a title, and nothing else that
    /// `check_bio_deflines_of` checks.
    const NAMES: Shown = Shown {
        names: true,
        local_ids: false,
        outfmt0_titles: false,
    };

    // The nucleotide 20-letter check of the title reads the trailing white space.
    #[test]
    fn deflines_ending_with_nucleotides_and_white_space_are_rejected() {
        assert!(check_deflines_of(
            format!(">id {} \nAC\n", "ACGT".repeat(6)).as_bytes(),
            "query",
            "BLASTN"
        )
        .is_err());
    }

    fn records(text: &str) -> Vec<fasta::Record> {
        fasta::Reader::new(text.as_bytes())
            .records()
            .collect::<std::result::Result<_, _>>()
            .unwrap()
    }

    #[test]
    fn only_inputs_read_as_ncbi_reads_them_are_accepted() {
        let text = ">q1 a b\nACGTNRY\nacgtn\n>q2\t\nU\r\n";
        assert!(check_deflines(text.as_bytes(), "query").is_ok());
        assert!(check_residues_of(&records(text), "query", "BLASTN").is_ok());
        for (text, problem) in [
            (">q1\tb\nACGT\n", "control character 0x09"),
            (
                ">q0\nA\n> q1\nACGT\n",
                "record 2 has a defline that starts with white space",
            ),
            (">q1 \u{e9}\nACGT\n", "non-ASCII"),
            (">\nACGT\n", "record 1 has a defline that is empty"),
            (
                ">q0\nA\n>?100\nACGT\n",
                "record 2 has a defline that starts with '?'",
            ),
            (">?_q1\nACGT\n", "starts with '?'"),
            (
                ">q1 ACGTACGTACGTACGTACGTA \nACGT\n",
                "ends with white space after 20",
            ),
            (
                ">q0\nA\n>   \nACGT\n",
                "record 2 has a defline that is empty",
            ),
        ] {
            let error = check_deflines(text.as_bytes(), "query")
                .unwrap_err()
                .to_string();
            assert!(error.contains(problem), "{text:?}: {error}");
            assert!(error.contains("not supported by LOSAT"), "{error}");
        }
        for (text, problem) in [
            (">q1\nACXGT\n", "'X' at residue 3"),
            (">q1\nAC-GT\n", "'-' at residue 3"),
            (">q1\nAC GT\n", "0x20 at residue 3"),
        ] {
            let error = check_residues_of(&records(text), "query", "BLASTN")
                .unwrap_err()
                .to_string();
            assert!(error.contains(problem), "{text:?}: {error}");
            assert!(error.contains("not supported by LOSAT"), "{error}");
        }
        let error =
            check_records_have_residues_of(&records(">q1 empty\n>q2\nACGT\n"), "query", "BLASTN")
                .unwrap_err()
                .to_string();
        assert!(
            error.contains("query record 1 (q1) has no residues"),
            "{error}"
        );
        // Record by record, as ABI v1: record 1's letter before record 2's emptiness.
        let both = records(">q1\nAXG\n>q2 empty\n");
        assert!(check_records(&both, "query")
            .unwrap_err()
            .to_string()
            .contains("query record 1 (q1) has 'X' at residue 2"));
        assert!(check_records_have_residues_of(&both, "query", "BLASTN")
            .unwrap_err()
            .to_string()
            .contains("query record 2 (q2) has no residues"));
    }

    #[test]
    fn titles_that_end_with_nucleotides_are_accepted() {
        // NCBI 2.17.0: the text after `>` longer than 20 bytes, its last 20 ACGT, gets the
        // reader's title warning (`fasta_reader`); the defline itself is read alike.
        let text = ">q1 ACGTACGTACGTACGTACGTA\nACGT\n>ACGTACGTACGTACGTACGT\nACGT\n>ACGTACGTACGTACGTACGTA\nACGT\n>q2 ACGTNACGTACGTACGTACGTA\nACGT\n";
        assert!(check_deflines(text.as_bytes(), "query").is_ok());
        assert!(check_deflines(b">q1 ACGTACGTACGTACGTACGTA\r\nACGT\r\n", "query").is_ok());
        let warned = text
            .lines()
            .filter_map(|line| line.strip_prefix('>'))
            .filter(|defline| ends_with_nucleotides(defline.as_bytes()))
            .count();
        assert_eq!(warned, 2);
    }

    #[test]
    fn non_ascii_bytes_in_sequence_lines_are_rejected() {
        assert!(check_sequence_lines(b">q1 d\nACGT\r\nacgt\n\n>q2\nNN\n", "query").is_ok());
        // A byte order mark before the first defline is left to the reader's rejection.
        assert!(check_sequence_lines(b"\xef\xbb\xbf>q1\nACGT\n", "query").is_ok());
        for (bytes, record) in [
            (&b">q1\nACGT\xc2\xa0\n"[..], "query record 1"),
            (b">q1\nACGT\n>q2\nAC\n\xe3\x80\x80\n", "query record 2"),
            (b">q1\n\xef\xbb\xbfACGT\n", "query record 1"),
        ] {
            let error = check_sequence_lines(bytes, "query")
                .unwrap_err()
                .to_string();
            assert!(error.contains(record), "{error}");
            assert!(error.contains("not supported by LOSAT"), "{error}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447, 1616-1679,
    // 2038-2043 and fasta_reader_utils.cpp:466-483 (`from_bio` gives the reader's record for
    // the inputs that today's checks accept: the title, `U` stored as `T` for nucleotides,
    // the lowercase letters, the local ID and the title warning).
    #[test]
    fn from_bio_gives_the_readers_record_for_accepted_inputs() {
        let fifty = "A".repeat(50);
        let inputs: [(&[u8], bool); 6] = [
            (b">q1 a description\nACGTacgtNNRY\n", false),
            (b">q1\nACGUacguTT\nAC\n", false),
            (b">q1 x ACGTACGTACGTACGTACGTA\nACGT\n", false),
            (b">q1  two  spaces\nAC\n>q2 second\nGG\n", false),
            (b">p1 desc\nMKVLUuX*\n", true),
            (format!(">p1 q{fifty}\nMKV\n").leak().as_bytes(), true),
        ];
        for (bytes, protein) in inputs {
            let config = ReaderConfig::query("BLASTN", protein, false);
            let reader = read_all(
                &mut FastaInputSource::from_bytes(bytes, config),
                &mut |_| Ok(()),
            )
            .unwrap();
            let bio: Vec<FastaRecord> = fasta::Reader::new(bytes)
                .records()
                .enumerate()
                .map(|(index, record)| from_bio(&record.unwrap(), index + 1, "Query_", protein))
                .collect();
            assert_eq!(bio, reader, "{:?}", String::from_utf8_lossy(bytes));
        }
        let record = from_bio(
            &fasta::Record::with_attrs("s", None, b"ACGU"),
            7,
            "Subject_",
            false,
        );
        assert_eq!(
            (record.local_id.as_str(), record.sequence.as_slice()),
            ("Subject_7", &b"ACGT"[..])
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-213
    // ```c++
    //     // trim leading whitespace from title (is this appropriate?)
    //     while (title_start < len
    //         &&  isspace((unsigned char)defline[title_start])) {
    //         ++title_start;
    //     }
    // ```
    // ABI v1's `bio` bridge (`from_bio`): the title is the ID,
    // a space and the description, without the white space at its start (a defline that starts
    // with white space has an empty `bio` ID), and a nucleotide `U` is `T`.
    #[test]
    fn the_bio_bridge_keeps_ncbis_titles() {
        let bio_records = [
            fasta::Record::with_attrs("q1", Some("first query"), b"ACGTACGTAC"),
            fasta::Record::with_attrs("q2", None, b"AC"),
            fasta::Record::with_attrs("", Some("q3 x"), b"acgu"),
        ];
        let records: Vec<FastaRecord> = bio_records
            .iter()
            .enumerate()
            .map(|(index, record)| from_bio(record, index + 1, "Query_", false))
            .collect();
        assert_eq!(
            records,
            vec![
                FastaRecord::new("Query_1", b"q1 first query", b"ACGTACGTAC"),
                FastaRecord::new("Query_2", b"q2", b"AC"),
                FastaRecord::new("Query_3", b"q3 x", b"acgt"),
            ]
        );
        assert_eq!(records[2].shown_id(), b"q3");
        let protein = from_bio(&bio_records[2], 1, "Subject_", true);
        assert_eq!(protein, FastaRecord::new("Subject_1", b"q3 x", b"acgu"));
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-213
    // ```c++
    //     // trim leading whitespace from title (is this appropriate?)
    //     while (title_start < len
    //         &&  isspace((unsigned char)defline[title_start])) {
    //         ++title_start;
    //     }
    // ```
    // SFd re-audit D-1 (round 1): BLASTP and TBLASTX reject the deflines with an empty
    // `bio` ID that NCBI reads as another title: a non-ASCII white space character at the
    // start (or as the whole defline, or at the end after white space at the start), a
    // carriage return inside the line, a control character after white space at the start.
    #[test]
    fn empty_id_deflines_that_ncbi_reads_otherwise_are_rejected() {
        for (defline, problem) in [
            ("\u{a0}x y", "the non-ASCII white space character U+00A0"),
            ("\u{a0}", "the non-ASCII white space character U+00A0"),
            (
                "\u{3000}\u{3000}x",
                "the non-ASCII white space character U+3000",
            ),
            ("\u{85}x", "the non-ASCII white space character U+0085"),
            ("\u{2028}", "the non-ASCII white space character U+2028"),
            (" \u{a0}", "the non-ASCII white space character U+00A0"),
            ("\x0c\u{1680}", "the non-ASCII white space character U+1680"),
            (" x\u{a0}", "the non-ASCII white space character U+00A0"),
            (" x\u{202f}\r", "the non-ASCII white space character U+202F"),
            ("\rx y", "a carriage return before the end of its line"),
            (" \rx", "a carriage return before the end of its line"),
            ("\r\x0b x", "a carriage return before the end of its line"),
            (" x\ry", "a carriage return before the end of its line"),
            ("\r\u{a0}", "a carriage return before the end of its line"),
            (
                " x \r\u{a0}",
                "a carriage return before the end of its line",
            ),
            ("\r>x", "a carriage return before the end of its line"),
            (" x\ty", "the control character 0x09"),
            ("\tx\x01y", "the control character 0x01"),
            ("  x\x0b y", "the control character 0x0b"),
        ] {
            let bytes = format!(">q0\nMK\n>{defline}\nMKV\n");
            let error = check_bio_deflines_of(bytes.as_bytes(), "subject", "BLASTP", NAMES)
                .unwrap_err()
                .to_string();
            assert_eq!(
                error,
                format!("subject record 2 has a defline that starts with white space and has {problem}; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTP (begin the defline with a character that is not white space)"),
                "{defline:?}"
            );
        }
        // The first record that NCBI reads otherwise is named, in a CRLF file too.
        let error = check_bio_deflines_of(
            b">\xc2\xa0a\r\nMK\r\n>\xe3\x80\x80\r\nMK\r\n",
            "query",
            "TBLASTX",
            NAMES,
        )
        .unwrap_err()
        .to_string();
        assert!(
            error.starts_with("query record 1 has a defline that starts with white space and has the non-ASCII white space character U+00A0;"),
            "{error}"
        );
        assert!(error.contains("LOSAT's TBLASTX"), "{error}");
        // A record left without a title is rejected only where its name shows.
        let untitled = Shown::default();
        assert!(check_bio_deflines_of(b">\xc2\xa0\nMK\n", "query", "BLASTP", untitled).is_ok());
        assert!(check_bio_deflines_of(b">\xc2\xa0\nMK\n", "query", "BLASTP", NAMES).is_err());
        assert!(check_bio_deflines_of(b">\xc2\xa0x\nMK\n", "query", "BLASTP", untitled).is_err());
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:376-380
    // ```c++
    //         CTempString line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    //
    //         if (line.empty()) {
    //             continue; // ignore lines containing only whitespace
    //         }
    // ```
    // The deflines whose v1 records are NCBI's are accepted (README gate item 12): empty,
    // C white space (space, tab, VT, FF; carriage returns included), a title after C white
    // space (also when a non-ASCII white space character follows it), a CRLF line end; and
    // every defline with an ID, which `bio` reads as ABI v1 did.
    #[test]
    fn empty_id_deflines_that_ncbi_reads_alike_are_accepted() {
        for defline in [
            "",
            " ",
            "   ",
            "\t",
            "\x0b",
            "\x0c",
            "\x0c ",
            " \t\x0b ",
            "\r",
            "\r\r",
            " \r ",
            "  x y",
            " x",
            "\tx y",
            "\x0bx y",
            "\x0cx y",
            " \t x",
            " \u{a0}x",
            "\x0b\u{3000}x y",
            " x\r",
            " x \r ",
            " x\t",
            " x\t\u{a0}",
            " \x01x",
            " x\u{7f}",
            " \u{e9} x \u{e9}",
            "x\ry z",
            "x\u{a0}",
            "x\u{a0}y z",
            "\u{feff}x y",
            "x\ty",
        ] {
            let bytes = format!(">{defline}\nMKV\n");
            assert!(
                check_bio_deflines_of(bytes.as_bytes(), "query", "BLASTP", NAMES).is_ok(),
                "{defline:?}"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038
    // ```c++
    //     NStr::TruncateSpacesInPlace(processed_title);
    // ```
    // The check accepts a defline with an empty `bio` ID exactly when `from_bio` gives the
    // local IDs, titles and residues of NCBI's reader (`fasta_reader`), for every defline of
    // up to four characters of white space (ASCII and not), carriage returns, a letter, a
    // control character, `>` and `;`, in a protein and a nucleotide input. (The reader's
    // messages go to standard error, which ABI v1 does not return.)
    #[test]
    fn the_empty_id_check_accepts_exactly_the_records_of_ncbis_reader() {
        let alphabet = [
            " ", "\t", "\x0b", "\x0c", "\r", "\u{a0}", "\u{3000}", "\u{85}", "x", "\x01", ">", ";",
        ];
        let mut deflines = vec![String::new()];
        let mut last = vec![String::new()];
        for _ in 0..4 {
            last = last
                .iter()
                .flat_map(|prefix| alphabet.iter().map(move |c| format!("{prefix}{c}")))
                .collect();
            deflines.extend(last.iter().cloned());
        }
        let (mut accepted, mut rejected) = (0, 0);
        for defline in &deflines {
            if defline.chars().next().is_some_and(|c| !c.is_whitespace()) {
                continue;
            }
            for (protein, residues) in [(true, "MKV"), (false, "ACGT")] {
                let text = format!(">{defline}\n{residues}\n");
                let bio: Vec<(String, Vec<u8>, Vec<u8>)> = fasta::Reader::new(text.as_bytes())
                    .records()
                    .enumerate()
                    .map(|(index, record)| {
                        let record = from_bio(&record.unwrap(), index + 1, "Query_", protein);
                        (record.local_id, record.title, record.sequence)
                    })
                    .collect();
                let config = ReaderConfig::query("BLASTP", protein, false);
                let ncbi = read_all(
                    &mut FastaInputSource::from_bytes(text.as_bytes(), config),
                    &mut |_| Ok(()),
                )
                .map(|records| {
                    records
                        .into_iter()
                        .map(|record| (record.local_id, record.title, record.sequence))
                        .collect::<Vec<_>>()
                });
                let alike = ncbi.as_ref().is_ok_and(|ncbi| *ncbi == bio);
                let checked = check_bio_deflines_of(text.as_bytes(), "query", "BLASTP", NAMES);
                assert_eq!(
                    checked.is_ok(),
                    alike,
                    "{defline:?} protein {protein}: {checked:?}"
                );
                if alike {
                    accepted += 1;
                } else {
                    rejected += 1;
                }
            }
        }
        assert!(accepted > 1000 && rejected > 1000, "{accepted} {rejected}");
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-362
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
    // ```
    // A record without a title after a `>?` line (not the first line) is numbered lower
    // by NCBI's reader than by `bio`; where the search shows the local IDs, it is rejected.
    #[test]
    fn untitled_records_after_gap_lines_are_rejected_where_their_numbers_show() {
        let ids = Shown {
            names: true,
            local_ids: true,
            outfmt0_titles: false,
        };
        let error =
            check_bio_deflines_of(b">s1\nMK\n>?100\nMK\n>\nMK\n", "subject", "TBLASTX", ids)
                .unwrap_err()
                .to_string();
        assert_eq!(
            error,
            "subject record 3 has a defline without a title after a line that starts with '>?', which NCBI BLAST+ reads as a gap, not as a record (NCBI BLAST+ numbers this record 2); this is not supported by LOSAT's TBLASTX (give the record a title)"
        );
        for (text, rejected) in [
            (">s1\nMK\n>?100\nMK\n>\nMK\n", true),
            (">s1\nMK\n>?100\nMK\n> \t\r\nMK\n", true),
            (">s1\nMK\n>?_?5\nMK\n>\nMK\n", true),
            (">s1\nMK\n>?\nMK\n>s2\nMK\n>\nMK\n", true),
            // A first `>?` line opens a record without a title for NCBI: no shift.
            (">?5\nMK\n>\nMK\n", false),
            // `>?_x` is a defline for NCBI.
            (">s1\nMK\n>?_x\nMK\n>\nMK\n", false),
            // A record with a title is shown by its first word.
            (">s1\nMK\n>?5\nMK\n> x\nMK\n", false),
            (">\nMK\n>?5\nMK\n>s2\nMK\n", false),
        ] {
            assert_eq!(
                check_bio_deflines_of(text.as_bytes(), "query", "BLASTP", ids).is_err(),
                rejected,
                "{text:?}"
            );
            assert!(
                check_bio_deflines_of(text.as_bytes(), "query", "BLASTP", Shown::default()).is_ok()
            );
        }
    }

    // The gap rule against NCBI's reader (`fasta_reader`): for every input of one to four
    // records with the deflines below, the check rejects exactly when some record without
    // a title has a local ID that no record without a title of NCBI's reader has.
    #[test]
    fn the_gap_rule_follows_the_local_ids_of_ncbis_reader() {
        let deflines = ["s1", "", " ", "?5", "?", "?_x", "?_?3", " x", "?unk100"];
        let ids = Shown {
            names: true,
            local_ids: true,
            outfmt0_titles: false,
        };
        let mut inputs = vec![Vec::<&str>::new()];
        let mut last = inputs.clone();
        for _ in 0..4 {
            last = last
                .iter()
                .flat_map(|prefix| {
                    deflines.iter().map(move |defline| {
                        let mut next = prefix.clone();
                        next.push(*defline);
                        next
                    })
                })
                .collect();
            inputs.extend(last.iter().cloned());
        }
        let (mut accepted, mut rejected) = (0, 0);
        for records in inputs.iter().filter(|records| !records.is_empty()) {
            let text: String = records
                .iter()
                .map(|defline| format!(">{defline}\nMKV\n"))
                .collect();
            let untitled = |records: &[FastaRecord]| -> Vec<String> {
                records
                    .iter()
                    .filter(|record| record.title.is_empty())
                    .map(|record| record.local_id.clone())
                    .collect()
            };
            let bio: Vec<FastaRecord> = fasta::Reader::new(text.as_bytes())
                .records()
                .enumerate()
                .map(|(index, record)| from_bio(&record.unwrap(), index + 1, "Query_", true))
                .collect();
            let config = ReaderConfig::query("BLASTP", true, false);
            let Ok(ncbi) = read_all(
                &mut FastaInputSource::from_bytes(text.as_bytes(), config),
                &mut |_| Ok(()),
            ) else {
                continue;
            };
            let ncbi_untitled = untitled(&ncbi);
            let shifted = untitled(&bio).iter().any(|id| !ncbi_untitled.contains(id));
            let checked = check_bio_deflines_of(text.as_bytes(), "query", "BLASTP", ids);
            assert_eq!(checked.is_err(), shifted, "{text:?}: {checked:?}");
            if shifted {
                rejected += 1;
            } else {
                accepted += 1;
            }
        }
        assert!(accepted > 1000 && rejected > 1000, "{accepted} {rejected}");
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:304-306
    // ```c++
    //     if (curr > 0 && curr != ' ') {
    //         *out++ = curr;
    //     }
    // ```
    // BLASTP's outfmt 0 title of a defline with an ID that `bio` reads otherwise is rejected
    // when it ends with a non-ASCII character (the shared report drops its last byte, as
    // NCBI does for its own title, and ABI v1's report kept it).
    #[test]
    fn outfmt0_titles_that_bio_reads_otherwise_and_that_end_with_non_ascii_are_rejected() {
        let titles = Shown {
            names: true,
            local_ids: false,
            outfmt0_titles: true,
        };
        for (defline, problem) in [
            ("s1\tx \u{e9}", "it has the control character 0x09"),
            ("s1 x\x01\u{e9}", "it has the control character 0x01"),
            (
                "s1\u{a0}x \u{e9}",
                "it has the non-ASCII white space character U+00A0",
            ),
            (
                "s1 \u{e9}\u{a0}",
                "it has the non-ASCII white space character U+00A0",
            ),
            ("s1\tx \u{e9}.,; ", "it has the control character 0x09"),
            (
                "s1 \u{e9}\r\u{e9}",
                "it has a carriage return before the end of its line",
            ),
            ("s1\tx\u{3000}y\u{e9}", "it has the control character 0x09"),
            ("?;\u{e9}", "it starts with '?', a gap in the sequence"),
            ("?_x\u{e9}", "it starts with '?_', which NCBI BLAST+ drops"),
        ] {
            let text = format!(">q0\nMK\n>{defline}\nMKV\n");
            let error = check_bio_deflines_of(text.as_bytes(), "subject", "BLASTP", titles)
                .unwrap_err()
                .to_string();
            assert_eq!(
                error,
                format!("subject record 2 has a defline that NCBI BLAST+ reads with another title ({problem}) and that ends with a non-ASCII character, whose last byte NCBI BLAST+'s outfmt 0 drops; this is not supported by LOSAT's BLASTP outfmt 0 (use ASCII deflines without control characters)"),
                "{defline:?}"
            );
            assert!(
                check_bio_deflines_of(text.as_bytes(), "subject", "BLASTP", Shown::default())
                    .is_ok()
            );
        }
        // Read alike, or ending with an ASCII character (ABI v1's bytes, as before).
        for defline in [
            "s1 x \u{e9}",
            "s1 \u{e9}\t",
            "s1\tx",
            "s1\tx \u{e9} z",
            "\u{e9}\ty",
            "?x",
            "?_x y",
            "s1\u{a0}x",
            "x\ry z",
        ] {
            let text = format!(">{defline}\nMKV\n");
            assert!(
                check_bio_deflines_of(text.as_bytes(), "subject", "BLASTP", titles).is_ok(),
                "{defline:?}"
            );
        }
        // The shared report drops the last byte of a non-ASCII title end, and of nothing else.
        let shown = |title: &str| {
            crate::report::defline::generate_defline(title.as_bytes(), true, true)
                .unwrap()
                .text
        };
        assert_eq!(shown("s1 x \u{e9}"), b"s1 x \xc3");
        assert_eq!(shown("s1 x \u{e9}.;"), b"s1 x \xc3");
        assert_eq!(shown("s1 x \u{e9} z"), "s1 x \u{e9} z".as_bytes());
    }
}
