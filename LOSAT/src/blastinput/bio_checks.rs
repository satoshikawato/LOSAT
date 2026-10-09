//! The checks of `bio` FASTA records that ABI v1 (plan TD-1, which freezes its accepted
//! inputs and messages) and the Web adapter's `register` (until step S10 of the port plan)
//! make before they search records made from them (`FastaRecord::from_bio`). LOSAT's
//! programs read their inputs with NCBI's reader (`fasta_reader`) and make none of them.
//!
//! LOSAT's ABI v1 parses FASTA with `bio` (and the Web adapter's index scan reproduces it),
//! while NCBI BLAST+ reads it with `CFastaReader`. For the inputs accepted here the two
//! agree, except that `U` is read as `T`, which `FastaRecord::from_bio` applies. The others
//! are rejected, because the reader of NCBI would change them:
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

use anyhow::{bail, Result};
use bio::io::fasta;

use crate::blastinput::input_files::is_blank;

/// Why a FASTA file that `bio` cannot parse is not read: `bio` fails on text before the
/// first defline (blank lines, `;` comments, a byte order mark) and on bytes that are not
/// UTF-8, which NCBI reads.
pub fn unreadable_fasta(program: &str) -> String {
    format!("FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-8), which NCBI BLAST+ may read, is not supported by LOSAT's {program}")
}

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
    check_deflines_with(bytes, role, program, false)
}

/// `check_deflines_of` for protein FASTA input: NCBI's protein reader checks the end of a
/// title for 50 amino-acid letters instead of 20 nucleotides (the title warning of
/// `fasta_reader`), so the white space at the end of a defline matters after 50 letters.
pub fn check_protein_deflines_of(bytes: &[u8], role: &str, program: &str) -> Result<()> {
    check_deflines_with(bytes, role, program, true)
}

fn check_deflines_with(bytes: &[u8], role: &str, program: &str, protein: bool) -> Result<()> {
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
        } else if raw.len() != defline.len() && !protein && ends_with_nucleotides(defline) {
            "ends with white space after 20 nucleotide letters (NCBI BLAST+'s warning about the letters depends on that white space, which LOSAT's reader drops)".to_string()
        } else if raw.len() != defline.len() && protein && ends_with_amino_acids(defline) {
            "ends with white space after 50 amino-acid letters (NCBI BLAST+'s warning about the letters depends on that white space, which LOSAT's reader drops)".to_string()
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

/// Whether the text of a defline (after `>`) makes NCBI's protein reader warn that the title
/// ends with amino acids: it is longer than 50 bytes and its last 50 are ASCII letters
/// (the title warning of `fasta_reader`).
fn ends_with_amino_acids(text: &[u8]) -> bool {
    text.len() > 50 && text[text.len() - 50..].iter().all(u8::is_ascii_alphabetic)
}

/// Rejects a protein sequence line with a byte that NCBI's reader removes with a warning
/// ("FASTA-Reader: Ignoring invalid residues", or "CFastaReader: Hyphens are invalid" for
/// `-`): any byte but a letter, `*` and the white space and `;` comment that NCBI skips
/// without a message. NCBI writes the warning when it reads the file, so LOSAT rejects the
/// file there (before `Query is Empty!` and the option checks); the bytes NCBI skips
/// silently are rejected later, with the other deferred checks (`check_protein_input_of`).
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:920-931
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
/// ```
/// `role` is `query` or `subject`.
pub fn check_protein_sequence_lines_of(bytes: &[u8], role: &str, program: &str) -> Result<()> {
    let mut record = 0;
    for line in bytes.split(|&byte| byte == b'\n') {
        if line.starts_with(b">") {
            record += 1;
            continue;
        }
        // NCBI reference: c++/src/objtools/readers/fasta.cpp:376-385
        // ```c
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
        // ```
        let is_space = |byte: &u8| matches!(byte, b'\t' | b'\r' | 0x0b | 0x0c | b' ');
        let start = line.iter().position(|byte| !is_space(byte));
        if record == 0 || start.is_none_or(|start| matches!(line[start], b'!' | b'#' | b';')) {
            continue;
        }
        for &byte in line {
            match byte {
                b';' => break,
                b'\t' | b'\r' | 0x0b | 0x0c | b' ' | b'*' => {}
                byte if byte.is_ascii_alphabetic() => {}
                byte => {
                    let shown = if byte.is_ascii_graphic() {
                        format!("'{}'", byte as char)
                    } else {
                        format!("byte 0x{byte:02x}")
                    };
                    bail!(
                        "{role} record {record} has {shown} in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by LOSAT's {program} (use letters and '*')"
                    );
                }
            }
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

/// Rejects the protein records whose deflines or residues NCBI reads differently: the
/// deflines of `check_deflines_of`, and a residue that is not an ASCII letter or `*`
/// (NCBI's protein reader removes the other bytes, with a warning for each line).
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:967-979
/// ```c
///         case eCharType_HyphenToIgnoreAndWarn:
///             bIgnorableHyphenSeen = true;
///             break;
///         case eCharType_Comment:
///             // artificially advance pos to the end to break the pos loop
///             pos = s_len;
///             break;
///         case eCharType_Bad:
///             if( bad_pos_line_num < 0 ) {
///                 bad_pos_line_num = LineNumber();
///             }
///             bad_pos_vec.push_back(pos);
///             break;
/// ```
/// `role` is `query` or `subject`.
pub fn check_protein_input_of(
    bytes: &[u8],
    records: &[fasta::Record],
    role: &str,
    program: &str,
) -> Result<()> {
    if is_blank(bytes) {
        return Ok(());
    }
    check_protein_deflines_of(bytes, role, program)?;
    check_protein_residues_of(records, role, program)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_lines_follow_ncbis_protein_reader() {
        // NCBI fasta.cpp:376-385,920-931: white space and `;` comments are skipped without
        // a message; other non-letters are removed with a warning (LOSAT rejects them).
        for text in [
            ">p1\nACDE FGH\t*\n; comment 123\n#x 1\n!y\nKLM;9\n",
            ">p\n\n  \n",
        ] {
            assert!(
                check_protein_sequence_lines_of(text.as_bytes(), "query", "BLASTP").is_ok(),
                "{text:?}"
            );
        }
        for text in [
            ">p1\nAC1DE\n",
            ">p1\nAC-DE\n",
            ">p1\nACDE\u{a0}\n",
            ">p1\nAC.DE\n",
        ] {
            let error =
                check_protein_sequence_lines_of(text.as_bytes(), "subject", "TBLASTN").unwrap_err();
            assert!(
                format!("{error:#}").contains("not supported by LOSAT's TBLASTN"),
                "{text:?}"
            );
        }
        // NCBI fasta.cpp:1651-1673: the 50-letter check reads the trailing white space.
        let letters = "A".repeat(50);
        assert!(check_protein_deflines_of(
            format!(">id {letters}\nAC\n").as_bytes(),
            "query",
            "BLASTP"
        )
        .is_ok());
        assert!(check_protein_deflines_of(
            format!(">id {letters} \nAC\n").as_bytes(),
            "query",
            "BLASTP"
        )
        .is_err());
        // The nucleotide 20-letter check does not apply to protein deflines.
        assert!(check_protein_deflines_of(
            format!(">id {} \nAC\n", "ACGT".repeat(6)).as_bytes(),
            "query",
            "BLASTP"
        )
        .is_ok());
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
}
