//! The FASTA records that BLASTN reads as NCBI BLAST+ does.
//!
//! LOSAT parses FASTA with `bio` (and the Web adapter's index scan reproduces it), while
//! NCBI BLAST+ reads it with `CFastaReader`. For the inputs accepted here the two agree,
//! except that `U` is read as `T`, which `with_u_as_t` applies. The others are rejected,
//! because the reader of NCBI would change them:
//!
//! - an empty defline, or one that starts with white space or has a control character or
//!   a non-ASCII byte (NCBI names an empty one `Query_1`, skips the white space, ends the
//!   title at the control character, and
//!   splits the ID and wraps the title by bytes, where `bio` splits by Unicode white space
//!   and the report wraps by characters): `check_deflines`, on the bytes of the file,
//!   because `bio` drops the character that ends the ID;
//! - a record without residues, or a residue that is not an IUPAC nucleotide letter (NCBI
//!   warns that the sequence contains no data, ignores white space and hyphens, ends the
//!   line at `;`, and removes other characters with a warning): `check_residues`.

use anyhow::{bail, Result};
use bio::io::fasta;

/// Why a FASTA file that `bio` cannot parse is not read: `bio` fails on text before the
/// first defline (blank lines, `;` comments, a byte order mark) and on bytes that are not
/// UTF-8, which NCBI reads.
pub const UNREADABLE_FASTA: &str =
    "FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-8), which NCBI BLAST+ may read, is not supported by LOSAT's BLASTN";

/// Whether a FASTA file has no character but white space: NCBI reads it as a file without
/// records (an empty query, or no subject).
///
/// NCBI reference: c++/src/app/blast/blast_app_util.cpp:862-866
/// ```c
/// 	IOS_BASE::iostate orig_state = in.rdstate();
/// 	IOS_BASE::fmtflags orig_flags = in.setf(ios::skipws);
///
/// 	if(! (in >> c))
/// 		return true;
/// ```
/// `in >> c` skips the characters of C's `isspace`.
pub fn is_blank(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .all(|byte| matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
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
            "{role} record {record} has a defline that {problem}; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTN (use ASCII deflines without control characters)"
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
    let mut record = 0;
    for line in bytes.split(|&byte| byte == b'\n') {
        if line.starts_with(b">") {
            record += 1;
        } else if record > 0 && !line.is_ascii() {
            bail!(
                "{role} record {record} has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's BLASTN (use IUPAC nucleotide letters)"
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

/// NCBI's warning for a record whose defline ends with nucleotides (`ends_with_nucleotides`),
/// written when NCBI reads the record: the subjects when they are read, the queries after
/// `Query is Empty!`. The deflines are those of `bio` (the ID, a space and the rest), which
/// are NCBI's text for the deflines that `check_deflines` accepts.
pub fn write_title_warnings(
    records: &[fasta::Record],
    diagnostics: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    for record in records {
        let text = match record.desc() {
            Some(desc) => format!("{} {desc}", record.id()),
            None => record.id().to_string(),
        };
        if ends_with_nucleotides(text.as_bytes()) {
            diagnostics.write_all(b"FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n")?;
        }
    }
    Ok(())
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
/// `role` is `query` or `subject`.
pub fn check_residues(records: &[fasta::Record], role: &str) -> Result<()> {
    first_problem(records, |index, record| {
        invalid_residue(index, record, role)
    })
}

/// Rejects a record without residues: NCBI reads it without a message and reports it when
/// it sets up the search ("Sequence contains no data"), which LOSAT does not reproduce.
/// `role` is `query` or `subject`.
pub fn check_records_have_residues(records: &[fasta::Record], role: &str) -> Result<()> {
    first_problem(records, |index, record| no_residues(index, record, role))
}

/// Both checks, record by record (the order of ABI v1, which plan TD-1 freezes).
pub fn check_records(records: &[fasta::Record], role: &str) -> Result<()> {
    first_problem(records, |index, record| {
        no_residues(index, record, role).or_else(|| invalid_residue(index, record, role))
    })
}

fn first_problem(
    records: &[fasta::Record],
    problem: impl Fn(usize, &fasta::Record) -> Option<String>,
) -> Result<()> {
    match records
        .iter()
        .enumerate()
        .find_map(|(index, record)| problem(index, record))
    {
        Some(message) => bail!("{message}"),
        None => Ok(()),
    }
}

fn no_residues(index: usize, record: &fasta::Record, role: &str) -> Option<String> {
    record.seq().is_empty().then(|| {
        format!(
            "{role} record {} ({}) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTN",
            index + 1,
            record.id()
        )
    })
}

fn invalid_residue(index: usize, record: &fasta::Record, role: &str) -> Option<String> {
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
        "{role} record {} ({}) has {shown} at residue {}, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not supported by LOSAT's BLASTN (use IUPAC nucleotide letters)",
        index + 1,
        record.id(),
        position + 1
    ))
}

/// The records with `U` read as `T` (and `u` as `t`), or `None` when none has a `U`.
///
/// NCBI reference: c++/src/objtools/readers/fasta.cpp:856-873
/// ```c
///         case 'A': case 'B': case 'C': case 'D':
///         case 'G': case 'H':
///         case 'K':
///         case 'M':
///         case 'R': case 'S': case 'T': case 'U': case 'V': case 'W':
///         case 'Y':
///             CloseGap(pos == 0);
///             m_SeqData[m_CurrentPos] = c;
///     ...
///         case 'r': case 's': case 't': case 'u': case 'v': case 'w':
///         case 'y':
///             char_type = eCharType_MaskedNonGap;
/// ```
/// `CFastaReader` keeps `U` as a residue; NCBI BLAST+ 2.17.0 searches it as `T` and shows
/// it as `T`, with or without other `T` and in either case (the oracle runs of
/// docs/evidence/losat_web_e2c/).
pub fn with_u_as_t(records: &[fasta::Record]) -> Option<Vec<fasta::Record>> {
    if !records
        .iter()
        .any(|record| record.seq().iter().any(|&byte| matches!(byte, b'U' | b'u')))
    {
        return None;
    }
    Some(
        records
            .iter()
            .map(|record| {
                let seq: Vec<u8> = record
                    .seq()
                    .iter()
                    .map(|&byte| match byte {
                        b'U' => b'T',
                        b'u' => b't',
                        other => other,
                    })
                    .collect();
                fasta::Record::with_attrs(record.id(), record.desc(), &seq)
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!(check_residues(&records(text), "query").is_ok());
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
            let error = check_residues(&records(text), "query")
                .unwrap_err()
                .to_string();
            assert!(error.contains(problem), "{text:?}: {error}");
            assert!(error.contains("not supported by LOSAT"), "{error}");
        }
        let error = check_records_have_residues(&records(">q1 empty\n>q2\nACGT\n"), "query")
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
        assert!(check_records_have_residues(&both, "query")
            .unwrap_err()
            .to_string()
            .contains("query record 2 (q2) has no residues"));
    }

    #[test]
    fn titles_that_end_with_nucleotides_get_ncbi_warning() {
        // NCBI 2.17.0: the text after `>` longer than 20 bytes, its last 20 ACGT.
        let text = ">q1 ACGTACGTACGTACGTACGTA\nACGT\n>ACGTACGTACGTACGTACGT\nACGT\n>ACGTACGTACGTACGTACGTA\nACGT\n>q2 ACGTNACGTACGTACGTACGTA\nACGT\n";
        assert!(check_deflines(text.as_bytes(), "query").is_ok());
        assert!(check_deflines(b">q1 ACGTACGTACGTACGTACGTA\r\nACGT\r\n", "query").is_ok());
        let mut out = Vec::new();
        write_title_warnings(&records(text), &mut out).unwrap();
        let warning = "FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n";
        assert_eq!(String::from_utf8(out).unwrap(), warning.repeat(2));
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

    #[test]
    fn u_is_read_as_t() {
        assert!(with_u_as_t(&records(">q\nACGT\n")).is_none());
        let read = with_u_as_t(&records(">q d\nACGUu\n")).unwrap();
        assert_eq!(read[0].seq(), b"ACGTt");
        assert_eq!(read[0].desc(), Some("d"));
    }
}
