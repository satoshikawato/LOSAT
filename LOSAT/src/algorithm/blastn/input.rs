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
        let defline = defline.trim_ascii_end();
        let problem = if defline.is_empty() {
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
    for (index, record) in records.iter().enumerate() {
        if record.seq().is_empty() {
            bail!(
                "{role} record {} ({}) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTN",
                index + 1,
                record.id()
            );
        }
        if let Some(position) = record
            .seq()
            .iter()
            .position(|&byte| !IUPAC_NUCLEOTIDE[byte as usize])
        {
            let byte = record.seq()[position];
            let shown = if byte.is_ascii_graphic() {
                format!("'{}'", byte as char)
            } else {
                format!("0x{byte:02x}")
            };
            bail!(
                "{role} record {} ({}) has {shown} at residue {}, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not supported by LOSAT's BLASTN (use IUPAC nucleotide letters)",
                index + 1,
                record.id(),
                position + 1
            );
        }
    }
    Ok(())
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
            (
                ">q1 empty\n>q2\nACGT\n",
                "query record 1 (q1) has no residues",
            ),
        ] {
            let error = check_residues(&records(text), "query")
                .unwrap_err()
                .to_string();
            assert!(error.contains(problem), "{text:?}: {error}");
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
