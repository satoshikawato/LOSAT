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
    let mut title = Vec::with_capacity(record.id().len() + 1);
    title.extend_from_slice(record.id().as_bytes());
    if let Some(desc) = record.desc() {
        title.push(b' ');
        title.extend_from_slice(desc.as_bytes());
    }
    // `isspace` of the C locale.
    let start = title
        .iter()
        .position(|&byte| !matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
        .unwrap_or(title.len());
    title.drain(..start);
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
}
