//! The query sets of probe `F/f2/tools/qprobe.py` (NCBI BLAST+ 2.17.0 on ws, 2026-11-20): each
//! case's Bioseq-set as NCBI's strategy held it (`<case>.asn`), and what NCBI did with it: the
//! FASTA text of its temporary file (`<case>.fasta`) or the exception it reported
//! (`<case>.exception`, the report's lines in the probe build's tree). `cases.tsv` says which
//! cases LOSAT refuses and why.

use super::*;
use crate::report::serial::asn_types::{read_document, table::named};

const FIXTURES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/strategy_queries"
);

/// The tree that the probe build's exception reports name.
const PROBE_BUILD_ROOT: &str = "/opt/conda/conda-bld/blast_1787399212531/work/c++";

fn read(name: &str) -> Option<Vec<u8>> {
    std::fs::read(format!("{FIXTURES}/{name}")).ok()
}

/// The text of the query set `text` held whole, after the check that the text written to a file a
/// block at a time is the same, or the same error.
fn render(text: &[u8]) -> Result<Vec<u8>, RenderError> {
    let (value, posts) = read_document(text, named::BIOSEQ_SET);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.unwrap_or_else(|error| panic!("{error:?}"));
    let set = V::new(named::BIOSEQ_SET, &value);
    let held = render_queries(set, Text::held()).map(bytes);
    let file = crate::blastinput::input_files::temporary_output("query").expect("a file");
    let written = render_queries(set, Text::file(file)).map(bytes);
    assert_eq!(format!("{written:?}"), format!("{held:?}"));
    held
}

/// The bytes of a text written.
fn bytes(rendered: Rendered) -> Vec<u8> {
    match rendered {
        Rendered::Bytes(bytes) => bytes,
        Rendered::File(mut file) => {
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(&mut file, &mut bytes).expect("the file read");
            bytes
        }
    }
}

#[test]
fn every_probed_query_set_is_written_as_ncbi_wrote_it() {
    let manifest = String::from_utf8(read("cases.tsv").expect("cases.tsv")).unwrap();
    let mut checked = 0;
    for line in manifest.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let case = fields.next().unwrap();
        let losat = fields.next().unwrap();
        let text = read(&format!("{case}.asn")).unwrap_or_else(|| panic!("{case}.asn"));
        let got = render(&text);
        match (
            losat,
            read(&format!("{case}.fasta")),
            read(&format!("{case}.exception")),
        ) {
            ("refused", _, _) => match got {
                Err(RenderError::Unsupported(part)) => assert!(!part.is_empty(), "{case}"),
                other => panic!("{case}: {other:?}"),
            },
            ("same", Some(fasta), None) => match got {
                Ok(bytes) => assert!(
                    bytes == fasta,
                    "{case}:\n{}\n{}",
                    String::from_utf8_lossy(&bytes),
                    String::from_utf8_lossy(&fasta)
                ),
                Err(error) => panic!("{case}: {error:?}"),
            },
            ("same", None, Some(exception)) => match got {
                Err(RenderError::Exception(e)) => assert_eq!(
                    String::from_utf8_lossy(&e.what_in(PROBE_BUILD_ROOT)),
                    crate::wording::losat_text(&String::from_utf8(exception).unwrap()),
                    "{case}"
                ),
                other => panic!("{case}: {other:?}"),
            },
            _ => panic!("{case}: no expectation"),
        }
        checked += 1;
    }
    assert!(checked >= 100, "{checked} cases");
}

#[test]
fn fills_of_a_long_seq_data_end_at_1024_letter_steps() {
    // `reference out of range` names the end of the first fill past the data: 2000 letters of
    // which 1500 are held end their fills at 1024 and 2000; 3000 at 1024, 2048.
    for (length, end) in [(2000, 2000), (3000, 2048), (1501, 1501)] {
        let text = format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr raw, mol na, length {length}, seq-data iupacna \"{}\" }} }} }} }}",
            "ACGT".repeat(375)
        );
        match render(text.as_bytes()) {
            Err(RenderError::Exception(e)) => assert_eq!(
                e.get_msg(),
                format!("reference out of range of Seq-inst data: {end}").as_bytes()
            ),
            other => panic!("{length}: {other:?}"),
        }
    }
}

#[test]
fn a_segment_that_passes_the_length_is_filled_up_to_it() {
    // A literal of 3000 letters that holds 1500, in a Bioseq of `length` letters: its fills end at
    // 1024 and at `length` (`dl_cut_short`); only a fill past the 1500 letters throws.
    let text = |length: u32| {
        format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr delta, mol na, length {length}, ext delta {{ literal {{ length 3000, seq-data iupacna \"{}\" }} }} }} }} }} }}",
            "ACGT".repeat(375)
        )
    };
    for length in [1024, 1500] {
        let fasta = render(text(length).as_bytes()).unwrap();
        let letters = fasta.iter().filter(|b| b"ACGT".contains(b)).count();
        assert_eq!(letters, length as usize);
    }
    match render(text(1600).as_bytes()) {
        Err(RenderError::Exception(e)) => {
            assert_eq!(
                e.get_msg(),
                b"reference out of range of Seq-inst data: 1600"
            )
        }
        other => panic!("{other:?}"),
    }
}

/// Sequences of more letters than the blocks in which they are converted and written are the
/// lines of the letters written whole: a gap, letters copied and letters converted, held or
/// written to a file (`render` checks the two).
#[test]
fn long_sequences_are_written_in_blocks_as_whole() {
    let gap = 3 * LETTER_BLOCK + 5;
    let packed = FILE_BLOCK + 7;
    let copied = 2 * LETTER_BLOCK + 3;
    let set = format!(
        "Bioseq-set ::= {{ seq-set {{ \
         seq {{ id {{ local str \"Query_1\" }}, descr {{ title \"v\" }}, inst {{ repr virtual, mol na, length {gap} }} }}, \
         seq {{ id {{ local str \"Query_2\" }}, descr {{ title \"p\" }}, inst {{ repr raw, mol na, length {}, seq-data ncbi2na '{}'H }} }}, \
         seq {{ id {{ local str \"Query_3\" }}, descr {{ title \"c\" }}, inst {{ repr raw, mol aa, length {copied}, seq-data ncbieaa \"{}\" }} }} }} }}",
        4 * packed,
        "1B".repeat(packed),
        "M".repeat(copied),
    );
    let text = render(set.as_bytes()).unwrap_or_else(|error| panic!("{error:?}"));
    let mut expected = b">lcl|Query_1 v\n".to_vec();
    crate::sequence::write_fasta_lines(&vec![b'N'; gap], WIDTH, &mut expected);
    expected.extend_from_slice(b">lcl|Query_2 p\n");
    crate::sequence::write_fasta_lines(&b"ACGT".repeat(packed), WIDTH, &mut expected);
    expected.extend_from_slice(b">lcl|Query_3 c\n");
    crate::sequence::write_fasta_lines(&vec![b'M'; copied], WIDTH, &mut expected);
    assert!(expected.len() > 4 * FILE_BLOCK);
    assert!(text == expected, "the long sequences' lines");
}
